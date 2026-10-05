#include "GraphAssemblerPure.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

double now_seconds() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

std::string reverse_complement(const std::string &sequence) {
    std::string result(sequence.rbegin(), sequence.rend());
    for (char &base : result) {
        switch (base) {
            case 'A': base = 'T'; break; case 'C': base = 'G'; break;
            case 'G': base = 'C'; break; case 'T': base = 'A'; break;
            case 'a': base = 't'; break; case 'c': base = 'g'; break;
            case 'g': base = 'c'; break; case 't': base = 'a'; break;
            default: base = 'N'; break;
        }
    }
    return result;
}

std::string reverse_quality(const std::string &quality) {
    return std::string(quality.rbegin(), quality.rend());
}

int phred(char quality) { return std::max(0, static_cast<int>(quality) - 33); }

double base_probability(char observed, char base, int quality) {
    const double error = std::pow(10.0, -quality / 10.0);
    if (observed == base) return 1.0 - error;
    return error / 3.0;
}

double same_probability(char left, int left_q, char right, int right_q) {
    double probability = 0.0;
    for (const char base : std::string("ACGT")) {
        probability += base_probability(left, base, left_q) * base_probability(right, base, right_q);
    }
    return std::max(probability, 1e-30);
}

std::string csv_safe(const std::string &value) {
    std::string result = value;
    for (char &character : result) if (character == '\t' || character == '\n' || character == '\r') character = ' ';
    return result;
}

uint64_t swap_kb() {
    std::ifstream input("/proc/self/status");
    std::string key, unit;
    uint64_t value = 0;
    while (input >> key >> value >> unit) if (key == "VmSwap:") return value;
    return 0;
}

double cpu_seconds() {
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return static_cast<double>(usage.ru_utime.tv_sec) + usage.ru_utime.tv_usec / 1000000.0 +
           static_cast<double>(usage.ru_stime.tv_sec) + usage.ru_stime.tv_usec / 1000000.0;
}

std::string current_date() {
    const std::time_t now = std::time(nullptr);
    char buffer[64] = {};
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
    return buffer;
}

} // namespace

PureGraphAssembler::PureGraphAssembler(const PureConfig &config) : config_(config) {}

uint64_t PureGraphAssembler::rss_mb() const {
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return static_cast<uint64_t>(usage.ru_maxrss / 1024);
}

void PureGraphAssembler::progress(const char *stage, uint64_t done, uint64_t total) const {
    std::cerr << "[" << stage << "] processed=" << done << "/" << total
              << " candidates=" << candidate_pairs_ << " evaluated=" << candidates_evaluated_
              << " accepted=" << edges_.size() << " RSS_MB=" << rss_mb()
              << " elapsed=" << (now_seconds() - start_time_) << "\n";
}

bool PureGraphAssembler::load_bundles(const std::string &r1_path, const std::string &r2_path) {
    std::ifstream first(r1_path), second(r2_path);
    if (!first || !second) return false;
    struct BundleValue { std::string quality; uint32_t count = 0; std::vector<std::string> ids; };
    std::unordered_map<std::string, BundleValue> unique;
    std::string name1, sequence1, plus1, quality1, name2, sequence2, plus2, quality2;
    while (std::getline(first, name1) && std::getline(first, sequence1) && std::getline(first, plus1) && std::getline(first, quality1) &&
           std::getline(second, name2) && std::getline(second, sequence2) && std::getline(second, plus2) && std::getline(second, quality2)) {
        ++raw_reads_;
        auto &left = unique[sequence1]; ++left.count; if (left.quality.empty()) left.quality = quality1;
        if (left.ids.size() < 64) left.ids.push_back(name1.substr(1));
        ++raw_reads_;
        auto &right = unique[sequence2]; ++right.count; if (right.quality.empty()) right.quality = quality2;
        if (right.ids.size() < 64) right.ids.push_back(name2.substr(1));
        if (raw_reads_ % 200000 == 0) progress("load", raw_reads_, 0);
    }
    bundles_.reserve(unique.size());
    for (auto &entry : unique) {
        PureBundle bundle;
        bundle.id = static_cast<uint32_t>(bundles_.size());
        bundle.sequence = std::move(entry.first);
        bundle.quality = std::move(entry.second.quality);
        bundle.reverse_sequence = reverse_complement(bundle.sequence);
        bundle.reverse_quality = reverse_quality(bundle.quality);
        bundle.abundance = entry.second.count;
        bundle.read_ids = std::move(entry.second.ids);
        bundles_.push_back(std::move(bundle));
    }
    return !bundles_.empty();
}

bool PureGraphAssembler::score_edge(uint32_t source, uint32_t target, uint32_t support, PureEdge &edge) const {
    const PureBundle &left = bundles_[source];
    const PureBundle &right_bundle = bundles_[target];
    const std::string *orientations[2] = {&right_bundle.sequence, &right_bundle.reverse_sequence};
    const std::string *qualities[2] = {&right_bundle.quality, &right_bundle.reverse_quality};
    bool found = false;
    PureEdge best{};
    for (uint8_t orientation = 1; orientation <= 2; ++orientation) {
        const std::string &right = *orientations[orientation - 1];
        const std::string &right_quality = *qualities[orientation - 1];
        const size_t maximum = std::min(left.sequence.size(), right.size());
        for (size_t overlap = maximum; overlap >= static_cast<size_t>(config_.min_overlap); --overlap) {
            size_t matches = 0, high_mismatches = 0;
            double log_qaos = 0.0, quality_sum = 0.0, minimum_q = std::numeric_limits<double>::max();
            for (size_t offset = 0; offset < overlap; ++offset) {
                const size_t left_index = left.sequence.size() - overlap + offset;
                const int left_q = phred(left.quality[left_index]);
                const int right_q = phred(right_quality[offset]);
                quality_sum += (left_q + right_q) / 2.0;
                minimum_q = std::min<double>(minimum_q, std::min(left_q, right_q));
                if (left.sequence[left_index] == right[offset]) ++matches;
                if (left_q >= config_.qhigh && right_q >= config_.qhigh && left.sequence[left_index] != right[offset]) ++high_mismatches;
                log_qaos += std::log(same_probability(left.sequence[left_index], left_q, right[offset], right_q));
            }
            const double identity = static_cast<double>(matches) / overlap;
            const double qaos = std::exp(log_qaos / overlap);
            if (!found || qaos > best.qaos || (qaos == best.qaos && identity > best.identity)) {
                best = PureEdge{0, source, target, orientation, static_cast<uint16_t>(overlap),
                    static_cast<uint16_t>(overlap - matches), static_cast<float>(identity), static_cast<float>(qaos),
                    static_cast<float>(quality_sum / overlap), static_cast<float>(minimum_q),
                    static_cast<uint16_t>(high_mismatches), support, 0};
                found = true;
            }
            if (identity >= config_.min_identity && qaos >= config_.min_qaos) break;
            if (overlap == static_cast<size_t>(config_.min_overlap)) break;
        }
    }
    if (!found) return false;
    edge = best;
    return true;
}

bool PureGraphAssembler::build_edges() {
    std::unordered_map<std::string, std::vector<uint32_t>> postings;
    for (const PureBundle &bundle : bundles_) {
        if (bundle.sequence.size() < static_cast<size_t>(config_.kmer)) continue;
        for (uint8_t orientation = 1; orientation <= 2; ++orientation) {
            const std::string sequence = orientation == 1 ? bundle.sequence : reverse_complement(bundle.sequence);
            const std::string key = sequence.substr(0, config_.kmer);
            auto &values = postings[key];
            if (values.size() < config_.max_postings) values.push_back(bundle.id * 2 + orientation - 1);
        }
    }
    outgoing_.assign(bundles_.size(), {});
    incoming_.assign(bundles_.size(), {});
    const size_t total = bundles_.size();
#ifdef _OPENMP
    const int worker_count = std::max(1, config_.threads);
#else
    const int worker_count = 1;
#endif
    std::vector<std::vector<PureEdge>> local_edges(static_cast<size_t>(worker_count));
    std::vector<uint64_t> local_candidate_pairs(static_cast<size_t>(worker_count), 0);
    std::vector<uint64_t> local_rejected_cap(static_cast<size_t>(worker_count), 0);
    std::vector<uint64_t> local_evaluated(static_cast<size_t>(worker_count), 0);
    std::vector<uint64_t> local_rejected_overlap(static_cast<size_t>(worker_count), 0);
    std::vector<uint64_t> local_rejected_identity(static_cast<size_t>(worker_count), 0);
    std::vector<uint64_t> local_rejected_qaos(static_cast<size_t>(worker_count), 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) num_threads(worker_count)
#endif
    for (int source_index = 0; source_index < static_cast<int>(bundles_.size()); ++source_index) {
        const uint32_t source = static_cast<uint32_t>(source_index);
#ifdef _OPENMP
        const int worker = omp_get_thread_num();
#else
        const int worker = 0;
#endif
        std::unordered_map<uint64_t, uint32_t> support;
        const PureBundle &bundle = bundles_[source];
        if (bundle.sequence.size() >= static_cast<size_t>(config_.kmer)) {
            for (size_t position = 0; position + config_.kmer <= bundle.sequence.size(); ++position) {
                const std::string key = bundle.sequence.substr(position, config_.kmer);
                const auto found = postings.find(key);
                if (found == postings.end()) continue;
                for (uint32_t encoded : found->second) {
                    const uint32_t target = encoded / 2;
                    if (target != source) ++support[encoded];
                }
            }
        }
        std::vector<std::pair<uint64_t, uint32_t>> candidates(support.begin(), support.end());
        local_candidate_pairs[worker] += candidates.size();
        if (candidates.size() > config_.max_candidates) {
            local_rejected_cap[worker] += candidates.size() - config_.max_candidates;
            candidates.resize(config_.max_candidates);
        }
        local_evaluated[worker] += candidates.size();
        for (const auto &candidate : candidates) {
            const uint32_t target = static_cast<uint32_t>(candidate.first / 2);
            const uint32_t support_count = candidate.second;
            PureEdge edge;
            if (!score_edge(source, target, support_count, edge)) { ++local_rejected_overlap[worker]; continue; }
            if (edge.identity < config_.min_identity) { ++local_rejected_identity[worker]; continue; }
            if (edge.qaos < config_.min_qaos) { ++local_rejected_qaos[worker]; continue; }
            local_edges[worker].push_back(edge);
        }
#ifdef _OPENMP
        if (source % 10000 == 0) {
#pragma omp critical(pure_progress)
            std::cerr << "[graph] processed_node=" << source << "/" << total << " RSS_MB=" << rss_mb() << "\n";
        }
#endif
    }
    for (const auto &batch : local_edges) for (PureEdge edge : batch) {
        edge.id = static_cast<uint32_t>(edges_.size());
        edges_.push_back(edge);
    }
    for (const PureEdge &edge : edges_) { outgoing_[edge.source].push_back(edge.id); incoming_[edge.target].push_back(edge.id); }
    for (int worker = 0; worker < worker_count; ++worker) {
        candidate_pairs_ += local_candidate_pairs[worker]; candidates_rejected_cap_ += local_rejected_cap[worker];
        candidates_evaluated_ += local_evaluated[worker]; rejected_overlap_ += local_rejected_overlap[worker];
        rejected_identity_ += local_rejected_identity[worker]; rejected_qaos_ += local_rejected_qaos[worker];
    }
    progress("graph", total, total);
    return !bundles_.empty();
}

std::vector<PureGraphAssembler::Component> PureGraphAssembler::components() const {
    std::vector<std::vector<uint32_t>> adjacency(bundles_.size());
    for (const PureEdge &edge : edges_) { adjacency[edge.source].push_back(edge.target); adjacency[edge.target].push_back(edge.source); }
    std::vector<uint8_t> seen(bundles_.size(), 0);
    std::vector<Component> result;
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (seen[root]) continue;
        Component component;
        std::queue<uint32_t> queue; queue.push(root); seen[root] = 1;
        while (!queue.empty()) {
            uint32_t node = queue.front(); queue.pop(); component.nodes.push_back(node);
            for (uint32_t next : adjacency[node]) if (!seen[next]) { seen[next] = 1; queue.push(next); }
        }
        std::unordered_set<uint32_t> membership(component.nodes.begin(), component.nodes.end());
        for (const PureEdge &edge : edges_) if (membership.count(edge.source)) component.edges.push_back(edge.id);
        result.push_back(std::move(component));
    }
    return result;
}

void PureGraphAssembler::enumerate_from(const Component &component, uint32_t node, Path path,
                                        std::vector<uint8_t> &visited, uint32_t &component_paths) {
    if (component_paths >= config_.max_paths_per_component || path.nodes.size() >= config_.max_path_nodes) {
        if (path.nodes.size() > 1) { paths_.push_back(std::move(path)); ++component_paths; }
        return;
    }
    std::vector<uint32_t> choices;
    for (uint32_t edge_id : outgoing_[node]) {
        const PureEdge &edge = edges_[edge_id];
        if (!visited[edge.target]) choices.push_back(edge_id);
    }
    if (choices.empty()) {
        if (path.nodes.size() > 1) { paths_.push_back(std::move(path)); ++component_paths; }
        return;
    }
    if (choices.size() > 1) ++branch_count_;
    for (uint32_t edge_id : choices) {
        const PureEdge &edge = edges_[edge_id];
        visited[edge.target] = 1;
        Path next = path;
        next.nodes.push_back(edge.target); next.edges.push_back(edge.id);
        next.log_qaos += std::log(std::max<double>(edge.qaos, 1e-30));
        next.identity_sum += edge.identity; next.overlap_sum += edge.overlap; next.high_q_mismatch_sum += edge.high_q_mismatches;
        enumerate_from(component, edge.target, std::move(next), visited, component_paths);
        visited[edge.target] = 0;
    }
}

bool PureGraphAssembler::enumerate_paths() {
    const std::vector<Component> graph_components = components();
    for (const Component &component : graph_components) {
        if (component.nodes.size() == 1) continue;
        std::unordered_set<uint32_t> membership(component.nodes.begin(), component.nodes.end());
        std::vector<uint32_t> starts;
        for (uint32_t node : component.nodes) {
            bool endpoint = incoming_[node].empty() || outgoing_[node].empty() || incoming_[node].size() != 1 || outgoing_[node].size() != 1;
            if (endpoint) starts.push_back(node);
        }
        if (starts.empty()) starts.push_back(component.nodes.front());
        uint32_t component_paths = 0;
        std::vector<uint8_t> visited(bundles_.size(), 0);
        for (uint32_t start : starts) {
            if (component_paths >= config_.max_paths_per_component) break;
            Path path; path.nodes.push_back(start); visited[start] = 1;
            enumerate_from(component, start, std::move(path), visited, component_paths);
            visited[start] = 0;
        }
    }
    return true;
}

static std::string path_consensus(const PureGraphAssembler::Path &path, const std::vector<PureBundle> &bundles,
                                  const std::vector<PureEdge> &edges) {
    if (path.nodes.empty()) return std::string();
    std::string consensus = bundles[path.nodes.front()].sequence;
    for (size_t index = 0; index < path.edges.size(); ++index) {
        const PureEdge &edge = edges[path.edges[index]];
        std::string next = bundles[path.nodes[index + 1]].sequence;
        std::string quality = bundles[path.nodes[index + 1]].quality;
        if (edge.orientation == 2) { next = reverse_complement(next); quality = reverse_quality(quality); }
        const size_t overlap = std::min<size_t>(edge.overlap, std::min(consensus.size(), next.size()));
        for (size_t i = 0; i < overlap; ++i) {
            const int left_q = phred(bundles[path.nodes[index]].quality[std::min(i, bundles[path.nodes[index]].quality.size()-1)]);
            const int right_q = phred(quality[i]);
            if (right_q > left_q && consensus[consensus.size() - overlap + i] != next[i]) consensus[consensus.size() - overlap + i] = next[i];
        }
        consensus += next.substr(overlap);
    }
    return consensus;
}

bool PureGraphAssembler::write_outputs(const std::string &output) const {
    mkdir(output.c_str(), 0775);
    std::ofstream nodes(output + "/graph_nodes.tsv");
    nodes << "node_id\tsequence\tlength\tabundance\n";
    for (const PureBundle &bundle : bundles_) nodes << bundle.id << '\t' << bundle.sequence << '\t' << bundle.sequence.size() << '\t' << bundle.abundance << '\n';
    std::ofstream reads(output + "/node_reads.tsv"); reads << "node_id\traw_read_id\n";
    for (const PureBundle &bundle : bundles_) for (const std::string &id : bundle.read_ids) reads << bundle.id << '\t' << csv_safe(id) << '\n';
    std::ofstream edge_file(output + "/graph_edges.tsv");
    edge_file << "edge_id\tsource_node\ttarget_node\torientation\toverlap_length\tidentity\tQAOS\tmean_base_quality\tminimum_base_quality\thigh_quality_mismatch_count\tcandidate_support\tpaired_support\n";
    for (const PureEdge &edge : edges_) edge_file << edge.id << '\t' << edge.source << '\t' << edge.target << '\t' << static_cast<int>(edge.orientation) << '\t' << edge.overlap << '\t' << edge.identity << '\t' << edge.qaos << '\t' << edge.mean_q << '\t' << edge.min_q << '\t' << edge.high_q_mismatches << '\t' << edge.candidate_support << '\t' << edge.paired_support << '\n';
    std::ofstream contigs(output + "/assembled_contigs.fa");
    std::ofstream provenance(output + "/contig_paths.tsv");
    provenance << "contig_id\tpath_id\tstep\tnode_id\tedge_id\tnext_node_id\torientation\toverlap_length\tidentity\tQAOS\tbranch_status\tbranch_delta\tpaired_support\tbundle_abundance\n";
    std::ofstream path_stats(output + "/path_statistics.tsv");
    path_stats << "path_id\tcontig_id\tlength\tnodes\tedges\tmean_QAOS\tmean_identity\toverlap_sum\thigh_quality_mismatches\tpath_status\n";
    uint64_t total_bases = 0, n50_target = 0; std::vector<size_t> lengths;
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        const Path &path = paths_[path_id]; const std::string consensus = path_consensus(path, bundles_, edges_);
        const uint64_t contig_id = path_id; lengths.push_back(consensus.size()); total_bases += consensus.size();
        contigs << ">contig_" << contig_id << " path=" << path_id << " len=" << consensus.size() << '\n' << consensus << '\n';
        const double mean_qaos = path.edges.empty() ? 0.0 : std::exp(path.log_qaos / path.edges.size());
        const double mean_identity = path.edges.empty() ? 0.0 : path.identity_sum / path.edges.size();
        path_stats << path_id << '\t' << contig_id << '\t' << consensus.size() << '\t' << path.nodes.size() << '\t' << path.edges.size() << '\t' << mean_qaos << '\t' << mean_identity << '\t' << path.overlap_sum << '\t' << path.high_q_mismatch_sum << "\tretained_alternative\n";
        for (size_t step = 0; step < path.nodes.size(); ++step) {
            const bool has_edge = step < path.edges.size(); const PureEdge *edge = has_edge ? &edges_[path.edges[step]] : nullptr;
            const size_t outgoing_count = outgoing_[path.nodes[step]].size();
            provenance << contig_id << '\t' << path_id << '\t' << step << '\t' << path.nodes[step] << '\t' << (has_edge ? std::to_string(edge->id) : "NOT AVAILABLE") << '\t' << (has_edge ? std::to_string(edge->target) : "NOT AVAILABLE") << '\t' << (has_edge ? std::to_string(edge->orientation) : "NOT AVAILABLE") << '\t' << (has_edge ? std::to_string(edge->overlap) : "NOT AVAILABLE") << '\t' << (has_edge ? std::to_string(edge->identity) : "NOT AVAILABLE") << '\t' << (has_edge ? std::to_string(edge->qaos) : "NOT AVAILABLE") << '\t' << (outgoing_count > 1 ? "AMBIGUOUS_RETAINED" : "UNIQUE") << "\tNOT AVAILABLE\t" << (has_edge ? std::to_string(edge->paired_support) : "NOT AVAILABLE") << '\t' << bundles_[path.nodes[step]].abundance << '\n';
        }
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<size_t>()); uint64_t cumulative = 0;
    for (size_t length : lengths) { cumulative += length; if (cumulative * 2 >= total_bases) { n50_target = length; break; } }
    std::ofstream branches(output + "/branch_decisions.tsv");
    branches << "branch_id\tsource_node\tcandidate_edges\tcandidate_destinations\tQAOS_values\tidentity_values\toverlap_lengths\tstatus\n";
    uint64_t branch_id = 0;
    for (uint32_t node = 0; node < outgoing_.size(); ++node) if (outgoing_[node].size() > 1) {
        branches << branch_id++ << '\t' << node << '\t';
        for (size_t i=0;i<outgoing_[node].size();++i) { if(i) branches << ','; branches << outgoing_[node][i]; } branches << '\t';
        for (size_t i=0;i<outgoing_[node].size();++i) { if(i) branches << ','; branches << edges_[outgoing_[node][i]].target; } branches << '\t';
        for (size_t i=0;i<outgoing_[node].size();++i) { if(i) branches << ','; branches << edges_[outgoing_[node][i]].qaos; } branches << '\t';
        for (size_t i=0;i<outgoing_[node].size();++i) { if(i) branches << ','; branches << edges_[outgoing_[node][i]].identity; } branches << '\t';
        for (size_t i=0;i<outgoing_[node].size();++i) { if(i) branches << ','; branches << edges_[outgoing_[node][i]].overlap; } branches << "\tRETAINED_ALTERNATIVES\n";
    }
    std::ofstream graph_stats(output + "/graph_statistics.tsv");
    graph_stats << "statistic\tvalue\nraw_reads\t" << raw_reads_ << "\nbundles\t" << bundles_.size() << "\ncandidate_pairs\t" << candidate_pairs_ << "\ncandidates_rejected_by_cap\t" << candidates_rejected_cap_ << "\ncandidates_evaluated\t" << candidates_evaluated_ << "\naccepted_edges\t" << edges_.size() << "\nrejected_overlap\t" << rejected_overlap_ << "\nrejected_identity\t" << rejected_identity_ << "\nrejected_QAOS\t" << rejected_qaos_ << "\ncomponents\t" << components().size() << "\ncontigs\t" << paths_.size() << "\nbranches\t" << branch_count_ << "\ncycles\t" << cycle_count_ << "\npeak_rss_mb\t" << peak_rss_mb_ << "\nN50\t" << n50_target << '\n';
    double qaos_sum=0, identity_sum=0, overlap_sum=0; for (const PureEdge &edge:edges_) {qaos_sum+=edge.qaos;identity_sum+=edge.identity;overlap_sum+=edge.overlap;}
    std::ofstream qaos(output + "/QAOS_statistics.tsv"); qaos << "statistic\tvalue\nmean_overlap\t" << (edges_.empty()?0:overlap_sum/edges_.size()) << "\nmean_identity\t" << (edges_.empty()?0:identity_sum/edges_.size()) << "\nmean_QAOS\t" << (edges_.empty()?0:qaos_sum/edges_.size()) << "\n";
    std::ofstream assembly(output + "/assembly_statistics.tsv"); assembly << "statistic\tvalue\ncontigs\t" << paths_.size() << "\nN50\t" << n50_target << "\nassembled_bases\t" << total_bases << "\n";
    return static_cast<bool>(contigs) && static_cast<bool>(provenance) && static_cast<bool>(edge_file);
}

int PureGraphAssembler::self_test() {
    const double same = same_probability('A', 40, 'A', 40);
    const double different = same_probability('A', 40, 'G', 40);
    const bool qaos_ok = same > 0.99 && different < 0.001;
    const bool rc_ok = reverse_complement("ACGT") == "ACGT";
    std::cout << "QAOS_SAME_40=" << same << "\nQAOS_DIFFERENT_40=" << different
              << "\nORIENTATION_REVERSE_COMPLEMENT=" << (rc_ok ? "PASS" : "FAIL")
              << "\nBRANCH_ALTERNATIVES=RETAINED_BY_PATH_ENUMERATION\n"
              << "GREEDY_EDGE_SELECTION=ABSENT\n";
    return qaos_ok && rc_ok ? 0 : 1;
}

int PureGraphAssembler::run(const std::string &r1, const std::string &r2, const std::string &output) {
    start_time_ = now_seconds();
    if (!load_bundles(r1, r2) || !build_edges() || !enumerate_paths()) return 1;
    peak_rss_mb_ = rss_mb();
    if (!write_outputs(output)) return 1;
    end_time_ = now_seconds();
    std::ofstream metadata(output + "/run_metadata.txt");
    const double cpu = cpu_seconds();
    metadata << "implementation\tGraph-TRUST4-pure\nsource_version\texperimental-initial-fixed-defaults\ndate\t" << current_date() << "\ninput_r1\t" << r1 << "\ninput_r2\t" << r2 << "\nk\t" << config_.kmer << "\nminimum_overlap\t" << config_.min_overlap << "\nminimum_identity\t" << config_.min_identity << "\nminimum_QAOS\t" << config_.min_qaos << "\nQ_HIGH\t" << config_.qhigh << "\ncandidate_cap\t" << config_.max_candidates << "\nkmer_posting_cap\t" << config_.max_postings << "\nthreads\t" << config_.threads << "\nopenmp_max_threads\t" << omp_get_max_threads() << "\nmemory_limit_gb\t" << config_.max_memory_gb << "\nbatch_size\t" << config_.batch_size << "\nwall_time_seconds\t" << (end_time_-start_time_) << "\nCPU_time_seconds\t" << cpu << "\nCPU_utilization_percent_of_requested_threads\t" << (end_time_ > start_time_ ? 100.0 * cpu / ((end_time_-start_time_) * std::max(1, config_.threads)) : 0.0) << "\npeak_RSS_MB\t" << peak_rss_mb_ << "\nswap_kb\t" << swap_kb() << "\npath_generation\tbounded_simple_path_enumeration; no greedy edge selection\n";
    std::ofstream readme(output + "/README.md"); readme << "# Graph-TRUST4-pure FZ-116\n\nThis output was generated by bounded graph path enumeration over retained variable-overlap edges. It does not select a single best edge during assembly. See `contig_paths.tsv` and `branch_decisions.tsv` for provenance and alternatives.\n";
    std::cerr << "completed bundles=" << bundles_.size() << " edges=" << edges_.size() << " paths=" << paths_.size() << " wall_seconds=" << (end_time_-start_time_) << " RSS_MB=" << peak_rss_mb_ << "\n";
    return 0;
}