#include "clean_graph.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;
double seconds() { return std::chrono::duration<double>(Clock::now().time_since_epoch()).count(); }
int phred(char value) { return std::max(0, static_cast<int>(value) - 33); }

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

uint64_t hash_kmer(const std::string &sequence, size_t start, size_t k) {
    uint64_t hash = 1469598103934665603ULL;
    for (size_t i = start; i < start + k; ++i) {
        const char base = sequence[i];
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') return 0;
        hash ^= static_cast<uint64_t>(base);
        hash *= 1099511628211ULL;
    }
    return hash;
}

double base_probability(char observed, char truth, int quality) {
    const double error = std::pow(10.0, -quality / 10.0);
    return observed == truth ? 1.0 - error : error / 3.0;
}

double agreement_probability(char left, int left_q, char right, int right_q) {
    double result = 0;
    for (char truth : std::string("ACGT"))
        result += base_probability(left, truth, left_q) * base_probability(right, truth, right_q);
    return std::max(result, 1e-30);
}

uint64_t read_count(const std::string &path) {
    std::ifstream input(path);
    if (!input) return 0;
    uint64_t lines = 0;
    std::string line;
    while (std::getline(input, line)) ++lines;
    return lines / 4;
}

struct OrientedPosting { uint32_t read; uint8_t orientation; };
uint32_t oriented_vertex(uint32_t read, uint8_t orientation) { return read * 2U + orientation; }

uint64_t peak_rss_mb() {
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return static_cast<uint64_t>(usage.ru_maxrss / 1024);
}

double cpu_seconds() {
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1e6 +
           usage.ru_stime.tv_sec + usage.ru_stime.tv_usec / 1e6;
}

} // namespace

CleanGraphAssembler::CleanGraphAssembler(CleanConfig config) : config_(config) {}

bool CleanGraphAssembler::load_candidates(const std::string &r1, const std::string &r2,
                                           const std::string &candidate_stats) {
    std::unordered_map<std::string, uint64_t> existing_stats;
    {
        std::ifstream prior(candidate_stats);
        std::string line;
        while (std::getline(prior, line)) {
            const size_t tab = line.find('\t');
            if (tab != std::string::npos && line.substr(0, tab) != "metric")
                existing_stats[line.substr(0, tab)] = std::strtoull(line.c_str() + tab + 1, nullptr, 10);
        }
    }
    std::ifstream first(r1), second(r2);
    if (!first || !second) return false;
    std::string n1, s1, p1, q1, n2, s2, p2, q2;
    while (std::getline(first, n1)) {
        if (!std::getline(first, s1) || !std::getline(first, p1) || !std::getline(first, q1) ||
            !std::getline(second, n2) || !std::getline(second, s2) ||
            !std::getline(second, p2) || !std::getline(second, q2)) return false;
        if (s1.size() != q1.size() || s2.size() != q2.size()) return false;
        const uint32_t left = static_cast<uint32_t>(reads_.size());
        reads_.push_back({left, n1.substr(1), n2.substr(1), std::move(s1), std::move(q1), left + 1});
        const uint32_t right = static_cast<uint32_t>(reads_.size());
        reads_.push_back({right, n2.substr(1), n1.substr(1), std::move(s2), std::move(q2), left});
        reads_[left].mate_id = right;
    }
    const uint64_t candidate_pairs = reads_.size() / 2;
    const uint64_t raw_pairs = existing_stats.count("raw_read_pairs")
        ? existing_stats["raw_read_pairs"]
        : std::min(read_count(r1), read_count(r2));
    raw_reads_ = raw_pairs * 2;
    std::ofstream stats(candidate_stats);
    if (!stats) return false;
    stats << "metric\tvalue\nraw_read_pairs\t" << raw_pairs
          << "\ncandidate_read_pairs\t" << candidate_pairs
          << "\ncandidate_fraction\t" << (raw_pairs ? static_cast<double>(candidate_pairs) / raw_pairs : 0)
          << "\nraw_reads\t" << raw_reads_ << "\ncandidate_reads\t" << reads_.size() << '\n';
    return !reads_.empty();
}

bool CleanGraphAssembler::construct_graph() {
    std::unordered_map<uint64_t, std::vector<OrientedPosting>> prefix_index;
    prefix_index.reserve(reads_.size() * 2);
    for (const CleanRead &read : reads_) {
        const std::string rc = reverse_complement(read.sequence);
        for (uint8_t orientation = 0; orientation < 2; ++orientation) {
            const std::string &seq = orientation ? rc : read.sequence;
            if (seq.size() < static_cast<size_t>(config_.kmer)) continue;
            const uint64_t key = hash_kmer(seq, 0, config_.kmer);
            auto &posting = prefix_index[key];
            if (posting.size() < config_.max_postings)
                posting.push_back({read.id, orientation});
        }
    }

    outgoing_.assign(reads_.size() * 2, {});
    incoming_.assign(reads_.size() * 2, {});
    std::vector<std::vector<CleanEdge>> local_edges(config_.threads);
    std::vector<uint64_t> local_candidate_hits(config_.threads, 0);
    std::vector<uint64_t> local_scored_pairs(config_.threads, 0);
    std::vector<uint64_t> local_capped_hits(config_.threads, 0);
    std::vector<uint64_t> local_rejected(config_.threads, 0);
#pragma omp parallel for schedule(dynamic, 16) num_threads(config_.threads)
    for (int64_t read_index = 0; read_index < static_cast<int64_t>(reads_.size()); ++read_index) {
        const unsigned worker = static_cast<unsigned>(omp_get_thread_num());
        const CleanRead &source_read = reads_[static_cast<size_t>(read_index)];
        std::unordered_set<uint64_t> candidates;
        candidates.reserve(config_.max_candidates_per_read * 2U);
        const size_t max_overlap = std::min<size_t>(config_.max_overlap, source_read.sequence.size());
        for (uint8_t source_orientation = 0; source_orientation < 2; ++source_orientation) {
            const std::string source = source_orientation ? reverse_complement(source_read.sequence) : source_read.sequence;
            if (source.size() < static_cast<size_t>(config_.min_overlap)) continue;
            const size_t max_ov = std::min(max_overlap, source.size());
            for (size_t overlap = static_cast<size_t>(config_.min_overlap); overlap <= max_ov; ++overlap) {
                const size_t position = source.size() - overlap;
                const uint64_t key = hash_kmer(source, position, config_.kmer);
                auto found = prefix_index.find(key);
                if (found == prefix_index.end()) continue;
                for (const OrientedPosting &posting : found->second) {
                    if (posting.read == source_read.id) continue;
                    ++local_candidate_hits[worker];
                    const uint64_t encoded = (static_cast<uint64_t>(source_orientation) << 63) |
                                             (static_cast<uint64_t>(posting.orientation) << 62) |
                                             posting.read;
                    if (candidates.find(encoded) != candidates.end()) continue;
                    if (candidates.size() < config_.max_candidates_per_read) candidates.insert(encoded);
                    else ++local_capped_hits[worker];
                }
            }
        }
        std::vector<uint64_t> ordered(candidates.begin(), candidates.end());
        std::sort(ordered.begin(), ordered.end());
        local_scored_pairs[worker] += ordered.size();
        for (uint64_t code : ordered) {
            const uint8_t source_orientation = static_cast<uint8_t>((code >> 63) & 1);
            const uint8_t target_orientation = static_cast<uint8_t>((code >> 62) & 1);
            const uint32_t target = static_cast<uint32_t>(code & 0xffffffffU);
            const CleanRead &target_read = reads_[target];
            const std::string source_seq = source_orientation ? reverse_complement(source_read.sequence) : source_read.sequence;
            const std::string source_qual = source_orientation ? reverse_quality(source_read.quality) : source_read.quality;
            const std::string target_seq = target_orientation ? reverse_complement(target_read.sequence) : target_read.sequence;
            const std::string target_qual = target_orientation ? reverse_quality(target_read.quality) : target_read.quality;
            const size_t maximum = std::min({static_cast<size_t>(config_.max_overlap), source_seq.size(), target_seq.size()});
            std::vector<CleanEdge> accepted_overlaps;
            for (size_t overlap = maximum; overlap >= static_cast<size_t>(config_.min_overlap); --overlap) {
                size_t matches = 0;
                double log_qaos = 0, quality_total = 0, evidence_total = 0;
                for (size_t i = 0; i < overlap; ++i) {
                    const size_t si = source_seq.size() - overlap + i;
                    const int q1 = phred(source_qual[si]);
                    const int q2 = phred(target_qual[i]);
                    matches += source_seq[si] == target_seq[i];
                    quality_total += (q1 + q2) / 2.0;
                    evidence_total += 1.0 - std::pow(10.0, -std::min(q1, q2) / 10.0);
                    log_qaos += std::log(agreement_probability(source_seq[si], q1, target_seq[i], q2));
                }
                const double identity = static_cast<double>(matches) / overlap;
                const double qaos = std::exp(log_qaos / overlap);
                if (identity >= config_.min_identity && qaos >= config_.min_qaos) {
                    CleanEdge edge;
                    edge.source = source_read.id;
                    edge.target = target;
                    edge.source_orientation = source_orientation;
                    edge.target_orientation = target_orientation;
                    edge.overlap = static_cast<uint16_t>(overlap);
                    edge.matches = static_cast<uint16_t>(matches);
                    edge.mismatches = static_cast<uint16_t>(overlap - matches);
                    edge.identity = identity;
                    edge.qaos = qaos;
                    edge.mean_quality = quality_total / overlap;
                    edge.quality_support = evidence_total / overlap;
                    accepted_overlaps.push_back(edge);
                }
                if (overlap == static_cast<size_t>(config_.min_overlap)) break;
            }
            if (!accepted_overlaps.empty()) {
                for (const CleanEdge &edge : accepted_overlaps) local_edges[worker].push_back(edge);
            } else {
                ++local_rejected[worker];
            }
        }
    }
    for (const auto &batch : local_edges) {
        for (CleanEdge edge : batch) {
            edge.id = static_cast<uint32_t>(edges_.size());
            edges_.push_back(edge);
            outgoing_[oriented_vertex(edge.source, edge.source_orientation)].push_back(edge.id);
            incoming_[oriented_vertex(edge.target, edge.target_orientation)].push_back(edge.id);
        }
    }
    candidate_edges_ = std::accumulate(local_candidate_hits.begin(), local_candidate_hits.end(), uint64_t{0});
    candidate_pairs_scored_ = std::accumulate(local_scored_pairs.begin(), local_scored_pairs.end(), uint64_t{0});
    rejected_cap_ = std::accumulate(local_capped_hits.begin(), local_capped_hits.end(), uint64_t{0});
    rejected_alignment_ = std::accumulate(local_rejected.begin(), local_rejected.end(), uint64_t{0});
    return true;
}

void CleanGraphAssembler::traverse_paths() {
    std::vector<std::vector<uint32_t>> weak(reads_.size());
    for (const CleanEdge &edge : edges_) {
        weak[edge.source].push_back(edge.target);
        weak[edge.target].push_back(edge.source);
    }
    std::vector<int64_t> component(reads_.size(), -1);
    std::vector<std::vector<uint32_t>> members;
    for (uint32_t read = 0; read < reads_.size(); ++read) {
        if (component[read] >= 0) continue;
        const int64_t cid = static_cast<int64_t>(members.size());
        members.emplace_back();
        std::queue<uint32_t> queue;
        queue.push(read);
        component[read] = cid;
        while (!queue.empty()) {
            const uint32_t current = queue.front(); queue.pop();
            members.back().push_back(current);
            for (uint32_t next : weak[current]) if (component[next] < 0) {
                component[next] = cid;
                queue.push(next);
            }
        }
    }

    std::vector<uint32_t> path_count(members.size(), 0);
    std::vector<uint8_t> visited(reads_.size() * 2, 0);
    std::vector<CleanPath> path;
    auto emit = [&](CleanPath &&candidate, const std::string &why) {
        candidate.termination = why;
        paths_.push_back(std::move(candidate));
        ++explored_paths_;
    };
    auto walk = [&](auto &&self, CleanPath current, uint32_t cid) -> void {
        const uint32_t state = oriented_vertex(current.nodes.back(), current.orientations.back());
        if (path_count[cid] >= config_.max_paths_per_component) return;
        const auto &next_edges = outgoing_[state];
        if (current.nodes.size() >= config_.max_nodes_per_path) {
            current.length_limited = true;
            emit(std::move(current), "node_limit");
            ++path_count[cid];
            return;
        }
        bool extended = false;
        for (uint32_t edge_id : next_edges) {
            const CleanEdge &edge = edges_[edge_id];
            const uint32_t next_state = oriented_vertex(edge.target, edge.target_orientation);
            if (visited[next_state]) {
                CleanPath cycle_path = current;
                cycle_path.cycle_terminated = true;
                cycle_path.edges.push_back(edge_id);
                cycle_path.cumulative_qaos += edge.qaos;
                cycle_path.cumulative_identity += edge.identity;
                emit(std::move(cycle_path), "cycle");
                ++cycles_;
                ++path_count[cid];
                continue;
            }
            if (path_count[cid] >= config_.max_paths_per_component) break;
            extended = true;
            CleanPath next = current;
            next.nodes.push_back(edge.target);
            next.orientations.push_back(edge.target_orientation);
            next.edges.push_back(edge_id);
            next.cumulative_qaos += edge.qaos;
            next.cumulative_identity += edge.identity;
            visited[next_state] = 1;
            self(self, std::move(next), cid);
            visited[next_state] = 0;
        }
        if (!extended && path_count[cid] < config_.max_paths_per_component) {
            emit(std::move(current), next_edges.empty() ? "dead_end" : "cycle_blocked");
            ++path_count[cid];
        }
    };

    for (uint32_t cid = 0; cid < members.size(); ++cid) {
        std::vector<uint32_t> starts;
        for (uint32_t read : members[cid]) {
            for (uint8_t orientation = 0; orientation < 2; ++orientation) {
                const uint32_t state = oriented_vertex(read, orientation);
                if (outgoing_[state].empty() || incoming_[state].empty() ||
                    outgoing_[state].size() != 1 || incoming_[state].size() != 1)
                    starts.push_back(state);
            }
        }
        if (members[cid].size() == 1 && weak[members[cid][0]].empty()) {
            CleanPath singleton;
            singleton.component = cid;
            singleton.nodes.push_back(members[cid][0]);
            singleton.orientations.push_back(0);
            emit(std::move(singleton), "isolated_node");
            ++path_count[cid];
            continue;
        }
        if (starts.empty() && !members[cid].empty()) starts.push_back(oriented_vertex(members[cid][0], 0));
        std::sort(starts.begin(), starts.end());
        starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
        for (uint32_t start : starts) {
            if (path_count[cid] >= config_.max_paths_per_component) break;
            CleanPath initial;
            initial.component = cid;
            initial.nodes.push_back(start / 2);
            initial.orientations.push_back(static_cast<uint8_t>(start % 2));
            visited[start] = 1;
            walk(walk, std::move(initial), cid);
            visited[start] = 0;
        }
        if (path_count[cid] >= config_.max_paths_per_component) ++capped_components_;
    }
}

CleanPath CleanGraphAssembler::consensus_path(CleanPath path, std::string &sequence,
                                               std::vector<std::array<double, 4>> &support) const {
    if (path.nodes.empty()) return path;
    std::vector<size_t> starts(path.nodes.size(), 0);
    for (size_t i = 1; i < path.nodes.size() && i - 1 < path.edges.size(); ++i) {
        starts[i] = starts[i - 1] +
            reads_[path.nodes[i - 1]].sequence.size() - edges_[path.edges[i - 1]].overlap;
    }
    size_t length = 0;
    for (size_t i = 0; i < path.nodes.size(); ++i)
        length = std::max(length, starts[i] + reads_[path.nodes[i]].sequence.size());
    support.assign(length, {0.0, 0.0, 0.0, 0.0});
    for (size_t i = 0; i < path.nodes.size(); ++i) {
        const CleanRead &read = reads_[path.nodes[i]];
        const std::string seq = path.orientations[i] ? reverse_complement(read.sequence) : read.sequence;
        const std::string qual = path.orientations[i] ? reverse_quality(read.quality) : read.quality;
        for (size_t pos = 0; pos < seq.size(); ++pos) {
            const char observed = seq[pos];
            const int q = phred(qual[pos]);
            const int observed_index = observed == 'A' ? 0 : observed == 'C' ? 1 : observed == 'G' ? 2 : observed == 'T' ? 3 : -1;
            if (observed_index < 0) continue;
            for (int base = 0; base < 4; ++base)
                support[starts[i] + pos][base] += base_probability(observed, "ACGT"[base], q);
        }
    }
    sequence.reserve(length);
    constexpr char bases[] = "ACGT";
    for (const auto &position : support) {
        const auto best = std::max_element(position.begin(), position.end());
        sequence.push_back(bases[std::distance(position.begin(), best)]);
    }
    return path;
}

void CleanGraphAssembler::write_outputs(const std::string &output) const {
    mkdir(output.c_str(), 0775);
    std::ofstream nodes(output + "/graph_nodes.tsv");
    nodes << "node_id\tread_id\tmate_id\tlength\tsequence\tquality\n";
    std::ofstream node_bed(output + "/graph_nodes.bed");
    std::ofstream edges_out(output + "/graph_edges.tsv");
    edges_out << "edge_id\tsource_read\ttarget_read\tsource_orientation\ttarget_orientation\toverlap_length\tmatches\tmismatches\tidentity\tQAOS\tmean_base_quality\tquality_support\n";
    std::ofstream edges_bed(output + "/graph_edges.bed");
    edges_bed << "# BEDPE-style overlap edges; node-relative coordinates, not genomic coordinates\n";
    for (const CleanRead &read : reads_) {
        nodes << read.id << '\t' << read.name << '\t' << read.mate_id << '\t'
              << read.sequence.size() << '\t' << read.sequence << '\t' << read.quality << '\n';
        node_bed << "read_" << read.id << "\t0\t" << read.sequence.size() << "\t"
                 << read.name << "\t0\t+\n";
    }
    for (const CleanEdge &edge : edges_) {
        edges_out << edge.id << '\t' << reads_[edge.source].name << '\t' << reads_[edge.target].name
                  << '\t' << static_cast<int>(edge.source_orientation) << '\t'
                  << static_cast<int>(edge.target_orientation) << '\t' << edge.overlap << '\t'
                  << edge.matches << '\t' << edge.mismatches << '\t' << edge.identity << '\t'
                  << edge.qaos << '\t' << edge.mean_quality << '\t' << edge.quality_support << '\n';
        edges_bed << "read_" << edge.source << '\t'
                  << reads_[edge.source].sequence.size() - edge.overlap << '\t'
                  << reads_[edge.source].sequence.size() << "\tread_" << edge.target << "\t0\t"
                  << edge.overlap << '\t' << edge.identity << '\t'
                  << static_cast<int>(edge.source_orientation) << '\t'
                  << static_cast<int>(edge.target_orientation) << '\t' << edge.qaos << '\n';
    }

    std::ofstream fasta(output + "/assembled_contigs.fa");
    std::ofstream path_stats(output + "/path_statistics.tsv");
    path_stats << "contig_id\tcomponent_id\tlength\tnodes\tedges\tmean_identity\tmean_QAOS\tcumulative_QAOS\ttermination\tcycle_terminated\tnode_limit_reached\n";
    std::ofstream provenance(output + "/contig_paths.tsv");
    provenance << "contig_id\tcomponent_id\tstep\tread_id\tmate_id\tedge_id\torientation\tpath_start\n";
    std::ofstream path_bed(output + "/graph_paths.bed");
    std::ofstream consensus_out(output + "/consensus_support.tsv");
    consensus_out << "contig_id\tposition\tA_support\tC_support\tG_support\tT_support\tselected_base\tconfidence\n";
    std::vector<size_t> lengths;
    uint64_t total_bases = 0;
    for (size_t path_id = 0; path_id < paths_.size(); ++path_id) {
        std::string sequence;
        std::vector<std::array<double, 4>> support;
        CleanPath path = consensus_path(paths_[path_id], sequence, support);
        if (sequence.empty()) continue;
        lengths.push_back(sequence.size());
        total_bases += sequence.size();
        fasta << ">contig_" << path_id << " component=" << path.component
              << " nodes=" << path.nodes.size() << " termination=" << path.termination << '\n'
              << sequence << '\n';
        const double mean_identity = path.edges.empty() ? 0 :
            path.cumulative_identity / path.edges.size();
        const double mean_qaos = path.edges.empty() ? 0 :
            path.cumulative_qaos / path.edges.size();
        path_stats << path_id << '\t' << path.component << '\t' << sequence.size() << '\t'
                   << path.nodes.size() << '\t' << path.edges.size() << '\t'
                   << mean_identity << '\t' << mean_qaos << '\t' << path.cumulative_qaos
                   << '\t' << path.termination << '\t' << path.cycle_terminated << '\t'
                   << path.length_limited << '\n';
        size_t offset = 0;
        for (size_t step = 0; step < path.nodes.size(); ++step) {
            if (step > 0 && step - 1 < path.edges.size())
                offset += reads_[path.nodes[step - 1]].sequence.size() - edges_[path.edges[step - 1]].overlap;
            const CleanRead &read = reads_[path.nodes[step]];
            provenance << path_id << '\t' << path.component << '\t' << step << '\t'
                       << read.name << '\t' << read.mate_name << '\t'
                       << (step < path.edges.size() ? std::to_string(path.edges[step]) : "NA")
                       << '\t' << static_cast<int>(path.orientations[step]) << '\t' << offset << '\n';
            path_bed << "contig_" << path_id << '\t' << offset << '\t'
                     << offset + read.sequence.size() << "\tread_" << read.id << '\t'
                     << (step < path.edges.size() ? edges_[path.edges[step]].qaos : 0) << '\t'
                     << (path.orientations[step] ? '-' : '+') << '\n';
        }
        for (size_t pos = 0; pos < support.size(); ++pos) {
            const auto best = std::max_element(support[pos].begin(), support[pos].end());
            const size_t selected = static_cast<size_t>(std::distance(support[pos].begin(), best));
            const double total = std::accumulate(support[pos].begin(), support[pos].end(), 0.0);
            consensus_out << path_id << '\t' << pos << '\t' << support[pos][0] << '\t'
                          << support[pos][1] << '\t' << support[pos][2] << '\t'
                          << support[pos][3] << "\tACGT"[selected] << '\t'
                          << (total ? *best / total : 0) << '\n';
        }
    }

    std::vector<std::vector<uint32_t>> weak(reads_.size());
    for (const CleanEdge &edge : edges_) {
        weak[edge.source].push_back(edge.target);
        weak[edge.target].push_back(edge.source);
    }
    std::vector<uint8_t> seen(reads_.size(), 0);
    std::vector<std::vector<uint32_t>> components;
    for (uint32_t root = 0; root < reads_.size(); ++root) {
        if (seen[root]) continue;
        components.emplace_back();
        std::queue<uint32_t> q; q.push(root); seen[root] = 1;
        while (!q.empty()) {
            uint32_t node = q.front(); q.pop();
            components.back().push_back(node);
            for (uint32_t next : weak[node]) if (!seen[next]) { seen[next] = 1; q.push(next); }
        }
    }
    uint64_t isolated = 0, branch_nodes = 0;
    std::vector<uint32_t> degree(reads_.size(), 0);
    for (const CleanEdge &edge : edges_) { ++degree[edge.source]; ++degree[edge.target]; }
    for (uint32_t d : degree) { isolated += d == 0; branch_nodes += d > 2; }

    std::vector<double> overlaps, identities, qaos;
    for (const CleanEdge &edge : edges_) {
        overlaps.push_back(edge.overlap);
        identities.push_back(edge.identity);
        qaos.push_back(edge.qaos);
    }
    auto mean = [](const std::vector<double> &v) {
        return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };
    auto median = [](std::vector<double> v) {
        if (v.empty()) return 0.0;
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        return v[mid];
    };
    std::sort(lengths.begin(), lengths.end(), std::greater<size_t>());
    uint64_t half = 0, n50 = 0;
    for (size_t length : lengths) { half += length; if (half * 2 >= total_bases) { n50 = length; break; } }
    const double elapsed = seconds() - start_seconds_;
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::ofstream graph_stats(output + "/graph_statistics.tsv");
    graph_stats << "statistic\tvalue\nraw_reads\t" << raw_reads_
                << "\ncandidate_reads\t" << reads_.size()
                << "\ncandidate_pairs\t" << reads_.size() / 2
                << "\ngraph_nodes\t" << reads_.size()
                << "\ncandidate_kmer_hits\t" << candidate_edges_
                << "\ncandidate_pairs_scored\t" << candidate_pairs_scored_
                << "\ngraph_edges\t" << edges_.size()
                << "\nunretained_candidate_hits_at_cap\t" << rejected_cap_
                << "\nrejected_by_alignment_thresholds\t" << rejected_alignment_
                << "\nconnected_components\t" << components.size()
                << "\nisolated_nodes\t" << isolated
                << "\nbranch_nodes\t" << branch_nodes
                << "\ncycles\t" << cycles_
                << "\npaths_generated\t" << explored_paths_
                << "\npaths_emitted\t" << paths_.size()
                << "\ncomponents_at_path_limit\t" << capped_components_
                << "\nmax_paths_per_component\t" << config_.max_paths_per_component
                << "\nmax_nodes_per_path\t" << config_.max_nodes_per_path
                << "\nmean_overlap\t" << mean(overlaps) << "\nmedian_overlap\t" << median(overlaps)
                << "\nmean_identity\t" << mean(identities) << "\nmedian_identity\t" << median(identities)
                << "\nmean_QAOS\t" << mean(qaos) << "\nmedian_QAOS\t" << median(qaos)
                << "\nN50\t" << n50 << "\nassembled_bases\t" << total_bases
                << "\nrequested_threads\t" << config_.threads
                << "\nopenmp_threads\t" << omp_get_max_threads()
                << "\nCPU_utilization_percent\t" << (elapsed ? 100 * cpu_seconds() / (elapsed * config_.threads) : 0)
                << "\nelapsed_seconds\t" << elapsed
                << "\nuser_CPU_seconds\t" << usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1e6
                << "\npeak_RSS_MB\t" << peak_rss_mb()
                << "\nmajor_page_faults\t" << usage.ru_majflt
                << "\nminor_page_faults\t" << usage.ru_minflt << '\n';

    mkdir((output + "/visualization").c_str(), 0775);
    std::ofstream visualization(output + "/visualization/README.txt");
    visualization << "DOT files show selected diagnostic overlap-graph components.\n"
                  << "A simple linear component, a branched component, an alternative-path component,\n"
                  << "and the largest long-read component are selected when available. The latter is\n"
                  << "only a long-path candidate; it is not asserted to be a biological BCR.\n";
    std::vector<uint32_t> selected;
    if (!components.empty()) {
        auto largest = static_cast<uint32_t>(std::max_element(components.begin(), components.end(),
            [](const auto &a, const auto &b) { return a.size() < b.size(); }) - components.begin());
        selected.push_back(largest);
        for (uint32_t cid = 0; cid < components.size() && selected.size() < 4; ++cid)
            if (cid != largest) selected.push_back(cid);
    }
    for (uint32_t cid : selected) {
        const auto &component = components[cid];
        std::unordered_set<uint32_t> members(component.begin(), component.end());
        const std::string dot_path = output + "/visualization/component_" + std::to_string(cid) + ".dot";
        std::ofstream dot(dot_path);
        dot << "digraph component_" << cid << " {\nrankdir=LR;\n";
        const size_t node_limit = std::min<size_t>(component.size(), 100);
        for (size_t i = 0; i < node_limit; ++i) {
            uint32_t n = component[i];
            dot << "r" << n << " [label=\"" << n << "\\n" << reads_[n].sequence.size() << "bp\"];\n";
        }
        unsigned emitted = 0;
        for (const CleanEdge &edge : edges_) {
            if (emitted >= 300) break;
            if (members.count(edge.source) && members.count(edge.target) &&
                edge.source < component[node_limit - 1] + 1 && edge.target < component[node_limit - 1] + 1) {
                dot << "r" << edge.source << " -> r" << edge.target << " [label=\"ov="
                    << edge.overlap << " id=" << edge.identity << " q=" << edge.qaos << "\"];\n";
                ++emitted;
            }
        }
        dot << "}\n";
    }
}

int CleanGraphAssembler::run(const std::string &r1, const std::string &r2,
                             const std::string &output, const std::string &candidate_stats) {
    start_seconds_ = seconds();
    if (!load_candidates(r1, r2, candidate_stats)) {
        std::cerr << "failed to read paired candidate FASTQs\n";
        return 1;
    }
    if (!construct_graph()) {
        std::cerr << "graph construction failed\n";
        return 1;
    }
    traverse_paths();
    write_outputs(output);
    std::cerr << "candidate_reads=" << reads_.size() << " graph_edges=" << edges_.size()
              << " paths=" << paths_.size() << " peak_RSS_MB=" << peak_rss_mb() << '\n';
    return 0;
}
