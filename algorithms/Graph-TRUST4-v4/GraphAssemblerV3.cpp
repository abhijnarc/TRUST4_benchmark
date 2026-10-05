#include "GraphAssemblerV3.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <numeric>
#include <queue>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

double now_seconds() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int phred(char c) {
    return std::max(0, std::min(60, static_cast<int>(c) - 33));
}

double error_probability(int q) {
    return std::pow(10.0, -q / 10.0);
}

float compatibility(char a, int qa, char b, int qb) {
    const double e1 = error_probability(qa);
    const double e2 = error_probability(qb);
    if (a == b) return static_cast<float>((1.0 - e1) * (1.0 - e2));
    return static_cast<float>(e1 / 3.0 + e2 / 3.0 - e1 * e2 / 9.0);
}

bool read_fastq_record(std::ifstream &file, std::string &sequence, std::string &quality) {
    std::string name, plus;
    return static_cast<bool>(std::getline(file, name) && std::getline(file, sequence) &&
                             std::getline(file, plus) && std::getline(file, quality));
}

int cheap_overlap(const std::string &left, const std::string &right, int minimum) {
    const int maximum = static_cast<int>(std::min(left.size(), right.size()));
    for (int length = maximum; length >= minimum; --length) {
        const int start = static_cast<int>(left.size()) - length;
        int matches = 0;
        for (int i = 0; i < length; ++i) matches += left[start + i] == right[i];
        if (static_cast<float>(matches) / length >= 0.80f) return length;
    }
    return 0;
}

} // namespace

GraphAssemblerV4::GraphAssemblerV4(const V4Config &config) : config_(config), start_seconds_(now_seconds()) {}

uint64_t GraphAssemblerV4::rss_mb() const {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return static_cast<uint64_t>(usage.ru_maxrss / 1024);
}

bool GraphAssemblerV4::memory_ok(const char *stage) {
    const uint64_t rss = rss_mb();
    stats_.peak_rss_mb = std::max(stats_.peak_rss_mb, rss);
    std::fprintf(stderr, "[Memory] stage=%s rss=%llu MB limit=%d GB\n", stage,
                 static_cast<unsigned long long>(rss), config_.max_memory_gb);
    if (rss >= static_cast<uint64_t>(config_.max_memory_gb) * 1024ULL) {
        std::fprintf(stderr, "[Memory] limit reached; stopping before next batch\n");
        return false;
    }
    return true;
}

void GraphAssemblerV4::progress(const char *stage, uint64_t done, uint64_t total,
                                uint64_t candidates, uint64_t accepted) const {
    std::fprintf(stderr, "[%s] bundles=%llu/%llu threads=%d candidate_pairs=%llu accepted_edges=%llu RSS=%llu MB elapsed=%.1fs\n",
                 stage, static_cast<unsigned long long>(done), static_cast<unsigned long long>(total),
                 config_.threads, static_cast<unsigned long long>(candidates),
                 static_cast<unsigned long long>(accepted), static_cast<unsigned long long>(rss_mb()),
                 now_seconds() - start_seconds_);
}

bool GraphAssemblerV4::load_and_bundle(const std::string &r1, const std::string &r2) {
    std::ifstream first(r1), second(r2);
    if (!first || !second) {
        std::fprintf(stderr, "ERROR: cannot open FASTQ inputs\n");
        return false;
    }
    struct BundleValue { std::string quality; uint32_t count = 0; };
    std::unordered_map<std::string, BundleValue> unique;
    std::string sequence1, quality1, sequence2, quality2;
    while (read_fastq_record(first, sequence1, quality1) && read_fastq_record(second, sequence2, quality2)) {
        ++stats_.raw_reads;
        ++unique[sequence1].count;
        if (unique[sequence1].quality.empty()) unique[sequence1].quality = quality1;
        ++stats_.raw_reads;
        ++unique[sequence2].count;
        if (unique[sequence2].quality.empty()) unique[sequence2].quality = quality2;
        if (stats_.raw_reads % 100000 == 0) progress("Load", stats_.raw_reads, 0, 0, 0);
    }
    bundles_.reserve(unique.size());
    for (auto &entry : unique) {
        V4Bundle bundle;
        bundle.id = static_cast<uint32_t>(bundles_.size());
        bundle.sequence = std::move(entry.first);
        bundle.quality = std::move(entry.second.quality);
        bundle.abundance = entry.second.count;
        bundles_.push_back(std::move(bundle));
    }
    stats_.bundles = bundles_.size();
    unique.clear();
    unique.rehash(0);
    std::fprintf(stderr, "[Bundle] raw_reads=%llu bundles=%llu compression=%.2fx\n",
                 static_cast<unsigned long long>(stats_.raw_reads),
                 static_cast<unsigned long long>(stats_.bundles),
                 stats_.bundles ? static_cast<double>(stats_.raw_reads) / stats_.bundles : 0.0);
    return memory_ok("bundling");
}

bool GraphAssemblerV4::build_index() {
    std::fprintf(stderr, "[Index] building compact k-mer index k=%d max_postings=%u\n",
                 config_.kmer, config_.max_postings);
    index_.reserve(bundles_.size() * 2);
    for (const V4Bundle &bundle : bundles_) {
        if (bundle.sequence.size() < static_cast<size_t>(config_.kmer)) continue;
        std::vector<std::string> seen;
        seen.reserve(bundle.sequence.size());
        for (size_t pos = 0; pos + config_.kmer <= bundle.sequence.size(); ++pos) {
            std::string kmer = bundle.sequence.substr(pos, config_.kmer);
            if (std::find(seen.begin(), seen.end(), kmer) != seen.end()) continue;
            seen.push_back(kmer);
            std::vector<uint32_t> &posting = index_[kmer];
            if (posting.size() < config_.max_postings) posting.push_back(bundle.id);
        }
        if ((bundle.id + 1) % 5000 == 0) progress("Index", bundle.id + 1, bundles_.size(), 0, 0);
    }
    return memory_ok("index");
}

bool GraphAssemblerV4::score_overlap(uint32_t a, uint32_t b, V4Edge &edge) const {
    const V4Bundle &left = bundles_[a];
    const V4Bundle &right = bundles_[b];
    const int maximum = static_cast<int>(std::min(left.sequence.size(), right.sequence.size()));
    int best_length = 0, best_mismatches = 0, best_high = 0;
    float best_qaos = 0.0f, best_expected = 0.0f, best_mean_q = 0.0f, best_min_q = 0.0f;
    for (int length = config_.min_overlap; length <= maximum; ++length) {
        const int start = static_cast<int>(left.sequence.size()) - length;
        double log_probability = 0.0;
        float expected = 0.0f;
        float quality_sum = 0.0f;
        float minimum_q = 60.0f;
        int mismatches = 0, high = 0;
        for (int i = 0; i < length; ++i) {
            const char base1 = left.sequence[start + i], base2 = right.sequence[i];
            const int q1 = phred(left.quality[start + i]), q2 = phred(right.quality[i]);
            const double probability = std::max(1e-30, static_cast<double>(compatibility(base1, q1, base2, q2)));
            log_probability += std::log(probability);
            quality_sum += (q1 + q2) / 2.0f;
            minimum_q = std::min(minimum_q, static_cast<float>((q1 + q2) / 2.0f));
            if (base1 != base2) {
                ++mismatches;
                const double e1 = error_probability(q1), e2 = error_probability(q2);
                expected += static_cast<float>(e1 + e2 - e1 * e2);
                if (q1 >= config_.qhigh && q2 >= config_.qhigh) ++high;
            }
        }
        const float qaos = static_cast<float>(std::exp(log_probability / length));
        const float identity = 1.0f - static_cast<float>(mismatches) / length;
        if (identity >= config_.min_identity && qaos >= config_.min_qaos &&
            (qaos > best_qaos || (qaos == best_qaos && length > best_length))) {
            best_length = length;
            best_mismatches = mismatches;
            best_high = high;
            best_qaos = qaos;
            best_expected = expected;
            best_mean_q = quality_sum / length;
            best_min_q = minimum_q;
        }
    }
    if (!best_length) return false;
    edge.source = a;
    edge.target = b;
    edge.overlap = static_cast<uint16_t>(best_length);
    edge.mismatches = static_cast<uint16_t>(best_mismatches);
    edge.identity = 1.0f - static_cast<float>(best_mismatches) / best_length;
    edge.qaos = best_qaos;
    edge.expected_mismatches = best_expected;
    edge.high_quality_mismatches = static_cast<uint16_t>(best_high);
    edge.orientation = 1;
    edge.paired_support = 0;
    edge.mean_q = best_mean_q;
    edge.min_q = best_min_q;
    return true;
}

bool GraphAssemblerV4::generate_and_score() {
    const uint64_t total = bundles_.size();
    const uint32_t batch_size = std::max<uint32_t>(1, config_.batch_size);
    for (uint64_t begin = 0; begin < total; begin += batch_size) {
        const uint64_t end = std::min(total, begin + batch_size);
        const int workers = std::max(1, config_.threads);
        std::vector<std::vector<V4Edge>> local_edges(static_cast<size_t>(workers));
        std::vector<uint64_t> local_candidates(static_cast<size_t>(workers), 0);
        std::vector<uint64_t> local_accepted(static_cast<size_t>(workers), 0);
#pragma omp parallel num_threads(workers)
        {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif
            std::vector<uint32_t> candidates;
            candidates.reserve(config_.max_candidates);
#pragma omp for schedule(dynamic)
            for (long long id = static_cast<long long>(begin); id < static_cast<long long>(end); ++id) {
                candidates.clear();
                const V4Bundle &bundle = bundles_[static_cast<size_t>(id)];
                std::unordered_map<uint32_t, uint16_t> support_map;
                if (bundle.sequence.size() >= static_cast<size_t>(config_.kmer)) {
                    const size_t first = bundle.sequence.size() > 100 ? bundle.sequence.size() - 100 : 0;
                    for (size_t pos = first; pos + config_.kmer <= bundle.sequence.size(); ++pos) {
                        const std::string kmer = bundle.sequence.substr(pos, config_.kmer);
                        auto found = index_.find(kmer);
                        if (found == index_.end()) continue;
                        for (uint32_t candidate : found->second) {
                            if (candidate != static_cast<uint32_t>(id)) ++support_map[candidate];
                        }
                    }
                }
                candidates.reserve(support_map.size());
                for (const auto &entry : support_map) candidates.push_back(entry.first);
                struct RankedCandidate { uint32_t id; uint16_t support; uint16_t overlap; float identity; };
                std::vector<RankedCandidate> ranked;
                ranked.reserve(candidates.size());
                for (uint32_t candidate : candidates) {
                    const uint16_t support = support_map[candidate];
                    const int overlap = cheap_overlap(bundle.sequence, bundles_[candidate].sequence, config_.min_overlap);
                    const int start = overlap ? static_cast<int>(bundle.sequence.size()) - overlap : 0;
                    int matches = 0;
                    for (int i = 0; overlap && i < overlap; ++i) matches += bundle.sequence[start + i] == bundles_[candidate].sequence[i];
                    ranked.push_back({candidate, support, static_cast<uint16_t>(overlap), overlap ? static_cast<float>(matches) / overlap : 0.0f});
                }
                std::sort(ranked.begin(), ranked.end(), [](const RankedCandidate &a, const RankedCandidate &b) {
                    if (a.support != b.support) return a.support > b.support;
                    if (a.overlap != b.overlap) return a.overlap > b.overlap;
                    return a.identity > b.identity;
                });
                if (ranked.size() > config_.max_candidates) ranked.resize(config_.max_candidates);
                local_candidates[thread] += ranked.size();
                for (const RankedCandidate &candidate : ranked) {
                    V4Edge edge;
                    if (score_overlap(static_cast<uint32_t>(id), candidate.id, edge)) {
                        local_edges[thread].push_back(edge);
                        ++local_accepted[thread];
                    }
                }
            }
        }
        uint64_t batch_candidates = 0, batch_accepted = 0;
        for (int thread = 0; thread < workers; ++thread) {
            batch_candidates += local_candidates[thread];
            batch_accepted += local_accepted[thread];
            edges_.insert(edges_.end(), local_edges[thread].begin(), local_edges[thread].end());
        }
        stats_.candidate_pairs += batch_candidates;
        stats_.accepted_edges += batch_accepted;
        progress("Candidates", end, total, stats_.candidate_pairs, stats_.accepted_edges);
        if (!memory_ok("candidate-batch")) return false;
    }
    return true;
}

void GraphAssemblerV4::build_contigs_v4() {
    struct PathState {
        std::vector<uint32_t> nodes;
        std::vector<size_t> edge_ids;
        bool cycle;
    };
    std::vector<std::vector<size_t>> outgoing(bundles_.size());
    adjacency_.assign(bundles_.size(), std::vector<uint32_t>());
    for (size_t edge_id = 0; edge_id < edges_.size(); ++edge_id) {
        const V4Edge &edge = edges_[edge_id];
        outgoing[edge.source].push_back(edge_id);
        adjacency_[edge.source].push_back(edge.target);
        adjacency_[edge.target].push_back(edge.source);
    }
    std::vector<uint32_t> indegree(bundles_.size(), 0);
    for (const V4Edge &edge : edges_) ++indegree[edge.target];
    std::vector<uint8_t> claimed(bundles_.size(), 0);

    auto reverse_complement = [](const std::string &sequence) {
        std::string result;
        result.reserve(sequence.size());
        for (std::string::const_reverse_iterator it = sequence.rbegin(); it != sequence.rend(); ++it) {
            switch (*it) {
                case 'A': result.push_back('T'); break; case 'C': result.push_back('G'); break;
                case 'G': result.push_back('C'); break; case 'T': result.push_back('A'); break;
                case 'a': result.push_back('t'); break; case 'c': result.push_back('g'); break;
                case 'g': result.push_back('c'); break; case 't': result.push_back('a'); break;
                default: result.push_back('N'); break;
            }
        }
        return result;
    };
    auto ranked_edges = [&](uint32_t node, const std::unordered_set<uint32_t> &path_nodes) {
        std::vector<size_t> result;
        for (size_t edge_id : outgoing[node]) {
            if (!path_nodes.count(edges_[edge_id].target)) result.push_back(edge_id);
        }
        std::sort(result.begin(), result.end(), [&](size_t left, size_t right) {
            const V4Edge &a = edges_[left], &b = edges_[right];
            if (a.qaos != b.qaos) return a.qaos > b.qaos;
            if (a.overlap != b.overlap) return a.overlap > b.overlap;
            if (a.high_quality_mismatches != b.high_quality_mismatches) return a.high_quality_mismatches < b.high_quality_mismatches;
            if (a.identity != b.identity) return a.identity > b.identity;
            if (a.paired_support != b.paired_support) return a.paired_support > b.paired_support;
            if (bundles_[a.target].abundance != bundles_[b.target].abundance) return bundles_[a.target].abundance > bundles_[b.target].abundance;
            return a.target < b.target;
        });
        return result;
    };
    auto emit_path = [&](const PathState &state) {
        if (state.nodes.size() < 2) return;
        std::vector<std::string> sequences;
        std::vector<std::string> qualities;
        std::vector<size_t> offsets(state.nodes.size(), 0);
        sequences.push_back(bundles_[state.nodes.front()].sequence);
        qualities.push_back(bundles_[state.nodes.front()].quality);
        size_t end = sequences.front().size();
        uint64_t total_overlap = 0;
        for (size_t i = 0; i < state.edge_ids.size(); ++i) {
            const V4Edge &edge = edges_[state.edge_ids[i]];
            std::string sequence = bundles_[state.nodes[i + 1]].sequence;
            std::string quality = bundles_[state.nodes[i + 1]].quality;
            if (edge.orientation != 1) { sequence = reverse_complement(sequence); std::reverse(quality.begin(), quality.end()); }
            offsets[i + 1] = offsets[i] + sequences[i].size() - edge.overlap;
            end = std::max(end, offsets[i + 1] + sequence.size());
            total_overlap += edge.overlap;
            sequences.push_back(std::move(sequence));
            qualities.push_back(std::move(quality));
        }
        std::vector<std::array<double, 5>> scores(end);
        for (size_t i = 0; i < state.nodes.size(); ++i) {
            for (size_t j = 0; j < sequences[i].size(); ++j) {
                const char observation = sequences[i][j];
                const int q = phred(qualities[i].size() > j ? qualities[i][j] : 'I');
                const double error = error_probability(q);
                for (int base = 0; base < 5; ++base) {
                    const char expected = "ACGTN"[base];
                    const bool match = observation == expected || (expected == 'N' && observation != 'A' && observation != 'C' && observation != 'G' && observation != 'T');
                    const double probability = match ? 1.0 - error : error / 3.0;
                    scores[offsets[i] + j][base] += std::log(std::max(1e-30, probability));
                }
            }
        }
        std::string consensus;
        consensus.reserve(end);
        double confidence_sum = 0.0;
        for (const std::array<double, 5> &score : scores) {
            int best = 0, second = 1;
            for (int base = 1; base < 5; ++base) {
                if (score[base] > score[best]) { second = best; best = base; }
                else if (base != best && score[base] > score[second]) second = base;
            }
            consensus.push_back("ACGTN"[best]);
            confidence_sum += score[best] - score[second];
        }
        const uint64_t expected_length = end;
        V4AssemblyDebug debug;
        debug.contig_id = contigs_.size();
        debug.contig_length = consensus.size();
        debug.num_nodes = state.nodes.size();
        debug.unique_nodes = std::unordered_set<uint32_t>(state.nodes.begin(), state.nodes.end()).size();
        debug.repeated_nodes = debug.num_nodes - debug.unique_nodes;
        debug.total_overlap = total_overlap;
        debug.expected_length = expected_length;
        debug.max_node_length = 0;
        for (uint32_t node : state.nodes) debug.max_node_length = std::max<uint64_t>(debug.max_node_length, bundles_[node].sequence.size());
        debug.cycle_detected = state.cycle;
        debug.valid_path = !state.cycle && debug.repeated_nodes == 0 && debug.contig_length == debug.expected_length;
        debug.consensus_confidence = consensus.empty() ? 0.0f : static_cast<float>(confidence_sum / consensus.size());
        if (state.cycle) ++stats_.cycles;
        contigs_.push_back(std::move(consensus));
        assembly_debug_.push_back(debug);
        consensus_confidences_.push_back(debug.consensus_confidence);
        for (uint32_t node : state.nodes) claimed[node] = 1;
    };

    auto extend = [&](uint32_t root) {
        std::vector<PathState> frontier(1, PathState{{root}, {}, false});
        bool active = true;
        while (active && !frontier.empty()) {
            active = false;
            std::vector<PathState> next;
            for (const PathState &state : frontier) {
                const uint32_t node = state.nodes.back();
                std::unordered_set<uint32_t> path_nodes(state.nodes.begin(), state.nodes.end());
                std::vector<size_t> choices = ranked_edges(node, path_nodes);
                if (choices.empty()) { emit_path(state); continue; }
                if (choices.size() > 1) {
                    ++stats_.branches;
                    const V4Edge &best = edges_[choices[0]], &second = edges_[choices[1]];
                    const float delta = best.qaos - second.qaos;
                    V4BranchDecision decision;
                    decision.branch_node = node; decision.best_edge = choices[0]; decision.second_edge = choices[1];
                    decision.best_qaos = best.qaos; decision.second_qaos = second.qaos; decision.qaos_delta = delta;
                    decision.overlap_best = best.overlap; decision.overlap_second = second.overlap;
                    decision.high_q_mismatch_best = best.high_quality_mismatches; decision.high_q_mismatch_second = second.high_quality_mismatches;
                    decision.paired_support_best = best.paired_support; decision.paired_support_second = second.paired_support;
                    if (delta >= config_.min_branch_qaos_delta) {
                        ++stats_.resolved_branches; decision.decision = "best_only"; decision.reason = "QAOS_delta_threshold";
                        choices.resize(1);
                    } else {
                        ++stats_.ambiguous_branches; decision.decision = "retain_alternatives"; decision.reason = "QAOS_delta_below_threshold";
                        if (choices.size() > config_.max_alternative_paths_per_branch) choices.resize(config_.max_alternative_paths_per_branch);
                    }
                    branch_decisions_.push_back(decision);
                }
                for (size_t edge_id : choices) {
                    const V4Edge &edge = edges_[edge_id];
                    PathState child = state;
                    if (path_nodes.count(edge.target)) { child.cycle = true; next.push_back(child); continue; }
                    child.nodes.push_back(edge.target); child.edge_ids.push_back(edge_id); next.push_back(std::move(child));
                }
                active = true;
            }
            if (next.size() > config_.max_alternative_paths_per_branch) next.resize(config_.max_alternative_paths_per_branch);
            frontier.swap(next);
        }
    };
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (!claimed[root] && !outgoing[root].empty() && (indegree[root] != 1 || outgoing[root].size() > 1)) extend(root);
    }
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (!claimed[root] && !outgoing[root].empty()) extend(root);
    }
    std::vector<uint8_t> component_seen(bundles_.size(), 0);
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (component_seen[root]) continue;
        std::queue<uint32_t> queue; queue.push(root); component_seen[root] = 1; size_t size = 0;
        while (!queue.empty()) { uint32_t node = queue.front(); queue.pop(); ++size; for (uint32_t neighbor : adjacency_[node]) if (!component_seen[neighbor]) { component_seen[neighbor] = 1; queue.push(neighbor); } }
        ++stats_.components; if (size == 1) ++stats_.isolated;
    }
    stats_.contigs = contigs_.size();
    std::vector<size_t> lengths; for (const std::string &contig : contigs_) lengths.push_back(contig.size());
    std::sort(lengths.begin(), lengths.end(), std::greater<size_t>());
    size_t total = std::accumulate(lengths.begin(), lengths.end(), static_cast<size_t>(0)), cumulative = 0;
    for (size_t length : lengths) { cumulative += length; if (cumulative >= total / 2) { stats_.n50 = length; break; } }
}

void GraphAssemblerV4::build_contigs() {
    std::vector<std::vector<size_t>> outgoing(bundles_.size());
    adjacency_.assign(bundles_.size(), std::vector<uint32_t>());
    for (size_t edge_id = 0; edge_id < edges_.size(); ++edge_id) {
        const V4Edge &edge = edges_[edge_id];
        outgoing[edge.source].push_back(edge_id);
        adjacency_[edge.source].push_back(edge.target);
        adjacency_[edge.target].push_back(edge.source);
    }

    std::vector<uint32_t> indegree(bundles_.size(), 0);
    for (const V4Edge &edge : edges_) ++indegree[edge.target];
    std::vector<uint8_t> used_node(bundles_.size(), 0);
    std::vector<uint8_t> used_edge(edges_.size(), 0);

    auto choose_edge = [&](uint32_t node, const std::unordered_set<uint32_t> &path_nodes) -> size_t {
        size_t selected = static_cast<size_t>(-1);
        for (size_t edge_id : outgoing[node]) {
            const V4Edge &candidate = edges_[edge_id];
            if (used_edge[edge_id] || (used_node[candidate.target] && !path_nodes.count(candidate.target))) continue;
            if (selected == static_cast<size_t>(-1) ||
                candidate.qaos > edges_[selected].qaos ||
                (candidate.qaos == edges_[selected].qaos && candidate.overlap > edges_[selected].overlap)) {
                selected = edge_id;
            }
        }
        return selected;
    };

    auto emit_path = [&](const std::vector<uint32_t> &path, const std::vector<size_t> &path_edges,
                         bool cycle_detected) {
        if (path.size() < 2) return;
        std::unordered_set<uint32_t> unique_nodes(path.begin(), path.end());
        uint64_t total_overlap = 0;
        uint64_t expected_length = bundles_[path.front()].sequence.size();
        uint64_t max_node_length = bundles_[path.front()].sequence.size();
        std::string contig = bundles_[path.front()].sequence;
        bool valid_path = true;
        for (size_t i = 0; i < path_edges.size(); ++i) {
            const V4Edge &edge = edges_[path_edges[i]];
            const std::string &node_sequence = bundles_[path[i + 1]].sequence;
            const size_t overlap = std::min<size_t>(edge.overlap, node_sequence.size());
            if (edge.overlap > node_sequence.size() || edge.overlap > contig.size()) valid_path = false;
            total_overlap += edge.overlap;
            expected_length += node_sequence.size() >= edge.overlap ? node_sequence.size() - edge.overlap : 0;
            max_node_length = std::max<uint64_t>(max_node_length, node_sequence.size());
            if (edge.orientation == 1) {
                contig += node_sequence.substr(overlap);
            } else {
                std::string reverse_complement;
                reverse_complement.reserve(node_sequence.size());
                for (std::string::const_reverse_iterator it = node_sequence.rbegin(); it != node_sequence.rend(); ++it) {
                    switch (*it) {
                        case 'A': reverse_complement.push_back('T'); break;
                        case 'C': reverse_complement.push_back('G'); break;
                        case 'G': reverse_complement.push_back('C'); break;
                        case 'T': reverse_complement.push_back('A'); break;
                        case 'a': reverse_complement.push_back('t'); break;
                        case 'c': reverse_complement.push_back('g'); break;
                        case 'g': reverse_complement.push_back('c'); break;
                        case 't': reverse_complement.push_back('a'); break;
                        default: reverse_complement.push_back('N'); break;
                    }
                }
                const size_t reverse_overlap = std::min<size_t>(edge.overlap, reverse_complement.size());
                contig += reverse_complement.substr(reverse_overlap);
            }
        }
        V4AssemblyDebug debug;
        debug.contig_id = contigs_.size();
        debug.contig_length = contig.size();
        debug.num_nodes = path.size();
        debug.unique_nodes = unique_nodes.size();
        debug.repeated_nodes = path.size() - unique_nodes.size();
        debug.total_overlap = total_overlap;
        debug.expected_length = expected_length;
        debug.max_node_length = max_node_length;
        debug.cycle_detected = cycle_detected;
        debug.valid_path = valid_path && !cycle_detected && debug.repeated_nodes == 0 &&
                           debug.contig_length == debug.expected_length;
        contigs_.push_back(std::move(contig));
        assembly_debug_.push_back(debug);
    };

    auto walk = [&](uint32_t root) {
        std::vector<uint32_t> path(1, root);
        std::vector<size_t> path_edges;
        std::unordered_set<uint32_t> path_nodes;
        path_nodes.insert(root);
        used_node[root] = 1;
        uint32_t current = root;
        bool cycle_detected = false;
        while (true) {
            const size_t edge_id = choose_edge(current, path_nodes);
            if (edge_id == static_cast<size_t>(-1)) break;
            const V4Edge &edge = edges_[edge_id];
            used_edge[edge_id] = 1;
            if (path_nodes.count(edge.target)) { cycle_detected = true; break; }
            path_edges.push_back(edge_id);
            path.push_back(edge.target);
            path_nodes.insert(edge.target);
            used_node[edge.target] = 1;
            current = edge.target;
            if (outgoing[current].size() != 1) break;
        }
        emit_path(path, path_edges, cycle_detected);
    };

    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (!used_node[root] && !outgoing[root].empty() &&
            (indegree[root] != 1 || outgoing[root].size() != 1)) walk(root);
    }
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (!used_node[root] && !outgoing[root].empty()) walk(root);
    }

    std::vector<uint8_t> component_seen(bundles_.size(), 0);
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (component_seen[root]) continue;
        std::queue<uint32_t> queue;
        queue.push(root);
        component_seen[root] = 1;
        size_t component_size = 0;
        while (!queue.empty()) {
            const uint32_t node = queue.front(); queue.pop();
            ++component_size;
            for (uint32_t neighbor : adjacency_[node]) {
                if (!component_seen[neighbor]) { component_seen[neighbor] = 1; queue.push(neighbor); }
            }
        }
        ++stats_.components;
        if (component_size == 1) ++stats_.isolated;
    }
    stats_.contigs = contigs_.size();
    std::vector<size_t> lengths;
    for (const std::string &contig : contigs_) lengths.push_back(contig.size());
    std::sort(lengths.begin(), lengths.end(), std::greater<size_t>());
    size_t total = std::accumulate(lengths.begin(), lengths.end(), static_cast<size_t>(0));
    size_t cumulative = 0;
    for (size_t length : lengths) { cumulative += length; if (cumulative >= total / 2) { stats_.n50 = length; break; } }
}

bool GraphAssemblerV4::write_outputs() const {
    if (mkdir(output_dir_.c_str(), 0775) != 0 && errno != EEXIST) return false;
    std::ofstream bundles(output_dir_ + "/bundles.tsv");
    bundles << "bundle_id\tconsensus_seq\tabundance\tlength\n";
    for (const V4Bundle &bundle : bundles_) bundles << bundle.id << '\t' << bundle.sequence << '\t' << bundle.abundance << '\t' << bundle.sequence.size() << '\n';
    std::ofstream edges(output_dir_ + "/graph_edges.tsv");
    edges << "source\ttarget\toverlap_length\tidentity\tqaos\torientation\tmismatches\texpected_mismatches\thigh_quality_mismatches\tpaired_support\tmean_q\tmin_q\n";
    double identity = 0, qaos = 0, overlap = 0;
    for (const V4Edge &edge : edges_) {
        edges << edge.source << '\t' << edge.target << '\t' << edge.overlap << '\t' << edge.identity << '\t' << edge.qaos << '\t' << static_cast<int>(edge.orientation) << '\t' << edge.mismatches << '\t' << edge.expected_mismatches << '\t' << edge.high_quality_mismatches << '\t' << edge.paired_support << '\t' << edge.mean_q << '\t' << edge.min_q << '\n';
        identity += edge.identity; qaos += edge.qaos; overlap += edge.overlap;
    }
    std::ofstream contigs(output_dir_ + "/assembled_contigs.fa");
    for (size_t i = 0; i < contigs_.size(); ++i) contigs << ">contig_" << i << " len=" << contigs_[i].size() << '\n' << contigs_[i] << '\n';
    std::ofstream graph_stats(output_dir_ + "/graph_statistics.tsv");
    graph_stats << "statistic\tvalue\nraw_reads\t" << stats_.raw_reads << "\nbundles\t" << stats_.bundles << "\ncandidate_pairs\t" << stats_.candidate_pairs << "\naccepted_edges\t" << stats_.accepted_edges << "\ncomponents\t" << stats_.components << "\nisolated_nodes\t" << stats_.isolated << "\nbranches\t" << stats_.branches << "\nresolved_branches\t" << stats_.resolved_branches << "\nambiguous_branches\t" << stats_.ambiguous_branches << "\ncycles\t" << stats_.cycles << "\npeak_rss_mb\t" << stats_.peak_rss_mb << '\n';
    std::ofstream assembly_stats(output_dir_ + "/assembly_statistics.tsv");
    assembly_stats << "statistic\tvalue\ncontigs\t" << stats_.contigs << "\nn50\t" << stats_.n50 << "\n";
    std::ofstream qaos_stats(output_dir_ + "/QAOS_statistics.tsv");
    qaos_stats << "statistic\tvalue\nmean_overlap\t" << (edges_.empty() ? 0 : overlap / edges_.size()) << "\nmean_identity\t" << (edges_.empty() ? 0 : identity / edges_.size()) << "\nmean_QAOS\t" << (edges_.empty() ? 0 : qaos / edges_.size()) << "\n";
    std::ofstream metadata(output_dir_ + "/run_metadata.txt");
    metadata << "implementation\tGraph-TRUST4-v4\nthreads\t" << config_.threads << "\nkmer\t" << config_.kmer << "\nmin_overlap\t" << config_.min_overlap << "\nmin_identity\t" << config_.min_identity << "\nmin_QAOS\t" << config_.min_qaos << "\nQ_HIGH\t" << config_.qhigh << "\nmax_candidates_per_bundle\t" << config_.max_candidates << "\nmax_kmer_postings\t" << config_.max_postings << "\nbatch_size\t" << config_.batch_size << "\nmax_memory_gb\t" << config_.max_memory_gb << "\nmin_branch_QAOS_delta\t" << config_.min_branch_qaos_delta << "\nmax_alternative_paths_per_branch\t" << config_.max_alternative_paths_per_branch << "\npaired_end_support\tunavailable_input_loader\npeak_rss_mb\t" << stats_.peak_rss_mb << "\n";
    std::ofstream debug(output_dir_ + "/assembly_debug.tsv");
    debug << "contig_id\tcontig_length\tnum_nodes\tunique_nodes\trepeated_nodes\ttotal_overlap\texpected_length\tmax_node_length\tcycle_detected\tvalid_path\n";
    for (const V4AssemblyDebug &record : assembly_debug_) {
        debug << record.contig_id << '\t' << record.contig_length << '\t' << record.num_nodes << '\t'
              << record.unique_nodes << '\t' << record.repeated_nodes << '\t' << record.total_overlap << '\t'
              << record.expected_length << '\t' << record.max_node_length << '\t'
              << (record.cycle_detected ? 1 : 0) << '\t' << (record.valid_path ? 1 : 0) << '\n';
    }
    std::ofstream branches(output_dir_ + "/branch_decisions.tsv");
    branches << "branch_node\tbest_edge\tsecond_edge\tbest_QAOS\tsecond_QAOS\tQAOS_delta\toverlap_best\toverlap_second\thighQ_mismatch_best\thighQ_mismatch_second\tpaired_support_best\tpaired_support_second\tdecision\treason\n";
    for (const V4BranchDecision &record : branch_decisions_) {
        branches << record.branch_node << '\t' << record.best_edge << '\t' << record.second_edge << '\t'
                 << record.best_qaos << '\t' << record.second_qaos << '\t' << record.qaos_delta << '\t'
                 << record.overlap_best << '\t' << record.overlap_second << '\t' << record.high_q_mismatch_best << '\t'
                 << record.high_q_mismatch_second << '\t' << record.paired_support_best << '\t'
                 << record.paired_support_second << '\t' << record.decision << '\t' << record.reason << '\n';
    }
    std::ofstream paths(output_dir_ + "/path_statistics.tsv");
    paths << "contig_id\tcontig_length\tnum_nodes\tunique_nodes\trepeated_nodes\ttotal_overlap\texpected_length\tcycle_detected\tvalid_path\tconsensus_confidence\n";
    for (size_t i = 0; i < assembly_debug_.size(); ++i) {
        const V4AssemblyDebug &record = assembly_debug_[i];
        paths << record.contig_id << '\t' << record.contig_length << '\t' << record.num_nodes << '\t' << record.unique_nodes << '\t'
              << record.repeated_nodes << '\t' << record.total_overlap << '\t' << record.expected_length << '\t'
              << (record.cycle_detected ? 1 : 0) << '\t' << (record.valid_path ? 1 : 0) << '\t' << record.consensus_confidence << '\n';
    }
    std::ofstream readme(output_dir_ + "/README.md");
    readme << "# Graph-TRUST4-v4\n\nThis isolated implementation preserves v3-debug memory limits, k=9 compact indexing, variable QAOS overlaps, and 64-thread batched candidate generation. v4 adds cheap candidate ranking, branch-aware bounded alternatives, lexicographic edge ranking, and Phred-weighted consensus.\n\nPaired-end support is recorded as zero/unavailable because the inherited loader does not retain mate identifiers. No biological evaluation is run automatically.\n";
    return bundles.good() && edges.good() && contigs.good() && graph_stats.good() && assembly_stats.good() && qaos_stats.good() && metadata.good() && debug.good() && branches.good() && paths.good() && readme.good();
}

int GraphAssemblerV4::run(const std::string &r1, const std::string &r2, const std::string &output_dir) {
    output_dir_ = output_dir;
    if (!load_and_bundle(r1, r2) || !build_index() || !generate_and_score()) return 1;
    build_contigs_v4();
    if (!edges_.empty()) {
        for (const V4Edge &edge : edges_) { stats_.mean_identity += edge.identity; stats_.mean_qaos += edge.qaos; stats_.mean_overlap += edge.overlap; }
        stats_.mean_identity /= edges_.size(); stats_.mean_qaos /= edges_.size(); stats_.mean_overlap /= edges_.size();
    }
    if (!write_outputs()) return 1;
    progress("Complete", stats_.bundles, stats_.bundles, stats_.candidate_pairs, stats_.accepted_edges);
    return 0;
}
