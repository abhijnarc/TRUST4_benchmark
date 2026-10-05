#include "GraphAssemblerV3.hpp"

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <numeric>
#include <queue>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unordered_map>

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

} // namespace

GraphAssemblerV3::GraphAssemblerV3(const V3Config &config) : config_(config), start_seconds_(now_seconds()) {}

uint64_t GraphAssemblerV3::rss_mb() const {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return static_cast<uint64_t>(usage.ru_maxrss / 1024);
}

bool GraphAssemblerV3::memory_ok(const char *stage) {
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

void GraphAssemblerV3::progress(const char *stage, uint64_t done, uint64_t total,
                                uint64_t candidates, uint64_t accepted) const {
    std::fprintf(stderr, "[%s] bundles=%llu/%llu threads=%d candidate_pairs=%llu accepted_edges=%llu RSS=%llu MB elapsed=%.1fs\n",
                 stage, static_cast<unsigned long long>(done), static_cast<unsigned long long>(total),
                 config_.threads, static_cast<unsigned long long>(candidates),
                 static_cast<unsigned long long>(accepted), static_cast<unsigned long long>(rss_mb()),
                 now_seconds() - start_seconds_);
}

bool GraphAssemblerV3::load_and_bundle(const std::string &r1, const std::string &r2) {
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
        V3Bundle bundle;
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

bool GraphAssemblerV3::build_index() {
    std::fprintf(stderr, "[Index] building compact k-mer index k=%d max_postings=%u\n",
                 config_.kmer, config_.max_postings);
    index_.reserve(bundles_.size() * 2);
    for (const V3Bundle &bundle : bundles_) {
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

bool GraphAssemblerV3::score_overlap(uint32_t a, uint32_t b, V3Edge &edge) const {
    const V3Bundle &left = bundles_[a];
    const V3Bundle &right = bundles_[b];
    const int maximum = static_cast<int>(std::min(left.sequence.size(), right.sequence.size()));
    int best_length = 0, best_mismatches = 0, best_high = 0;
    float best_qaos = 0.0f, best_expected = 0.0f;
    for (int length = config_.min_overlap; length <= maximum; ++length) {
        const int start = static_cast<int>(left.sequence.size()) - length;
        double log_probability = 0.0;
        float expected = 0.0f;
        float quality_sum = 0.0f;
        int mismatches = 0, high = 0;
        for (int i = 0; i < length; ++i) {
            const char base1 = left.sequence[start + i], base2 = right.sequence[i];
            const int q1 = phred(left.quality[start + i]), q2 = phred(right.quality[i]);
            const double probability = std::max(1e-30, static_cast<double>(compatibility(base1, q1, base2, q2)));
            log_probability += std::log(probability);
            quality_sum += (q1 + q2) / 2.0f;
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
    return true;
}

bool GraphAssemblerV3::generate_and_score() {
    const uint64_t total = bundles_.size();
    const uint32_t batch_size = std::max<uint32_t>(1, config_.batch_size);
    for (uint64_t begin = 0; begin < total; begin += batch_size) {
        const uint64_t end = std::min(total, begin + batch_size);
        const int workers = std::max(1, config_.threads);
        std::vector<std::vector<V3Edge>> local_edges(static_cast<size_t>(workers));
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
                const V3Bundle &bundle = bundles_[static_cast<size_t>(id)];
                if (bundle.sequence.size() >= static_cast<size_t>(config_.kmer)) {
                    const size_t first = bundle.sequence.size() > 100 ? bundle.sequence.size() - 100 : 0;
                    for (size_t pos = first; pos + config_.kmer <= bundle.sequence.size(); ++pos) {
                        const std::string kmer = bundle.sequence.substr(pos, config_.kmer);
                        auto found = index_.find(kmer);
                        if (found == index_.end()) continue;
                        candidates.insert(candidates.end(), found->second.begin(), found->second.end());
                    }
                }
                std::sort(candidates.begin(), candidates.end());
                candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
                candidates.erase(std::remove(candidates.begin(), candidates.end(), static_cast<uint32_t>(id)), candidates.end());
                if (candidates.size() > config_.max_candidates) candidates.resize(config_.max_candidates);
                local_candidates[thread] += candidates.size();
                for (uint32_t candidate : candidates) {
                    V3Edge edge;
                    if (score_overlap(static_cast<uint32_t>(id), candidate, edge)) {
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

void GraphAssemblerV3::build_contigs() {
    adjacency_.assign(bundles_.size(), std::vector<uint32_t>());
    for (const V3Edge &edge : edges_) {
        adjacency_[edge.source].push_back(edge.target);
        adjacency_[edge.target].push_back(edge.source);
    }
    std::vector<uint8_t> visited(bundles_.size(), 0);
    for (uint32_t root = 0; root < bundles_.size(); ++root) {
        if (visited[root]) continue;
        std::queue<uint32_t> queue;
        std::vector<uint32_t> component;
        queue.push(root);
        visited[root] = 1;
        while (!queue.empty()) {
            uint32_t node = queue.front(); queue.pop();
            component.push_back(node);
            for (uint32_t neighbor : adjacency_[node]) {
                if (!visited[neighbor]) { visited[neighbor] = 1; queue.push(neighbor); }
            }
        }
        ++stats_.components;
        if (component.size() == 1) { ++stats_.isolated; continue; }
        std::sort(component.begin(), component.end());
        std::string contig;
        for (uint32_t id : component) contig += bundles_[id].sequence;
        contigs_.push_back(std::move(contig));
    }
    stats_.contigs = contigs_.size();
    std::vector<size_t> lengths;
    for (const std::string &contig : contigs_) lengths.push_back(contig.size());
    std::sort(lengths.begin(), lengths.end(), std::greater<size_t>());
    size_t total = std::accumulate(lengths.begin(), lengths.end(), static_cast<size_t>(0));
    size_t cumulative = 0;
    for (size_t length : lengths) { cumulative += length; if (cumulative >= total / 2) { stats_.n50 = length; break; } }
}

bool GraphAssemblerV3::write_outputs() const {
    if (mkdir(output_dir_.c_str(), 0775) != 0 && errno != EEXIST) return false;
    std::ofstream bundles(output_dir_ + "/bundles.tsv");
    bundles << "bundle_id\tconsensus_seq\tabundance\tlength\n";
    for (const V3Bundle &bundle : bundles_) bundles << bundle.id << '\t' << bundle.sequence << '\t' << bundle.abundance << '\t' << bundle.sequence.size() << '\n';
    std::ofstream edges(output_dir_ + "/graph_edges.tsv");
    edges << "source\ttarget\toverlap_length\tidentity\tqaos\torientation\tmismatches\texpected_mismatches\thigh_quality_mismatches\n";
    double identity = 0, qaos = 0, overlap = 0;
    for (const V3Edge &edge : edges_) {
        edges << edge.source << '\t' << edge.target << '\t' << edge.overlap << '\t' << edge.identity << '\t' << edge.qaos << '\t' << static_cast<int>(edge.orientation) << '\t' << edge.mismatches << '\t' << edge.expected_mismatches << '\t' << edge.high_quality_mismatches << '\n';
        identity += edge.identity; qaos += edge.qaos; overlap += edge.overlap;
    }
    std::ofstream contigs(output_dir_ + "/assembled_contigs.fa");
    for (size_t i = 0; i < contigs_.size(); ++i) contigs << ">contig_" << i << " len=" << contigs_[i].size() << '\n' << contigs_[i] << '\n';
    std::ofstream graph_stats(output_dir_ + "/graph_statistics.tsv");
    graph_stats << "statistic\tvalue\nraw_reads\t" << stats_.raw_reads << "\nbundles\t" << stats_.bundles << "\ncandidate_pairs\t" << stats_.candidate_pairs << "\naccepted_edges\t" << stats_.accepted_edges << "\ncomponents\t" << stats_.components << "\nisolated_nodes\t" << stats_.isolated << "\npeak_rss_mb\t" << stats_.peak_rss_mb << '\n';
    std::ofstream assembly_stats(output_dir_ + "/assembly_statistics.tsv");
    assembly_stats << "statistic\tvalue\ncontigs\t" << stats_.contigs << "\nn50\t" << stats_.n50 << "\n";
    std::ofstream qaos_stats(output_dir_ + "/QAOS_statistics.tsv");
    qaos_stats << "statistic\tvalue\nmean_overlap\t" << (edges_.empty() ? 0 : overlap / edges_.size()) << "\nmean_identity\t" << (edges_.empty() ? 0 : identity / edges_.size()) << "\nmean_QAOS\t" << (edges_.empty() ? 0 : qaos / edges_.size()) << "\n";
    std::ofstream metadata(output_dir_ + "/run_metadata.txt");
    metadata << "implementation\tGraph-TRUST4-v3\nthreads\t" << config_.threads << "\nkmer\t" << config_.kmer << "\nmin_overlap\t" << config_.min_overlap << "\nmax_candidates_per_bundle\t" << config_.max_candidates << "\nmax_kmer_postings\t" << config_.max_postings << "\nbatch_size\t" << config_.batch_size << "\nmax_memory_gb\t" << config_.max_memory_gb << "\npeak_rss_mb\t" << stats_.peak_rss_mb << "\n";
    return bundles.good() && edges.good() && contigs.good() && graph_stats.good() && assembly_stats.good() && qaos_stats.good() && metadata.good();
}

int GraphAssemblerV3::run(const std::string &r1, const std::string &r2, const std::string &output_dir) {
    output_dir_ = output_dir;
    if (!load_and_bundle(r1, r2) || !build_index() || !generate_and_score()) return 1;
    build_contigs();
    if (!edges_.empty()) {
        for (const V3Edge &edge : edges_) { stats_.mean_identity += edge.identity; stats_.mean_qaos += edge.qaos; stats_.mean_overlap += edge.overlap; }
        stats_.mean_identity /= edges_.size(); stats_.mean_qaos /= edges_.size(); stats_.mean_overlap /= edges_.size();
    }
    if (!write_outputs()) return 1;
    progress("Complete", stats_.bundles, stats_.bundles, stats_.candidate_pairs, stats_.accepted_edges);
    return 0;
}
