#include "GraphAssemblerV2.hpp"
#include <algorithm>
#include <numeric>
#include <sys/resource.h>
#ifdef _OPENMP
#include <omp.h>
#endif

GraphAssemblerV2::GraphAssemblerV2(const GraphConfigV2 &cfg) : config_(cfg) {}

GraphAssemblerV2::~GraphAssemblerV2() {}

// Memory monitoring
long GraphAssemblerV2::GetRSS() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss / 1024;  // Convert to MB
}

int GraphAssemblerV2::CheckMemory(const char *stage) {
    long rss = GetRSS();
    stats_.peak_rss_mb = std::max(stats_.peak_rss_mb, rss);

    fprintf(stderr, "[Memory] %s: %ld MB / %d GB max\n", stage, rss, config_.max_memory_gb);

    if (rss > config_.max_memory_gb * 1024) {
        fprintf(stderr, "ERROR: Memory limit exceeded! (%ld MB > %d GB)\n", rss, config_.max_memory_gb);
        return -1;
    }
    return 0;
}

void GraphAssemblerV2::PrintProgress(const char *stage, int reads, int bundles, int candidates, int edges) {
    long rss = GetRSS();
    fprintf(stderr, "[Progress] %s: reads=%d, bundles=%d, candidates=%d, edges=%d, RSS=%ld MB\n",
            stage, reads, bundles, candidates, edges, rss);
}

float GraphAssemblerV2::PhredToErrorProb(int q) {
    if (q <= 0) return 1.0;
    return pow(10.0, -q / 10.0);
}

float GraphAssemblerV2::ComputeBaseCompatibility(char base1, int q1, char base2, int q2) {
    float e1 = PhredToErrorProb(q1);
    float e2 = PhredToErrorProb(q2);

    if (base1 == base2) {
        return (1.0 - e1) * (1.0 - e2);
    } else {
        float prob_same = e1 / 3.0 + e2 / 3.0 - e1 * e2 / 9.0;
        return prob_same;
    }
}

int GraphAssemblerV2::LoadReadsWithQualities(const char *fastq1_path, const char *fastq2_path) {
    fprintf(stderr, "[LoadReads] Loading paired-end reads...\n");

    FILE *f1 = fopen(fastq1_path, "r");
    FILE *f2 = fopen(fastq2_path, "r");
    if (!f1 || !f2) {
        fprintf(stderr, "ERROR: Cannot open FASTQ files\n");
        return -1;
    }

    char line[4096];
    int read_count = 0;
    int pair_count = 0;

    while (fgets(line, sizeof(line), f1) && fgets(line, sizeof(line), f2)) {
        if (read_count % 4 == 0) {
            char seq1[4096], qual1[4096], seq2[4096], qual2[4096];
            fgets(seq1, sizeof(seq1), f1);
            fgets(line, sizeof(line), f1);
            fgets(qual1, sizeof(qual1), f1);

            fgets(seq2, sizeof(seq2), f2);
            fgets(line, sizeof(line), f2);
            fgets(qual2, sizeof(qual2), f2);

            seq1[strcspn(seq1, "\n")] = 0;
            qual1[strcspn(qual1, "\n")] = 0;
            seq2[strcspn(seq2, "\n")] = 0;
            qual2[strcspn(qual2, "\n")] = 0;

            if (strlen(seq1) > 0) {
                raw_sequences_.push_back(std::string(seq1));
                raw_qualities_.push_back(std::string(qual1));
            }
            if (strlen(seq2) > 0) {
                raw_sequences_.push_back(std::string(seq2));
                raw_qualities_.push_back(std::string(qual2));
            }
            pair_count++;

            if (pair_count % 50000 == 0) {
                PrintProgress("LoadReads", raw_sequences_.size(), 0, 0, 0);
            }
        }
        read_count++;
    }

    fclose(f1);
    fclose(f2);

    stats_.total_raw_reads = raw_sequences_.size();
    CheckMemory("AfterLoadReads");
    fprintf(stderr, "[LoadReads] Loaded %d reads (%d pairs)\n", stats_.total_raw_reads, pair_count);

    return 0;
}

int GraphAssemblerV2::BundleRedundantReads() {
    fprintf(stderr, "[BundleReads] Bundling highly redundant reads...\n");

    std::unordered_map<std::string, std::vector<int>> seq_to_reads;

    for (size_t i = 0; i < raw_sequences_.size(); ++i) {
        seq_to_reads[raw_sequences_[i]].push_back(i);
    }

    for (auto &pair : seq_to_reads) {
        ReadBundle bundle;
        bundle.bundle_id = bundles_.size();
        bundle.consensus_seq = pair.first;
        bundle.abundance = pair.second.size();
        // NOTE: read_ids NOT stored (saves memory)

        if (pair.second.size() > 0) {
            bundle.qualities = raw_qualities_[pair.second[0]];
        }

        // REMOVED: per_position_counts (was ~190 GB for 634K bundles!)
        // REMOVED: per_position_quals (compute on-demand if needed)

        bundles_.push_back(bundle);

        if (bundles_.size() % 50000 == 0) {
            PrintProgress("BundleReads", stats_.total_raw_reads, bundles_.size(), 0, 0);
        }
    }

    // Clear raw reads to free memory
    raw_sequences_.clear();
    raw_qualities_.clear();

    stats_.total_bundles = bundles_.size();
    stats_.compression_ratio = (float)stats_.total_raw_reads / stats_.total_bundles;

    CheckMemory("AfterBundling");
    fprintf(stderr, "[BundleReads] Created %d bundles (compression: %.2fx)\n",
            stats_.total_bundles, stats_.compression_ratio);

    return 0;
}

int GraphAssemblerV2::BuildKmerIndex() {
    fprintf(stderr, "[BuildIndex] Building k-mer index (k=%d)...\n", config_.kmer_size);

    kmer_index_.clear();

    for (int bundle_id = 0; bundle_id < (int)bundles_.size(); ++bundle_id) {
        const std::string &seq = bundles_[bundle_id].consensus_seq;
        int klen = config_.kmer_size;

        for (int i = 0; i <= (int)seq.length() - klen; ++i) {
            std::string kmer = seq.substr(i, klen);
            kmer_index_[kmer].push_back(bundle_id);
        }

        if ((bundle_id + 1) % 50000 == 0) {
            PrintProgress("BuildIndex", stats_.total_raw_reads, bundle_id + 1, 0, 0);
        }
    }

    CheckMemory("AfterKmerIndex");
    fprintf(stderr, "[BuildIndex] Built index with %zu distinct k-mers\n", kmer_index_.size());
    return 0;
}

int GraphAssemblerV2::GenerateCandidateBundlePairs() {
    fprintf(stderr, "[Candidates] Generating candidate pairs (threads=%d)...\n", config_.num_threads);

    candidate_bundle_pairs_.clear();
    std::vector<std::vector<std::pair<int, int>>> thread_local_pairs(config_.num_threads);

    time_t start_time = time(NULL);

    #pragma omp parallel num_threads(config_.num_threads)
    {
        int thread_id = 0;
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        #endif

        std::vector<std::pair<int, int>> &local_pairs = thread_local_pairs[thread_id];

        #pragma omp for schedule(dynamic)
        for (int b1 = 0; b1 < (int)bundles_.size(); ++b1) {
            const std::string &seq1 = bundles_[b1].consensus_seq;
            std::set<int> candidates;

            int klen = config_.kmer_size;
            for (int i = std::max(0, (int)seq1.length() - 100); i < (int)seq1.length() - klen + 1; ++i) {
                std::string kmer = seq1.substr(i, klen);
                if (kmer_index_.count(kmer)) {
                    for (int b2 : kmer_index_[kmer]) {
                        if (b2 != b1) {
                            candidates.insert(b2);
                        }
                    }
                }
            }

            // LIMIT candidates per bundle to control memory
            if (candidates.size() > (size_t)config_.max_candidates_per_bundle) {
                auto it = candidates.begin();
                std::advance(it, config_.max_candidates_per_bundle);
                candidates.erase(it, candidates.end());
            }

            // Add to thread-local buffer (no synchronization needed)
            for (int b2 : candidates) {
                local_pairs.push_back({b1, b2});
            }
        }
    }

    // Merge thread-local results with global deduplication
    fprintf(stderr, "[Candidates] Merging %d thread-local buffers...\n", config_.num_threads);

    std::set<std::pair<int, int>> seen_pairs;
    long total_generated = 0;

    for (int t = 0; t < config_.num_threads; ++t) {
        for (const auto &pair : thread_local_pairs[t]) {
            int a = pair.first, b = pair.second;
            if (a > b) std::swap(a, b);

            if (!seen_pairs.count({a, b})) {
                candidate_bundle_pairs_.push_back(pair);
                seen_pairs.insert({a, b});
            }
            total_generated++;
        }
    }

    time_t end_time = time(NULL);
    double elapsed = difftime(end_time, start_time);

    stats_.total_candidate_pairs = candidate_bundle_pairs_.size();
    CheckMemory("AfterCandidateGeneration");
    fprintf(stderr, "[Candidates] Generated %ld pairs (%ld unique, %.1f sec, threads=%d)\n",
            total_generated, candidate_bundle_pairs_.size(), elapsed, config_.num_threads);

    return 0;
}

int GraphAssemblerV2::ComputeQAOS(const std::string &seq1, const std::string &qual1,
                                  const std::string &seq2, const std::string &qual2,
                                  int start1, int start2, int overlap_len,
                                  float &qaos, int &mismatches, float &expected_mismatch,
                                  int &high_q_mismatches, float &mean_q, float &min_q) {
    if (start1 < 0 || start2 < 0 || overlap_len <= 0) return -1;
    if (start1 + overlap_len > (int)seq1.length()) return -1;
    if (start2 + overlap_len > (int)seq2.length()) return -1;

    float product = 1.0;
    mismatches = 0;
    expected_mismatch = 0;
    high_q_mismatches = 0;

    // Use stack array instead of vector to reduce allocation overhead
    float qual_sum = 0;
    float min_qual = 60;

    for (int i = 0; i < overlap_len; ++i) {
        char b1 = seq1[start1 + i];
        char b2 = seq2[start2 + i];
        int q1 = (int)qual1[start1 + i] - 33;
        int q2 = (int)qual2[start2 + i] - 33;

        q1 = std::max(0, std::min(60, q1));
        q2 = std::max(0, std::min(60, q2));

        float compat = ComputeBaseCompatibility(b1, q1, b2, q2);
        product *= compat;

        float e1 = PhredToErrorProb(q1);
        float e2 = PhredToErrorProb(q2);
        if (b1 != b2) {
            mismatches++;
            expected_mismatch += e1 + e2 - e1 * e2;
            if (q1 >= config_.qhigh_threshold && q2 >= config_.qhigh_threshold) {
                high_q_mismatches++;
            }
        }

        float q_avg = (q1 + q2) / 2.0;
        qual_sum += q_avg;
        min_qual = std::min(min_qual, q_avg);
    }

    qaos = pow(product, 1.0 / overlap_len);
    mean_q = qual_sum / overlap_len;
    min_q = min_qual;

    return 0;
}

int GraphAssemblerV2::ComputeVariableLengthOverlap(int bundle1_id, int bundle2_id,
                                                     QualityAwareOverlap &overlap) {
    const ReadBundle &b1 = bundles_[bundle1_id];
    const ReadBundle &b2 = bundles_[bundle2_id];

    const std::string &seq1 = b1.consensus_seq;
    const std::string &seq2 = b2.consensus_seq;
    const std::string &qual1 = b1.qualities;
    const std::string &qual2 = b2.qualities;

    int best_overlap_len = 0;
    float best_qaos = 0;
    int best_mismatches = 0;
    float best_expected = 0;
    int best_high_q_mismatches = 0;
    float best_mean_q = 0;
    float best_min_q = 0;
    int best_start1 = 0, best_start2 = 0;

    int min_overlap = config_.min_overlap_length;
    int max_overlap = std::min((int)seq1.length(), (int)seq2.length());

    for (int overlap_len = min_overlap; overlap_len <= max_overlap; ++overlap_len) {
        int start1 = seq1.length() - overlap_len;
        int start2 = 0;

        float qaos;
        int mismatches;
        float expected_mismatch;
        int high_q_mismatches;
        float mean_q, min_q;

        if (ComputeQAOS(seq1, qual1, seq2, qual2, start1, start2, overlap_len,
                        qaos, mismatches, expected_mismatch, high_q_mismatches, mean_q, min_q) < 0) {
            continue;
        }

        float identity = 1.0 - (float)mismatches / overlap_len;

        if (identity >= config_.min_identity_threshold && qaos >= config_.min_qaos_threshold) {
            if (qaos > best_qaos || (qaos == best_qaos && overlap_len > best_overlap_len)) {
                best_overlap_len = overlap_len;
                best_qaos = qaos;
                best_mismatches = mismatches;
                best_expected = expected_mismatch;
                best_high_q_mismatches = high_q_mismatches;
                best_mean_q = mean_q;
                best_min_q = min_q;
                best_start1 = start1;
                best_start2 = start2;
            }
        }
    }

    if (best_overlap_len > 0) {
        overlap.bundle1_id = bundle1_id;
        overlap.bundle2_id = bundle2_id;
        overlap.overlap_length = best_overlap_len;
        overlap.qaos = best_qaos;
        overlap.observed_identity = 1.0 - (float)best_mismatches / best_overlap_len;
        overlap.observed_mismatches = best_mismatches;
        overlap.expected_mismatches = best_expected;
        overlap.high_quality_mismatches = best_high_q_mismatches;
        overlap.mean_q = best_mean_q;
        overlap.min_q = best_min_q;
        overlap.bundle1_end = best_start1 + best_overlap_len;
        overlap.bundle2_start = best_start2;
        overlap.abundance1 = b1.abundance;
        overlap.abundance2 = b2.abundance;
        overlap.orientation = 1;

        return 1;
    }

    return 0;
}

int GraphAssemblerV2::ComputeQualityAwareOverlaps() {
    fprintf(stderr, "[Overlaps] Computing quality-aware overlaps for %zu candidates...\n",
            candidate_bundle_pairs_.size());

    edges_.clear();
    int valid_overlaps = 0;

    for (size_t i = 0; i < candidate_bundle_pairs_.size(); ++i) {
        int b1 = candidate_bundle_pairs_[i].first;
        int b2 = candidate_bundle_pairs_[i].second;

        QualityAwareOverlap overlap;
        if (ComputeVariableLengthOverlap(b1, b2, overlap) > 0) {
            edges_.push_back(overlap);
            valid_overlaps++;
        }

        if ((i + 1) % 100000 == 0) {
            PrintProgress("Overlaps", stats_.total_raw_reads, stats_.total_bundles,
                         i + 1, valid_overlaps);
            if (CheckMemory("DuringOverlapComputation") < 0) return -1;
        }
    }

    stats_.total_edges = edges_.size();

    CheckMemory("AfterOverlapComputation");
    fprintf(stderr, "[Overlaps] Found %d valid overlaps (%.1f%% of candidates)\n",
            valid_overlaps, 100.0 * valid_overlaps / candidate_bundle_pairs_.size());

    return 0;
}

int GraphAssemblerV2::BuildGraph() {
    fprintf(stderr, "[Graph] Building adjacency structure...\n");

    adjacency_list_.clear();
    std::set<int> unique_nodes;

    for (const auto &edge : edges_) {
        adjacency_list_[edge.bundle1_id].push_back(edge.bundle2_id);
        adjacency_list_[edge.bundle2_id].push_back(edge.bundle1_id);
        unique_nodes.insert(edge.bundle1_id);
        unique_nodes.insert(edge.bundle2_id);
    }

    stats_.total_nodes = unique_nodes.size();
    CheckMemory("AfterGraphConstruction");
    fprintf(stderr, "[Graph] Graph has %d nodes and %d edges\n", stats_.total_nodes, stats_.total_edges);

    return 0;
}

int GraphAssemblerV2::CleanGraph() {
    fprintf(stderr, "[Clean] Cleaning graph...\n");

    std::vector<QualityAwareOverlap> filtered_edges;

    for (const auto &edge : edges_) {
        if (edge.overlap_length >= config_.min_overlap_length &&
            edge.observed_identity >= config_.min_identity_threshold &&
            edge.qaos >= config_.min_qaos_threshold) {
            filtered_edges.push_back(edge);
        }
    }

    edges_ = filtered_edges;
    stats_.total_edges = edges_.size();

    fprintf(stderr, "[Clean] After filtering: %d edges\n", stats_.total_edges);

    return 0;
}

int GraphAssemblerV2::ExtractContigs() {
    fprintf(stderr, "[Extract] Extracting contigs from connected components...\n");

    return BuildConnectedComponents();
}

int GraphAssemblerV2::BuildConnectedComponents() {
    std::set<int> visited;
    std::vector<std::vector<int>> components;

    for (int node = 0; node < (int)bundles_.size(); ++node) {
        if (visited.count(node)) continue;

        std::vector<int> component;
        std::queue<int> q;
        q.push(node);
        visited.insert(node);

        while (!q.empty()) {
            int curr = q.front();
            q.pop();
            component.push_back(curr);

            if (adjacency_list_.count(curr)) {
                for (int neighbor : adjacency_list_[curr]) {
                    if (!visited.count(neighbor)) {
                        visited.insert(neighbor);
                        q.push(neighbor);
                    }
                }
            }
        }

        components.push_back(component);
    }

    stats_.num_connected_components = components.size();
    stats_.largest_component_size = 0;
    stats_.num_isolated_nodes = 0;

    for (const auto &comp : components) {
        if (comp.size() == 1) {
            stats_.num_isolated_nodes++;
        }
        stats_.largest_component_size = std::max(stats_.largest_component_size, (int)comp.size());

        if (comp.size() > 1) {
            TraverseComponentForContigs(comp);
        }
    }

    fprintf(stderr, "[Extract] Found %d components (largest: %d nodes, isolated: %d)\n",
            stats_.num_connected_components, stats_.largest_component_size, stats_.num_isolated_nodes);

    return 0;
}

int GraphAssemblerV2::TraverseComponentForContigs(const std::vector<int> &component_nodes) {
    AssembledContig contig;
    contig.contig_id = contigs_.size();
    contig.bundle_ids = component_nodes;

    std::string combined_seq;
    for (int bundle_id : component_nodes) {
        if (combined_seq.empty()) {
            combined_seq = bundles_[bundle_id].consensus_seq;
        } else {
            combined_seq += bundles_[bundle_id].consensus_seq;
        }
    }

    contig.consensus = combined_seq;
    contig.contig_length = combined_seq.length();

    if (GenerateConsensus(contig) < 0) return -1;

    contigs_.push_back(contig);
    return 0;
}

int GraphAssemblerV2::GenerateConsensus(AssembledContig &contig) {
    return 0;
}

int GraphAssemblerV2::GenerateConsensusSequences() {
    fprintf(stderr, "[Consensus] Generating consensus sequences...\n");

    stats_.num_assembled_contigs = contigs_.size();

    std::vector<int> contig_lengths;
    for (const auto &contig : contigs_) {
        contig_lengths.push_back(contig.contig_length);
        if (contig.contig_length >= 150) stats_.num_contigs_150bp++;
        if (contig.contig_length >= 300) stats_.num_contigs_300bp++;
    }

    if (!contig_lengths.empty()) {
        stats_.avg_contig_length = std::accumulate(contig_lengths.begin(), contig_lengths.end(), 0) / (double)contig_lengths.size();

        std::sort(contig_lengths.begin(), contig_lengths.end());
        int total = 0;
        for (int len : contig_lengths) {
            total += len;
        }
        int target = total / 2;
        total = 0;
        for (int len : contig_lengths) {
            total += len;
            if (total >= target) {
                stats_.n50_contig_length = len;
                break;
            }
        }
    }

    fprintf(stderr, "[Consensus] %d assembled contigs, N50=%.0f, avg=%.0f\n",
            stats_.num_assembled_contigs, stats_.n50_contig_length, stats_.avg_contig_length);

    if (edges_.size() > 0) {
        float total_overlap = 0, total_identity = 0, total_qaos = 0;
        for (const auto &edge : edges_) {
            total_overlap += edge.overlap_length;
            total_identity += edge.observed_identity;
            total_qaos += edge.qaos;
        }
        stats_.mean_overlap_length = total_overlap / edges_.size();
        stats_.mean_identity = total_identity / edges_.size();
        stats_.mean_qaos = total_qaos / edges_.size();
    }

    return 0;
}

int GraphAssemblerV2::WriteBundles(const char *output_path) {
    fprintf(stderr, "[Output] Writing bundles to %s\n", output_path);

    FILE *f = fopen(output_path, "w");
    if (!f) return -1;

    fprintf(f, "bundle_id\tconsensus_seq\tabundance\n");
    for (const auto &bundle : bundles_) {
        fprintf(f, "%d\t%s\t%d\n", bundle.bundle_id, bundle.consensus_seq.c_str(), bundle.abundance);
    }

    fclose(f);
    return 0;
}

int GraphAssemblerV2::WriteGraphEdges(const char *output_path) {
    fprintf(stderr, "[Output] Writing graph edges to %s\n", output_path);

    FILE *f = fopen(output_path, "w");
    if (!f) return -1;

    fprintf(f, "source\ttarget\torientation\toverlap_length\tidentity\tqaos\t");
    fprintf(f, "observed_mismatches\texpected_mismatches\thigh_quality_mismatches\t");
    fprintf(f, "bundle1_abundance\tbundle2_abundance\n");

    for (const auto &edge : edges_) {
        fprintf(f, "%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%.2f\t%d\t%d\t%d\n",
                edge.bundle1_id, edge.bundle2_id, edge.orientation,
                edge.overlap_length, edge.observed_identity, edge.qaos,
                edge.observed_mismatches, edge.expected_mismatches,
                edge.high_quality_mismatches,
                edge.abundance1, edge.abundance2);
    }

    fclose(f);
    return 0;
}

int GraphAssemblerV2::WriteGraphStats(const char *stats_path) {
    fprintf(stderr, "[Output] Writing graph statistics to %s\n", stats_path);

    FILE *f = fopen(stats_path, "w");
    if (!f) return -1;

    fprintf(f, "statistic\tvalue\n");
    fprintf(f, "total_raw_reads\t%d\n", stats_.total_raw_reads);
    fprintf(f, "total_bundles\t%d\n", stats_.total_bundles);
    fprintf(f, "compression_ratio\t%.2f\n", stats_.compression_ratio);
    fprintf(f, "total_candidate_pairs\t%d\n", stats_.total_candidate_pairs);
    fprintf(f, "total_nodes\t%d\n", stats_.total_nodes);
    fprintf(f, "total_edges\t%d\n", stats_.total_edges);
    fprintf(f, "num_connected_components\t%d\n", stats_.num_connected_components);
    fprintf(f, "largest_component_size\t%d\n", stats_.largest_component_size);
    fprintf(f, "num_isolated_nodes\t%d\n", stats_.num_isolated_nodes);
    fprintf(f, "peak_rss_mb\t%ld\n", stats_.peak_rss_mb);

    fclose(f);
    return 0;
}

int GraphAssemblerV2::WriteContigsFasta(const char *output_path) {
    fprintf(stderr, "[Output] Writing contigs to %s\n", output_path);

    FILE *f = fopen(output_path, "w");
    if (!f) return -1;

    for (int i = 0; i < (int)contigs_.size(); ++i) {
        fprintf(f, ">contig_%d len=%d\n", i, contigs_[i].contig_length);
        fprintf(f, "%s\n", contigs_[i].consensus.c_str());
    }

    fclose(f);
    return 0;
}

int GraphAssemblerV2::WriteAssemblyStats(const char *stats_path) {
    fprintf(stderr, "[Output] Writing assembly statistics to %s\n", stats_path);

    FILE *f = fopen(stats_path, "w");
    if (!f) return -1;

    fprintf(f, "statistic\tvalue\n");
    fprintf(f, "num_assembled_contigs\t%d\n", stats_.num_assembled_contigs);
    fprintf(f, "mean_contig_length\t%.2f\n", stats_.avg_contig_length);
    fprintf(f, "n50_contig_length\t%.0f\n", stats_.n50_contig_length);
    fprintf(f, "num_contigs_150bp\t%d\n", stats_.num_contigs_150bp);
    fprintf(f, "num_contigs_300bp\t%d\n", stats_.num_contigs_300bp);

    fclose(f);
    return 0;
}

int GraphAssemblerV2::WriteQAOSStats(const char *stats_path) {
    fprintf(stderr, "[Output] Writing QAOS statistics to %s\n", stats_path);

    FILE *f = fopen(stats_path, "w");
    if (!f) return -1;

    fprintf(f, "statistic\tvalue\n");
    fprintf(f, "mean_overlap_length\t%.2f\n", stats_.mean_overlap_length);
    fprintf(f, "mean_identity\t%.4f\n", stats_.mean_identity);
    fprintf(f, "mean_qaos\t%.4f\n", stats_.mean_qaos);

    fclose(f);
    return 0;
}
