#include "GraphAssembler.hpp"
#include <stdarg.h>
#include <set>
#include <map>

// ==================================================
// Logging Helper
// ==================================================

void GraphAssembler::PrintLog(const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);

    time_t mytime = time(NULL);
    struct tm *localT = localtime(&mytime);
    char stime[500];
    strftime(stime, sizeof(stime), "%Y-%m-%d %H:%M:%S", localT);

    fprintf(stderr, "[%s] ", stime);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    fflush(stderr);

    va_end(args);
}

// ==================================================
// Constructor & Destructor
// ==================================================

GraphAssembler::GraphAssembler(const GraphConfig &cfg) : config_(cfg) {
}

GraphAssembler::~GraphAssembler() {
}

// ==================================================
// Load FASTQ Paired-End Reads
// ==================================================

int GraphAssembler::LoadReads(const char *fastq1_path, const char *fastq2_path) {
    PrintLog("Loading reads from %s and %s...", fastq1_path, fastq2_path);

    FILE *f1 = fopen(fastq1_path, "r");
    FILE *f2 = fopen(fastq2_path, "r");
    if (!f1 || !f2) {
        PrintLog("ERROR: Cannot open FASTQ files");
        return -1;
    }

    char line1[100001];
    char line2[100001];
    int read_count = 0;

    while (true) {
        // Read header lines
        if (!fgets(line1, sizeof(line1), f1) || !fgets(line2, sizeof(line2), f2)) break;

        // Read sequence lines
        if (!fgets(line1, sizeof(line1), f1)) break;
        if (!fgets(line2, sizeof(line2), f2)) break;

        // Remove newlines
        int len1 = strlen(line1);
        while (len1 > 0 && (line1[len1-1] == '\n' || line1[len1-1] == '\r')) len1--;
        line1[len1] = '\0';

        // Skip quality lines
        if (!fgets(line2, sizeof(line2), f1)) break;
        if (!fgets(line2, sizeof(line2), f2)) break;

        // Store read
        reads_.emplace_back(read_count, std::string(line1));
        read_count++;

        if (read_count % 100000 == 0) {
            PrintLog("Loaded %d reads...", read_count);
        }
    }

    fclose(f1);
    fclose(f2);

    PrintLog("Total reads loaded: %d", read_count);
    stats_.total_reads = read_count;
    stats_.total_nodes = read_count;

    return 0;
}

// ==================================================
// Deduplicate Identical Reads
// ==================================================

int GraphAssembler::DeduplicateReads() {
    PrintLog("Deduplicating reads...");

    std::map<std::string, int> seq_to_id;
    std::vector<Read> deduplicated;
    int dedup_count = 0;

    for (auto &r : reads_) {
        if (seq_to_id.count(r.seq)) {
            deduplicated[seq_to_id[r.seq]].frequency++;
            dedup_count++;
        } else {
            r.id = deduplicated.size();
            seq_to_id[r.seq] = r.id;
            deduplicated.push_back(r);
        }
    }

    reads_ = std::move(deduplicated);
    stats_.total_reads = reads_.size();
    stats_.total_nodes = reads_.size();

    PrintLog("After deduplication: %d unique reads (removed %d duplicates)",
            (int)reads_.size(), dedup_count);

    return 0;
}

// ==================================================
// Build K-mer Index
// ==================================================

int GraphAssembler::BuildKmerIndex() {
    PrintLog("Building k-mer index (k=%d)...", config_.kmer_size);

    int indexed_reads = 0;
    for (const auto &r : reads_) {
        if ((int)r.seq.length() < config_.kmer_size) continue;

        for (int i = 0; i <= (int)r.seq.length() - config_.kmer_size; ++i) {
            std::string kmer = r.seq.substr(i, config_.kmer_size);
            kmer_index_[kmer].push_back({r.id, i});
        }

        indexed_reads++;
        if (indexed_reads % 100000 == 0) {
            PrintLog("Indexed %d reads...", indexed_reads);
        }
    }

    PrintLog("K-mer index built: %zu unique kmers", kmer_index_.size());
    return 0;
}

// ==================================================
// Generate Candidate Pairs
// ==================================================

int GraphAssembler::GenerateCandidatePairs() {
    PrintLog("Generating candidate pairs using k-mer index...");

    std::set<std::pair<int, int>> unique_pairs;

    for (auto &r : reads_) {
        std::set<int> candidates;

        if ((int)r.seq.length() < config_.kmer_size) continue;

        for (int i = 0; i <= (int)r.seq.length() - config_.kmer_size; ++i) {
            std::string kmer = r.seq.substr(i, config_.kmer_size);

            auto &hits = kmer_index_[kmer];
            for (const auto &hit : hits) {
                int other_id = hit.first;
                if (other_id != r.id) {
                    candidates.insert(other_id);
                }
            }

            if ((int)candidates.size() > config_.max_candidates_per_read) break;
        }

        for (int cand_id : candidates) {
            int a = r.id;
            int b = cand_id;
            if (a > b) std::swap(a, b);
            unique_pairs.insert({a, b});

            if (unique_pairs.size() > 5000000) break;
        }
    }

    candidate_pairs_.clear();
    for (const auto &p : unique_pairs) {
        candidate_pairs_.push_back(p);
    }

    stats_.total_candidates_tested = candidate_pairs_.size();
    PrintLog("Generated %zu candidate pairs", candidate_pairs_.size());

    return 0;
}

// ==================================================
// Compute Overlaps
// ==================================================

int GraphAssembler::ComputeOverlap(int read1_id, int read2_id, OverlapEdge &edge) {
    if (read1_id >= (int)reads_.size() || read2_id >= (int)reads_.size()) {
        return -1;
    }

    const Read &r1 = reads_[read1_id];
    const Read &r2 = reads_[read2_id];

    int best_overlap_len = 0;
    int best_mismatches = 1000000;
    float best_identity = 0;
    int best_pos1_start = -1;

    int max_overlap = std::min(r1.seq.length(), r2.seq.length());
    for (int overlap_len = config_.min_overlap_length; overlap_len <= max_overlap; ++overlap_len) {
        int pos1_start = r1.seq.length() - overlap_len;
        int mismatches = 0;
        for (int i = 0; i < overlap_len; ++i) {
            if (r1.seq[pos1_start + i] != r2.seq[i]) {
                mismatches++;
            }
        }

        float identity_local = 1.0 - (float)mismatches / overlap_len;
        if (identity_local >= config_.min_identity_threshold && mismatches < best_mismatches) {
            best_overlap_len = overlap_len;
            best_mismatches = mismatches;
            best_identity = identity_local;
            best_pos1_start = pos1_start;
        }
    }

    if (best_overlap_len >= config_.min_overlap_length) {
        edge.read1_id = read1_id;
        edge.read2_id = read2_id;
        edge.overlap_len = best_overlap_len;
        edge.overlap_pos1_start = best_pos1_start;
        edge.overlap_pos1_end = r1.seq.length() - 1;
        edge.overlap_pos2_start = 0;
        edge.overlap_pos2_end = best_overlap_len - 1;
        edge.mismatches = best_mismatches;
        edge.identity = best_identity;
        edge.orientation = 1;
        return 0;
    }

    return -1;
}

int GraphAssembler::ComputeOverlaps() {
    PrintLog("Computing overlaps for %zu candidate pairs...", candidate_pairs_.size());

    int overlap_found = 0;

    for (size_t i = 0; i < candidate_pairs_.size(); ++i) {
        int id1 = candidate_pairs_[i].first;
        int id2 = candidate_pairs_[i].second;

        OverlapEdge edge;
        if (ComputeOverlap(id1, id2, edge) == 0) {
            edges_.push_back(edge);
            overlap_found++;
        }

        if ((i + 1) % 100000 == 0) {
            PrintLog("Processed %zu pairs, %d overlaps found...", i + 1, overlap_found);
        }
    }

    stats_.total_candidates_with_overlap = overlap_found;
    stats_.total_edges = edges_.size();

    PrintLog("Overlap computation complete: %d edges found", overlap_found);

    return 0;
}

// ==================================================
// Build Graph from Edges
// ==================================================

int GraphAssembler::BuildGraph() {
    PrintLog("Building graph from %d edges...", (int)edges_.size());

    adjacency_list_.clear();

    for (const auto &e : edges_) {
        adjacency_list_[e.read1_id].push_back(e.read2_id);
        adjacency_list_[e.read2_id].push_back(e.read1_id);
    }

    stats_.degree_histogram.resize(100, 0);
    for (const auto &pair : adjacency_list_) {
        int degree = pair.second.size();
        if (degree < 100) {
            stats_.degree_histogram[degree]++;
        } else {
            stats_.degree_histogram[99]++;
        }
    }

    int isolated_nodes = 0;
    for (const auto &r : reads_) {
        if (adjacency_list_.find(r.id) == adjacency_list_.end()) {
            isolated_nodes++;
        }
    }

    stats_.num_isolated_nodes = isolated_nodes;

    PrintLog("Graph built: %zu nodes with edges, %d isolated nodes",
            adjacency_list_.size(), isolated_nodes);

    return 0;
}

// ==================================================
// Clean Graph
// ==================================================

int GraphAssembler::CleanGraph() {
    PrintLog("Cleaning graph (minimal)...");

    int removed_edges = 0;
    std::vector<OverlapEdge> cleaned_edges;

    for (const auto &e : edges_) {
        if (e.identity >= 0.85) {
            cleaned_edges.push_back(e);
        } else {
            removed_edges++;
        }
    }

    edges_ = std::move(cleaned_edges);
    PrintLog("Removed %d low-identity edges", removed_edges);

    BuildGraph();

    return 0;
}

// ==================================================
// Build Connected Components
// ==================================================

int GraphAssembler::BuildConnectedComponents() {
    std::set<int> visited;
    int components = 0;
    int largest = 0;

    for (const auto &r : reads_) {
        if (visited.count(r.id)) continue;

        std::vector<int> component;
        std::queue<int> q;
        q.push(r.id);
        visited.insert(r.id);

        while (!q.empty()) {
            int node = q.front();
            q.pop();
            component.push_back(node);

            if (adjacency_list_.count(node)) {
                for (int neighbor : adjacency_list_[node]) {
                    if (!visited.count(neighbor)) {
                        visited.insert(neighbor);
                        q.push(neighbor);
                    }
                }
            }
        }

        components++;
        if ((int)component.size() > largest) {
            largest = component.size();
        }
    }

    stats_.num_connected_components = components;
    stats_.largest_component_size = largest;

    PrintLog("Found %d connected components, largest: %d nodes",
            components, largest);

    return 0;
}

// ==================================================
// Extract Contigs
// ==================================================

int GraphAssembler::ExtractContigs() {
    PrintLog("Extracting contigs from graph...");

    BuildConnectedComponents();

    std::set<int> used_reads;
    int contig_id = 0;

    for (const auto &r : reads_) {
        if (used_reads.count(r.id)) continue;

        Contig contig;
        contig.contig_id = contig_id++;
        contig.read_ids.push_back(r.id);
        used_reads.insert(r.id);

        int current_node = r.id;
        bool extended = true;

        while (extended) {
            extended = false;

            if (adjacency_list_.count(current_node)) {
                int best_neighbor = -1;
                int best_freq = 0;

                for (int neighbor : adjacency_list_[current_node]) {
                    if (!used_reads.count(neighbor) && reads_[neighbor].frequency > best_freq) {
                        best_neighbor = neighbor;
                        best_freq = reads_[neighbor].frequency;
                    }
                }

                if (best_neighbor != -1) {
                    contig.read_ids.push_back(best_neighbor);
                    used_reads.insert(best_neighbor);
                    current_node = best_neighbor;
                    extended = true;
                }
            }
        }

        contigs_.push_back(contig);
    }

    stats_.num_assembled_contigs = contigs_.size();
    PrintLog("Extracted %d contigs", (int)contigs_.size());

    return 0;
}

// ==================================================
// Generate Consensus
// ==================================================

int GraphAssembler::ComputeConsensus(Contig &contig) {
    if (contig.read_ids.empty()) return -1;

    int max_len = 0;
    for (int read_id : contig.read_ids) {
        if ((int)reads_[read_id].seq.length() > max_len) {
            max_len = reads_[read_id].seq.length();
        }
    }

    std::vector<std::vector<int>> base_counts(max_len, std::vector<int>(4, 0));

    for (int read_id : contig.read_ids) {
        const std::string &seq = reads_[read_id].seq;
        for (int j = 0; j < (int)seq.length() && j < max_len; ++j) {
            char c = seq[j];
            if (c == 'A') base_counts[j][0]++;
            else if (c == 'C') base_counts[j][1]++;
            else if (c == 'G') base_counts[j][2]++;
            else if (c == 'T') base_counts[j][3]++;
        }
    }

    char consensus_bases[] = {'A', 'C', 'G', 'T'};
    contig.consensus.clear();

    for (int i = 0; i < max_len; ++i) {
        int best_base = 0;
        int best_count = base_counts[i][0];

        for (int j = 1; j < 4; ++j) {
            if (base_counts[i][j] > best_count) {
                best_base = j;
                best_count = base_counts[i][j];
            }
        }

        contig.consensus += consensus_bases[best_base];
    }

    return 0;
}

int GraphAssembler::GenerateConsensus() {
    PrintLog("Generating consensus sequences for %d contigs...", (int)contigs_.size());

    for (auto &contig : contigs_) {
        ComputeConsensus(contig);
    }

    std::vector<int> lengths;
    double total_len = 0;
    for (const auto &c : contigs_) {
        lengths.push_back(c.consensus.length());
        total_len += c.consensus.length();
        stats_.avg_contig_length += c.consensus.length();
    }

    if (!contigs_.empty()) {
        stats_.avg_contig_length /= contigs_.size();
    }

    std::sort(lengths.begin(), lengths.end(), std::greater<int>());
    double cumsum = 0;
    for (int len : lengths) {
        cumsum += len;
        if (cumsum >= total_len / 2) {
            stats_.n50_contig_length = len;
            break;
        }
    }

    PrintLog("Consensus generated. Avg length: %.0f, N50: %.0f",
            stats_.avg_contig_length, stats_.n50_contig_length);

    return 0;
}

// ==================================================
// Write Output
// ==================================================

int GraphAssembler::WriteContigsFasta(const char *output_path) {
    PrintLog("Writing contigs to %s...", output_path);

    FILE *f = fopen(output_path, "w");
    if (!f) {
        PrintLog("ERROR: Cannot open output file");
        return -1;
    }

    for (const auto &c : contigs_) {
        fprintf(f, ">contig_%d\n%s\n", c.contig_id, c.consensus.c_str());
    }

    fclose(f);
    PrintLog("Wrote %d contigs", (int)contigs_.size());

    return 0;
}

int GraphAssembler::WriteGraphStats(const char *stats_path) {
    PrintLog("Writing graph statistics to %s...", stats_path);

    FILE *f = fopen(stats_path, "w");
    if (!f) {
        PrintLog("ERROR: Cannot open stats file");
        return -1;
    }

    fprintf(f, "GRAPH_STATISTICS\n");
    fprintf(f, "total_reads\t%d\n", stats_.total_reads);
    fprintf(f, "total_nodes\t%d\n", stats_.total_nodes);
    fprintf(f, "total_edges\t%d\n", stats_.total_edges);
    fprintf(f, "candidates_tested\t%d\n", stats_.total_candidates_tested);
    fprintf(f, "candidates_with_overlap\t%d\n", stats_.total_candidates_with_overlap);
    fprintf(f, "num_connected_components\t%d\n", stats_.num_connected_components);
    fprintf(f, "largest_component_size\t%d\n", stats_.largest_component_size);
    fprintf(f, "num_isolated_nodes\t%d\n", stats_.num_isolated_nodes);
    fprintf(f, "num_assembled_contigs\t%d\n", stats_.num_assembled_contigs);
    fprintf(f, "avg_contig_length\t%.1f\n", stats_.avg_contig_length);
    fprintf(f, "n50_contig_length\t%.0f\n", stats_.n50_contig_length);

    fprintf(f, "\nDEGREE_DISTRIBUTION\n");
    for (int i = 0; i < (int)stats_.degree_histogram.size(); ++i) {
        if (stats_.degree_histogram[i] > 0) {
            fprintf(f, "degree_%d\t%d\n", i, stats_.degree_histogram[i]);
        }
    }

    fclose(f);
    PrintLog("Statistics written");

    return 0;
}
