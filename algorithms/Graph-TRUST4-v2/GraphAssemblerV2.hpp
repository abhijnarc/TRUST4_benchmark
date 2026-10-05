#ifndef _GRAPH_ASSEMBLER_V2_H
#define _GRAPH_ASSEMBLER_V2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <memory>
#include <set>

// ==================================================
// Data Structures for v2 Assembly
// ==================================================

struct ReadBundle {
    int bundle_id;
    std::string consensus_seq;
    int abundance;
    std::vector<int> read_ids;
    std::string qualities;
    std::vector<std::vector<int>> per_position_counts;
    std::vector<float> per_position_quals;
    bool is_rc;

    ReadBundle() : bundle_id(-1), abundance(0), is_rc(false) {}
};

struct QualityAwareOverlap {
    int bundle1_id;
    int bundle2_id;
    int overlap_length;
    float observed_identity;
    int observed_mismatches;
    float expected_mismatches;
    float qaos;
    float mean_q;
    float min_q;
    int high_quality_mismatches;
    int orientation;
    int bundle1_end;
    int bundle2_start;
    int abundance1;
    int abundance2;

    QualityAwareOverlap() : bundle1_id(-1), bundle2_id(-1), overlap_length(0),
                            observed_identity(0), observed_mismatches(0),
                            expected_mismatches(0), qaos(0), mean_q(0), min_q(0),
                            high_quality_mismatches(0), orientation(1),
                            bundle1_end(0), bundle2_start(0),
                            abundance1(0), abundance2(0) {}
};

struct AssembledContig {
    int contig_id;
    std::string consensus;
    std::vector<int> bundle_ids;
    int contig_length;
    std::vector<float> position_qualities;

    AssembledContig() : contig_id(-1), contig_length(0) {}
};

// ==================================================
// Configuration
// ==================================================

struct GraphConfigV2 {
    int kmer_size;
    int min_overlap_length;
    float min_identity_threshold;
    float min_qaos_threshold;
    int qhigh_threshold;
    float bundle_identity_threshold;
    int max_candidates_per_bundle;
    bool allow_rc_overlaps;

    GraphConfigV2() : kmer_size(9),
                     min_overlap_length(31),
                     min_identity_threshold(0.90),
                     min_qaos_threshold(0.90),
                     qhigh_threshold(20),
                     bundle_identity_threshold(0.95),
                     max_candidates_per_bundle(1000),
                     allow_rc_overlaps(true) {}
};

// ==================================================
// Graph Statistics
// ==================================================

struct GraphStatsV2 {
    int total_raw_reads;
    int total_bundles;
    float compression_ratio;
    int total_nodes;
    int total_edges;
    int num_connected_components;
    int largest_component_size;
    int num_isolated_nodes;
    int num_assembled_contigs;
    double avg_contig_length;
    double n50_contig_length;
    int num_contigs_150bp;
    int num_contigs_300bp;
    float mean_overlap_length;
    float mean_identity;
    float mean_qaos;

    GraphStatsV2() : total_raw_reads(0), total_bundles(0), compression_ratio(0),
                    total_nodes(0), total_edges(0),
                    num_connected_components(0), largest_component_size(0),
                    num_isolated_nodes(0), num_assembled_contigs(0),
                    avg_contig_length(0), n50_contig_length(0),
                    num_contigs_150bp(0), num_contigs_300bp(0),
                    mean_overlap_length(0), mean_identity(0), mean_qaos(0) {}
};

// ==================================================
// Main Assembler Class
// ==================================================

class GraphAssemblerV2 {
public:
    GraphAssemblerV2(const GraphConfigV2 &cfg);
    ~GraphAssemblerV2();

    int LoadReadsWithQualities(const char *fastq1_path, const char *fastq2_path);
    int BundleRedundantReads();
    int BuildKmerIndex();
    int GenerateCandidateBundlePairs();
    int ComputeQualityAwareOverlaps();
    int BuildGraph();
    int CleanGraph();
    int ExtractContigs();
    int GenerateConsensusSequences();
    int WriteBundles(const char *output_path);
    int WriteGraphEdges(const char *output_path);
    int WriteGraphStats(const char *stats_path);
    int WriteContigsFasta(const char *output_path);
    int WriteAssemblyStats(const char *stats_path);
    int WriteQAOSStats(const char *stats_path);

    const GraphStatsV2& GetStats() const { return stats_; }

private:
    GraphConfigV2 config_;
    std::vector<std::string> raw_sequences_;
    std::vector<std::string> raw_qualities_;
    std::vector<ReadBundle> bundles_;
    std::vector<QualityAwareOverlap> edges_;
    std::unordered_map<int, std::vector<int>> adjacency_list_;
    std::vector<AssembledContig> contigs_;
    GraphStatsV2 stats_;
    std::unordered_map<std::string, std::vector<int>> kmer_index_;
    std::vector<std::pair<int, int>> candidate_bundle_pairs_;

    float PhredToErrorProb(int q);
    float ComputeBaseCompatibility(char base1, int q1, char base2, int q2);
    int ComputeVariableLengthOverlap(int bundle1_id, int bundle2_id,
                                      QualityAwareOverlap &overlap);
    int ComputeQAOS(const std::string &seq1, const std::string &qual1,
                    const std::string &seq2, const std::string &qual2,
                    int start1, int start2, int overlap_len,
                    float &qaos, int &mismatches, float &expected_mismatch,
                    int &high_q_mismatches, float &mean_q, float &min_q);
    int BuildConnectedComponents();
    int TraverseComponentForContigs(const std::vector<int> &component_nodes);
    int GenerateConsensus(AssembledContig &contig);
};

#endif
