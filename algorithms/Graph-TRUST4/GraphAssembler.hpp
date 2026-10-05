#ifndef _GRAPH_ASSEMBLER_H
#define _GRAPH_ASSEMBLER_H

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

// ==================================================
// Data Structures for Graph Assembly
// ==================================================

struct Read {
    int id;
    std::string seq;
    int len;
    int frequency; // multiplicity if deduplicated

    Read() : id(-1), len(0), frequency(1) {}
    Read(int id_, const std::string &seq_) : id(id_), seq(seq_), len(seq_.length()), frequency(1) {}
};

struct OverlapEdge {
    int read1_id;
    int read2_id;
    int overlap_len;
    int overlap_pos1_start, overlap_pos1_end;
    int overlap_pos2_start, overlap_pos2_end;
    int mismatches;
    float identity;
    int orientation; // 1: same strand, -1: RC

    OverlapEdge() : read1_id(-1), read2_id(-1), overlap_len(0),
                   overlap_pos1_start(0), overlap_pos1_end(0),
                   overlap_pos2_start(0), overlap_pos2_end(0),
                   mismatches(0), identity(0), orientation(1) {}
};

struct Contig {
    int contig_id;
    std::string consensus;
    std::vector<int> read_ids;

    Contig() : contig_id(-1) {}
};

// ==================================================
// Configuration
// ==================================================

struct GraphConfig {
    int kmer_size;
    int min_kmer_freq_for_candidate;
    int min_overlap_length;
    float min_identity_threshold;
    int max_candidates_per_read;
    bool allow_rc_overlaps;

    GraphConfig() : kmer_size(9),
                   min_kmer_freq_for_candidate(2),
                   min_overlap_length(20),
                   min_identity_threshold(0.90),
                   max_candidates_per_read(1000),
                   allow_rc_overlaps(true) {}
};

// ==================================================
// Graph Statistics
// ==================================================

struct GraphStats {
    int total_reads;
    int total_nodes;
    int total_edges;
    int total_candidates_tested;
    int total_candidates_with_overlap;
    std::vector<int> degree_histogram;
    int num_connected_components;
    int largest_component_size;
    int num_isolated_nodes;
    int num_assembled_contigs;
    double avg_contig_length;
    double n50_contig_length;

    GraphStats() : total_reads(0), total_nodes(0), total_edges(0),
                  total_candidates_tested(0), total_candidates_with_overlap(0),
                  num_connected_components(0), largest_component_size(0),
                  num_isolated_nodes(0), num_assembled_contigs(0),
                  avg_contig_length(0), n50_contig_length(0) {}
};

// ==================================================
// Graph Assembler Main Class
// ==================================================

class GraphAssembler {
public:
    GraphAssembler(const GraphConfig &cfg);
    ~GraphAssembler();

    int LoadReads(const char *fastq1_path, const char *fastq2_path);
    int DeduplicateReads();
    int BuildKmerIndex();
    int GenerateCandidatePairs();
    int ComputeOverlaps();
    int BuildGraph();
    int CleanGraph();
    int ExtractContigs();
    int GenerateConsensus();
    int WriteContigsFasta(const char *output_path);
    int WriteGraphStats(const char *stats_path);

    const GraphStats& GetStats() const { return stats_; }

private:
    GraphConfig config_;
    std::vector<Read> reads_;
    std::vector<OverlapEdge> edges_;
    std::unordered_map<int, std::vector<int>> adjacency_list_;
    std::vector<Contig> contigs_;
    GraphStats stats_;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> kmer_index_;
    std::vector<std::pair<int, int>> candidate_pairs_;

    int ComputeOverlap(int read1_id, int read2_id, OverlapEdge &edge);
    int BuildConnectedComponents();
    int TraverseComponentForContigs(const std::vector<int> &component_nodes);
    int ComputeConsensus(Contig &contig);
    void PrintLog(const char *fmt, ...);
};

#endif
