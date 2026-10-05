#ifndef GRAPH_ASSEMBLER_V4_HPP
#define GRAPH_ASSEMBLER_V4_HPP

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct V4Bundle {
    uint32_t id;
    std::string sequence;
    std::string quality;
    uint32_t abundance;
};

struct V4Edge {
    uint32_t source;
    uint32_t target;
    uint16_t overlap;
    uint16_t mismatches;
    float identity;
    float qaos;
    float expected_mismatches;
    uint16_t high_quality_mismatches;
    uint8_t orientation;
    uint16_t paired_support;
    float mean_q;
    float min_q;
};

struct V4Config {
    int kmer = 9;
    int min_overlap = 31;
    float min_identity = 0.90f;
    float min_qaos = 0.90f;
    int qhigh = 20;
    uint32_t max_candidates = 500;
    uint32_t max_postings = 1000;
    uint32_t batch_size = 10000;
    int max_memory_gb = 32;
    int threads = 64;
    float min_branch_qaos_delta = 0.01f;
    uint32_t max_alternative_paths_per_branch = 2;
};

struct V4Stats {
    uint64_t raw_reads = 0;
    uint64_t bundles = 0;
    uint64_t candidate_pairs = 0;
    uint64_t accepted_edges = 0;
    uint64_t components = 0;
    uint64_t isolated = 0;
    uint64_t contigs = 0;
    uint64_t n50 = 0;
    uint64_t peak_rss_mb = 0;
    double mean_identity = 0;
    double mean_qaos = 0;
    double mean_overlap = 0;
    uint64_t branches = 0;
    uint64_t resolved_branches = 0;
    uint64_t ambiguous_branches = 0;
    uint64_t cycles = 0;
};

struct V4AssemblyDebug {
    uint64_t contig_id;
    uint64_t contig_length;
    uint64_t num_nodes;
    uint64_t unique_nodes;
    uint64_t repeated_nodes;
    uint64_t total_overlap;
    uint64_t expected_length;
    uint64_t max_node_length;
    bool cycle_detected;
    bool valid_path;
    float consensus_confidence;
};

struct V4BranchDecision {
    uint32_t branch_node;
    uint32_t best_edge;
    uint32_t second_edge;
    float best_qaos;
    float second_qaos;
    float qaos_delta;
    uint16_t overlap_best;
    uint16_t overlap_second;
    uint16_t high_q_mismatch_best;
    uint16_t high_q_mismatch_second;
    uint16_t paired_support_best;
    uint16_t paired_support_second;
    std::string decision;
    std::string reason;
};

class GraphAssemblerV4 {
public:
    explicit GraphAssemblerV4(const V4Config &config);
    int run(const std::string &r1, const std::string &r2, const std::string &output_dir);

private:
    V4Config config_;
    V4Stats stats_;
    std::vector<V4Bundle> bundles_;
    std::vector<V4Edge> edges_;
    std::unordered_map<std::string, std::vector<uint32_t>> index_;
    std::vector<std::vector<uint32_t>> adjacency_;
    std::vector<std::string> contigs_;
    std::vector<V4AssemblyDebug> assembly_debug_;
    std::vector<V4BranchDecision> branch_decisions_;
    std::vector<float> consensus_confidences_;
    std::string output_dir_;
    double start_seconds_;

    bool load_and_bundle(const std::string &r1, const std::string &r2);
    bool build_index();
    bool generate_and_score();
    bool score_overlap(uint32_t a, uint32_t b, V4Edge &edge) const;
    void build_contigs();
    void build_contigs_v4();
    bool memory_ok(const char *stage);
    uint64_t rss_mb() const;
    void progress(const char *stage, uint64_t done, uint64_t total, uint64_t candidates, uint64_t accepted) const;
    bool write_outputs() const;
};

#endif
