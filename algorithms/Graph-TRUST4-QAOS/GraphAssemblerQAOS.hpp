#ifndef GRAPH_ASSEMBLER_QAOS_HPP
#define GRAPH_ASSEMBLER_QAOS_HPP

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

struct QAOSBundle {
    uint32_t id;
    std::string sequence;
    std::string quality;
    std::string reverse_sequence;
    std::string reverse_quality;
    uint32_t abundance;
    uint8_t bcr_evidence = 0;
    std::vector<std::string> read_ids;
};

struct QAOSEdge {
    uint32_t id;
    uint32_t source;
    uint32_t target;
    uint8_t orientation;
    uint16_t overlap;
    uint16_t mismatches;
    float identity;
    float qaos;
    float edge_score;
    float mean_q;
    float min_q;
    uint16_t high_q_mismatches;
    uint32_t candidate_support;
    uint16_t paired_support;
};

struct QAOSConfig {
    int kmer = 9;
    int min_overlap = 31;
    float min_identity = 0.90f;
    float min_qaos = 0.90f;
    float min_edge_score = 0.75f;
    float weight_identity = 0.35f;
    float weight_qaos = 0.40f;
    float weight_overlap = 0.15f;
    float weight_quality = 0.10f;
    bool identity_only = false;
    int qhigh = 20;
    uint32_t max_candidates = 500;
    uint32_t max_postings = 1000;
    uint32_t batch_size = 10000;
    int max_memory_gb = 32;
    int threads = 64;
    uint32_t max_paths_per_component = 1000;
    uint32_t max_path_nodes = 100;
    std::string reference_path = "reference/TRUST4/human_IMGT+C.fa";
    bool v3_resolution = false;
};

class QAOSGraphAssembler {
public:
    struct Path {
        std::vector<uint32_t> nodes;
        std::vector<uint32_t> edges;
        double log_qaos = 0.0;
        double log_edge_score = 0.0;
        double abundance_score = 0.0;
        double identity_sum = 0.0;
        uint64_t overlap_sum = 0;
        uint64_t high_q_mismatch_sum = 0;
        uint64_t paired_support = 0;
            uint8_t bcr_state = 0;
            int biological_order_score = 0;
    };
    explicit QAOSGraphAssembler(const QAOSConfig &config);
    int run(const std::string &r1, const std::string &r2, const std::string &output);
    static int self_test();

private:
    struct Candidate {
        uint32_t target;
        uint32_t support;
    };
    struct Component {
        std::vector<uint32_t> nodes;
        std::vector<uint32_t> edges;
    };

    QAOSConfig config_;
    std::vector<QAOSBundle> bundles_;
    uint64_t paths_before_pruning_ = 0;
    uint64_t paths_removed_by_dominance_ = 0;
    uint64_t paths_removed_biologically_inconsistent_ = 0;
    uint64_t paths_retained_ambiguous_ = 0;
    std::vector<QAOSEdge> edges_;
    std::vector<std::vector<uint32_t>> outgoing_;
    std::vector<std::vector<uint32_t>> incoming_;
    std::unordered_set<uint64_t> paired_links_;
    std::vector<Path> paths_;
    uint64_t raw_reads_ = 0;
    uint64_t candidate_pairs_ = 0;
    uint64_t candidates_rejected_cap_ = 0;
    uint64_t candidates_evaluated_ = 0;
    uint64_t rejected_overlap_ = 0;
    uint64_t rejected_identity_ = 0;
    uint64_t rejected_qaos_ = 0;
    uint64_t branch_count_ = 0;
    uint64_t cycle_count_ = 0;
    double start_time_ = 0.0;
    double end_time_ = 0.0;
    uint64_t peak_rss_mb_ = 0;

    bool load_bundles(const std::string &r1, const std::string &r2);
    bool load_bcr_reference();
    bool build_edges();
    bool enumerate_paths();
    bool write_outputs(const std::string &output) const;
    void prune_paths();
    bool score_edge(uint32_t source, uint32_t target, uint32_t support, QAOSEdge &edge) const;
    void enumerate_from(const Component &component, uint32_t node, Path path,
                        std::vector<uint8_t> &visited, uint32_t &component_paths);
    std::vector<Component> components() const;
    uint64_t rss_mb() const;
    void progress(const char *stage, uint64_t done, uint64_t total) const;
};

#endif