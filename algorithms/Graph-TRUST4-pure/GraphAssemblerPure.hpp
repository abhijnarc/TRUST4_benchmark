#ifndef GRAPH_ASSEMBLER_PURE_HPP
#define GRAPH_ASSEMBLER_PURE_HPP

#include <cstdint>
#include <string>
#include <vector>

struct PureBundle {
    uint32_t id;
    std::string sequence;
    std::string quality;
    std::string reverse_sequence;
    std::string reverse_quality;
    uint32_t abundance;
    std::vector<std::string> read_ids;
};

struct PureEdge {
    uint32_t id;
    uint32_t source;
    uint32_t target;
    uint8_t orientation;
    uint16_t overlap;
    uint16_t mismatches;
    float identity;
    float qaos;
    float mean_q;
    float min_q;
    uint16_t high_q_mismatches;
    uint32_t candidate_support;
    uint16_t paired_support;
};

struct PureConfig {
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
    uint32_t max_paths_per_component = 1000;
    uint32_t max_path_nodes = 100;
};

class PureGraphAssembler {
public:
    struct Path {
        std::vector<uint32_t> nodes;
        std::vector<uint32_t> edges;
        double log_qaos = 0.0;
        double identity_sum = 0.0;
        uint64_t overlap_sum = 0;
        uint64_t high_q_mismatch_sum = 0;
    };
    explicit PureGraphAssembler(const PureConfig &config);
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

    PureConfig config_;
    std::vector<PureBundle> bundles_;
    std::vector<PureEdge> edges_;
    std::vector<std::vector<uint32_t>> outgoing_;
    std::vector<std::vector<uint32_t>> incoming_;
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
    bool build_edges();
    bool enumerate_paths();
    bool write_outputs(const std::string &output) const;
    bool score_edge(uint32_t source, uint32_t target, uint32_t support, PureEdge &edge) const;
    void enumerate_from(const Component &component, uint32_t node, Path path,
                        std::vector<uint8_t> &visited, uint32_t &component_paths);
    std::vector<Component> components() const;
    uint64_t rss_mb() const;
    void progress(const char *stage, uint64_t done, uint64_t total) const;
};

#endif