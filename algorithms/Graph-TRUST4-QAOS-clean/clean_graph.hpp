#ifndef GRAPH_TRUST4_QAOS_CLEAN_HPP
#define GRAPH_TRUST4_QAOS_CLEAN_HPP

#include <cstdint>
#include <array>
#include <string>
#include <vector>

struct CleanConfig {
    int kmer = 9;
    int min_overlap = 31;
    int max_overlap = 150;
    double min_identity = 0.90;
    double min_qaos = 0.90;
    unsigned max_postings = 1000;
    unsigned max_candidates_per_read = 50;
    unsigned max_nodes_per_path = 64;
    unsigned max_paths_per_component = 100;
    unsigned threads = 64;
};

struct CleanRead {
    uint32_t id = 0;
    std::string name;
    std::string mate_name;
    std::string sequence;
    std::string quality;
    uint32_t mate_id = UINT32_MAX;
};

struct CleanEdge {
    uint32_t id = 0;
    uint32_t source = 0;
    uint32_t target = 0;
    uint8_t source_orientation = 0;
    uint8_t target_orientation = 0;
    uint8_t orientation = 0;
    uint16_t overlap = 0;
    uint16_t matches = 0;
    uint16_t mismatches = 0;
    double identity = 0.0;
    double qaos = 0.0;
    double mean_quality = 0.0;
    double quality_support = 0.0;
};

struct CleanPath {
    std::vector<uint32_t> nodes;
    std::vector<uint32_t> edges;
    std::vector<uint8_t> orientations;
    uint32_t component = 0;
    bool cycle_terminated = false;
    bool length_limited = false;
    double cumulative_qaos = 0.0;
    double cumulative_identity = 0.0;
    std::string termination;
};

class CleanGraphAssembler {
public:
    explicit CleanGraphAssembler(CleanConfig config);
    int run(const std::string &r1, const std::string &r2,
            const std::string &output, const std::string &candidate_stats);

private:
    CleanConfig config_;
    std::vector<CleanRead> reads_;
    std::vector<CleanEdge> edges_;
    std::vector<std::vector<uint32_t>> outgoing_;
    std::vector<std::vector<uint32_t>> incoming_;
    std::vector<CleanPath> paths_;
    uint64_t candidate_edges_ = 0;
    uint64_t candidate_pairs_scored_ = 0;
    uint64_t raw_reads_ = 0;
    uint64_t rejected_cap_ = 0;
    uint64_t rejected_alignment_ = 0;
    uint64_t explored_paths_ = 0;
    uint64_t capped_components_ = 0;
    uint64_t cycles_ = 0;
    double start_seconds_ = 0;

    bool load_candidates(const std::string &r1, const std::string &r2,
                         const std::string &candidate_stats);
    bool construct_graph();
    void traverse_paths();
    void write_outputs(const std::string &output) const;
    CleanPath consensus_path(CleanPath path, std::string &sequence,
                             std::vector<std::array<double, 4>> &support) const;
};

#endif
