#ifndef GRAPH_ASSEMBLER_V3_HPP
#define GRAPH_ASSEMBLER_V3_HPP

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct V3Bundle {
    uint32_t id;
    std::string sequence;
    std::string quality;
    uint32_t abundance;
};

struct V3Edge {
    uint32_t source;
    uint32_t target;
    uint16_t overlap;
    uint16_t mismatches;
    float identity;
    float qaos;
    float expected_mismatches;
    uint16_t high_quality_mismatches;
    uint8_t orientation;
};

struct V3Config {
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
};

struct V3Stats {
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
};

class GraphAssemblerV3 {
public:
    explicit GraphAssemblerV3(const V3Config &config);
    int run(const std::string &r1, const std::string &r2, const std::string &output_dir);

private:
    V3Config config_;
    V3Stats stats_;
    std::vector<V3Bundle> bundles_;
    std::vector<V3Edge> edges_;
    std::unordered_map<std::string, std::vector<uint32_t>> index_;
    std::vector<std::vector<uint32_t>> adjacency_;
    std::vector<std::string> contigs_;
    std::string output_dir_;
    double start_seconds_;

    bool load_and_bundle(const std::string &r1, const std::string &r2);
    bool build_index();
    bool generate_and_score();
    bool score_overlap(uint32_t a, uint32_t b, V3Edge &edge) const;
    void build_contigs();
    bool memory_ok(const char *stage);
    uint64_t rss_mb() const;
    void progress(const char *stage, uint64_t done, uint64_t total, uint64_t candidates, uint64_t accepted) const;
    bool write_outputs() const;
};

#endif
