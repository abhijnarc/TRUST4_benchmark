#include "clean_graph.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
void usage(const char *name) {
    std::cerr << "Usage: " << name << " -1 candidate_R1.fq -2 candidate_R2.fq "
              << "-o OUTPUT -s candidate_statistics.tsv [options]\n"
              << "  -k INT k-mer seed (default 9)\n"
              << "  -m INT minimum overlap (default 31)\n"
              << "  --max-overlap INT (default 150)\n"
              << "  -i FLOAT minimum identity (default 0.90)\n"
              << "  -q FLOAT minimum QAOS (default 0.90)\n"
              << "  --max-postings INT (default 1000)\n"
              << "  --max-candidates INT (default 50 per read)\n"
              << "  --max-path-nodes INT (default 64)\n"
              << "  --max-paths INT (default 100 per component)\n"
              << "  -t INT threads (default 64)\n";
}
}

int main(int argc, char **argv) {
    CleanConfig config;
    std::string r1, r2, output, candidate_stats;
    for (int i = 1; i < argc; ++i) {
        if (i + 1 >= argc) { usage(argv[0]); return 2; }
        if (!std::strcmp(argv[i], "-1")) r1 = argv[++i];
        else if (!std::strcmp(argv[i], "-2")) r2 = argv[++i];
        else if (!std::strcmp(argv[i], "-o")) output = argv[++i];
        else if (!std::strcmp(argv[i], "-s")) candidate_stats = argv[++i];
        else if (!std::strcmp(argv[i], "-k")) config.kmer = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-m")) config.min_overlap = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-overlap")) config.max_overlap = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-i")) config.min_identity = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-q")) config.min_qaos = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-postings")) config.max_postings = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-candidates")) config.max_candidates_per_read = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-path-nodes")) config.max_nodes_per_path = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-paths")) config.max_paths_per_component = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-t")) config.threads = std::atoi(argv[++i]);
        else { usage(argv[0]); return 2; }
    }
    if (r1.empty() || r2.empty() || output.empty() || candidate_stats.empty() ||
        config.threads == 0 || config.kmer < 1 || config.min_overlap < config.kmer) {
        usage(argv[0]);
        return 2;
    }
#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(config.threads));
    std::cerr << "requested_threads=" << config.threads
              << " omp_max_threads=" << omp_get_max_threads() << "\n";
#else
    std::cerr << "OpenMP support is required\n";
    return 1;
#endif
    return CleanGraphAssembler(config).run(r1, r2, output, candidate_stats);
}
