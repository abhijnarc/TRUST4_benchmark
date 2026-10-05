#include "GraphAssemblerPure.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(const char *program) {
    std::cerr << "Usage: " << program << " -1 R1 -2 R2 -o OUTPUT_DIR [options]\n"
              << "  -k INT k-mer (default 9)  -m INT minimum overlap (default 31)\n"
              << "  -i FLOAT minimum identity (default 0.90)  -q FLOAT minimum QAOS (default 0.90)\n"
              << "  -H INT Q_HIGH (default 20)  -c INT candidate cap (default 500)\n"
              << "  -p INT posting cap (default 1000)  -b INT batch size (default 10000)\n"
              << "  -M INT memory GB (default 32)  -t INT threads (default 64)\n"
              << "  --max-paths INT (default 1000)  --max-path-nodes INT (default 100)\n"
              << "  --self-test run QAOS/orientation/branch sanity tests\n";
}

int main(int argc, char **argv) {
    PureConfig config;
    std::string r1, r2, output;
    bool self_test = false;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--self-test")) { self_test = true; continue; }
        if (i + 1 >= argc) { usage(argv[0]); return 2; }
        if (!std::strcmp(argv[i], "-1")) r1 = argv[++i];
        else if (!std::strcmp(argv[i], "-2")) r2 = argv[++i];
        else if (!std::strcmp(argv[i], "-o")) output = argv[++i];
        else if (!std::strcmp(argv[i], "-k")) config.kmer = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-m")) config.min_overlap = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-i")) config.min_identity = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-q")) config.min_qaos = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "-H")) config.qhigh = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-c")) config.max_candidates = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-p")) config.max_postings = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-b")) config.batch_size = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-M")) config.max_memory_gb = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-t")) config.threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-paths")) config.max_paths_per_component = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--max-path-nodes")) config.max_path_nodes = std::atoi(argv[++i]);
        else { usage(argv[0]); return 2; }
    }
    if (self_test) return PureGraphAssembler::self_test();
    if (r1.empty() || r2.empty() || output.empty() || config.threads < 1) { usage(argv[0]); return 2; }
#ifdef _OPENMP
    omp_set_num_threads(config.threads);
    std::cerr << "requested_threads=" << config.threads
              << " omp_max_threads=" << omp_get_max_threads() << " cpu_count=" << omp_get_num_procs() << "\n";
#else
    std::cerr << "ERROR: OpenMP is required\n";
    return 1;
#endif
    return PureGraphAssembler(config).run(r1, r2, output);
}