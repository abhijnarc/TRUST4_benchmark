#include "GraphAssemblerV3.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(const char *program) {
    std::fprintf(stderr, "Usage: %s -1 R1 -2 R2 -o OUTPUT_DIR [options]\n", program);
    std::fprintf(stderr, "  -k INT  k-mer size (default 9)\n  -m INT  minimum overlap (default 31)\n");
    std::fprintf(stderr, "  -i FLOAT minimum identity (default 0.90)\n  -q FLOAT minimum QAOS (default 0.90)\n");
    std::fprintf(stderr, "  -H INT  high-quality threshold (default 20)\n  -c INT  max candidates (default 500)\n");
    std::fprintf(stderr, "  -M INT  RSS limit GB (default 32)\n  -t INT  OpenMP threads (default 64)\n");
    std::fprintf(stderr, "  -b INT  batch size (default 10000)\n  -p INT  max postings per k-mer (default 1000)\n");
}

int main(int argc, char **argv) {
    std::string r1, r2, output;
    V4Config config;
    for (int i = 1; i < argc; ++i) {
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
        else if (!std::strcmp(argv[i], "-M")) config.max_memory_gb = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-t")) config.threads = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-b")) config.batch_size = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-p")) config.max_postings = std::atoi(argv[++i]);
        else { usage(argv[0]); return 2; }
    }
    if (r1.empty() || r2.empty() || output.empty() || config.threads < 1) { usage(argv[0]); return 2; }
#ifdef _OPENMP
    omp_set_num_threads(config.threads);
    std::fprintf(stderr, "Graph-TRUST4-v4 OpenMP enabled; requested threads=%d\n", config.threads);
#else
    std::fprintf(stderr, "ERROR: this binary was not compiled with OpenMP\n");
    return 1;
#endif
    std::fprintf(stderr, "[Config] k=%d min_overlap=%d identity=%.3f QAOS=%.3f qhigh=%d candidates=%u postings=%u batch=%u memory=%dGB threads=%d\n",
                 config.kmer, config.min_overlap, config.min_identity, config.min_qaos, config.qhigh,
                 config.max_candidates, config.max_postings, config.batch_size, config.max_memory_gb, config.threads);
    return GraphAssemblerV4(config).run(r1, r2, output);
}
