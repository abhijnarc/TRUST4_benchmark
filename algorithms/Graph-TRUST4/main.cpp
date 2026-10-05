#include "GraphAssembler.hpp"
#include <sys/time.h>
#include <sys/resource.h>

void PrintUsage(const char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS]\n", prog);
    fprintf(stderr, "  -1 FILE    : R1 FASTQ file (required)\n");
    fprintf(stderr, "  -2 FILE    : R2 FASTQ file (required)\n");
    fprintf(stderr, "  -o PREFIX  : output prefix (default: graph_out)\n");
    fprintf(stderr, "  -k INT     : k-mer size (default: 9)\n");
    fprintf(stderr, "  -m INT     : min overlap length (default: 20)\n");
    fprintf(stderr, "  -i FLOAT   : min identity threshold (default: 0.90)\n");
}

int main(int argc, char *argv[]) {
    fprintf(stderr, "Graph-TRUST4 Assembly (Pilot v1.0)\n");
    fprintf(stderr, "==================================\n\n");

    const char *fastq1 = NULL;
    const char *fastq2 = NULL;
    const char *output_prefix = "graph_out";
    GraphConfig config;

    // Parse command line
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-1") == 0 && i + 1 < argc) {
            fastq1 = argv[++i];
        } else if (strcmp(argv[i], "-2") == 0 && i + 1 < argc) {
            fastq2 = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_prefix = argv[++i];
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            config.kmer_size = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            config.min_overlap_length = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            config.min_identity_threshold = atof(argv[++i]);
        } else {
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (!fastq1 || !fastq2) {
        PrintUsage(argv[0]);
        return 1;
    }

    fprintf(stderr, "Input R1: %s\n", fastq1);
    fprintf(stderr, "Input R2: %s\n", fastq2);
    fprintf(stderr, "Output prefix: %s\n", output_prefix);
    fprintf(stderr, "Config: k=%d, min_overlap=%d, min_identity=%.2f\n\n",
            config.kmer_size, config.min_overlap_length, config.min_identity_threshold);

    time_t start_time = time(NULL);

    // Create assembler
    GraphAssembler assembler(config);

    // Run pipeline
    if (assembler.LoadReads(fastq1, fastq2) < 0) return 1;
    if (assembler.DeduplicateReads() < 0) return 1;
    if (assembler.BuildKmerIndex() < 0) return 1;
    if (assembler.GenerateCandidatePairs() < 0) return 1;
    if (assembler.ComputeOverlaps() < 0) return 1;
    if (assembler.BuildGraph() < 0) return 1;
    if (assembler.CleanGraph() < 0) return 1;
    if (assembler.ExtractContigs() < 0) return 1;
    if (assembler.GenerateConsensus() < 0) return 1;

    // Write output
    char output_file[1000];
    snprintf(output_file, sizeof(output_file), "%s_contigs.fa", output_prefix);
    if (assembler.WriteContigsFasta(output_file) < 0) return 1;

    snprintf(output_file, sizeof(output_file), "%s_stats.tsv", output_prefix);
    if (assembler.WriteGraphStats(output_file) < 0) return 1;

    time_t end_time = time(NULL);
    double elapsed = difftime(end_time, start_time);

    fprintf(stderr, "\n==================================\n");
    fprintf(stderr, "Assembly Complete\n");
    fprintf(stderr, "Total time: %.0f seconds (%.1f minutes)\n", elapsed, elapsed / 60);

    const auto &stats = assembler.GetStats();
    fprintf(stderr, "\nFinal Statistics:\n");
    fprintf(stderr, "  Reads: %d\n", stats.total_reads);
    fprintf(stderr, "  Nodes: %d\n", stats.total_nodes);
    fprintf(stderr, "  Edges: %d\n", stats.total_edges);
    fprintf(stderr, "  Contigs: %d\n", stats.num_assembled_contigs);
    fprintf(stderr, "  Avg contig length: %.0f\n", stats.avg_contig_length);
    fprintf(stderr, "  N50 contig length: %.0f\n", stats.n50_contig_length);
    fprintf(stderr, "  Connected components: %d\n", stats.num_connected_components);
    fprintf(stderr, "  Largest component: %d nodes\n", stats.largest_component_size);

    fprintf(stderr, "\nOutput files:\n");
    snprintf(output_file, sizeof(output_file), "%s_contigs.fa", output_prefix);
    fprintf(stderr, "  %s\n", output_file);
    snprintf(output_file, sizeof(output_file), "%s_stats.tsv", output_prefix);
    fprintf(stderr, "  %s\n", output_file);

    return 0;
}
