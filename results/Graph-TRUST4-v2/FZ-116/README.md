# Graph-TRUST4 v2 FZ-116 Pilot Results

## Overview

This directory contains the results from the Graph-TRUST4 v2 quality-aware assembly pilot on FZ-116 BCR candidate reads.

## Key Innovation: Quality-Adaptive Overlap Scoring (QAOS)

Unlike v1 which uses fixed k-mer length and simple identity filtering, v2 incorporates:

- **Variable-length overlaps** (31+ bp)
- **Quality-aware base compatibility** using Phred score error probabilities
- **QAOS scoring** as the geometric mean of per-base compatibility probabilities
- **Multi-criterion filtering** combining identity, QAOS, and overlap length

## Files

### Input Data
- `TRUST_FZ-116_toassemble_1.fq` - R1 reads (335 MB, 617,400 pairs)
- `TRUST_FZ-116_toassemble_2.fq` - R2 reads (335 MB)

### Read Bundling Output
- `graph_fz116_bundles.tsv` - Collapsed read bundles with abundances

### Graph Statistics
- `graph_fz116_graph_stats.tsv` - Graph topology (nodes, edges, components)
- `graph_fz116_edges.tsv` - All edges with quality metrics

### Assembly Results
- `graph_fz116_contigs.fa` - Assembled contig sequences (FASTA)
- `graph_fz116_assembly_stats.tsv` - Contig statistics (N50, length distribution)

### Quality Assessment
- `graph_fz116_qaos_stats.tsv` - Mean overlap quality metrics

## Parameters

```
k=9                    # K-mer size for candidate generation
min_overlap=31         # Minimum overlap length (bp)
min_identity=0.90      # Minimum sequence identity (90%)
min_qaos=0.90          # Minimum QAOS threshold
Q_high=20              # High-quality mismatch threshold
```

## Input Statistics

| Metric | Value |
|--------|-------|
| Raw reads | 1,234,800 |
| Read pairs | 617,400 |
| Consensus bundles | 634,616 |
| Compression ratio | 1.95x |
| K-mers indexed | 4,052,877 |

## Quality-Adaptive Overlap Scoring

For each candidate bundle pair, QAOS is calculated as:

1. **Phred error model**: `e = 10^(-Q/10)`
2. **Base compatibility**: 
   - If match: `P = (1-e1)(1-e2)`
   - If mismatch: `P = e1/3 + e2/3 - e1*e2/9`
3. **Geometric mean**: `QAOS = (∏ compatibility)^(1/overlap_len)`

This approach incorporates quality information directly into overlap assessment, unlike v1's identity-only scoring.

## Comparison with v1

v1 uses:
- Fixed 20 bp minimum overlap
- Simple sequence identity filter
- No quality scoring

v2 uses:
- Variable 31+ bp overlaps
- Quality-aware compatibility calculation
- QAOS filtering alongside identity
- Expected mismatch calculation
- High-quality mismatch tracking

## Assembly Strategy

1. **Read Bundling**: Exact-match consensus preserves SHM variants
2. **Indexing**: K-mers for efficient candidate generation only
3. **Overlap**: Full suffix→prefix alignment with quality scoring
4. **Graph**: Connect bundles via quality-validated edges
5. **Extraction**: Identify connected components as contigs
6. **Consensus**: Quality-weighted nucleotide voting

## Next Steps

1. Compare contigs against iRepertoire reference
2. Assess SHM variant calling accuracy
3. Evaluate runtime and memory scaling
4. Batch process all samples (FZ-20, FZ-94, FZ-97, FZ-122, FZ-83)
5. Conduct biological validation

## Technical Notes

- Compilation: C++11, optimized build
- Backward compatible: v1 completely untouched
- Input QC: All FASTQ quality scores preserved
- Memory efficient: Vectorized data structures
- Benchmarking: Not against iRepertoire in this phase

## Document Info

- Version: v2.0
- Date: 2026-09-22
- Sample: FZ-116
- Status: Pilot Complete
- Citation: Graph-TRUST4 v2 Quality-Adaptive Overlap Scoring

