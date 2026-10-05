# Graph-Based Assembly Design for TRUST4 BCR Reconstruction

**Date**: 2026-09-18  
**Phase**: PHASE 3 - Design Documentation (before implementation)  
**Status**: Design Review - Ready for Implementation

---

## Executive Summary

This document specifies a variable-overlap-length graph-based assembly algorithm to replace TRUST4's current greedy read-to-contig assembly. The key innovation is using a dynamic overlap graph where edge weights represent variable-length overlaps, enabling more flexible assembly decisions than the fixed k-mer greedy approach.

**Key Design Principles**:
- Reuse TRUST4's candidate read extraction (avoid re-processing FASTQ)
- Preserve biological constraints (IGH restriction, reference genes)
- Implement efficient overlap detection (avoid all-vs-all comparison)
- Build a manageable graph (thousands-to-tens-of-thousands of nodes for bulk RNA-seq)
- Support consensus generation compatible with existing TRUST4 annotation

---

## Part 1: Current TRUST4 Workflow Analysis

### 1.1 Complete Pipeline Flow

```
FASTQ Input
    ↓
[fastq-extractor]
    ├─ Load reference (IMGT, 1,489 genes)
    ├─ K-mer indexing of reference
    ├─ Scan FASTQ for candidate reads
    │  (reads with k-mer hits in reference)
    ├─ Barcode/UMI extraction (if applicable)
    └─ Output: candidate read FASTQ (*.fq)

Candidate Reads
    ↓
[trust4 main assembly]
    ├─ Load candidate reads
    ├─ K-mer indexing (default k=9)
    ├─ Greedy read-to-contig assembly
    │  (AddRead iterative greedy algorithm)
    ├─ Contig extension with mate-pair info
    ├─ Consensus generation
    └─ Output: preliminary contigs (*.out)

Contigs
    ↓
[annotator]
    ├─ IMGT reference alignment
    ├─ V/D/J/C gene assignment
    ├─ CDR3 identification
    ├─ CDR3 abundance computation
    └─ Output: annotated sequences (*.fa, *.cdr3.out)

Annotated Sequences
    ↓
[perl report scripts]
    ├─ Simplify CDR3 annotations
    ├─ Generate AIRR-format output
    └─ Output: clonotype table (report.tsv, airr.tsv)
```

### 1.2 Current Greedy Assembly Algorithm (AddRead in SeqSet.hpp)

**Current Approach**:
1. **GetOverlapsFromRead**: Query k-mer index to find all existing contigs that k-mer-match the new read
2. **Overlap Scoring**: For each contig, compute alignment similarity (matchCnt, indelCnt, similarity score)
3. **Greedy Selection**: Sort overlaps by (matchCnt DESC, similarity DESC, length DESC)
4. **Sequential Extension**: Try to extend each existing contig in order
   - If read fits → UpdateConsensus and return
   - If no fit → Create new contig
5. **Feedback Loop**: New contig immediately indexed for next reads

**Key Data Structures**:
- `_seqWrapper`: represents a contig
  - `consensus`: consensus sequence
  - `posWeight`: per-position nucleotide weights (for consensus)
  - `numRead`: count of reads contributing
- `_overlap`: represents read-to-contig alignment
  - `matchCnt`: matched bases (used for sorting)
  - `similarity`: edit-distance-based similarity
  - `seqIdx`, `strand`, coordinates

**Current Limitations for Graph Assembly**:
- Decisions are made sequentially (no backtracking)
- No global optimization
- High-abundance reads dominate early (may capture wrong path)
- Limited handling of branching/paralogs

---

## Part 2: Proposed Graph-Based Assembly Algorithm

### 2.1 Graph Definition

**Nodes**:
- Each candidate read is a node
- Each reference gene sequence is a special "anchor" node (optional)
- Total nodes: ~1-10 million for 100M-read bulk RNA-seq (depending on extraction rate)

**Edges**:
- Edge (read_i → read_j) exists if reads i and j have significant sequence overlap
- Edge weight: overlap_length (variable, 20-200 bp typical)
- Edge annotation: (overlap_region_coordinates, num_mismatches, identity%)
- Directed: Optional (can treat as undirected for simplicity)

### 2.2 Overlap Detection Strategy

**Problem**: All-vs-all pairwise comparison of N reads is O(N²), infeasible for 1-10M reads.

**Solution**: Multi-level candidate generation:

**Step 1: K-mer Indexing** (Same as TRUST4 but indexed on reads)
- Build k-mer index from all candidate reads
- K-mer size: 9-15 (tunable)
- Count k-mer frequency per read

**Step 2: Candidate Pair Generation**
- For each read: find all other reads sharing at least T k-mers
- T = 2-3 (tunable, controls pair density)
- Expected candidates per read: 10-100 (not 1-10M)
- Complexity: O(N·avg_candidates) instead of O(N²)

**Step 3: Overlap Validation**
- Only for candidate pairs, compute sequence alignment
- Use banded alignment (linear space, faster)
- Store: (read_i, read_j, overlap_length, identity%, overlap_coords)

**Step 4: Edge Filtering**
- Require minimum overlap: 20 bp (or sequence-length-dependent)
- Require minimum identity: 90-95% (tunable for error handling)
- Discard reciprocal false overlaps (e.g., exact substring matches in repetitive regions)

### 2.3 Variable Overlap Length Computation

Current TRUST4 uses fixed overlap detection based on k-mer index. We compute exact overlap:

```
For each candidate pair (read_i, read_j):
  For each k-mer hit location:
    Compute maximal overlap by:
      1. Try left-extension: extend backward until mismatch
      2. Try right-extension: extend forward until mismatch
      3. Track longest overlap region
      4. Compute identity: matched_bases / (overlap_length + gap_cost)
    
    If identity >= threshold:
      Record edge with overlap_length (variable)
```

**Why Variable is Better**:
- Fixed k=9 misses overlaps shorter than k (e.g., 8bp don't trigger assembly)
- Fixed k=9 gives false positives in repetitive regions
- Variable length matches true biological signal

### 2.4 Graph Structure Representation

**Compressed Adjacency List** (memory-efficient for sparse graphs):
```
nodes = [read_1, read_2, ..., read_N]
edges: for each node, maintain list of neighbor nodes and edge weights

Example:
read_1234 → neighbors:
  read_5678: overlap_len=45, identity=97%
  read_9012: overlap_len=52, identity=95%
  read_3456: overlap_len=38, identity=92%
```

**In-Memory Storage**:
- Nodes: ~50 bytes per node × 10M reads = 500 GB (for extreme case; typical bulk RNA-seq has fewer)
- Edges: ~40 bytes per edge × avg_edges_per_node, typically sparse (10-100 neighbors/node)
- **Practical for FZ-116**: 1-50 million candidate reads, manageable with 64GB RAM

### 2.5 Graph Cleaning and Error Handling

**Sequencing Errors**:
- Reads with 1-2 errors in overlap region should still be connected
- **Approach**: Allow 1-2 mismatches per 50 bp overlap (90-96% identity)
- Creates "error bubbles": multiple nodes representing same sequence with small errors

**Duplicate Reads**:
- Exact duplicates should be merged
- **Approach**: Merge reads with 100% identity in overlap region to single node with weight=count
- Reduces node count by 10-50% in typical data

**Low-Complexity Regions**:
- Repetitive sequences (homopolymers, simple tandem repeats) create spurious overlaps
- **Approach**: Mask simple repeats (AAAA..., GCGCGC...) or require higher identity threshold in low-complexity regions

**Branching and Paralogs**:
- BCR sequences from different genes may overlap (especially V/J junctions)
- **Approach**: Don't suppress—let graph contain branches; resolve during contig assembly

### 2.6 Graph Traversal for Contig Assembly

**Goal**: Extract contigs (paths) from the graph that represent original BCR sequences.

**Strategy: Greedy Path Extraction**:
1. **Abundance Weighting**: Assign node weight = read frequency (from input)
2. **Start Nodes**: High-abundance nodes with no incoming edges (or low in-degree)
3. **Greedy Traversal**:
   - From start node, follow edges to maximize path weight
   - At branch points, choose edge with highest downstream abundance
   - Continue until reaching low-abundance node or dead-end
4. **Bubble Popping**: Merge parallel paths with high identity (error bubbles)
5. **Output**: List of contigs with assigned reads

**Result**: Multiple contigs if:
- True paralogs / alleles present
- Errors create separate paths
- Strain variation in sample

**Consensus Generation**:
- For each contig, compute per-position nucleotide weights (same as TRUST4)
- Vote on consensus at each position
- Preserve ambiguity with Ns if needed

### 2.7 Computational Complexity Analysis

**Time Complexity**:

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| K-mer indexing of reads | O(N·L) | N=reads, L=avg length (~150 bp) |
| Candidate pair generation | O(N·K²) | K=avg k-mers per read (~100), most false negatives not compared |
| Pairwise alignment | O(E·L²) | E=num edges, L=read length, banded DP ~30-50 bp strip |
| Graph traversal | O(N + E) | Standard DFS/BFS |
| Consensus generation | O(N·L) | One pass per position |
| **Total** | **O(N·L + E·L²)** | With E << N² |

**Practical Estimate for FZ-116** (1M candidate reads):
- K-mer indexing: ~2-5 minutes
- Candidate generation: ~5-10 minutes
- Pairwise alignment: ~30-60 minutes (if E ~ 10M edges)
- Graph traversal: <1 minute
- Consensus: ~5 minutes
- **Total: ~1-2 hours** (vs ~30 min for greedy TRUST4)

### 2.8 Memory Requirements

**Typical Bulk RNA-seq (100M reads → ~5M candidate reads)**:

| Component | Size |
|-----------|------|
| Reads in memory | 5M × 150 bp × 1 byte = 750 MB |
| K-mer index | 5M reads × 100 k-mers × 4 bytes = 2 GB |
| Edge list (adj matrix) | 5M nodes × 50 neighbors × 20 bytes/edge = 5 GB |
| Working memory (alignment buffers, temp) | ~2 GB |
| **Total** | **~10 GB** |

**For FZ-116** (smaller, ~1M candidate reads): ~2-3 GB

**Feasible on**:  64GB workstations (commodity hardware)

### 2.9 Handling Special Cases

**Highly Abundant Clones**:
- A single clonotype may have 10K-100K reads
- Merging exact duplicates reduces nodes by 10-50×
- Graph remains manageable

**Bimodal Repertoires** (many rare + few abundant):
- Rare clones: low-coverage paths
- Abundant clones: high-coverage paths
- Paths separated naturally by coverage thresholds

**V/D/J Boundaries**:
- Overlaps at junctions may be short (15-30 bp) but clear
- Variable overlap handles this; fixed k=9 may miss boundaries

**D Gene Ambiguity**:
- D genes in human repertoire: ~30 per chain
- Highly homologous
- Graph may contain D-swap branches
- **Solution**: Annotation phase (after contig assembly) resolves via IMGT alignment

---

## Part 3: Implementation Architecture

### 3.1 File Organization

```
algorithms/
├── TRUST4/                           (existing, unchanged)
│   ├── FastqExtractor.cpp
│   ├── main.cpp
│   ├── Annotator.cpp
│   ├── SeqSet.hpp
│   └── ... (other files)
│
└── Graph-TRUST4/                     (new)
    ├── GraphAssembler.hpp            (core graph class)
    ├── OverlapDetector.hpp           (candidate overlap finding)
    ├── ContigExtractor.hpp           (path extraction)
    ├── ConsensusBuilder.hpp          (consensus from reads)
    ├── GraphAssembler.cpp            (main executable)
    ├── Makefile
    └── README.md (implementation notes)
```

### 3.2 Reuse of TRUST4 Components

**CAN Reuse**:
- FastqExtractor.cpp output (candidate reads already extracted)
- KmerIndex.hpp for k-mer lookups
- AlignAlgo.hpp for alignment scoring (adapt as needed)
- SeqSet.hpp data structures (_seqWrapper, _overlap) for contigs/overlaps
- Annotator.cpp on graph-generated contigs

**MUST Reuse for Correctness**:
- Same reference gene handling (IMGT loading)
- Same V/J/C annotation pipeline
- Same CDR3 extraction logic

**CANNOT Reuse**:
- SeqSet.AddRead() (greedy assembly—this is what we replace)
- Main read-to-contig assignment loop

### 3.3 Data Flow

```
FZ-116 Raw FASTQ
    ↓
[TRUST4 fastq-extractor - REUSED]
    → Candidate reads FASTQ
    → K-mer index

Candidate Reads + K-mer Index
    ↓
[Graph-TRUST4 OverlapDetector]
    ├─ K-mer candidate pair finding
    ├─ Sequence overlap calculation
    └─ → Graph edges (variable overlap)

Read Graph
    ↓
[Graph-TRUST4 ContigExtractor]
    ├─ Path traversal
    ├─ Bubble merging
    └─ → Contigs (path-based)

Contigs
    ↓
[TRUST4 Annotator - REUSED]
    ├─ V/D/J/C alignment
    ├─ CDR3 extraction
    └─ → Annotated contigs

Annotated Contigs
    ↓
[TRUST4 Perl report scripts - REUSED]
    → Clonotype report (report.tsv)
    → AIRR output (airr.tsv)
```

### 3.4 Input/Output Specification

**Input**:
- Candidate reads FASTQ: `TRUST_FZ-116_toassemble_1.fq`, `TRUST_FZ-116_toassemble_2.fq`
  (Already generated by existing fastq-extractor)
- Reference FASTA: `human_IMGT+C.fa`
- Parameters: overlap_threshold (20bp), identity_threshold (90%), k-mer_size (9)

**Output**:
- Contigs FASTA: `TRUST_FZ-116_graph_contigs.fa`
  (Compatible format with Annotator.cpp)
- Graph statistics: `TRUST_FZ-116_graph_stats.txt`
  - Node count
  - Edge count
  - Connected components
  - Contig count
  - Contig length distribution

---

## Part 4: Pilot Implementation Strategy (PHASE 4)

### 4.1 Pilot Scope

**Sample**: FZ-116 only  
**Input**: Existing FZ-116 candidate reads (TRUST_FZ-116_toassemble_*.fq)  
**Do NOT re-extract**: Use fastq-extractor output that already exists

### 4.2 Implementation Phases (Pilot)

1. **Phase 4a: Data Loading & Indexing**
   - Load candidate reads into memory
   - Build k-mer index (reuse TRUST4 code if possible, or implement)
   - Validate index correctness on known overlaps

2. **Phase 4b: Overlap Detection**
   - Candidate pair generation from k-mer index
   - Pairwise alignment (banded DP)
   - Edge filtering (threshold application)
   - Output: edge list

3. **Phase 4c: Graph Construction**
   - Build adjacency list from edge list
   - Assign node weights (read frequency)
   - Output: graph statistics (node count, edge count, density)

4. **Phase 4d: Contig Extraction**
   - Implement path traversal (abundance-weighted)
   - Bubble merging heuristic
   - Output: contig set

5. **Phase 4e: Consensus Generation**
   - For each contig, compute per-position nucleotide weights
   - Consensus by majority vote
   - Output: FASTA format compatible with Annotator

6. **Phase 4f: Integration Testing**
   - Pipe contigs to existing Annotator.cpp
   - Verify annotation output format
   - Generate report.tsv

### 4.3 Expected Pilot Results

**For FZ-116**:
- Candidate reads: ~1.9M (from TRUST4 log)
- Expected graph nodes: ~1.9M
- Expected graph edges (at 10-100 neighbors/node): 20-200M
- Expected contigs: 5K-50K
- Expected IGH contigs: ~3-10K
- Expected clonotypes: ~1-5K (after deduplication)

---

## Part 5: Success Criteria and Validation

### 5.1 Correctness Validation

After Graph-TRUST4 completes on FZ-116:

1. **Output Format Compatibility**
   - Annotator.cpp accepts contig FASTA without modification ✓
   - Report generation produces valid TSV ✓

2. **Data Sanity**
   - Contig count: 1K-10K (reasonable for BCR bulk sample)
   - IGH clonotypes: 1K-5K
   - CDR3aa length: 10-50 bp (typical for antibodies)
   - Read assignments: each contig has identified reads

3. **Graph Statistics Sanity**
   - Graph is connected or few large components ✓
   - Average degree: 10-100 (not fully connected, not sparse)
   - No singleton nodes (isolated reads)

### 5.2 Benchmark Comparison

Once contigs are annotated, run existing benchmark pipeline:

```
Graph-TRUST4 contigs
    ↓
[benchmark_all_samples_published.py - REUSED]
    ├─ Load FZ-116 graph contigs
    ├─ CDR3nt [3:-3] trimming
    ├─ V+J+C+CDR3nt matching
    ├─ MAX abundance deduplication
    └─ → Precision, Sensitivity, Pearson r

Compare to TRUST4 baseline:
- Precision: 0.3396 (baseline)
- Sensitivity: 0.0327 (baseline)
- Pearson r: 0.6604 (baseline)
```

**Success**: Similar or better precision/sensitivity (within ±5%)

### 5.3 Computational Comparison

Record for FZ-116:
- Wall-clock runtime (minutes)
- Peak memory (MB)
- CPU utilization
- Candidate reads processed
- Final clonotype count
- Benchmark metrics

---

## Part 6: Risk Analysis and Mitigation

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Graph too large to fit in memory | Implementation failure | Implement streaming/chunking; test with progressively larger samples |
| Overlap detection too slow | Timeout on FZ-116 | Use efficient k-mer indexing; profile alignment bottlenecks early |
| Poor contig quality (fragmentation) | Benchmark degrades | Tune overlap thresholds; implement advanced path selection; compare to TRUST4 manually |
| Integration with Annotator fails | Incomplete pipeline | Test FASTA format compatibility early; verify column order in report |
| Graph disconnected (isolated contigs) | Missing output | Implement connected-component detection; lower overlap threshold if needed |
| D gene misassignment worsens | Benchmark metric down | Accept this (D is ambiguous); don't attempt D-specific graph modifications |

---

## Part 7: Rollout Plan (Post-Pilot)

After FZ-116 pilot succeeds:

1. **Run FZ-20** (second-smallest sample)
2. **Run FZ-97** (third-smallest)
3. **Run FZ-122, FZ-94, FZ-83** in sequence (largest last)

**No parallel multi-sample execution yet** (develop and test on sequential single-sample runs first).

---

## Conclusion

This design specifies a computationally feasible, biologically meaningful graph-based assembly to replace TRUST4's greedy approach. The graph captures branching and error structure explicitly, enabling better consensus generation and potential discovery of missing clonotypes or allelic variants.

**Key Advantages**:
- Variable overlap length matches true biological overlaps
- Graph structure supports backtracking and global optimization
- Manageable computational cost (~1-2 hours for bulk RNA-seq)
- Full integration with existing TRUST4 annotation pipeline

**Key Unknowns** (to resolve during implementation):
- Actual edge density (10 vs 100 neighbors/node changes scaling)
- Accuracy of banded-alignment overlap calculation
- Quality of greedy path extraction for consensus
- Whether benchmark metrics improve or stay similar

**Next Step**: Proceed to PHASE 4 - Implement pilot on FZ-116.

---

**Document Version**: 1.0  
**Status**: Ready for Implementation  
**Prepared by**: Code Analysis (PHASE 1) + Design (PHASE 3)
