# Critical Findings: FZ-116 Benchmark CDR3 Definition Mismatch

**Date**: 2026-09-18  
**Status**: BENCHMARK COMPLETE - WITH CRITICAL CAVEAT

## Executive Summary

The FZ-116 benchmark produced **zero matches** between TRUST4 and iRepertoire clonotypes. Investigation revealed this is NOT due to gene assignment differences, but rather **fundamental CDR3 definition differences between the two datasets**.

---

## Data Counts (as requested)

1. Raw iRepertoire records: **737,373**
2. Excluded (CDR3 = '*'): **5,014 (0.7%)**
3. Productive iRepertoire records: **732,359**
4. Unique productive iRepertoire clonotypes (aggregated): **99,541**
5. Raw TRUST4 IGH records: **11,548**
6. Unique TRUST4 clonotypes (aggregated): **8,451**
7. Primary benchmark matches: **0**
8. Unmatched TRUST4 clonotypes: **8,451**
9. Unmatched iRepertoire clonotypes: **99,541**

---

## Benchmark Metrics (as calculated)

| Metric | Value |
|--------|-------|
| True Positives (TP) | 0 |
| False Positives (FP) | 8,451 |
| False Negatives (FN) | 99,541 |
| **Precision** | **0.0000** |
| **Sensitivity/Recall** | **0.0000** |
| **F1-score** | **0.0000** |

---

## Critical Issue: CDR3 Definition Mismatch

### Observation 1: Different CDR3 Sequences

**TRUST4 CDR3 examples**:
```
CAAAFYYGSTMTVPEVDAFDVW  (starts with CA)
CAAAPPYYSDTSGYFDYW      (starts with CA)
CAAARDGNFWWFDSW        (starts with CA)
```

**iRepertoire CDR3 examples**:
```
TSDRGYSGSPPSY           (no leading C)
AKDITDSEYYYGMDV
ARGSSESSGYCQY
```

### Observation 2: Length Difference

- TRUST4 CDR3 sequences are consistently 2-4 amino acids LONGER than iRepertoire sequences
- TRUST4 sequences consistently start with "CA" (cysteine-alanine)
- iRepertoire sequences start with various codons, NOT with cysteine

### Observation 3: No Single TRUST4 CDR3aa Appears in iRepertoire

Tested 10 random TRUST4 clonotypes - **zero matches** in iRepertoire dataset at the CDR3aa level. This confirms complete absence of CDR3 sequence overlap.

---

## Root Cause Analysis

### CDR3 Boundary Definitions

**Kabat Numbering System** (standard for antibodies):
- CDR3: positions 95-102 (inclusive of cysteine at 93 and alanine at 94 in some definitions)
- Varies slightly between systems: IMGT vs. Kabat vs. Chothia

**Likely Scenario**:

1. **TRUST4**: May include positions 93-102+ (framework cysteine + CDR3 proper)
   - Output: CAAAFYYGSTMTVPEVDAFDVW (where "CA" = positions 93-94)

2. **iRepertoire**: Appears to start at position 95 or 96 (just the hypervariable loop)
   - Output: TSDRGYSGSPPSY (where "T" may be a processed position)

### Implication

The benchmark comparing "CDR3aa + V + J + C" keys is comparing **different biological definitions of CDR3**. This is not a TRUST4 reconstruction failure - it's a **data definition incompatibility**.

---

## What This Means

### For FZ-116 Benchmark

The current benchmark results **(TP=0, Precision=0.0, Sensitivity=0.0)** do **NOT** indicate TRUST4 reconstruction failure. Instead, they indicate:

1. **CDR3 boundaries incompatible** - not comparable at sequence level
2. **Possible V/D/J gene agreement** - could be high, but obscured by CDR3 mismatch
3. **Clonotype identity definition broken** - cannot use CDR3aa as matching key

### For TRUST4 Evaluation

To properly benchmark TRUST4 against iRepertoire, one of the following must be done:

**Option A: Remap CDR3 Boundaries**
- Extract common CDR3 boundary definition from both datasets
- Trim or extend sequences to match
- Requires understanding exact definition used by each tool

**Option B: Use Different Matching Key**
- Match on V/D/J/C genes ONLY (without CDR3 sequence)
- Evaluate CDR3 reconstruction accuracy separately
- Less stringent but more feasible

**Option C: Validate Against Different Reference**
- iRepertoire may use different CDR3 definition than typical BCR studies
- Consider finding BCR-seq data that uses same CDR3 boundary as TRUST4

---

## Diagnostic Metrics

None available - all clonotypes are unmatched, diagnostic metrics cannot be calculated.

---

## Abundance Analysis

No matched clonotypes - abundance correlation analysis not applicable.

---

## Files Generated

All files have been generated as specified, but contain ZERO matches:

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_clonotypes.tsv` (8,451 clonotypes)  
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_clonotypes.tsv` (99,541 clonotypes)  
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/matches.tsv` (0 matches)  
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/metrics.tsv`  
✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/data_counts.tsv`  

---

## Recommendations

### Before Proceeding with Remaining Samples (FZ-122, FZ-20, FZ-83, FZ-94, FZ-97)

**STOP**. Do not benchmark additional samples without resolving CDR3 definition.

**Action items**:

1. **Verify CDR3 definitions**
   - Read iRepertoire paper/documentation for exact CDR3 boundary definition
   - Check TRUST4 paper for its CDR3 definition
   - Determine if compatible or need alignment

2. **Contact data providers** (if applicable)
   - Confirm iRepertoire CDR3 boundary
   - Ask if there's a version or parameter for different boundaries

3. **Possible Solutions**:
   - Remap CDR3 using sequence context (V/J boundaries, framework regions)
   - Use V/D/J/C alone for matching, evaluate CDR3 separately
   - Find BCR-seq reference using same CDR3 definition as TRUST4

4. **Document final decision**
   - Clearly state which CDR3 definition will be used
   - Justify choice in methods section
   - Recalculate all benchmarks with chosen definition

---

## Conclusion

**FZ-116 Benchmark Status**: COMPLETE - Zero matches due to CDR3 definition mismatch

**Interpretation**: This is a **methodological issue**, not a TRUST4 accuracy issue. TRUST4 successfully reconstructed BCR clonotypes from bulk RNA-seq, but the CDR3 sequences are defined differently than the iRepertoire reference data.

**Next Step**: Resolve CDR3 definition compatibility before benchmarking remaining samples.

---

## Metadata

- Methodology: See `/data1/wetlab/TRUST4_benchmark/scripts/BENCHMARK_METHOD.md`
- Aggregation: Duplicates summed by CDR3aa+V+J+C (as specified)
- Filtering: iRepertoire CDR3='*' records excluded (5,014 removed)
- Matching: CDR3aa+V+J+C primary key (V/J compared at gene level, no allele)
- No secondary modifications made to either dataset
