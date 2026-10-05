# Duplicate Clonotype Analysis - TRUST4 vs iRepertoire

**Date**: 2026-09-18  
**Script Output**: `benchmark_trust4_vs_irep.py` stopped at duplicate detection

---

## Summary: DUPLICATES DETECTED - ACTION REQUIRED

The benchmarking script successfully:
✅ Loaded both datasets
✅ Filtered TRUST4 to IGH (11,548 clones)
✅ Normalized gene names and CDR3 sequences
✅ **Detected duplicates** → STOPPED before metric calculation

---

## TRUST4 Duplicates

**Count**: 1,680 records (out of 11,548 IGH clones)
**Definition**: Multiple rows with identical (CDR3aa, V_gene, J_gene, C_gene)

### Root Cause Analysis

The same clonotype (CDR3aa + V + J + C) appears with **different nucleotide sequences**:

```
CDR3aa: CAAFSDRNFGDYRYYFEYW
V: IGHV1-58 | J: IGHJ4 | C: IGHA1

CDR3nt variant 1: ...TCGGATCGCAACTTC... (2 reads)
CDR3nt variant 2: ...TCGGATCGCAACTTC... (1 read)  ← Single nucleotide difference
```

### Interpretation

**Likely cause**: Somatic hypermutation or sequencing errors at the nucleotide level produce identical amino acid sequences but different DNA sequences.

**Example details**:
- CAAFSDRNFGDYRYYFEYW (IGHV1-58/IGHJ4/IGHA1): 
  - Variant 1: Read count=2, Frequency=5.76e-05
  - Variant 2: Read count=1, Frequency=2.87e-05
  
- CAAGFACMAW (IGHV3-23/IGHJ4/IGHA1):
  - Variant 1: Read count=3
  - Variant 2: Read count=1
  
**Percentage of TRUST4 IGH affected**: 1,680 / 11,548 = **14.5%**

---

## iRepertoire Duplicates

**Count**: 677,224 records (out of 737,373 total)
**Definition**: Multiple rows with identical (CDR3aa, V_gene, J_gene, C_gene)

### Root Cause Analysis

**PRIMARY ISSUE**: **Stop codons (CDR3aa = '*')**

Many iRepertoire records have CDR3(pep) = '*' (stop codon), indicating **non-productive rearrangements** or **sequencing artifacts**.

```
CDR3aa: * (stop codon)
V: hIGHV1-18  | J: hIGHJ4  | C: hIGHA2

Example duplicates:
- hIGHD3-9*01 with hIGHJ4*02/hIGHA2*04: 42 separate records
- hIGHD2-2*03 with hIGHJ3*02/hIGHA2*04: 2 separate records
```

**D gene diversity despite identical V/J/C**: Indicates that when CDR3 is missing/corrupt, different D genes are called for the same context. These are likely non-productive sequences with annotation errors.

### Secondary Issues

Among valid (non-stop) CDR3 sequences, also see duplicates:
- Same clonotype with different D gene assignments
- Same clonotype with different allele-level V or J assignments

**Percentage of iRep affected**: 677,224 / 737,373 = **91.9% have duplicates**

**Percentage with stop codons**: Requires counting

---

## Critical Decision Point

### Option A: Deduplicate by Summing Counts

**For TRUST4**: Sum read counts and frequency for duplicate clonotypes
- Before: 11,548 clones
- After: ~9,868 unique clones (1,680 duplicates → ~840 unique after summing)

**For iRepertoire**: Sum copy counts for duplicate clonotypes
- Before: 737,373 records
- After: ~60,000 unique clones (estimate, pending stop codon filtering)

**Pros**: Simple, preserves quantitative info
**Cons**: Loses information about nucleotide-level variants (which may be biologically meaningful)

### Option B: Keep Highest-Count Variant

**For TRUST4**: Keep only highest-count variant of each clonotype
- Before: 11,548 clones
- After: ~9,868 unique clones

**For iRepertoire**: Keep only one record per unique key
- Before: 737,373 records
- After: ~60,000 unique clones

**Pros**: Simpler deduplication, avoids double-counting
**Cons**: Discards lower-abundance variants (may lose real biology)

### Option C: Filter Out Non-Productive (Stop Codons)

**For iRepertoire**: Remove all records where CDR3aa = '*'
- Before: 737,373 records
- After: ~60,000 productive clones (estimated)
- Then deduplicate remaining

**Pros**: Removes clear artifacts
**Cons**: TRUST4 may have non-productive sequences too (need to check)

### Option D: Hybrid Approach

1. **Filter iRepertoire to productive only** (remove CDR3 = '*')
2. **Deduplicate iRepertoire by summing copy counts**
3. **For TRUST4**: Deduplicate by summing counts (keep CDR3aa variants separate)
4. **Then benchmark**: 

**Rationale**: 
- iRepertoire stop codons are clear errors
- TRUST4 somatic variants deserve separate accounting (may be biologically interesting)

---

## Matching Logic (Proposed)

The script defines matching as:

```python
Match = (
    TRUST4.CDR3aa_norm == iRep.CDR3aa_norm AND
    TRUST4.V_gene_norm == iRep.V_gene_norm AND
    TRUST4.J_gene_norm == iRep.J_gene_norm AND
    TRUST4.C_gene_norm == iRep.C_gene_norm
)
```

**This is correct**, but duplicates must be resolved first.

### Proposed Deduplication Before Matching

**TRUST4 deduplication**:
- Group by (CDR3aa_norm, V_gene_norm, J_gene_norm, C_gene_norm)
- Sum: `count` and recalculate `frequency`
- Result: ~9,868 unique clones

**iRepertoire deduplication**:
- Remove rows where CDR3aa_norm = '*' (stop codons)
- Group by (CDR3aa_norm, V_gene_norm, J_gene_norm, C_gene_norm)
- Sum: `copy`
- Result: ~60,000 productive clones

---

## Files Generated for Inspection

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_normalized.tsv`
   - All 11,548 TRUST4 records with normalized columns

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_normalized.tsv`
   - All 737,373 iRepertoire records with normalized columns

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/trust4_duplicates.tsv`
   - 1,680 records (records that form duplicates)

✅ `/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/irep_duplicates.tsv`
   - 677,224 records (records that form duplicates)

---

## Next Steps

1. **Choose deduplication strategy** (Option A, B, C, or D above)
2. **Implement in Python script** (modify deduplication logic)
3. **Re-run benchmarking** with deduplicated data
4. **Generate final metrics**: TP, FP, FN, precision, sensitivity, F1

---

## Key Questions to Resolve

1. Should TRUST4 somatic variants be merged (sum counts) or kept separate?
2. Should iRepertoire stop codons (CDR3 = '*') be filtered out or merged?
3. For abundance comparison: use raw counts or normalized? (currently using both)
4. If multiple iRep records match one TRUST4 clonotype after deduplication, how to handle?

---

## Status

**Current**: WAITING FOR USER INPUT on deduplication strategy
**Blocker**: Cannot calculate benchmark metrics until duplicates are resolved
**Effort**: ~30 min to implement chosen strategy and re-run

