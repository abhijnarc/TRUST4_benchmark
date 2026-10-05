# TRUST4 vs iRepertoire Benchmarking Methodology

**Date**: 2026-09-18  
**Sample**: FZ-116 (SRR7882936)  
**Script**: `benchmark_trust4_vs_irep.py`

---

## 1. Overview

This document specifies the exact matching criteria, normalization rules, and metric definitions used to benchmark TRUST4 BCR reconstruction against iRepertoire BCR-seq reference data.

**Key Principle**: Match TRUST4 clonotypes to iRepertoire based on biological identity (CDR3 + V + J + C), not abundance equivalence.

---

## 2. Data Filtering

### 2.1 TRUST4 Filtering

**Input**: `TRUST_FZ-116_report.tsv` (29,803 clonotypes)

**Filter**: Restrict to IGH (heavy chain) only

**Criterion**: Extract first V gene, check if starts with "IGH"

**Result**: 11,548 IGH clonotypes selected for benchmarking

**Rationale**: iRepertoire reference contains IGH only; light chains (IGK, IGL) and TCR cannot be validated

### 2.2 iRepertoire Data

**Input**: `FZ-116.csv.gz` (737,373 records)

**Status**: All records are IGH (human BCR heavy chain)

**No filtering applied**

---

## 3. Normalization Rules

### 3.1 CDR3 Amino Acid (CDR3aa)

**Original format**: Variable (uppercase, mixed case, unknown)

**Normalization**: Convert to uppercase

**Rationale**: Simple representation normalization

**Example**:
- Input: `CANGPPGLSRGVTRRHYYFGMDVW` → Output: `CANGPPGLSRGVTRRHYYFGMDVW`
- Input: `cangppglsrgvtrrhyyfgmdvw` → Output: `CANGPPGLSRGVTRRHYYFGMDVW`

**Preservation**: All 20 standard amino acids and stop codon (_) preserved as-is

### 3.2 CDR3 Nucleotide (CDR3nt)

**Original format**: Variable (uppercase, mixed case)

**Normalization**: Convert to uppercase

**Rationale**: Simple representation normalization

**Note**: iRepertoire sequences may use mixed case to indicate masked/uncertain regions. Converted to uppercase for comparison; exact nucleotide matching used as a **diagnostic metric, not primary criterion**.

**Example**:
- iRep: `acgagtgatcgaggttatagtgggagcccgccctcctac` → Normalized: `ACGAGTGATCGAGGTTATAGTGGGAGCCCGCCCTCCTAC`
- TRUST4: `ACGAGTGATCGAGGTTATAGTGGGAGCCCGCCCTCCTAC` → Normalized: `ACGAGTGATCGAGGTTATAGTGGGAGCCCGCCCTCCTAC`

### 3.3 V Gene Name

**iRepertoire format**: `hIGHV3-23*01` (has "h" prefix, allele suffix)

**TRUST4 format**: `IGHV3-23*05` (no "h" prefix, allele suffix)

**Normalization steps**:
1. Remove "h" prefix if present
2. Remove allele suffix (`*01`, `*02`, etc.)
3. Compare at gene level (e.g., `IGHV3-23`)

**Rationale**: Allele differences are not informative for clonotype identity; PCR/sequencing depth drives allele calling. Gene segment identity is the relevant biological criterion.

**Result**: `IGHV3-23`

**Handling multiple candidates**: TRUST4 may report multiple V genes separated by `|`, ranked by similarity. Use first (highest-ranked) only.

**Example**:
- iRep: `hIGHV3-23*01` → Normalized: `IGHV3-23`
- TRUST4: `IGHV3-23*05` → Normalized: `IGHV3-23`
- TRUST4 multi: `IGHV3-23*01|IGHV3-23D*01` → First: `IGHV3-23*01` → Normalized: `IGHV3-23`

### 3.4 J Gene Name

**Identical to V gene normalization**:
1. Remove "h" prefix
2. Remove allele suffix

**Example**:
- iRep: `hIGHJ6*02` → Normalized: `IGHJ6`
- TRUST4: `IGHJ6*02` → Normalized: `IGHJ6`

### 3.5 D Gene Name

**Normalization**: Same as V/J (remove "h" prefix, remove allele suffix)

**Status in benchmarking**: **NOT used for primary TP/FP/FN classification**

**Rationale**: D gene is highly variable and harder to call; included as separate diagnostic metric

**Example**:
- iRep: `hIGHD3-10*01` → Normalized: `IGHD3-10`
- TRUST4: `IGHD3-10*01` → Normalized: `IGHD3-10`

### 3.6 C Gene (Constant Region / Isotype)

**iRepertoire format**: `hIGHA1`, `hIGHA2`, `hIGHG1`, etc. (with "h" prefix)

**TRUST4 format**: `IGHA1`, `IGHA2`, `IGHG1`, etc. (no "h" prefix)

**Normalization**:
1. Remove "h" prefix only
2. **Do NOT remove allele suffix** (e.g., preserve IGHA1 vs IGHA2 distinction)

**Rationale**: Isotype is biologically meaningful (IgA vs IgG vs IgM); alleles carry immunological significance. Keep full identity.

**Example**:
- iRep: `hIGHA1*03` → Normalized: `IGHA1*03`
- TRUST4: `IGHA1` → Normalized: `IGHA1`

---

## 4. Primary Matching Criterion

### 4.1 Definition of True Positive (TP)

A TRUST4 clonotype is a **True Positive** if ALL of the following match an iRepertoire record:

```
CDR3aa (normalized) AND
V gene (normalized, no allele) AND
J gene (normalized, no allele) AND
C gene (normalized, preserve isotype)
```

Both must match exactly after normalization.

### 4.2 Matching Algorithm

1. **Candidate lookup**: For each TRUST4 record, use CDR3aa as candidate key
2. **Validation**: Check V + J + C agreement
3. **Confirmation**: If all four fields match, classify as TP

### 4.3 Duplicate Handling

**STOP CONDITION**: If either dataset contains duplicate clonotypes (same CDR3aa + V + J + C), **STOP before metric calculation**.

**Action required**: Manually define aggregation strategy (sum counts, average, or treat as separate).

**Rationale**: Duplicates indicate data quality issues or representation ambiguities; cannot proceed without explicit decision.

---

## 5. Metrics Definitions

### 5.1 Primary Benchmark Metrics

| Metric | Definition | Calculation |
|--------|-----------|-------------|
| **True Positive (TP)** | TRUST4 clonotype matches iRep by CDR3+V+J+C | Count of matches |
| **False Positive (FP)** | TRUST4 clonotype NOT in iRep | n_TRUST4 - TP - ambiguous |
| **False Negative (FN)** | iRep clonotype NOT in TRUST4 | n_iRep - TP - ambiguous |
| **Precision** | Fraction of TRUST4 calls that are correct | TP / (TP + FP) |
| **Sensitivity / Recall** | Fraction of iRep clones reconstructed by TRUST4 | TP / (TP + FN) |
| **F1-score** | Harmonic mean of precision and recall | 2 × (P × R) / (P + R) |

### 5.2 Secondary Diagnostic Metrics

Calculated separately for investigation (not part of primary benchmark):

- **CDR3aa-only matches**: Count of TRUST4 records with matching CDR3aa in iRep (ignoring V/J/C)
- **CDR3nt exact matches**: Among TP, count with exact nucleotide sequence match
- **V gene agreement**: Among CDR3-matched pairs, % with matching V gene
- **J gene agreement**: Among CDR3-matched pairs, % with matching J gene
- **C gene agreement**: Among CDR3-matched pairs, % with matching C gene
- **D gene agreement**: Among CDR3-matched pairs, % with matching D gene
- **Complete V/J/C agreement**: Among CDR3-matched pairs, % with all three matching

### 5.3 Abundance Concordance Analysis

**NOT a primary benchmark metric** (TRUST4 read counts ≠ iRep copy counts)

**Purpose**: Investigate correlation structure among matched clonotypes

**Measures**:
- Spearman correlation: TRUST4 `#count` vs iRepertoire `copy`
- Spearman correlation: TRUST4 `frequency` vs iRepertoire `copy`
- P-value for each correlation

**Interpretation**: Positive correlation ≠ quantitative equivalence; indicates rank-order consistency if present

**Restriction**: Calculated only for TP matches (n ≥ 2)

---

## 6. Ambiguous Matches

### 6.1 Definition

A TRUST4 clonotype that matches the CDR3+V+J+C key but finds **multiple iRepertoire entries** with identical key.

### 6.2 Handling

**Status**: Classified as `AMBIGUOUS_MULTI_IREP`

**Action**: Excluded from TP/FP/FN counts

**Rationale**: Cannot decide which iRep entry is "correct"; indicates potential duplicate in iRep or insufficient data to distinguish

**Reporting**: Counted and reported separately

---

## 7. Output Files

### 7.1 trust4_normalized.tsv

TRUST4 report with additional normalized columns:
- `cdr3aa_norm`: Uppercase CDR3 amino acid
- `cdr3nt_norm`: Uppercase CDR3 nucleotide
- `v_gene_norm`: Gene name without allele suffix
- `j_gene_norm`: Gene name without allele suffix
- `c_gene_norm`: Isotype name without "h" prefix
- `d_gene_norm`: Gene name without allele suffix

### 7.2 irep_normalized.tsv

iRepertoire data with normalized columns (same as above)

### 7.3 matches.tsv

Detailed matching results:

| Column | Meaning |
|--------|---------|
| `cdr3aa` | Matched CDR3 (amino acid) |
| `cdr3nt_trust4`, `cdr3nt_irep` | Nucleotide sequences (for comparison) |
| `cdr3nt_match` | Boolean: nucleotide sequences identical |
| `v_gene`, `j_gene`, `c_gene` | Matched gene names (normalized) |
| `d_gene_trust4`, `d_gene_irep` | D genes (for diagnostic comparison) |
| `d_gene_match` | Boolean: D genes agree |
| `trust4_count`, `trust4_frequency` | TRUST4 abundance measures |
| `irep_copy` | iRepertoire copy count |
| `match_status` | `TP` (true positive) or `AMBIGUOUS_MULTI_IREP` |
| `mismatch_reason` | Explanation if not TP |

### 7.4 metrics.tsv

Summary metrics:
- n_trust4_igh
- n_irep_total
- n_irep_unique_keys
- tp, fp, fn, ambiguous
- precision, sensitivity, recall, f1

### 7.5 abundance_correlation.tsv

Spearman correlations (TP pairs only):
- n_tp_pairs
- trust4_count_vs_irep_copy_spearman
- trust4_count_vs_irep_copy_pvalue
- trust4_frequency_vs_irep_copy_spearman
- trust4_frequency_vs_irep_copy_pvalue

---

## 8. Data Integrity Checks

### 8.1 Duplicates

Before metric calculation, script checks for duplicate clonotypes:

**TRUST4**: Check for multiple rows with identical (CDR3aa_norm, v_gene_norm, j_gene_norm, c_gene_norm)

**iRepertoire**: Check for multiple rows with identical (CDR3aa_norm, v_gene_norm, j_gene_norm, c_gene_norm)

**Action**: If duplicates found, **STOP** and report to user

**Rationale**: Duplicates require explicit aggregation decision (sum, average, weighted mean, etc.)

### 8.2 Missing Values

Handled as:
- `None` or `np.nan` → Recorded; excluded from comparisons
- `.` (TRUST4/iRep notation for absent gene) → Treated as `None`

---

## 9. Normalization Justification

### 9.1 Why Remove Allele Suffixes from V/J?

**Observation**: TRUST4 reports `IGHV3-23*05`, iRep reports `hIGHV3-23*01` for same clonotype

**Reason**: Allele calling depends on:
- PCR bias
- Sequencing depth
- Reference database version
- Alignment algorithm parameters

**Implication**: Allele differences do NOT indicate biological mismatch; are technical artifacts

**Solution**: Compare at gene level (e.g., `IGHV3-23` ← `*01` or `*05`)

### 9.2 Why Preserve C Gene Isotype?

**Observation**: C gene distinguishes IgA, IgG, IgM, IgD, IgE; also has variants (IGHA1 vs IGHA2)

**Reason**: Isotype and C-region genes are:
- Biologically significant (immune function differs by isotype)
- Part of published BCR benchmarks
- Used in repertoire classification

**Solution**: Preserve full C gene identity (IGHA1 ≠ IGHA2)

### 9.3 Why Use CDR3+V+J+C for Matching?

**Alternatives considered**:
- CDR3 alone: Too permissive (ignores gene identity)
- CDR3+V+J: Misses isotype info (C region is part of BCR)
- Nucleotide matching: Too strict (minor sequencing errors cause mismatches)

**Choice**: CDR3+V+J+C provides biologically meaningful matching while allowing for minor sequencing variation

---

## 10. Caveats and Limitations

1. **Abundance measures are uncorrelated by design**
   - TRUST4 reads from unselected RNA-seq
   - iRep copies from selected BCR-seq
   - No quantitative equivalence assumed

2. **D gene not used in primary match**
   - D gene calling is less reliable than V/J
   - Included as separate diagnostic metric

3. **Mixed-case sequences in iRep**
   - Converted to uppercase for comparison
   - Original masking information lost
   - Affects CDR3nt_match interpretation

4. **Duplicates cause STOP**
   - Cannot silently aggregate without bias
   - User must decide aggregation strategy

5. **No sequence quality scores**
   - TRUST4 provides quality/confidence (CDR3_score, germline_similarity)
   - iRep does not
   - Quality not factored into matching

---

## 11. Running the Script

```bash
cd /data1/wetlab/TRUST4_benchmark
python3 scripts/benchmark_trust4_vs_irep.py
```

**Output location**: `results/benchmark/FZ-116/`

**Expected output files**:
- trust4_normalized.tsv
- irep_normalized.tsv
- matches.tsv
- metrics.tsv
- abundance_correlation.tsv (if no duplicates)

**If duplicates found**: Script stops with diagnostic files for manual review

---

## 12. References

- TRUST4 publication: Song et al. Nat Methods (2021)
- AIRR standard: https://docs.airr-community.org/
- iRepertoire: BCR-seq benchmarking reference
