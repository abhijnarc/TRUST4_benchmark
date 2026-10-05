# Published TRUST4 Benchmark - All Samples Summary

**Date**: 2026-09-18
**Methodology**: Published evaluation from TRUST4_manuscript_evaluation/bulk/bcrval.py
**Samples**: 6 samples (FZ-116, FZ-20, FZ-83, FZ-94, FZ-97, FZ-122)
**Processing**: Individual sample processing, no pooling

---

## Results by Sample

| Sample | TRUST4 Clonotypes | iRepertoire Clonotypes | Matches | Precision | Sensitivity | Pearson r | D Agreement | Isotype Exact | Isotype Collapsed |
|--------|-------------------|------------------------|---------|-----------|-------------|-----------|-------------|---------------|-------------------|
| FZ-116 | 11,233 | 116,632 | 3,815 | 0.3396 | 0.0327 | 0.6604 | 0.4178 | 0.2375 | 1.0000 |
| FZ-20 | 7,567 | 68,064 | 2,756 | 0.3642 | 0.0405 | 0.7587 | 0.4793 | 0.3552 | 1.0000 |
| FZ-83 | 21,177 | 141,843 | 7,896 | 0.3729 | 0.0557 | 0.5918 | 0.4350 | 0.2618 | 1.0000 |
| FZ-94 | 17,955 | 117,969 | 5,957 | 0.3318 | 0.0505 | 0.5538 | 0.4249 | 0.2419 | 1.0000 |
| FZ-97 | 7,692 | 71,375 | 2,333 | 0.3033 | 0.0327 | 0.3186 | 0.4561 | 0.3476 | 1.0000 |
| FZ-122 | 4,034 | 80,663 | 1,343 | 0.3329 | 0.0166 | 0.6805 | 0.4989 | 0.1854 | 1.0000 |

---

## Summary Statistics

| Metric | Mean | Std Dev | Min | Max |
|--------|------|---------|-----|-----|
| Precision | 0.3408 | 0.0250 | 0.3033 | 0.3729 |
| Sensitivity | 0.0381 | 0.0140 | 0.0166 | 0.0557 |
| Pearson r | 0.5940 | 0.1526 | 0.3186 | 0.7587 |
| D Agreement | 0.4520 | 0.0321 | 0.4178 | 0.4989 |
| Isotype Exact | 0.2716 | 0.0668 | 0.1854 | 0.3552 |

---

## Notes

- **Precision**: Percentage of TRUST4 clonotypes matching iRepertoire
- **Sensitivity**: Percentage of iRepertoire clonotypes matching TRUST4
- **Pearson r**: Abundance correlation (TRUST4 count vs iRepertoire copy)
- **D Agreement**: Fraction of matched clonotypes with identical D genes (expected <50%, D is ambiguous)
- **Isotype Exact**: Fraction with identical original isotypes (IGHA1, IGHG1, etc.)
- **Isotype Collapsed**: Always 100% (by definition of matching key)
- **Samples processed independently**: No pooling, no combining datasets

## Data Files

For each sample, individual results in `results/benchmark/{sample}/published_eval/`:
- `trust4_parsed.tsv` — TRUST4 clonotypes (trimmed)
- `irep_parsed.tsv` — iRepertoire clonotypes
- `matches.tsv` — Matched clonotype pairs
- `metrics.tsv` — Precision/sensitivity/correlation

---

**Status**: ✓ Baseline established. Ready for next analysis step.
