#!/usr/bin/env python3
"""
Batch TRUST4 Benchmarking - Published Methodology

Applies the FZ-116 published methodology to all six samples:
FZ-116, FZ-20, FZ-83, FZ-94, FZ-97, FZ-122

Produces individual reports + consolidated summary.
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import gzip
from pathlib import Path

def remove_allele(gene_name):
    """Remove allele suffix from gene name"""
    if pd.isna(gene_name):
        return None
    gene = str(gene_name)
    if '*' in gene:
        gene = gene.split('*')[0]
    return gene

def collapse_isotype(c_gene):
    """Collapse isotype subclasses"""
    if pd.isna(c_gene):
        return None
    c = str(c_gene)
    if c in ['IGHA1', 'IGHA2']:
        return 'IGHA'
    elif c in ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4']:
        return 'IGHG'
    elif c in ['IGHD1', 'IGHD2', 'IGHD3', 'IGHD4']:
        return 'IGHD'
    elif c in ['IGHE1', 'IGHE2']:
        return 'IGHE'
    elif c in ['IGHM1', 'IGHM2']:
        return 'IGHM'
    else:
        return c

def trim_cdr3(cdr3_nt):
    """Trim 3 nucleotides from start and 3 from end"""
    if pd.isna(cdr3_nt):
        return None
    seq = str(cdr3_nt)
    if len(seq) < 6:
        return None
    return seq[3:-3]

def load_trust4(sample):
    """Load and parse TRUST4 report"""
    file_path = f'/data1/wetlab/TRUST4_benchmark/results/TRUST4/{sample}/TRUST_{sample}_report.tsv'
    df = pd.read_csv(file_path, sep='\t')
    df.columns = [col.lstrip('#') for col in df.columns]

    # Filter to IGH
    df = df[df['V'].str.startswith('IGH', na=False)].copy()

    # Parse and normalize
    df['v_gene'] = df['V'].apply(remove_allele)
    df['j_gene'] = df['J'].apply(remove_allele)
    df['c_gene_original'] = df['C'].apply(remove_allele)
    df['c_gene_collapsed'] = df['c_gene_original'].apply(collapse_isotype)
    df['cdr3nt_original'] = df['CDR3nt'].apply(lambda x: str(x).upper() if pd.notna(x) else None)
    df['cdr3nt_trimmed'] = df['cdr3nt_original'].apply(trim_cdr3)
    df['d_gene'] = df['D'].apply(remove_allele)

    # Remove records with incomplete CDR3
    df = df[df['cdr3nt_trimmed'].notna()].copy()

    return df

def load_irep(sample):
    """Load and parse iRepertoire reference"""
    file_path = f'/data1/wetlab/TRUST4_benchmark/reference/iRep/ng-bcr-validate/iRep/{sample}.csv.gz'
    df = pd.read_csv(file_path, compression='gzip')

    # Parse and normalize
    df['v_gene'] = df['V'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['j_gene'] = df['J'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['c_gene_original'] = df['C'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['c_gene_collapsed'] = df['c_gene_original'].apply(collapse_isotype)
    df['cdr3nt'] = df['CDR3(nuc)'].apply(lambda x: str(x).upper() if pd.notna(x) else None)
    df['d_gene'] = df['D'].apply(lambda x: remove_allele(str(x).lstrip('h')))

    return df

def coalesce_duplicates(df, is_trust4=True):
    """Coalesce identical clonotypes, keeping maximum abundance"""
    if is_trust4:
        abundance_col = 'count'
        key_cols = ['v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt_trimmed']
    else:
        abundance_col = 'copy'
        key_cols = ['v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt']

    df_coalesced = df.loc[df.groupby(key_cols)[abundance_col].idxmax()].copy()
    df_coalesced = df_coalesced.reset_index(drop=True)

    return df_coalesced

def match_clonotypes(df_t4, df_ir):
    """Match TRUST4 and iRepertoire clonotypes"""
    ir_index = {}
    for idx, row in df_ir.iterrows():
        key = (row['v_gene'], row['j_gene'], row['c_gene_collapsed'], row['cdr3nt'])
        if pd.notna(key[0]) and pd.notna(key[3]):
            ir_index[key] = row

    matches = []
    for idx, t4_row in df_t4.iterrows():
        key = (t4_row['v_gene'], t4_row['j_gene'], t4_row['c_gene_collapsed'], t4_row['cdr3nt_trimmed'])

        if key in ir_index:
            ir_row = ir_index[key]
            matches.append({
                'v_gene': t4_row['v_gene'],
                'j_gene': t4_row['j_gene'],
                'c_gene_collapsed': t4_row['c_gene_collapsed'],
                'cdr3nt_trimmed': t4_row['cdr3nt_trimmed'],
                'cdr3nt_original_trust4': t4_row['cdr3nt_original'],
                'cdr3nt_irep': ir_row['cdr3nt'],
                'cdr3aa_trust4': t4_row['CDR3aa'],
                'cdr3aa_irep': ir_row['CDR3(pep)'],
                'trust4_count': t4_row['count'],
                'irep_copy': ir_row['copy'],
                'd_gene_trust4': t4_row['d_gene'],
                'd_gene_irep': ir_row['d_gene'],
                'c_gene_trust4_original': t4_row['c_gene_original'],
                'c_gene_irep_original': ir_row['c_gene_original'],
            })

    return pd.DataFrame(matches)

def calculate_metrics(df_t4, df_ir, df_matches):
    """Calculate published-style metrics"""
    n_trust4 = len(df_t4)
    n_irep = len(df_ir)
    n_matches = len(df_matches)

    precision = n_matches / n_trust4 if n_trust4 > 0 else 0
    sensitivity = n_matches / n_irep if n_irep > 0 else 0

    metrics = {
        'n_trust4': n_trust4,
        'n_irep': n_irep,
        'n_matches': n_matches,
        'precision': precision,
        'sensitivity': sensitivity,
    }

    # Diagnostic metrics
    if len(df_matches) > 0:
        d_agree = (df_matches['d_gene_trust4'] == df_matches['d_gene_irep']).sum()
        iso_exact_agree = (df_matches['c_gene_trust4_original'] == df_matches['c_gene_irep_original']).sum()
        metrics['d_agreement'] = d_agree / len(df_matches)
        metrics['iso_exact_agreement'] = iso_exact_agree / len(df_matches)
    else:
        metrics['d_agreement'] = 0
        metrics['iso_exact_agreement'] = 0

    # Abundance correlation
    if len(df_matches) >= 2:
        df_corr = df_matches[['trust4_count', 'irep_copy']].dropna()
        if len(df_corr) >= 2:
            r, p = pearsonr(df_corr['trust4_count'], df_corr['irep_copy'])
            metrics['pearson_r'] = r
        else:
            metrics['pearson_r'] = np.nan
    else:
        metrics['pearson_r'] = np.nan

    # Isotype collapsed always 100% by definition (matching key)
    metrics['iso_collapsed_agreement'] = 1.0

    return metrics

def benchmark_sample(sample):
    """Run benchmark for a single sample"""
    print(f"\n{'='*80}")
    print(f"Benchmarking {sample}")
    print(f"{'='*80}")

    try:
        # Load
        print(f"Loading TRUST4...")
        df_t4 = load_trust4(sample)
        print(f"  Loaded {len(df_t4)} IGH records")

        print(f"Loading iRepertoire...")
        df_ir = load_irep(sample)
        print(f"  Loaded {len(df_ir)} records")

        # Coalesce
        df_t4 = coalesce_duplicates(df_t4, is_trust4=True)
        df_ir = coalesce_duplicates(df_ir, is_trust4=False)
        print(f"  Coalesced: TRUST4 {len(df_t4)} clonotypes, iRep {len(df_ir)} clonotypes")

        # Match
        df_matches = match_clonotypes(df_t4, df_ir)
        print(f"  Matched: {len(df_matches)} clonotypes")

        # Metrics
        metrics = calculate_metrics(df_t4, df_ir, df_matches)

        # Save results
        output_dir = f'/data1/wetlab/TRUST4_benchmark/results/benchmark/{sample}/published_eval'
        os.makedirs(output_dir, exist_ok=True)

        # Save parsed data
        t4_out = df_t4[['CDR3aa', 'v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt_trimmed', 'count']].copy()
        t4_out.columns = ['CDR3aa', 'V', 'J', 'C_collapsed', 'CDR3nt_trimmed', 'read_count']
        t4_out.to_csv(f'{output_dir}/trust4_parsed.tsv', sep='\t', index=False)

        ir_out = df_ir[['CDR3(pep)', 'v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt', 'copy']].copy()
        ir_out.columns = ['CDR3aa', 'V', 'J', 'C_collapsed', 'CDR3nt', 'copy_count']
        ir_out.to_csv(f'{output_dir}/irep_parsed.tsv', sep='\t', index=False)

        # Save matches
        if len(df_matches) > 0:
            df_matches.to_csv(f'{output_dir}/matches.tsv', sep='\t', index=False)

        # Save metrics
        metrics_df = pd.DataFrame([metrics])
        metrics_df.to_csv(f'{output_dir}/metrics.tsv', sep='\t', index=False)

        print(f"Results saved to {output_dir}/")

        return sample, metrics

    except Exception as e:
        print(f"ERROR: {e}")
        return sample, None

def main():
    samples = ['FZ-116', 'FZ-20', 'FZ-83', 'FZ-94', 'FZ-97', 'FZ-122']

    print("="*80)
    print("TRUST4 Published Methodology - All Samples")
    print("="*80)

    results = {}
    for sample in samples:
        sample_name, metrics = benchmark_sample(sample)
        if metrics is not None:
            results[sample_name] = metrics

    # Create summary TSV
    print(f"\n{'='*80}")
    print("Creating summary")
    print(f"{'='*80}")

    summary_rows = []
    for sample in samples:
        if sample in results:
            m = results[sample]
            summary_rows.append({
                'sample': sample,
                'trust4_clonotypes': m['n_trust4'],
                'irep_clonotypes': m['n_irep'],
                'matches': m['n_matches'],
                'precision': m['precision'],
                'sensitivity': m['sensitivity'],
                'pearson_r': m['pearson_r'],
                'd_agreement': m['d_agreement'],
                'isotype_exact_agreement': m['iso_exact_agreement'],
                'isotype_collapsed_agreement': m['iso_collapsed_agreement'],
            })

    summary_df = pd.DataFrame(summary_rows)

    # Save TSV
    summary_tsv_path = '/data1/wetlab/TRUST4_benchmark/results/benchmark/all_samples_published_summary.tsv'
    summary_df.to_csv(summary_tsv_path, sep='\t', index=False)
    print(f"\nSummary TSV saved: {summary_tsv_path}")

    # Calculate statistics
    precision_mean = summary_df['precision'].mean()
    precision_std = summary_df['precision'].std()
    sensitivity_mean = summary_df['sensitivity'].mean()
    sensitivity_std = summary_df['sensitivity'].std()
    pearson_mean = summary_df['pearson_r'].mean()
    pearson_std = summary_df['pearson_r'].std()

    # Create summary MD
    md_content = f"""# Published TRUST4 Benchmark - All Samples Summary

**Date**: 2026-09-18
**Methodology**: Published evaluation from TRUST4_manuscript_evaluation/bulk/bcrval.py
**Samples**: {len(summary_df)} samples (FZ-116, FZ-20, FZ-83, FZ-94, FZ-97, FZ-122)
**Processing**: Individual sample processing, no pooling

---

## Results by Sample

| Sample | TRUST4 Clonotypes | iRepertoire Clonotypes | Matches | Precision | Sensitivity | Pearson r | D Agreement | Isotype Exact | Isotype Collapsed |
|--------|-------------------|------------------------|---------|-----------|-------------|-----------|-------------|---------------|-------------------|
"""

    for _, row in summary_df.iterrows():
        md_content += f"| {row['sample']} | {row['trust4_clonotypes']:,} | {row['irep_clonotypes']:,} | {row['matches']:,} | {row['precision']:.4f} | {row['sensitivity']:.4f} | {row['pearson_r']:.4f} | {row['d_agreement']:.4f} | {row['isotype_exact_agreement']:.4f} | {row['isotype_collapsed_agreement']:.4f} |\n"

    md_content += f"""
---

## Summary Statistics

| Metric | Mean | Std Dev | Min | Max |
|--------|------|---------|-----|-----|
| Precision | {precision_mean:.4f} | {precision_std:.4f} | {summary_df['precision'].min():.4f} | {summary_df['precision'].max():.4f} |
| Sensitivity | {sensitivity_mean:.4f} | {sensitivity_std:.4f} | {summary_df['sensitivity'].min():.4f} | {summary_df['sensitivity'].max():.4f} |
| Pearson r | {pearson_mean:.4f} | {pearson_std:.4f} | {summary_df['pearson_r'].min():.4f} | {summary_df['pearson_r'].max():.4f} |
| D Agreement | {summary_df['d_agreement'].mean():.4f} | {summary_df['d_agreement'].std():.4f} | {summary_df['d_agreement'].min():.4f} | {summary_df['d_agreement'].max():.4f} |
| Isotype Exact | {summary_df['isotype_exact_agreement'].mean():.4f} | {summary_df['isotype_exact_agreement'].std():.4f} | {summary_df['isotype_exact_agreement'].min():.4f} | {summary_df['isotype_exact_agreement'].max():.4f} |

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

For each sample, individual results in `results/benchmark/{{sample}}/published_eval/`:
- `trust4_parsed.tsv` — TRUST4 clonotypes (trimmed)
- `irep_parsed.tsv` — iRepertoire clonotypes
- `matches.tsv` — Matched clonotype pairs
- `metrics.tsv` — Precision/sensitivity/correlation

---

**Status**: ✓ Baseline established. Ready for next analysis step.
"""

    # Save MD
    summary_md_path = '/data1/wetlab/TRUST4_benchmark/results/benchmark/all_samples_published_summary.md'
    with open(summary_md_path, 'w') as f:
        f.write(md_content)

    print(f"Summary MD saved: {summary_md_path}")

    # Print to console
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    print(summary_df.to_string(index=False))

    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nPrecision:    {precision_mean:.4f} ± {precision_std:.4f}")
    print(f"Sensitivity:  {sensitivity_mean:.4f} ± {sensitivity_std:.4f}")
    print(f"Pearson r:    {pearson_mean:.4f} ± {pearson_std:.4f}")

    print(f"\n{'='*80}")
    print("ALL SAMPLES COMPLETE")
    print(f"{'='*80}")

    return 0

if __name__ == '__main__':
    sys.exit(main())
