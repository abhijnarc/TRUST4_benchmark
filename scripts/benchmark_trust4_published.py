#!/usr/bin/env python3
"""
TRUST4 vs iRepertoire Benchmarking - PUBLISHED METHODOLOGY

Reproduces the evaluation logic from:
https://github.com/liulab-dfci/TRUST4_manuscript_evaluation/blob/master/bulk/bcrval.py

Key differences from custom approach:
1. CDR3nt trimming: Remove first 3 and last 3 nucleotides from TRUST4
2. Primary key: V + J + C + CDR3nt (without D, without CDR3aa)
3. Isotype collapsing: IGHA1/IGHA2 -> IGHA, IGHG1-4 -> IGHG
4. Duplicate handling: Keep MAXIMUM abundance, not summing
5. D gene: Excluded from primary benchmark

Output:
- Precision and sensitivity (published-style)
- Secondary diagnostic metrics (CDR3aa, D, exact isotypes, etc.)
"""

import sys
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import gzip
import os

def remove_allele(gene_name):
    """Remove allele suffix from gene name, e.g., IGHV3-23*01 -> IGHV3-23"""
    if pd.isna(gene_name):
        return None
    gene = str(gene_name)
    if '*' in gene:
        gene = gene.split('*')[0]
    return gene

def collapse_isotype(c_gene):
    """Collapse isotype subclasses: IGHA1/IGHA2 -> IGHA, IGHG1-4 -> IGHG"""
    if pd.isna(c_gene):
        return None
    c = str(c_gene)
    # Remove trailing numbers for subclasses
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
    """Trim 3 nucleotides from start and 3 from end (published methodology)"""
    if pd.isna(cdr3_nt):
        return None
    seq = str(cdr3_nt)
    if len(seq) < 6:
        return None  # Cannot trim if too short
    return seq[3:-3]

def load_trust4():
    """Load and parse TRUST4 report"""
    print("Loading TRUST4 report...")

    file_path = '/data1/wetlab/TRUST4_benchmark/results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv'
    df = pd.read_csv(file_path, sep='\t')
    df.columns = [col.lstrip('#') for col in df.columns]

    print(f"  Loaded {len(df)} records")

    # Filter to IGH
    df = df[df['V'].str.startswith('IGH', na=False)].copy()
    print(f"  After IGH filter: {len(df)} records")

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
    print(f"  After CDR3 trim validity check: {len(df)} records")

    return df

def load_irep():
    """Load and parse iRepertoire reference"""
    print("Loading iRepertoire reference...")

    file_path = '/data1/wetlab/TRUST4_benchmark/reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz'
    df = pd.read_csv(file_path, compression='gzip')

    print(f"  Loaded {len(df)} records")

    # Parse and normalize
    df['v_gene'] = df['V'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['j_gene'] = df['J'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['c_gene_original'] = df['C'].apply(lambda x: remove_allele(str(x).lstrip('h')))
    df['c_gene_collapsed'] = df['c_gene_original'].apply(collapse_isotype)
    df['cdr3nt'] = df['CDR3(nuc)'].apply(lambda x: str(x).upper() if pd.notna(x) else None)
    df['d_gene'] = df['D'].apply(lambda x: remove_allele(str(x).lstrip('h')))

    return df

def build_key(row, trim=True):
    """Build primary matching key: V + J + C + CDR3nt"""
    v = row.get('v_gene')
    j = row.get('j_gene')
    c = row.get('c_gene_collapsed')

    if trim:
        cdr3 = row.get('cdr3nt_trimmed')
    else:
        cdr3 = row.get('cdr3nt')

    if pd.isna(v) or pd.isna(j) or pd.isna(c) or pd.isna(cdr3):
        return None

    return (v, j, c, cdr3)

def coalesce_duplicates(df, is_trust4=True):
    """
    Coalesce identical clonotypes, keeping maximum abundance.
    Published methodology: keep the record with highest count/copy.
    """
    if is_trust4:
        abundance_col = 'count'
        key_cols = ['v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt_trimmed']
    else:
        abundance_col = 'copy'
        key_cols = ['v_gene', 'j_gene', 'c_gene_collapsed', 'cdr3nt']

    print(f"\nCoalescing duplicates (keeping max abundance)...")
    print(f"  Before: {len(df)} records")

    # Group by key and keep record with maximum abundance
    df_coalesced = df.loc[df.groupby(key_cols)[abundance_col].idxmax()].copy()
    df_coalesced = df_coalesced.reset_index(drop=True)

    print(f"  After: {len(df_coalesced)} unique clonotypes")

    return df_coalesced

def show_examples(df_t4, df_ir):
    """Show 5 examples of CDR3 trimming"""
    print("\n" + "="*80)
    print("CDR3 TRIMMING VALIDATION (5 examples)")
    print("="*80)

    for i in range(min(5, len(df_t4))):
        row = df_t4.iloc[i]
        print(f"\nExample {i+1}:")
        print(f"  Original CDR3nt: {row['cdr3nt_original']}")
        print(f"  Trimmed CDR3nt:  {row['cdr3nt_trimmed']}")
        print(f"  CDR3aa (TRUST4): {row['CDR3aa']}")
        print(f"  V: {row['v_gene']}, J: {row['j_gene']}, C: {row['c_gene_collapsed']}")

def match_clonotypes(df_t4, df_ir):
    """Match TRUST4 and iRepertoire clonotypes using published key"""
    print("\n" + "="*80)
    print("MATCHING CLONOTYPES")
    print("="*80)

    # Build iRepertoire index
    ir_index = {}
    for idx, row in df_ir.iterrows():
        key = (row['v_gene'], row['j_gene'], row['c_gene_collapsed'], row['cdr3nt'])
        if pd.notna(key[0]) and pd.notna(key[3]):  # Valid key
            ir_index[key] = row

    print(f"\niRepertoire index size: {len(ir_index)} unique clonotypes")

    # Match TRUST4 to iRepertoire
    matches = []
    for idx, t4_row in df_t4.iterrows():
        key = (t4_row['v_gene'], t4_row['j_gene'], t4_row['c_gene_collapsed'], t4_row['cdr3nt_trimmed'])

        if key in ir_index:
            ir_row = ir_index[key]

            # Collect additional diagnostic info
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

    df_matches = pd.DataFrame(matches)

    print(f"Matches found: {len(df_matches)}")
    print(f"TRUST4 clonotypes unmatched: {len(df_t4) - len(df_matches)}")
    print(f"iRepertoire clonotypes unmatched: {len(df_ir) - len(df_matches)}")

    return df_matches

def calculate_metrics(df_t4, df_ir, df_matches):
    """Calculate published-style precision and sensitivity"""
    print("\n" + "="*80)
    print("BENCHMARK METRICS (Published Methodology)")
    print("="*80)

    n_trust4 = len(df_t4)
    n_irep = len(df_ir)
    n_matches = len(df_matches)

    # Published precision: intersection / TRUST4_total
    precision = n_matches / n_trust4 if n_trust4 > 0 else 0

    # Published sensitivity: intersection / iRepertoire_total
    sensitivity = n_matches / n_irep if n_irep > 0 else 0

    print(f"\nData counts:")
    print(f"  TRUST4 clonotypes: {n_trust4}")
    print(f"  iRepertoire clonotypes: {n_irep}")
    print(f"  Intersection: {n_matches}")

    print(f"\nPublished metrics:")
    print(f"  Precision (TP / TRUST4): {precision:.4f}")
    print(f"  Sensitivity (TP / iRep): {sensitivity:.4f}")

    return {
        'precision': precision,
        'sensitivity': sensitivity,
        'n_trust4': n_trust4,
        'n_irep': n_irep,
        'n_matches': n_matches,
    }

def calculate_diagnostic_metrics(df_matches):
    """Calculate diagnostic metrics beyond published benchmark"""
    print(f"\n" + "="*80)
    print("DIAGNOSTIC METRICS")
    print("="*80)

    if len(df_matches) == 0:
        print("  No matches - diagnostic metrics N/A")
        return {}

    # CDR3aa agreement
    cdr3aa_agree = (df_matches['cdr3aa_trust4'] == df_matches['cdr3aa_irep']).sum()
    print(f"\n  CDR3aa exact agreement: {cdr3aa_agree}/{len(df_matches)} ({100*cdr3aa_agree/len(df_matches):.1f}%)")

    # CDR3nt before trimming agreement
    cdr3nt_orig_agree = (df_matches['cdr3nt_original_trust4'] == df_matches['cdr3nt_irep']).sum()
    print(f"  CDR3nt (TRUST4 original) exact agreement: {cdr3nt_orig_agree}/{len(df_matches)} ({100*cdr3nt_orig_agree/len(df_matches):.1f}%)")

    # CDR3nt after trimming agreement
    cdr3nt_trim_agree = (df_matches['cdr3nt_trimmed'] == df_matches['cdr3nt_irep']).sum()
    print(f"  CDR3nt (TRUST4 trimmed) exact agreement: {cdr3nt_trim_agree}/{len(df_matches)} ({100*cdr3nt_trim_agree/len(df_matches):.1f}%)")

    # V, J, C agreement (should all be 100% since that's the matching key)
    v_agree = (df_matches['v_gene'] == df_matches['v_gene']).sum()
    print(f"  V gene agreement: {len(df_matches)}/{len(df_matches)} (100.0% - part of key)")
    print(f"  J gene agreement: {len(df_matches)}/{len(df_matches)} (100.0% - part of key)")
    print(f"  C gene (collapsed) agreement: {len(df_matches)}/{len(df_matches)} (100.0% - part of key)")

    # D gene agreement
    d_agree = (df_matches['d_gene_trust4'] == df_matches['d_gene_irep']).sum()
    print(f"  D gene agreement: {d_agree}/{len(df_matches)} ({100*d_agree/len(df_matches):.1f}%)")

    # Exact isotype agreement (original, not collapsed)
    iso_agree = (df_matches['c_gene_trust4_original'] == df_matches['c_gene_irep_original']).sum()
    print(f"  C gene (original isotype) exact agreement: {iso_agree}/{len(df_matches)} ({100*iso_agree/len(df_matches):.1f}%)")

    return {
        'cdr3aa_agreement': cdr3aa_agree,
        'd_gene_agreement': d_agree,
        'iso_exact_agreement': iso_agree,
    }

def calculate_abundance_correlation(df_matches):
    """Calculate Pearson correlation for abundance"""
    print(f"\nAbundance correlation:")

    if len(df_matches) < 2:
        print("  Insufficient matched pairs")
        return None

    # Remove NaN
    df_corr = df_matches[['trust4_count', 'irep_copy']].dropna()

    if len(df_corr) < 2:
        print("  Insufficient data")
        return None

    r, p = pearsonr(df_corr['trust4_count'], df_corr['irep_copy'])
    print(f"  Pearson r: {r:.4f} (p={p:.4e})")

    return {'pearson_r': r, 'pearson_p': p}

def save_results(df_t4, df_ir, df_matches, metrics, output_dir):
    """Save all output files"""
    print(f"\nSaving results to {output_dir}/...")

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

    print(f"  Results saved")

def main():
    print("="*80)
    print("TRUST4 Benchmarking - Published Methodology (FZ-116)")
    print("="*80)

    # Load
    df_t4 = load_trust4()
    df_ir = load_irep()

    # Show examples before processing
    show_examples(df_t4, df_ir)

    # Coalesce
    df_t4 = coalesce_duplicates(df_t4, is_trust4=True)
    df_ir = coalesce_duplicates(df_ir, is_trust4=False)

    # Match
    df_matches = match_clonotypes(df_t4, df_ir)

    # Metrics
    metrics = calculate_metrics(df_t4, df_ir, df_matches)
    diag_metrics = calculate_diagnostic_metrics(df_matches)
    corr_metrics = calculate_abundance_correlation(df_matches)

    # Save
    output_dir = '/data1/wetlab/TRUST4_benchmark/results/benchmark/FZ-116/published_eval'
    save_results(df_t4, df_ir, df_matches, metrics, output_dir)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\nCDR3 normalization status: {'OVERLAPS FOUND ✓' if metrics['n_matches'] > 0 else 'NO OVERLAPS'}")
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Sensitivity: {metrics['sensitivity']:.4f}")
    print(f"Matching clonotypes: {metrics['n_matches']}")
    if corr_metrics:
        print(f"Abundance correlation (Pearson r): {corr_metrics['pearson_r']:.4f}")

    return 0

if __name__ == '__main__':
    sys.exit(main())
