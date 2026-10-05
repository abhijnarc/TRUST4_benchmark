#!/usr/bin/env python3
"""Evaluate existing Graph-TRUST4-v4 contig annotations against iRepertoire."""

import argparse
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ANNOTATION = HERE / "graph_v4_raw_airr.tsv"
FASTA = ROOT / "results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa"
REFERENCE = ROOT / "reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz"
TRUST4_REPORT = ROOT / "results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv"
METHOD_SCRIPT = ROOT / "scripts/benchmark_trust4_published.py"

EXPECTED = {
    "unique": 11233,
    "reference": 116632,
    "matches": 3815,
    "precision": 0.3396243212,
    "sensitivity": 0.03270971946,
    "pearson_r": 0.6604407834,
}


def load_method():
    spec = importlib.util.spec_from_file_location("benchmark_trust4_published", METHOD_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def value(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return None if text in ("", ".", "nan", "None") else text


def write_tsv(path, rows, columns=None):
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_csv(path, sep="\t", index=False, na_rep="NOT AVAILABLE")


def fasta_lengths(path):
    lengths, current_id, current_length = {}, None, 0
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_length
                current_id, current_length = line[1:].split()[0], 0
            else:
                current_length += len(line)
    if current_id is not None:
        lengths[current_id] = current_length
    return lengths


def nx(lengths, fraction=0.5):
    lengths = sorted((int(x) for x in lengths), reverse=True)
    if not lengths:
        return None
    target = sum(lengths) * fraction
    total = 0
    for length in lengths:
        total += length
        if total >= target:
            return length
    return lengths[-1]


def baseline_check(method):
    trust_raw = method.load_trust4()
    trust = method.coalesce_duplicates(trust_raw, is_trust4=True)
    reference = method.coalesce_duplicates(method.load_irep(), is_trust4=False)
    matches = method.match_clonotypes(trust, reference)
    metrics = method.calculate_metrics(trust, reference, matches)
    correlation = method.calculate_abundance_correlation(matches)
    observed = {
        "unique": len(trust),
        "reference": len(reference),
        "matches": metrics["n_matches"],
        "precision": metrics["precision"],
        "sensitivity": metrics["sensitivity"],
        "pearson_r": correlation["pearson_r"] if correlation else None,
    }
    valid = (
        observed["unique"] == EXPECTED["unique"]
        and observed["reference"] == EXPECTED["reference"]
        and observed["matches"] == EXPECTED["matches"]
        and all(abs(observed[k] - EXPECTED[k]) <= 1e-9 for k in ["precision", "sensitivity", "pearson_r"])
    )
    return valid, observed, trust_raw, trust, reference, matches


def normalize_annotations(method):
    raw = pd.read_csv(ANNOTATION, sep="\t", dtype=str, keep_default_na=False)
    lengths = fasta_lengths(FASTA)
    raw.rename(columns={column: "raw_" + column for column in raw.columns}, inplace=True)
    ann = pd.DataFrame(index=raw.index)
    ann["contig_id"] = raw["raw_sequence_id"]
    ann["contig_length"] = ann["contig_id"].map(lengths)
    ann["abundance"] = "NOT AVAILABLE"
    ann["V"] = raw["raw_v_call"].map(value)
    ann["D"] = raw["raw_d_call"].map(value)
    ann["J"] = raw["raw_j_call"].map(value)
    ann["C"] = raw["raw_c_call"].map(value)
    ann["CDR3_nt"] = raw["raw_junction"].str.upper().map(value)
    ann["CDR3_aa"] = raw["raw_junction_aa"].str.upper().map(value)
    ann["productive"] = raw["raw_productive"].map(value)
    ann["complete_vdj"] = raw["raw_complete_vdj"].map(value)
    ann["normalized_V"] = ann["V"].map(method.remove_allele)
    ann["normalized_D"] = ann["D"].map(method.remove_allele)
    ann["normalized_J"] = ann["J"].map(method.remove_allele)
    ann["normalized_C_original"] = ann["C"].map(method.remove_allele)
    ann["normalized_C_collapsed"] = ann["normalized_C_original"].map(method.collapse_isotype)
    ann["CDR3_nt_trimmed"] = ann["CDR3_nt"].map(method.trim_cdr3)
    ann["chain"] = ann.apply(lambda row: next((str(x)[:3] for x in [row["normalized_V"], row["normalized_J"], row["normalized_C_original"]] if value(x) and str(x).startswith(("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"))), "NOT AVAILABLE"), axis=1)

    def annotation_status(row):
        normalized_v = value(row["normalized_V"])
        if normalized_v and normalized_v.startswith("IGH") and row["productive"] == "T":
            return "productive_IGH"
        missing = [key for key, field in [("V", "normalized_V"), ("J", "normalized_J"), ("C", "normalized_C_original"), ("CDR3", "CDR3_nt")] if not value(row[field])]
        if missing:
            return "missing_" + "+".join(missing)
        return "annotated_nonproductive_or_partial"

    ann["annotation_status"] = ann.apply(annotation_status, axis=1)
    ann = pd.concat([ann, raw], axis=1)
    ann.to_csv(HERE / "graph_v4_annotated.tsv", sep="\t", index=False, na_rep="NOT AVAILABLE")
    normalized_columns = ["contig_id", "contig_length", "abundance", "chain", "V", "D", "J", "C", "CDR3_nt", "CDR3_aa",
                          "annotation_status", "productive", "complete_vdj", "normalized_V", "normalized_D", "normalized_J",
                          "normalized_C_original", "normalized_C_collapsed", "CDR3_nt_trimmed", "raw_v_call", "raw_d_call",
                          "raw_j_call", "raw_c_call", "raw_junction", "raw_junction_aa", "raw_v_cigar", "raw_d_cigar",
                          "raw_j_cigar", "raw_c_cigar", "raw_v_identity", "raw_j_identity"]
    ann[normalized_columns].to_csv(HERE / "graph_v4_normalized.tsv", sep="\t", index=False, na_rep="NOT AVAILABLE")
    return ann


def yield_and_length_tables(ann):
    total = len(ann)
    igh = ann["normalized_V"].fillna("").str.startswith("IGH")
    has_v = ann["normalized_V"].map(value).notna()
    has_j = ann["normalized_J"].map(value).notna()
    has_c = ann["normalized_C_original"].map(value).notna()
    has_cdr3 = ann["CDR3_nt"].map(value).notna()
    has_aa = ann["CDR3_aa"].map(value).notna()
    vj, vjcdr3, vjc, vjccdr3 = has_v & has_j, has_v & has_j & has_cdr3, has_v & has_j & has_c, has_v & has_j & has_c & has_cdr3
    productive = igh & (ann["productive"] == "T")
    flags = {
        "total_contigs": pd.Series(True, index=ann.index), "IGH_contigs": igh, "V_assigned": has_v, "J_assigned": has_j,
        "C_assigned": has_c, "CDR3_assigned": has_cdr3, "CDR3_amino_acid_assigned": has_aa,
        "V+J": vj, "V+J+CDR3": vjcdr3, "V+J+C": vjc, "V+J+C+CDR3": vjccdr3,
        "productive_IGH": productive, "missing_V": ~has_v, "missing_J": ~has_j, "missing_C": ~has_c, "missing_CDR3": ~has_cdr3,
    }
    yield_rows = []
    for name, mask in flags.items():
        all_count = int(mask.sum())
        igh_count = int((mask & igh).sum())
        yield_rows.append({"metric":name,"all_contigs_count":all_count,"all_contigs_percent":100*all_count/total if total else None,
                           "IGH_count":igh_count,"IGH_percent":100*igh_count/int(igh.sum()) if igh.any() else None})
    write_tsv(HERE / "annotation_yield.tsv", yield_rows)

    bins = [(150,199,"150-199"),(200,249,"200-249"),(250,299,"250-299"),(300,399,"300-399"),
            (400,499,"400-499"),(500,749,"500-749"),(750,999,"750-999"),(1000,1499,"1000-1499"),(1500,math.inf,">=1500")]
    length_rows = []
    for low, high, label in bins:
        in_bin = ann["contig_length"].between(low, high, inclusive="both") if math.isfinite(high) else ann["contig_length"] >= low
        denom, igh_denom = int(in_bin.sum()), int((in_bin & igh).sum())
        row = {"length_bin":label,"total_contigs":denom,"IGH_contigs":igh_denom}
        for label2, mask in [("V_assigned",has_v),("J_assigned",has_j),("CDR3_assigned",has_cdr3),
                             ("V_plus_J",vj),("V_plus_J_plus_CDR3",vjcdr3),("V_plus_J_plus_C_plus_CDR3",vjccdr3),
                             ("productive_IGH",productive)]:
            count=int((in_bin & mask).sum())
            igh_count=int((in_bin & mask & igh).sum())
            row[label2]=count
            row[label2+"_percent_bin"]=100*count/denom if denom else None
            row[label2+"_IGH"]=igh_count
            row[label2+"_percent_IGH_bin"]=100*igh_count/igh_denom if igh_denom else None
        length_rows.append(row)
    write_tsv(HERE / "length_vs_annotation.tsv", length_rows)

    completeness = [("all contigs",pd.Series(True,index=ann.index)),("IGH",igh),("V assigned",has_v & igh),
                    ("V+J",vj & igh),("V+J+CDR3",vjcdr3 & igh),("V+J+C+CDR3",vjccdr3 & igh),("productive IGH",productive)]
    complete_rows=[]
    for name, mask in completeness:
        vals=ann.loc[mask,"contig_length"].dropna().astype(int).tolist()
        complete_rows.append({"annotation_category":name,"count":int(mask.sum()),"percentage_all_contigs":100*int(mask.sum())/total if total else None,
                              "median_length":float(np.median(vals)) if vals else None,"mean_length":float(np.mean(vals)) if vals else None,"N50":nx(vals)})
    write_tsv(HERE / "completeness_summary.tsv",complete_rows)

    length_categories=[("V_assigned",has_v & igh),("J_assigned",has_j & igh),("CDR3_containing",has_cdr3 & igh),
                       ("V+J+CDR3",vjcdr3 & igh),("V+J+C+CDR3",vjccdr3 & igh)]
    length_summary={name:{"median":float(ann.loc[mask,"contig_length"].median()),"mean":float(ann.loc[mask,"contig_length"].mean()),"count":int(mask.sum())} for name,mask in length_categories if mask.any()}
    return yield_rows,length_rows,complete_rows,length_summary


def benchmark_tables(method, ann, reference, baseline_matches):
    # Published rules: uppercase CDR3nt, trim 3 nt from each end, remove V/J/D alleles,
    # collapse C subclasses, and retain only complete keys for primary evaluation.
    igh = ann[ann["normalized_V"].fillna("").str.startswith("IGH")].copy()
    candidate_rows={}
    for row in igh.to_dict("records"):
        key=(value(row["normalized_V"]),value(row["normalized_J"]),value(row["normalized_C_collapsed"]),value(row["CDR3_nt_trimmed"]))
        if all(key):
            candidate_rows.setdefault(key,row)
    ref_rows={}
    for row in reference.to_dict("records"):
        key=(value(row["v_gene"]),value(row["j_gene"]),value(row["c_gene_collapsed"]),value(row["cdr3nt"]))
        if all(key): ref_rows[key]=row
    candidate_keys=set(candidate_rows)
    reference_keys=set(ref_rows)
    matched_keys=candidate_keys & reference_keys
    n_candidate,n_reference,n_matches=len(candidate_keys),len(reference_keys),len(matched_keys)
    precision=n_matches/n_candidate if n_candidate else None
    sensitivity=n_matches/n_reference if n_reference else None
    rows=[]
    for key in sorted(matched_keys):
        v4=candidate_rows[key]; ir=ref_rows[key]
        rows.append({"v4_contig_id":v4["contig_id"],"v4_abundance":"NOT AVAILABLE","V":v4["normalized_V"],
                     "D":v4["normalized_D"],"J":v4["normalized_J"],"C":v4["normalized_C_original"],
                     "CDR3_nt":v4["CDR3_nt"],"CDR3_aa":v4["CDR3_aa"],"iRep_abundance":ir["copy"],"match_status":"matched",
                     "iRep_V":ir["v_gene"],"iRep_D":ir["d_gene"],"iRep_J":ir["j_gene"],"iRep_C":ir["c_gene_original"],"iRep_CDR3_nt":ir["cdr3nt"]})
    write_tsv(HERE/"graph_v4_vs_iRep_matches.tsv",rows,columns=["v4_contig_id","v4_abundance","V","D","J","C","CDR3_nt","CDR3_aa",
              "iRep_abundance","match_status","iRep_V","iRep_D","iRep_J","iRep_C","iRep_CDR3_nt"])

    metric_rows=[("v4_unique_clonotypes",n_candidate,"Unique complete V+J+collapsed-C+trimmed-CDR3nt keys; duplicate contigs collapsed."),
                 ("iRepertoire_reference_clonotypes",n_reference,"Published-method coalesced keys."),("matched_clonotypes",n_matches,"Intersection of primary keys."),
                 ("v4_only_clonotypes",n_candidate-n_matches,"Candidate keys absent from reference."),("iRepertoire_only_clonotypes",n_reference-n_matches,"Reference keys absent from v4."),
                 ("precision",precision,"matches / v4 unique clonotypes."),("sensitivity",sensitivity,"matches / reference clonotypes."),
                 ("pearson_r","NOT AVAILABLE","Per-contig/path abundance is absent from v4 artifacts; contig length is not substituted."),
                 ("abundance_definition","NOT AVAILABLE","Established graph evaluator sums bundle abundance across path nodes; v4 path node IDs are not stored."),
                 ("duplicate_key_policy","first contig retained for traceability","All v4 abundance values are unavailable; identical primary keys count once.")]
    write_tsv(HERE/"graph_v4_benchmark_metrics.tsv",[{"metric":a,"value":b,"notes":c} for a,b,c in metric_rows])

    criteria=[("A_CDR3_only",(3,)),("B_CDR3+V",(3,0)),("C_CDR3+J",(3,1)),("D_CDR3+V+J",(3,0,1)),("E_CDR3+V+J+C",(3,0,1,2))]
    breakdown=[]
    partial_candidate_keys=[(value(r["normalized_V"]),value(r["normalized_J"]),value(r["normalized_C_collapsed"]),value(r["CDR3_nt_trimmed"])) for r in igh.to_dict("records")]
    ref_fields=["v_gene","j_gene","c_gene_collapsed","cdr3nt"]
    for name,indices in criteria:
        cset=set(); rset=set()
        for k in partial_candidate_keys:
            reduced=tuple(k[i] for i in indices)
            if all(reduced): cset.add(reduced)
        for k in reference_keys:
            reduced=tuple(k[i] for i in indices); rset.add(reduced)
        nm=len(cset&rset)
        breakdown.append({"criterion":name,"candidate_count":len(cset),"reference_count":len(rset),"matches":nm,
                          "precision":nm/len(cset) if cset else None,"sensitivity":nm/len(rset) if rset else None})
    write_tsv(HERE/"matching_breakdown.tsv",breakdown)

    gene_rows=[]
    if matched_keys:
        d_pairs=[]; exact_c_pairs=[]; collapsed_c_pairs=[]
        for key in matched_keys:
            v4=candidate_rows[key]; ir=ref_rows[key]
            d_pairs.append((value(v4["normalized_D"]),value(ir["d_gene"])))
            exact_c_pairs.append((value(v4["normalized_C_original"]),value(ir["c_gene_original"])))
            collapsed_c_pairs.append((value(v4["normalized_C_collapsed"]),value(ir["c_gene_collapsed"])))
        for name,pairs in [("D_agreement",d_pairs),("exact_C_isotype_agreement",exact_c_pairs),("collapsed_C_isotype_agreement",collapsed_c_pairs)]:
            equal=sum(a==b for a,b in pairs)
            both=sum(a is not None and b is not None for a,b in pairs)
            gene_rows.append({"metric":name,"matched_clonotypes":len(pairs),"agreements":equal,"agreement_percent":100*equal/len(pairs),"both_assigned_pairs":both})
    else:
        for name in ["D_agreement","exact_C_isotype_agreement","collapsed_C_isotype_agreement"]:
            gene_rows.append({"metric":name,"matched_clonotypes":0,"agreements":0,"agreement_percent":"NOT AVAILABLE","both_assigned_pairs":0})
    write_tsv(HERE/"gene_agreement.tsv",gene_rows)

    ir_matched=np.asarray([number for number in [value(ref_rows[k]["copy"]) for k in matched_keys] if number is not None],dtype=float)
    abundance_rows=[{"metric":"matched_clonotypes","value":n_matches,"notes":"Primary matched keys."},
                    {"metric":"pairs_used_for_pearson","value":0,"notes":"V4 path abundance unavailable."},
                    {"metric":"pearson_r","value":"NOT AVAILABLE","notes":"No v4 abundance; do not substitute contig length."},
                    {"metric":"pearson_p_value","value":"NOT AVAILABLE","notes":"No v4 abundance."},
                    {"metric":"v4_abundance_minimum","value":"NOT AVAILABLE","notes":"No path node membership in v4 path tables."},
                    {"metric":"v4_abundance_median","value":"NOT AVAILABLE","notes":"No path node membership in v4 path tables."},
                    {"metric":"v4_abundance_mean","value":"NOT AVAILABLE","notes":"No path node membership in v4 path tables."},
                    {"metric":"v4_abundance_maximum","value":"NOT AVAILABLE","notes":"No path node membership in v4 path tables."}]
    for name,stat in [("iRep_abundance_minimum",np.min),("iRep_abundance_median",np.median),("iRep_abundance_mean",np.mean),("iRep_abundance_maximum",np.max)]:
        abundance_rows.append({"metric":name,"value":float(stat(ir_matched)) if len(ir_matched) else "NOT AVAILABLE","notes":"iRepertoire copy among primary matched clonotypes."})
    write_tsv(HERE/"abundance_correlation.tsv",abundance_rows)

    bins=[("copy = 1",1,1),("copy = 2-5",2,5),("copy = 6-10",6,10),("copy = 11-50",11,50),("copy = 51-100",51,100),("copy >100",101,math.inf)]
    stratified=[]
    matched_ref_keys={k for k in matched_keys}
    for label,low,high in bins:
        ref_bin={k for k,r in ref_rows.items() if value(r.get("copy")) is not None and low<=float(r["copy"])<=high}
        match_bin=len(ref_bin & matched_ref_keys)
        stratified.append({"reference_copy_bin":label,"reference_clonotypes":len(ref_bin),"v4_clonotypes":"NOT AVAILABLE",
                           "matches":match_bin,"precision":"NOT AVAILABLE","sensitivity":match_bin/len(ref_bin) if ref_bin else None,
                           "stratification_basis":"iRepertoire reference copy; v4 abundance unavailable, so v4 bucket count/precision cannot be calculated"})
    write_tsv(HERE/"abundance_stratified_benchmark.tsv",stratified)

    return {"unique":n_candidate,"reference":n_reference,"matches":n_matches,"precision":precision,"sensitivity":sensitivity,
            "pearson_r":None,"d_agreement":gene_rows[0]["agreement_percent"],"isotype_exact":gene_rows[1]["agreement_percent"],
            "isotype_collapsed":gene_rows[2]["agreement_percent"],"v4_only":n_candidate-n_matches,"irep_only":n_reference-n_matches,
            "breakdown":breakdown,"gene_rows":gene_rows,"stratified":stratified,"matched_rows":rows}


def comparison_table(method, trust_raw_report, trust_unique, trust_matches, trust_corr, trust_gene, v4, yield_rows):
    report=trust_raw_report.copy()
    report.columns=[c.lstrip("#") for c in report.columns]
    igh=report[report["V"].astype(str).str.startswith("IGH",na=False)].copy()
    def present(series): return series.astype(str).str.strip().ne("") & series.astype(str).str.strip().ne(".")
    def count(mask): return int(mask.sum())
    v=present(igh["V"]); j=present(igh["J"]); c=present(igh["C"]); cdr=present(igh["CDR3nt"])
    trust={"total_contigs":len(report),"IGH_contigs":len(igh),"V_assigned":count(v),"J_assigned":count(j),"C_assigned":count(c),
           "CDR3_assigned":count(cdr),"productive_IGH":"NOT AVAILABLE","V+J+CDR3":count(v&j&cdr),"V+J+C+CDR3":count(v&j&c&cdr),
           "unique_clonotypes":trust_unique,"reference_clonotypes":len(method.coalesce_duplicates(method.load_irep(),is_trust4=False)),
           "matched_clonotypes":trust_matches,"precision":trust_corr["precision"],"sensitivity":trust_corr["sensitivity"],
           "Pearson_r":trust_corr["pearson_r"],"D_agreement":trust_gene["D_agreement"],
           "exact_isotype_agreement":trust_gene["iso_exact_agreement"],"collapsed_isotype_agreement":100.0}
    yield_map={r["metric"]:r for r in yield_rows}
    v4={"total_contigs":yield_map["total_contigs"]["all_contigs_count"],"IGH_contigs":yield_map["IGH_contigs"]["all_contigs_count"],
        "V_assigned":yield_map["V_assigned"]["IGH_count"],"J_assigned":yield_map["J_assigned"]["IGH_count"],
        "C_assigned":yield_map["C_assigned"]["IGH_count"],"CDR3_assigned":yield_map["CDR3_assigned"]["IGH_count"],
        "productive_IGH":yield_map["productive_IGH"]["IGH_count"],"V+J+CDR3":yield_map["V+J+CDR3"]["IGH_count"],
        "V+J+C+CDR3":yield_map["V+J+C+CDR3"]["IGH_count"],"unique_clonotypes":v4["unique"],
        "reference_clonotypes":v4["reference"],"matched_clonotypes":v4["matches"],"precision":v4["precision"],"sensitivity":v4["sensitivity"],
        "Pearson_r":"NOT AVAILABLE","D_agreement":v4["d_agreement"],"exact_isotype_agreement":v4["isotype_exact"],"collapsed_isotype_agreement":v4["isotype_collapsed"]}
    metrics=["total_contigs","IGH_contigs","V_assigned","J_assigned","C_assigned","CDR3_assigned","productive_IGH","V+J+CDR3","V+J+C+CDR3",
             "unique_clonotypes","reference_clonotypes","matched_clonotypes","precision","sensitivity","Pearson_r","D_agreement","exact_isotype_agreement","collapsed_isotype_agreement"]
    labels={"TRUST4 FZ-116":"full FZ-116","Graph-TRUST4-v3-debug 100K":"100K subset","Graph-TRUST4-v4 100K":"100K subset"}
    rows=[]
    scale_note="TRUST4=full FZ-116 report (29,803 report records, not contigs); graph methods=100K subset contigs; raw differences are not input-size-equivalent"
    for metric in metrics:
        rows.append({"metric":metric,"TRUST4 FZ-116":trust.get(metric,"NOT AVAILABLE"),"Graph-TRUST4-v3-debug 100K":"NOT AVAILABLE",
                     "Graph-TRUST4-v4 100K":v4.get(metric,"NOT AVAILABLE"),"dataset_scales":scale_note})
    write_tsv(HERE/"TRUST4_V3_V4_COMPARISON.tsv",rows)
    return trust


def main():
    global ROOT, FASTA, REFERENCE, TRUST4_REPORT, METHOD_SCRIPT
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=ROOT)
    args=parser.parse_args()
    ROOT=args.root.resolve()
    FASTA=ROOT/"results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa"
    REFERENCE=ROOT/"reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz"
    TRUST4_REPORT=ROOT/"results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv"
    METHOD_SCRIPT=ROOT/"scripts/benchmark_trust4_published.py"
    method=load_method()
    valid,baseline,trust_raw,trust_unique,reference,trust_matches=baseline_check(method)
    print("TRUST4_BASELINE_VALIDATED="+str(valid).lower())
    if not valid:
        print("Baseline discrepancy; stopping before v4 evaluation outputs are written.")
        for key,val in baseline.items(): print(f"{key}={val}")
        raise SystemExit(2)

    ann=normalize_annotations(method)
    yields,length_rows,completeness,length_summary=yield_and_length_tables(ann)
    trust_corr={"precision":len(trust_matches)/len(trust_unique),"sensitivity":len(trust_matches)/len(reference),
                "pearson_r":method.calculate_abundance_correlation(trust_matches)["pearson_r"]}
    trust_gene={"D_agreement":float((trust_matches["d_gene_trust4"]==trust_matches["d_gene_irep"]).mean()),
                "iso_exact_agreement":float((trust_matches["c_gene_trust4_original"]==trust_matches["c_gene_irep_original"]).mean())}
    v4=benchmark_tables(method,ann,reference,trust_matches)
    fasta_length_by_id=dict(zip(ann["contig_id"],ann["contig_length"]))
    for row in length_rows:
        row["matched_primary_clonotypes"]=0
    for match in v4["matched_rows"]:
        contig_length=fasta_length_by_id.get(match["v4_contig_id"])
        for row,(low_bound,high_bound) in zip(length_rows,[(150,199),(200,249),(250,299),(300,399),(400,499),(500,749),(750,999),(1000,1499),(1500,math.inf)]):
            if contig_length is not None and low_bound<=contig_length<=high_bound:
                row["matched_primary_clonotypes"]+=1
                break
    write_tsv(HERE/"length_vs_annotation.tsv",length_rows)
    trust_raw_report=pd.read_csv(TRUST4_REPORT,sep="\t")
    baseline_summary=comparison_table(method,trust_raw_report,len(trust_unique),len(trust_matches),trust_corr,trust_gene,v4,yields)

    metric_file=pd.read_csv(HERE/"graph_v4_benchmark_metrics.tsv",sep="\t",dtype=str)
    values=dict(zip(metric_file["metric"],metric_file["value"]))
    strat=pd.read_csv(HERE/"abundance_stratified_benchmark.tsv",sep="\t",dtype=str).to_dict("records")
    breakdown=pd.read_csv(HERE/"matching_breakdown.tsv",sep="\t",dtype=str).to_dict("records")
    gene=pd.read_csv(HERE/"gene_agreement.tsv",sep="\t",dtype=str).to_dict("records")
    length_by_label={r["length_bin"]:r for r in length_rows}
    low,high=length_by_label["150-199"],length_by_label[">=1500"]
    igh_v_by_bin=", ".join(f"{r['length_bin']}: {r['V_assigned_IGH']}" for r in length_rows)
    igh_cdr3_by_bin=", ".join(f"{r['length_bin']}: {r['CDR3_assigned_IGH']}" for r in length_rows)
    validated_by_bin=", ".join(f"{r['length_bin']}: {r['matched_primary_clonotypes']}" for r in length_rows)
    summary_text=(
        "# Graph-TRUST4-v4 100K Biological Evaluation\n\n## Objective\n"
        "Determine whether branch-aware, quality-aware graph assembly improves recovery of biologically annotatable BCR sequence relative to the prior graph implementation.\n\n## Input\n"
        "- 100K paired-end subset; completed v4 contig FASTA.\n- FZ-116 iRepertoire reference.\n- TRUST4 baseline uses full FZ-116.\n\n## Annotation yield\n"
        + "| metric | all contigs | % all | IGH | % IGH |\n|---|---:|---:|---:|---:|\n"
        + "\n".join(f"| {r['metric']} | {r['all_contigs_count']} | {float(r['all_contigs_percent']):.3f} | {r['IGH_count']} | {float(r['IGH_percent']):.3f} |" for r in yields)
        + "\n\n## Contig length and completeness\n"
        + "Length-bin annotation counts are in `length_vs_annotation.tsv`; base columns count all chains, and corresponding `_IGH` columns/count percentages restrict to IGH V hits. Length summaries for IGH annotation categories are in `completeness_summary.tsv`.\n\n"
        + f"For 150-199 bp contigs (n={low['total_contigs']}, IGH n={low['IGH_contigs']}), IGH V assignment={low['V_assigned_IGH']}, CDR3={low['CDR3_assigned_IGH']}, V+J+CDR3={low['V_plus_J_plus_CDR3_IGH']}, V+J+C+CDR3={low['V_plus_J_plus_C_plus_CDR3_IGH']}.\n"
        + f"For >=1500 bp contigs (n={high['total_contigs']}, IGH n={high['IGH_contigs']}), IGH V assignment={high['V_assigned_IGH']}, CDR3={high['CDR3_assigned_IGH']}, V+J+CDR3={high['V_plus_J_plus_CDR3_IGH']}, V+J+C+CDR3={high['V_plus_J_plus_C_plus_CDR3_IGH']}.\n"
        + f"Across bins, IGH V assignments were {igh_v_by_bin}; IGH CDR3 assignments were {igh_cdr3_by_bin}; matched primary clonotypes were {validated_by_bin}. No IGH J or C assignments were present, so no length bin yielded a complete primary clonotype.\n"
        + "These are descriptive length strata; they do not establish that length caused annotation or biological recovery.\n\n## Primary benchmark\n"
        + f"- v4 unique clonotypes: {values['v4_unique_clonotypes']}\n- Reference clonotypes: {values['iRepertoire_reference_clonotypes']}\n- Matches: {values['matched_clonotypes']}\n"
        + f"- Precision: {values['precision']}\n- Sensitivity: {values['sensitivity']}\n- Pearson r: NOT AVAILABLE (v4 path abundance cannot be recovered from stored path tables).\n\n## Additional diagnostics\n"
        + "Relaxed matching results are in `matching_breakdown.tsv`. D and exact/collapsed isotype agreement are in `gene_agreement.tsv`.\n"
        + "Abundance-stratified reference counts are in `abundance_stratified_benchmark.tsv`; v4 abundance bins, precision, and Pearson correlation are unavailable because per-contig abundance is not stored. Matches are assigned to strata by iRepertoire reference copy.\n\n"
        + "## Comparison\nTRUST4 is full FZ-116; v3-debug/v4 are 100K subsets. The v3-debug biological output present in the project is full FZ-116, not the 100K subset; no v3 100K biological evaluation exists, so those values are NOT AVAILABLE. Raw count/metric differences across input scales are not directly equivalent. See `TRUST4_V3_V4_COMPARISON.tsv`.\n\n"
        + "## Interpretation\n### 1. Assembly-level observation\nThe v4 diagnostic assembly contains 131,032 contigs; its length statistics are in the adjacent computational diagnostics. This is an assembly-structure observation only.\n\n"
        + "### 2. Annotation-level observation\nThe yield and length-stratified tables show how V, J, CDR3, and productive IGH calls vary with contig length. Increased length is not itself evidence of improved biological recovery.\n\n"
        + f"### 3. Benchmark-level observation\nThe primary published-key evaluation yielded {values['matched_clonotypes']} validated clonotypes from {values['v4_unique_clonotypes']} v4 candidate clonotypes against {values['iRepertoire_reference_clonotypes']} reference clonotypes. Abundance concordance is unavailable.\n\n"
        + "The TRUST4 baseline was independently rerun in memory with `scripts/benchmark_trust4_published.py` and reproduced the specified values before v4 metrics were generated. No thresholds were tuned.\n\n"
        + "## Limitations\n- v4 was evaluated on a 100K-read subset.\n- iRepertoire is an external matched reference and is not necessarily molecule-identical to RNA-seq.\n- Low-abundance reference clonotypes can reduce sensitivity.\n- No threshold tuning was performed.\n- v4 path tables omit bundle/node membership and per-contig abundance; Pearson abundance correlation and v4 abundance-stratified precision cannot be calculated.\n- The available v3-debug biological results are full FZ-116, not a comparable 100K evaluation.\n")
    (HERE/"V4_100K_BIOLOGICAL_EVALUATION.md").write_text(summary_text)

    readme=(
        "# Graph-TRUST4-v4 100K Biological Evaluation\n\n"
        "## Inputs\n- FASTA: `results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa`\n"
        "- iRepertoire: `reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz`\n"
        "- TRUST4: `results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv`\n\n"
        "## Reproduction\nRun from project root. Annotation writes only into this evaluation folder:\n\n"
        "```sh\nalgorithms/TRUST4/annotator -f reference/TRUST4/human_IMGT+C.fa -a results/Graph-TRUST4-v4/test-100k/assembled_contigs.fa --fasta --needReverseComplement --noImpute -t 64 --outputFormat 1 > results/Graph-TRUST4-v4/test-100k/biological_evaluation/graph_v4_raw_airr.tsv\n"
        "python3 results/Graph-TRUST4-v4/test-100k/biological_evaluation/evaluate_graph_v4_biology.py\n```\n\n"
        "The evaluation script imports the functions in `scripts/benchmark_trust4_published.py` and first verifies the TRUST4 baseline. It stops before v4 output generation if any supplied baseline count/metric differs by more than 1e-9.\n\n"
        "## Normalization and matching\n"
        "Use uppercase CDR3 nucleotide sequence; trim 3 nt from both ends of TRUST4-style junctions; remove V/J/D/C allele suffixes; strip iRepertoire's leading `h`; collapse IGHA1/2 to IGHA, IGHG1-4 to IGHG, IGHD1-4 to IGHD, IGHE1/2 to IGHE, and IGHM1/2 to IGHM. Primary key is trimmed CDR3nt + normalized V + normalized J + collapsed C/isotype; D is diagnostic only. Incomplete primary keys are excluded. Identical v4 primary keys count once. Since no v4 abundance is stored, duplicates cannot be ranked by maximum abundance; the first contig is retained only for match-table traceability.\n\n"
        "TRUST4 and reference duplicate keys are coalesced by maximum count/copy using the published functions. Precision is matches / candidate unique clonotypes; sensitivity is matches / reference unique clonotypes. Baseline Pearson is the published Pearson correlation between TRUST4 count and iRepertoire copy among matches. V4 Pearson is unavailable because graph path node membership/abundance is not present; contig length is never used as a substitute.\n\n"
        "## Outputs\n`graph_v4_raw_airr.tsv` preserves raw annotator output; `graph_v4_annotated.tsv` includes standardized fields plus raw fields; `graph_v4_normalized.tsv` contains normalized calls. Other outputs: `annotation_yield.tsv`, `length_vs_annotation.tsv`, `graph_v4_benchmark_metrics.tsv`, `matching_breakdown.tsv`, `gene_agreement.tsv`, `abundance_correlation.tsv`, `graph_v4_vs_iRep_matches.tsv`, `abundance_stratified_benchmark.tsv`, `TRUST4_V3_V4_COMPARISON.tsv`, `completeness_summary.tsv`, and `V4_100K_BIOLOGICAL_EVALUATION.md`.\n\n"
        "No assembly, threshold tuning, or biological filtering based on iRepertoire was performed. The existing v3 biological evaluation is full FZ-116; no v3 100K biological output was found.\n")
    (HERE/"README.md").write_text(readme)

    print("=== V4 BIOLOGICAL EVALUATION ===\n")
    for key in ["total_contigs","IGH_contigs","V_assigned","J_assigned","C_assigned","CDR3_assigned","productive_IGH"]:
        yield_key={"total_contigs":"total_contigs","IGH_contigs":"IGH_contigs","V_assigned":"V_assigned","J_assigned":"J_assigned","C_assigned":"C_assigned","CDR3_assigned":"CDR3_assigned","productive_IGH":"productive_IGH"}[key]
        row=next(r for r in yields if r["metric"]==yield_key)
        print(f"{key.upper()}={row['all_contigs_count']}")
    for key,src in [("V_PLUS_J","V+J"),("V_PLUS_J_PLUS_CDR3","V+J+CDR3"),("V_PLUS_J_PLUS_C_PLUS_CDR3","V+J+C+CDR3")]:
        print(f"{key}={next(r['all_contigs_count'] for r in yields if r['metric']==src)}")
    for key,src in [("V4_UNIQUE_CLONOTYPES","unique"),("REFERENCE_CLONOTYPES","reference"),("MATCHES","matches"),("PRECISION","precision"),("SENSITIVITY","sensitivity"),("PEARSON_R","pearson_r")]:
        print(f"{key}={v4[src] if v4[src] is not None else 'NOT AVAILABLE'}")
    for key,src in [("D_AGREEMENT","d_agreement"),("ISOTYPE_EXACT_AGREEMENT","isotype_exact"),("ISOTYPE_COLLAPSED_AGREEMENT","isotype_collapsed")]:
        print(f"{key}={v4[src] if v4[src] is not None else 'NOT AVAILABLE'}")
    print("\n=== TRUST4 BASELINE VALIDATION ===\nTRUST4_BASELINE_VALIDATED=true")
    for label,key in [("TRUST4_UNIQUE","unique"),("TRUST4_REFERENCE","reference"),("TRUST4_MATCHES","matches"),("TRUST4_PRECISION","precision"),("TRUST4_SENSITIVITY","sensitivity"),("TRUST4_PEARSON_R","pearson_r")]:
        print(f"{label}={baseline[key]}")
    print("\n=== STATUS ===\nBIOLOGICAL_EVALUATION_COMPLETE=true")


if __name__ == "__main__":
    main()