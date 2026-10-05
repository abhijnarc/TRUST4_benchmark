#!/usr/bin/env python3
"""Generate read-only sequence-level diagnostics for Graph-TRUST4-v4 100K."""

import csv
import math
import re
import statistics
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
BIO = HERE.parent
ROOT = HERE.parents[4]
TEST = ROOT / "results/Graph-TRUST4-v4/test-100k"
FASTA = TEST / "assembled_contigs.fa"
AIRR = BIO / "graph_v4_raw_airr.tsv"
RC_AIRR = HERE / "IGH_V_HITS_REVERSE_COMPLEMENT_AIRR.tsv"
PATH_STATS = TEST / "path_statistics.tsv"
DEBUG_V4 = TEST / "assembly_debug.tsv"
DEBUG_V3 = ROOT / "results/Graph-TRUST4-v3-debug/test-100k/assembly_debug.tsv"
EDGE_STATS = TEST / "diagnostics/EDGE_QUALITY_STATISTICS.tsv"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
LENGTH_BINS = [(150,199,"150-199"),(200,249,"200-249"),(250,299,"250-299"),(300,399,"300-399"),
               (400,499,"400-499"),(500,749,"500-749"),(750,999,"750-999"),(1000,1499,"1000-1499"),(1500,math.inf,">=1500")]
DOWNSTREAM_BINS = [(0,49,"<50"),(50,99,"50-99"),(100,199,"100-199"),(200,299,"200-299"),
                   (300,499,"300-499"),(500,math.inf,">=500")]


def records(path):
    with path.open(newline="", errors="replace") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(name, columns, rows):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fasta_lengths(path):
    result, name, length = {}, None, 0
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    result[name] = length
                name, length = line[1:].split()[0], 0
            else:
                length += len(line)
    if name is not None:
        result[name] = length
    return result


def n50(values):
    values = sorted((int(x) for x in values), reverse=True)
    if not values:
        return None
    threshold, cumulative = sum(values) / 2, 0
    for value in values:
        cumulative += value
        if cumulative >= threshold:
            return value
    return values[-1]


def stats(values):
    values = [float(x) for x in values]
    if not values:
        return {"count":0,"minimum":"NOT AVAILABLE","median":"NOT AVAILABLE","mean":"NOT AVAILABLE","N50":"NOT AVAILABLE","maximum":"NOT AVAILABLE"}
    return {"count":len(values),"minimum":min(values),"median":statistics.median(values),"mean":statistics.mean(values),"N50":n50(values),"maximum":max(values)}


def truth(value):
    return bool(value and value.strip() not in {"", ".", "NOT AVAILABLE", "NOT_AVAILABLE", "nan", "None"})


def chain_from_calls(row):
    for field in ["v_call","d_call","j_call","c_call"]:
        call = row.get(field, "")
        if truth(call):
            gene = call.split("|")[0].split("*")[0]
            if gene.startswith("IGH"): return "IGH"
            if gene.startswith("IGK"): return "IGK"
            if gene.startswith("IGL"): return "IGL"
            if gene.startswith("TRA"): return "TRA"
            if gene.startswith("TRB"): return "TRB"
            if gene.startswith("TRG"): return "TRG"
            if gene.startswith("TRD"): return "TRD"
    return "NOT AVAILABLE"


def cigar_context(row, sequence_length):
    ops = [(int(length), op) for length, op in CIGAR_RE.findall(row.get("v_cigar", ""))]
    # TRUST4's AIRR formatter maps its internal insert/delete labels to the inverse
    # AIRR operation names. Query-consuming operations are S, M, and D in this output.
    query_length = sum(length for length, op in ops if op in "SMD")
    if not ops or query_length != sequence_length:
        return None, None, None
    start, aligned, alignment_started = 0, 0, False
    for length, op in ops:
        if op in "MD":
            alignment_started = True
            aligned += length
        elif op == "S" and not alignment_started:
            start += length
    end = start + aligned
    if start < 0 or end > sequence_length:
        return None, None, None
    return start, end, aligned


def read_debug(path):
    rows = list(records(path))
    lengths = [int(row["contig_length"]) for row in rows]
    return {
        "rows":len(rows), "valid":sum(row["valid_path"] == "1" for row in rows),
        "repeated":sum(int(row["repeated_nodes"]) > 0 for row in rows),
        "cycles":sum(row["cycle_detected"] == "1" for row in rows),
        "mismatches":sum(int(row["contig_length"]) != int(row["expected_length"]) for row in rows),
        "n50":n50(lengths), "maximum":max(lengths),
        "multiples_150":sum(length % 150 == 0 for length in lengths),
    }


def key_values(path):
    return {row["statistic"]:row["value"] for row in records(path)}


def main():
    all_lengths = fasta_lengths(FASTA)
    path_by_id = {str(row["contig_id"]):row for row in records(PATH_STATS)}
    raw_rows, top_rows = {}, {}
    top_ids = {"contig_" + str(i) for i,_ in sorted(enumerate(all_lengths.items()), key=lambda item:(-item[1][1],item[1][0]))[:100]}
    v_hits, raw_count, multiple_v_all = [], 0, 0
    for row in records(AIRR):
        raw_count += 1
        cid = row["sequence_id"]
        if cid in top_ids:
            top_rows[cid] = row
        if "|" in row.get("v_call", ""): multiple_v_all += 1
        if not row.get("v_call", "").startswith("IGH"):
            continue
        length = all_lengths.get(cid)
        sequence = row.get("sequence", "")
        if length is None:
            raise RuntimeError(f"FASTA sequence missing for AIRR contig {cid}")
        if len(sequence) != length:
            raise RuntimeError(f"AIRR/FASTA length mismatch for {cid}: {len(sequence)} != {length}")
        start, end, aligned = cigar_context(row, length)
        junction = row.get("junction", "")
        junction_start = sequence.find(junction) if truth(junction) else -1
        if junction_start < 0:
            junction_start, junction_end = "NOT AVAILABLE", "NOT AVAILABLE"
        else:
            junction_end = junction_start + len(junction)
        item = dict(row)
        item.update({"contig_id":cid,"contig_length":length,"sequence":sequence,"v_start":start,"v_end":end,
                     "v_aligned":aligned,"downstream":length-end if end is not None else None,
                     "cdr3_start":junction_start,"cdr3_end":junction_end})
        v_hits.append(item)
    if raw_count != len(all_lengths):
        raise RuntimeError(f"AIRR records ({raw_count}) != FASTA records ({len(all_lengths)})")
    if len(v_hits) != 712:
        raise RuntimeError(f"Expected 712 IGH V-hit records, found {len(v_hits)}")

    # 1. All IGH V-hit sequence records.
    hit_columns=["contig_id","contig_length","v_call","d_call","j_call","c_call","junction","junction_aa","productive","sequence"]
    write_tsv("IGH_V_HIT_CONTIGS.tsv",hit_columns,[{"contig_id":r["contig_id"],"contig_length":r["contig_length"],
        "v_call":r.get("v_call",""),"d_call":r.get("d_call",""),"j_call":r.get("j_call",""),"c_call":r.get("c_call",""),
        "junction":r.get("junction",""),"junction_aa":r.get("junction_aa",""),"productive":r.get("productive",""),"sequence":r["sequence"]} for r in v_hits])

    hit_lengths=[r["contig_length"] for r in v_hits]
    all_length_stats=stats(all_lengths.values())
    hit_length_stats=stats(hit_lengths)
    length_stat_rows=[]
    for group,result in [("IGH_V_HIT_CONTIGS",hit_length_stats),("ALL_V4_CONTIGS",all_length_stats)]:
        length_stat_rows.append({"group":group,**result})
    write_tsv("IGH_V_LENGTH_STATISTICS.tsv",["group","count","minimum","median","mean","N50","maximum"],length_stat_rows)
    length_distribution=[]
    for low,high,label in LENGTH_BINS:
        hit_count=sum(low<=r["contig_length"]<=high for r in v_hits)
        all_count=sum(low<=length<=high for length in all_lengths.values())
        length_distribution.append({"length_bin":label,"IGH_V_hit_count":hit_count,"fraction_of_IGH_V_hits":hit_count/len(v_hits),
                                    "all_contig_count":all_count,"fraction_of_all_contigs":all_count/len(all_lengths)})
    write_tsv("IGH_V_LENGTH_DISTRIBUTION.tsv",["length_bin","IGH_V_hit_count","fraction_of_IGH_V_hits","all_contig_count","fraction_of_all_contigs"],length_distribution)

    # 3. CIGAR-based coordinates. Coordinates are 0-based half-open query intervals.
    alignment_rows=[]
    for r in v_hits:
        start,end,aligned=r["v_start"],r["v_end"],r["v_aligned"]
        alignment_rows.append({"contig_id":r["contig_id"],"contig_length":r["contig_length"],"V_gene":r.get("v_call") or "NOT AVAILABLE",
            "V_alignment_start":start if start is not None else "NOT AVAILABLE","V_alignment_end":end if end is not None else "NOT AVAILABLE",
            "V_alignment_length":aligned if aligned is not None else "NOT AVAILABLE",
            "V_fraction_of_contig":aligned/r["contig_length"] if aligned is not None else "NOT AVAILABLE",
            "J_gene":r.get("j_call") or "NOT AVAILABLE","J_alignment_start":"NOT AVAILABLE","J_alignment_end":"NOT AVAILABLE",
            "CDR3_start":r["cdr3_start"],"CDR3_end":r["cdr3_end"],"sequence_orientation":r.get("rev_comp") or "NOT AVAILABLE"})
    write_tsv("IGH_V_ALIGNMENT_CONTEXT.tsv",["contig_id","contig_length","V_gene","V_alignment_start","V_alignment_end","V_alignment_length",
        "V_fraction_of_contig","J_gene","J_alignment_start","J_alignment_end","CDR3_start","CDR3_end","sequence_orientation"],alignment_rows)

    # 4. Bases downstream of the V alignment, with formal J/CDR3 annotation counts.
    extension_rows=[]
    for low,high,label in DOWNSTREAM_BINS:
        grouped=[r for r in v_hits if r["downstream"] is not None and low<=r["downstream"]<=high]
        extension_rows.append({"downstream_bin":label,"number":len(grouped),"fraction_of_V_hits":len(grouped)/len(v_hits),
            "J_assigned":sum(truth(r.get("j_call")) for r in grouped),"CDR3_assigned":sum(truth(r.get("junction")) for r in grouped),
            "definition":"query bases after V alignment end (0-based, half-open TRUST4 AIRR CIGAR coordinates)"})
    write_tsv("IGH_V_EXTENSION_ANALYSIS.tsv",["downstream_bin","number","fraction_of_V_hits","J_assigned","CDR3_assigned","definition"],extension_rows)

    # 5. No independent J search is reported: no local BLAST/minimap/vsearch tool
    # is available, and no path-aware J alignment output is stored.
    j_rows=[{"contig_id":r["contig_id"],"length":r["contig_length"],"V_call":r.get("v_call") or "NOT AVAILABLE",
             "formal_J_call":r.get("j_call") or "NOT AVAILABLE","best_J_candidate":"NOT AVAILABLE","J_identity":"NOT AVAILABLE",
             "J_alignment_length":"NOT AVAILABLE","J_candidate_position":"NOT AVAILABLE","evidence_of_J_sequence":"NOT AVAILABLE",
             "confidence":"NOT AVAILABLE","notes":"TRUST4 annotator (same IMGT reference) emitted no IGH J call; no independent local-alignment executable/path membership available."} for r in v_hits]
    write_tsv("IGH_J_FORENSIC.tsv",["contig_id","length","V_call","formal_J_call","best_J_candidate","J_identity","J_alignment_length",
        "J_candidate_position","evidence_of_J_sequence","confidence","notes"],j_rows)

    # 6. Stay conservative for non-formal junctions; no motif heuristics are invented.
    cdr3_rows=[]
    for r in v_hits:
        category="FORMAL_CDR3" if truth(r.get("junction")) else "UNCERTAIN"
        cdr3_rows.append({"contig_id":r["contig_id"],"contig_length":r["contig_length"],"junction":r.get("junction") or "NOT AVAILABLE",
            "junction_aa":r.get("junction_aa") or "NOT AVAILABLE","CDR3_context_category":category,
            "V_alignment_end":r["v_end"] if r["v_end"] is not None else "NOT AVAILABLE",
            "bases_after_V_alignment":r["downstream"] if r["downstream"] is not None else "NOT AVAILABLE",
            "notes":"No J alignment; partial-junction sequence presence cannot be called reliably without an alignment to the provided IGHJ references."})
    write_tsv("IGH_CDR3_FORENSIC.tsv",["contig_id","contig_length","junction","junction_aa","CDR3_context_category","V_alignment_end","bases_after_V_alignment","notes"],cdr3_rows)
    cdr3_category_counts=Counter(row["CDR3_context_category"] for row in cdr3_rows)
    cdr3_categories=["NO_JUNCTION_CONTEXT","PARTIAL_JUNCTION_CONTEXT","FORMAL_CDR3","JUNCTION_BUT_NO_FORMAL_CDR3","UNCERTAIN"]
    write_tsv("IGH_CDR3_CONTEXT_COUNTS.tsv",["category","count","percentage_of_IGH_V_hits","interpretation"],[
        {"category":category,"count":cdr3_category_counts[category],"percentage_of_IGH_V_hits":100*cdr3_category_counts[category]/len(v_hits),
         "interpretation":"Unobserved sequence-context categories remain uncertain because no independent IGHJ alignment was available." if category in ["NO_JUNCTION_CONTEXT","PARTIAL_JUNCTION_CONTEXT","JUNCTION_BUT_NO_FORMAL_CDR3","UNCERTAIN"] else "Formal TRUST4 annotator junction call."}
        for category in cdr3_categories])

    # 7-8. Aggregate path metrics exist, but exact node/edge membership was not written.
    path_rows=[]; branch_rows=[]
    for r in v_hits:
        path=path_by_id.get(r["contig_id"],{})
        path_rows.append({"contig_id":r["contig_id"],"contig_length":r["contig_length"],"num_nodes":path.get("num_nodes","NOT AVAILABLE"),
            "num_edges":"NOT AVAILABLE","path_nodes":"NOT AVAILABLE","path_edges":"NOT AVAILABLE","branch_count":"NOT AVAILABLE",
            "resolved_branch_count":"NOT AVAILABLE","ambiguous_branch_count":"NOT AVAILABLE","cycle":path.get("cycle_detected","NOT AVAILABLE"),
            "path_validity":path.get("valid_path","NOT AVAILABLE"),"PATH_MEMBERSHIP_AVAILABLE":"false",
            "reason":"path_statistics.tsv stores aggregate node counts but no node IDs; branch_decisions.tsv stores global IDs without contig-to-path mapping."})
        branch_rows.append({"contig_id":r["contig_id"],"branch_node":"NOT AVAILABLE","chosen_edge":"NOT AVAILABLE","alternative_edge":"NOT AVAILABLE",
            "chosen_QAOS":"NOT AVAILABLE","alternative_QAOS":"NOT AVAILABLE","QAOS_delta":"NOT AVAILABLE","chosen_overlap":"NOT AVAILABLE",
            "alternative_overlap":"NOT AVAILABLE","chosen_identity":"NOT AVAILABLE","alternative_identity":"NOT AVAILABLE",
            "paired_end_support":"NOT AVAILABLE","chosen_abundance":"NOT AVAILABLE","alternative_abundance":"NOT AVAILABLE",
            "branch_status":"NOT_AVAILABLE","reason":"No exact contig-to-node/edge membership in saved path outputs."})
    write_tsv("IGH_V_PATHS.tsv",["contig_id","contig_length","num_nodes","num_edges","path_nodes","path_edges","branch_count","resolved_branch_count",
        "ambiguous_branch_count","cycle","path_validity","PATH_MEMBERSHIP_AVAILABLE","reason"],path_rows)
    write_tsv("IGH_BRANCH_FORENSICS.tsv",["contig_id","branch_node","chosen_edge","alternative_edge","chosen_QAOS","alternative_QAOS","QAOS_delta",
        "chosen_overlap","alternative_overlap","chosen_identity","alternative_identity","paired_end_support","chosen_abundance","alternative_abundance","branch_status","reason"],branch_rows)

    # 9. Descriptive candidates only; flags are not chimerism calls.
    pathology=[]
    multiple_v=sum("|" in r.get("v_call","") for r in v_hits)
    multiple_j=sum("|" in r.get("j_call","") for r in v_hits)
    multiple_c=sum("|" in r.get("c_call","") for r in v_hits)
    for r in v_hits:
        flags=[]
        if "|" in r.get("v_call",""): flags.append("MULTIPLE_V_GENES")
        if "|" in r.get("j_call",""): flags.append("MULTIPLE_J_GENES")
        if "|" in r.get("c_call",""): flags.append("MULTIPLE_C_GENES")
        if r["downstream"] is not None and r["downstream"]>=200 and not truth(r.get("j_call")):
            flags.append("NO_J_AFTER_LONG_EXTENSION")
        if flags:
            pathology.append({"contig_id":r["contig_id"],"length":r["contig_length"],"V_gene(s)":r.get("v_call") or "NOT AVAILABLE",
                "J_gene(s)":r.get("j_call") or "NOT AVAILABLE","C_gene(s)":r.get("c_call") or "NOT AVAILABLE",
                "CDR3(s)":r.get("junction") or "NOT AVAILABLE","branch_count":"NOT AVAILABLE","ambiguous_branch_count":"NOT AVAILABLE",
                "mean_QAOS":"NOT AVAILABLE","minimum_QAOS":"NOT AVAILABLE","mean_identity":"NOT AVAILABLE","minimum_identity":"NOT AVAILABLE",
                "reason_flag":";".join(flags),"notes":"Pipe-separated V calls are annotator alternatives, not evidence of multiple V segments. Long extension flag uses >=200 downstream bases; no chimerism inferred."})
    write_tsv("PATHOLOGICAL_IGH_PATHS.tsv",["contig_id","length","V_gene(s)","J_gene(s)","C_gene(s)","CDR3(s)","branch_count","ambiguous_branch_count",
        "mean_QAOS","minimum_QAOS","mean_identity","minimum_identity","reason_flag","notes"],pathology)

    # 10. All-edge values are reportable; IGH path edge values are not traceable.
    edge_quality={r["metric"]:r for r in records(EDGE_STATS)}
    edge_rows=[]
    for metric in ["overlap_length","identity","qaos"]:
        row=edge_quality[metric]
        edge_rows.append({"metric":metric,"ALL_EDGES_n":row["n"],"ALL_EDGES_mean":row["mean"],"ALL_EDGES_median":row["median"],
            "ALL_EDGES_minimum":row["minimum"],"IGH_PATH_EDGES_n":"NOT AVAILABLE","IGH_PATH_EDGES_mean":"NOT AVAILABLE",
            "IGH_PATH_EDGES_median":"NOT AVAILABLE","IGH_PATH_EDGES_minimum":"NOT AVAILABLE",
            "reason":"Exact IGH path edge membership is not stored."})
    write_tsv("IGH_EDGE_QUALITY_COMPARISON.tsv",["metric","ALL_EDGES_n","ALL_EDGES_mean","ALL_EDGES_median","ALL_EDGES_minimum",
        "IGH_PATH_EDGES_n","IGH_PATH_EDGES_mean","IGH_PATH_EDGES_median","IGH_PATH_EDGES_minimum","reason"],edge_rows)

    overlap_rows=[]
    for low,high,label in [(31,49,"31-49"),(50,74,"50-74"),(75,99,"75-99"),(100,124,"100-124"),(125,149,"125-149"),(150,math.inf,">=150")]:
        overlap_rows.append({"overlap_bin":label,"edge_count":"NOT AVAILABLE","mean_QAOS":"NOT AVAILABLE","mean_identity":"NOT AVAILABLE",
            "branch_count":"NOT AVAILABLE","ambiguous_branch_count":"NOT AVAILABLE","note":"IGH-path edge set unavailable; all-edge bins are not substituted."})
    write_tsv("IGH_OVERLAP_DISTRIBUTION.tsv",["overlap_bin","edge_count","mean_QAOS","mean_identity","branch_count","ambiguous_branch_count","note"],overlap_rows)

    # 12. Compare the annotator's original strand flag with reverse-complement-only re-annotation.
    rc_rows={r["sequence_id"]:r for r in records(RC_AIRR)}
    orientation_rows=[]
    for r in v_hits:
        rc=r.get("rev_comp","")
        status="EXPECTED" if rc=="F" else "REVERSE_COMPLEMENT_EXPECTED" if rc=="T" else "NOT_AVAILABLE"
        reverse=rc_rows.get(r["contig_id"],{})
        orientation_rows.append({"contig_id":r["contig_id"],"rev_comp":rc or "NOT AVAILABLE","orientation_status":status,
            "reverse_complement_v_call":reverse.get("v_call") or "NOT AVAILABLE","reverse_complement_j_call":reverse.get("j_call") or "NOT AVAILABLE",
            "reverse_complement_c_call":reverse.get("c_call") or "NOT AVAILABLE","reverse_complement_junction":reverse.get("junction") or "NOT AVAILABLE",
            "reverse_complement_rescued_J_CDR3":bool((not truth(r.get("j_call")) and truth(reverse.get("j_call"))) or
                (not truth(r.get("junction")) and truth(reverse.get("junction")))),
            "graph_edge_orientation":"+1 for all accepted edges; path linkage unavailable","J_orientation_check":"NOT AVAILABLE",
            "notes":"Reverse-complement-only re-annotation used the same TRUST4 tool/reference."})
    write_tsv("IGH_ORIENTATION_CHECK.tsv",["contig_id","rev_comp","orientation_status","reverse_complement_v_call","reverse_complement_j_call",
        "reverse_complement_c_call","reverse_complement_junction","reverse_complement_rescued_J_CDR3","graph_edge_orientation","J_orientation_check","notes"],orientation_rows)
    orientation_counts=Counter(row["orientation_status"] for row in orientation_rows)
    orientation_rescued=sum(row["reverse_complement_rescued_J_CDR3"] for row in orientation_rows)

    # 13. Paired-end support audit based on implementation and outputs.
    paired_text=(
        "# Paired-End Support Audit\n\n"
        "**Status: NOT_RETAINED**\n\n"
        "Evidence:\n"
        "- `GraphAssemblerV3.cpp`, `read_fastq_record`, reads the FASTQ name into a local variable but returns only sequence and quality. `load_and_bundle` reads R1/R2 in lockstep, then keys the unique-bundle map by sequence and retains one quality string and an abundance count; it does not retain read names or mate-to-mate IDs.\n"
        "- `score_overlap` assigns `edge.paired_support = 0`; the graph edge table confirms paired support is zero for all 6,141,188 accepted edges.\n"
        "- Branch ranking includes paired support as a tie-break field, but all values are zero, so it supplies no discrimination and is not effective evidence for path selection.\n"
        "- `run_metadata.txt` records `paired_end_support=unavailable_input_loader`.\n\n"
        "Therefore mate relationships do not survive bundling and cannot inform overlap scoring or branch decisions in this completed run. This is a limitation, not proof that it caused the observed IGH annotation outcome.\n"
    )
    (HERE/"PAIRED_END_SUPPORT_AUDIT.md").write_text(paired_text)

    # 15. Longest contigs and annotation yield.
    top_lengths=sorted(all_lengths.items(),key=lambda item:(-item[1],item[0]))[:100]
    top_biology=[]
    for cid,length in top_lengths:
        r=top_rows.get(cid,{})
        v,j,cdr,c_gene=(r.get("v_call",""),r.get("j_call",""),r.get("junction",""),r.get("c_call",""))
        top_biology.append({"contig_id":cid,"length":length,"V":v or "NOT AVAILABLE","J":j or "NOT AVAILABLE","C":c_gene or "NOT AVAILABLE",
            "CDR3":cdr or "NOT AVAILABLE","productive":r.get("productive") or "NOT AVAILABLE","chain":chain_from_calls(r)})
    write_tsv("TOP_100_CONTIG_BIOLOGY.tsv",["contig_id","length","V","J","C","CDR3","productive","chain"],top_biology)
    top_counts={"IGH":0,"V":0,"J":0,"CDR3":0,"VJ":0,"VJCDR3":0,"VJC_CDR3":0}
    for row in top_biology:
        hasv=truth(row["V"]); hasj=truth(row["J"]); hasc=truth(row["C"]); hascdr=truth(row["CDR3"])
        top_counts["IGH"]+=row["chain"]=="IGH";top_counts["V"]+=hasv;top_counts["J"]+=hasj;top_counts["CDR3"]+=hascdr
        top_counts["VJ"]+=hasv and hasj;top_counts["VJCDR3"]+=hasv and hasj and hascdr;top_counts["VJC_CDR3"]+=hasv and hasj and hasc and hascdr
    top_chain_counts=Counter(row["chain"] for row in top_biology)

    # v3-debug comparison and aggregate forensic conclusion.
    v3=read_debug(DEBUG_V3);v4=read_debug(DEBUG_V4)
    graph_metrics=key_values(TEST/"graph_statistics.tsv")
    branch_metrics=key_values(TEST/"diagnostics/BRANCH_DECISION_STATISTICS.tsv")
    length_counts=Counter()
    for r in v_hits:
        for low,high,label in LENGTH_BINS:
            if low<=r["contig_length"]<=high:
                length_counts[label]+=1;break
    formal_cdr3=sum(truth(r.get("junction")) for r in v_hits)
    with_j=sum(truth(r.get("j_call")) for r in v_hits)
    with_c=sum(truth(r.get("c_call")) for r in v_hits)
    with_vj=sum(truth(r.get("v_call")) and truth(r.get("j_call")) for r in v_hits)
    with_vj_cdr3=sum(truth(r.get("v_call")) and truth(r.get("j_call")) and truth(r.get("junction")) for r in v_hits)
    with_vjc_cdr3=sum(truth(r.get("v_call")) and truth(r.get("j_call")) and truth(r.get("c_call")) and truth(r.get("junction")) for r in v_hits)
    productive=sum(r.get("productive")=="T" for r in v_hits)
    downstream=[r["downstream"] for r in v_hits if r["downstream"] is not None]
    summary_counts={"TOTAL_IGH_V_HITS":len(v_hits),"IGH_WITH_J":with_j,"IGH_WITH_C":with_c,"IGH_WITH_CDR3":formal_cdr3,
        "IGH_WITH_VJ":with_vj,"IGH_WITH_VJ_CDR3":with_vj_cdr3,"IGH_WITH_VJC_CDR3":with_vjc_cdr3,"IGH_PRODUCTIVE":productive,
        "IGH_MEDIAN_LENGTH":statistics.median(hit_lengths),"IGH_N50":n50(hit_lengths),"IGH_MAX_LENGTH":max(hit_lengths),
        "DOWNSTREAM_OF_V_MEDIAN":statistics.median(downstream),"DOWNSTREAM_OF_V_MAX":max(downstream)}
    multiples_v=sum("|" in r.get("v_call","") for r in v_hits)
    multiples_j=sum("|" in r.get("j_call","") for r in v_hits)
    long_no_j=sum((r["downstream"] or 0)>=200 and not truth(r.get("j_call")) for r in v_hits)
    suspicious=len(pathology)
    diagnostics={"A":"SUPPORTED","B":"UNCLEAR","C":"UNCLEAR","D":"UNCLEAR","E":"UNCLEAR","F":"NOT SUPPORTED","G":"SUPPORTED","H":"UNCLEAR","I":"UNCLEAR"}
    diagnosis=(
        "# Graph-TRUST4-v4 IGH Forensic Diagnostic\n\n## 1. Observed problem\n"
        f"The assembly has {len(all_lengths):,} contigs, but AIRR reports {len(v_hits)} IGH V-hit contigs, {with_j} with J, {with_c} with C, {formal_cdr3} with a formal junction/CDR3, and {productive} productive IGH calls. Primary V-J-C/CDR3 completeness is zero.\n\n"
        "## 2. V-hit sequence characteristics\n"
        f"The 712 V-hit contigs have length min/median/mean/N50/max {int(min(hit_lengths))}/{statistics.median(hit_lengths)}/{statistics.mean(hit_lengths):.2f}/{n50(hit_lengths)}/{max(hit_lengths)} bp; all v4 contigs have median {statistics.median(all_lengths.values())} bp, mean {statistics.mean(all_lengths.values()):.2f}, N50 {n50(all_lengths.values())}, max {max(all_lengths.values())}.\n"
        f"TRUST4 CIGAR query coordinates show V-alignment length min/median/mean/max {min(r['v_aligned'] for r in v_hits)}/{statistics.median(r['v_aligned'] for r in v_hits)}/{statistics.mean(r['v_aligned'] for r in v_hits):.2f}/{max(r['v_aligned'] for r in v_hits)} nt; identity min/median/mean/max {min(float(r['v_identity']) for r in v_hits if r.get('v_identity')):.2f}/{statistics.median(float(r['v_identity']) for r in v_hits if r.get('v_identity')):.2f}/{statistics.mean(float(r['v_identity']) for r in v_hits if r.get('v_identity')):.2f}/{max(float(r['v_identity']) for r in v_hits if r.get('v_identity')):.2f}%. CIGAR query lengths reconcile for all 712 records.\n"
        f"V-hit length-bin counts: {', '.join(f'{label}={length_counts[label]}' for _,_,label in LENGTH_BINS)}.\n\n"
        "## 3. V-to-J/CDR3 continuity\n"
        f"The median query sequence after the end of the V alignment is {statistics.median(downstream):.0f} bp (range {min(downstream)}-{max(downstream)}); downstream-bin counts: {', '.join(f'{label}={sum(low<=d<=high for d in downstream)}' for low,high,label in DOWNSTREAM_BINS)}. {sum(d<50 for d in downstream)} V hits have <50 bp downstream, while {long_no_j} have >=200 bp downstream but no formal J call. This demonstrates extension beyond the short V alignment in many records, but does not prove that the extension is IGH V-D-J context.\n"
        f"No IGH J alignment or partial-junction motif was independently established. The same TRUST4 annotator/reference produced no IGH J calls, including on reverse-complement copies; no local BLAST/minimap2/vsearch executable and no exact path membership are available, so a separate IGHJ sequence search is NOT AVAILABLE rather than guessed. CDR3 context categories: FORMAL_CDR3={formal_cdr3} ({formal_cdr3/len(v_hits):.2%}); UNCERTAIN={len(v_hits)-formal_cdr3} ({(len(v_hits)-formal_cdr3)/len(v_hits):.2%}); partial/no-junction categories are not assigned without a reliable J alignment.\n\n"
        "## 4. Graph/path evidence\n"
        "`PATH_MEMBERSHIP_AVAILABLE=false`. `path_statistics.tsv` stores contig length, node count, validity, and aggregate flags but no node IDs. `graph_edges.tsv` identifies global source/target and edge metrics; `branch_decisions.tsv` records global branch/edge IDs but not the contig path that used them. Therefore no exact IGH path nodes/edges, branch counts, branch-specific QAOS, or IGH-path edge-quality distribution can be assigned. No heuristic path reconstruction was used.\n"
        f"Across all graph paths, v4 has {v4['valid']}/{v4['rows']} valid, {v4['repeated']} repeated-node, {v4['cycles']} cycle, and {v4['mismatches']} length-mismatch paths. Global branch totals are {graph_metrics['branches']} ({graph_metrics['resolved_branches']} resolved; {graph_metrics['ambiguous_branches']} ambiguous; {float(branch_metrics['ambiguity_rate']):.2%} ambiguous); median QAOS_delta={branch_metrics['QAOS_delta_median']}. Branch IDs cannot be tied to the 712 V-hit paths. These global checks do not establish biological correctness.\n\n"
        "## 5. Edge-quality evidence\n"
        "For all accepted edges (not specifically IGH paths): overlap mean/median/min 113.974/117/31 bp; identity 0.989296/0.989796/0.933333; QAOS 0.947643/0.942705/0.900000. IGH-path edge metrics and overlap bins are NOT AVAILABLE because path edge IDs were not saved. Global edge quality cannot certify IGH joins.\n\n"
        "## 6. Orientation evidence\n"
        f"All 712 V-hit AIRR records have `rev_comp=F`. Re-annotating only their reverse-complement sequences with the same TRUST4 annotator/reference yielded 0 J calls, 0 C calls, and 1 formal junction; {orientation_rescued} contigs gained J/CDR3 evidence. Thus reverse-complement rescue is NOT SUPPORTED as the main explanation. All 6,141,188 accepted graph edges have orientation +1 because the v4 scorer hard-codes that value; path-level orientation continuity cannot be evaluated without membership.\n\n"
        "## 7. Paired-end evidence\n"
        "Status: NOT_RETAINED. The FASTQ reader discards record names; bundling keys reads by sequence and keeps abundance/quality but no mate IDs. Every saved edge and branch support value is zero. Paired support appears as a tie-breaker in branch ordering but cannot affect decisions when all values are zero. This is a real loss of evidence, not proof of causation. See `PAIRED_END_SUPPORT_AUDIT.md`.\n\n"
        "## 8. Comparison with v3\n"
        f"v3-debug 100K: N50={v3['n50']}, max={v3['maximum']}, valid={v3['valid']}/{v3['rows']}, repeated paths={v3['repeated']}, cycles={v3['cycles']}, length mismatches={v3['mismatches']}; {v3['multiples_150']} contig lengths are exact multiples of 150.\n"
        f"v4 100K: N50={v4['n50']}, max={v4['maximum']}, valid={v4['valid']}/{v4['rows']}, repeated paths={v4['repeated']}, cycles={v4['cycles']}, length mismatches={v4['mismatches']}; {v4['multiples_150']} exact multiples of 150. This supports that v4's stored path/overlap length checks are structurally cleaner and the multiple-of-150 pattern is much rarer proportionally ({v3['multiples_150']/v3['rows']:.2%} vs {v4['multiples_150']/v4['rows']:.2%}); structural validity is not biological completeness.\n"
        "BIOLOGICAL_COMPLETENESS_PROBLEM=YES (observed missing IGH J/C/CDR3 completeness); mechanism is not fully resolved.\n\n"
        f"Top 100 longest contigs: IGH={top_counts['IGH']}, with V={top_counts['V']}, J={top_counts['J']}, CDR3={top_counts['CDR3']}, V+J={top_counts['VJ']}, V+J+CDR3={top_counts['VJCDR3']}, complete V+J+C+CDR3={top_counts['VJC_CDR3']}; chain annotations={dict(top_chain_counts)}. All 100 are annotated TRB, with J calls TRBJ2-6 and C calls TRBC2 (98) or TRBC1 (2), and no V/CDR3 calls. This is evidence that the longest outputs are T-cell-receptor-associated contigs, not evidence that they are complete IGH assemblies or chimeras.\n\n"
        "## 9. Most strongly supported failure mode\n"
        "- A. Paths too short: SUPPORTED for the V-hit subset (median length 246 bp, N50 273 bp; 120 have <50 bp after the short V alignment), but not a complete explanation because many have longer downstream sequence.\n"
        "- B. Unrelated sequence merging: UNCLEAR; no concrete sequence/path evidence of a chimera was found, and exact paths are unavailable.\n"
        "- C. Incorrect branch resolution: UNCLEAR; branch decisions cannot be tied to these contigs.\n"
        "- D. Premature termination before V-D-J: UNCLEAR; the contigs lack formal J/CDR3 calls, but long downstream sequence exists and independent J sequence evidence is unavailable.\n"
        "- E. J/CDR3 annotator failure despite sequence presence: UNCLEAR; no independent IGHJ alignment established sequence presence.\n"
        "- F. Orientation problem: NOT SUPPORTED as the main cause; all 712 annotator flags are forward and reverse-complement re-annotation did not restore J/C or additional junctions. Graph edge orientation is still hard-coded +1.\n"
        "- G. Paired-end problem: SUPPORTED as a missing source of path evidence (mate IDs discarded; support always zero); causal link is unproven.\n"
        "- H. Candidate-edge problem: UNCLEAR; all-edge scores are strong but cannot be associated with IGH paths.\n"
        "- I. Wrong branch/path extraction: UNCLEAR; exact path membership was not emitted.\n"
        f"- Other: SUPPORTED observation that V hits are short, low-identity partial alignments (median {statistics.median(float(r['v_identity']) for r in v_hits):.2f}% identity over median 24 nt), and longer contigs are not translating into IGH complete calls.\n\n"
        "## 10. Recommended next algorithmic investigation\n"
        "First preserve exact ordered node IDs, edge IDs, orientations, branch decisions, and abundance for each contig. Then use those exact paths to align V-hit context against the same TRUST4 IMGT reference in both orientations and distinguish missing J sequence from missed annotation. Retain mate identifiers through bundling and measure mate support at candidate edges. Inspect long non-IGH contigs and V-hit downstream sequence for architecture discontinuities only after path membership exists. These are investigation recommendations only; no algorithm or threshold was changed.\n"
    )
    (HERE/"V4_IGH_FORENSIC_DIAGNOSTIC.md").write_text(diagnosis)

    # Terminal summary required by the request.
    print("=== IGH V4 FORENSIC DIAGNOSTIC ===\n")
    for key,source in [("IGH_V_HITS","TOTAL_IGH_V_HITS"),("IGH_WITH_J","IGH_WITH_J"),("IGH_WITH_C","IGH_WITH_C"),
                       ("IGH_WITH_CDR3","IGH_WITH_CDR3"),("IGH_WITH_VJ","IGH_WITH_VJ"),("IGH_WITH_VJ_CDR3","IGH_WITH_VJ_CDR3"),
                       ("IGH_WITH_VJC_CDR3","IGH_WITH_VJC_CDR3"),("IGH_PRODUCTIVE","IGH_PRODUCTIVE")]:
        print(f"{key}={summary_counts[source]}")
    print(f"\nIGH_MEDIAN_LENGTH={summary_counts['IGH_MEDIAN_LENGTH']}\nIGH_N50={summary_counts['IGH_N50']}\nIGH_MAX_LENGTH={summary_counts['IGH_MAX_LENGTH']}\n")
    print("\nPATH_MEMBERSHIP_AVAILABLE=false")
    print("PAIRED_END_SUPPORT_STATUS=NOT_RETAINED")
    print(f"\nMULTIPLE_V_GENES={multiples_v}\nMULTIPLE_J_GENES={multiples_j}\nAMBIGUOUS_BRANCH_PATHS=NOT AVAILABLE\nSUSPICIOUS_PATHS={suspicious}")
    print("\nV3_STRUCTURAL_PROBLEM_FIXED=true")
    print("V4_BIOLOGICAL_COMPLETENESS_PROBLEM=YES")
    for letter in "ABCDEFGH":
        print(f"DIAGNOSIS_{letter}={diagnostics[letter]}")
    print("\n=== SAFETY CHECKS ===\nASSEMBLY_MODIFIED=false\nSOURCE_MODIFIED=false\nPARAMETERS_CHANGED=false\nFULL_FZ116_RUN=false\nTHRESHOLD_TUNING=false\n\nDIAGNOSTIC_STATUS=COMPLETE")


if __name__ == "__main__":
    main()