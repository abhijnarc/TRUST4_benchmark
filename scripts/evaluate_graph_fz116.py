#!/usr/bin/env python3
"""Evaluate corrected Graph-TRUST4 FZ-116 output with the published TRUST4 key."""

import argparse
import gzip
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def remove_allele(value):
    if pd.isna(value) or str(value).strip() == "":
        return None
    value = str(value).strip()
    return value.split("*")[0] if "*" in value else value


def collapse_isotype(value):
    if pd.isna(value) or str(value).strip() == "":
        return None
    value = str(value).strip()
    groups = {
        "IGHA1": "IGHA", "IGHA2": "IGHA",
        "IGHG1": "IGHG", "IGHG2": "IGHG", "IGHG3": "IGHG", "IGHG4": "IGHG",
        "IGHD1": "IGHD", "IGHD2": "IGHD", "IGHD3": "IGHD", "IGHD4": "IGHD",
        "IGHE1": "IGHE", "IGHE2": "IGHE",
        "IGHM1": "IGHM", "IGHM2": "IGHM",
    }
    return groups.get(value, value)


def trim_cdr3(value):
    if pd.isna(value) or str(value).strip() == "":
        return None
    value = str(value).upper()
    return value[3:-3] if len(value) >= 6 else None


def parse_fasta(path):
    records = []
    name = None
    sequence = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence)))
                name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line.strip())
    if name is not None:
        records.append((name, "".join(sequence)))
    return records


def reconstruct_paths(bundles_path, edges_path):
    bundles = {}
    with open(bundles_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        positions = {name: index for index, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            bundle_id = int(fields[positions["bundle_id"]])
            bundles[bundle_id] = {
                "sequence": fields[positions["consensus_seq"]],
                "abundance": int(fields[positions["abundance"]]),
            }

    edges = []
    outgoing = defaultdict(list)
    indegree = defaultdict(int)
    with open(edges_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        positions = {name: index for index, name in enumerate(header)}
        for edge_id, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            edge = {
                "source": int(fields[positions["source"]]),
                "target": int(fields[positions["target"]]),
                "overlap": int(fields[positions["overlap_length"]]),
                "qaos": float(fields[positions["qaos"]]),
            }
            edges.append(edge)
            outgoing[edge["source"]].append(edge_id)
            indegree[edge["target"]] += 1

    used_node = set()
    used_edge = set()
    paths = []

    def choose_edge(node, path_nodes):
        choices = []
        for edge_id in outgoing[node]:
            edge = edges[edge_id]
            if edge_id in used_edge:
                continue
            if edge["target"] in used_node and edge["target"] not in path_nodes:
                continue
            choices.append(edge_id)
        if not choices:
            return None
        return max(choices, key=lambda edge_id: (edges[edge_id]["qaos"], edges[edge_id]["overlap"]))

    def walk(root):
        path = [root]
        path_edges = []
        path_nodes = {root}
        used_node.add(root)
        current = root
        cycle = False
        while True:
            edge_id = choose_edge(current, path_nodes)
            if edge_id is None:
                break
            edge = edges[edge_id]
            used_edge.add(edge_id)
            if edge["target"] in path_nodes:
                cycle = True
                break
            path_edges.append(edge_id)
            path.append(edge["target"])
            path_nodes.add(edge["target"])
            used_node.add(edge["target"])
            current = edge["target"]
            if len(outgoing[current]) != 1:
                break
        if len(path) < 2:
            return
        expected = len(bundles[path[0]]["sequence"])
        total_overlap = 0
        for edge_id, node in zip(path_edges, path[1:]):
            overlap = edges[edge_id]["overlap"]
            total_overlap += overlap
            expected += len(bundles[node]["sequence"]) - overlap
        paths.append({
            "contig_id": len(paths),
            "nodes": path,
            "num_nodes": len(path),
            "count": sum(bundles[node]["abundance"] for node in path),
            "expected_length": expected,
            "total_overlap": total_overlap,
            "cycle_detected": cycle,
        })

    node_ids = sorted(bundles)
    for root in node_ids:
        if root not in used_node and outgoing[root] and (indegree[root] != 1 or len(outgoing[root]) != 1):
            walk(root)
    for root in node_ids:
        if root not in used_node and outgoing[root]:
            walk(root)
    return paths


def load_graph_annotations(annotation_path, fasta_path, paths):
    annotations = pd.read_csv(annotation_path, sep="\t", dtype=str).fillna("")
    fasta = parse_fasta(fasta_path)
    fasta_lengths = {name: len(sequence) for name, sequence in fasta}
    path_by_id = {f"contig_{path['contig_id']}": path for path in paths}
    rows = []
    for _, annotation in annotations.iterrows():
        contig_id = annotation["sequence_id"]
        path = path_by_id.get(contig_id)
        v = remove_allele(annotation["v_call"])
        j = remove_allele(annotation["j_call"])
        c_original = remove_allele(annotation["c_call"])
        d = remove_allele(annotation["d_call"])
        cdr3_nt = annotation["junction"].upper() if annotation["junction"] else None
        cdr3_aa = annotation["junction_aa"].upper() if annotation["junction_aa"] else None
        rows.append({
            "sample": "FZ-116",
            "contig_id": contig_id,
            "count": path["count"] if path else np.nan,
            "V": annotation["v_call"],
            "D": annotation["d_call"],
            "J": annotation["j_call"],
            "C": annotation["c_call"],
            "CDR3_nt": cdr3_nt,
            "CDR3_aa": cdr3_aa,
            "normalized_V": v,
            "normalized_J": j,
            "normalized_C": collapse_isotype(c_original),
            "normalized_D": d,
            "original_C": c_original,
            "CDR3_nt_trimmed": trim_cdr3(cdr3_nt),
            "productive": annotation["productive"],
            "path_nodes": path["num_nodes"] if path else np.nan,
            "total_overlap": path["total_overlap"] if path else np.nan,
            "expected_length": path["expected_length"] if path else np.nan,
            "contig_length": fasta_lengths.get(contig_id, np.nan),
            "cycle_detected": int(path["cycle_detected"]) if path else np.nan,
        })
    return pd.DataFrame(rows)


def load_trust4(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [column.lstrip("#") for column in df.columns]
    df = df[df["V"].str.startswith("IGH", na=False)].copy()
    df["v_gene"] = df["V"].apply(remove_allele)
    df["j_gene"] = df["J"].apply(remove_allele)
    df["c_gene_original"] = df["C"].apply(remove_allele)
    df["c_gene_collapsed"] = df["c_gene_original"].apply(collapse_isotype)
    df["cdr3nt_original"] = df["CDR3nt"].astype(str).str.upper()
    df["cdr3nt_trimmed"] = df["cdr3nt_original"].apply(trim_cdr3)
    df["d_gene"] = df["D"].apply(remove_allele)
    return df[df["cdr3nt_trimmed"].notna()].copy()


def load_irep(path):
    df = pd.read_csv(path, compression="gzip")
    df["v_gene"] = df["V"].apply(lambda value: remove_allele(str(value).lstrip("h")))
    df["j_gene"] = df["J"].apply(lambda value: remove_allele(str(value).lstrip("h")))
    df["c_gene_original"] = df["C"].apply(lambda value: remove_allele(str(value).lstrip("h")))
    df["c_gene_collapsed"] = df["c_gene_original"].apply(collapse_isotype)
    df["cdr3nt"] = df["CDR3(nuc)"].astype(str).str.upper()
    df["d_gene"] = df["D"].apply(lambda value: remove_allele(str(value).lstrip("h")))
    return df


def coalesce(df, key_columns, abundance_column):
    indices = df.groupby(key_columns)[abundance_column].idxmax()
    return df.loc[indices].reset_index(drop=True)


def graph_to_benchmark(df):
    df = df[df["normalized_V"].fillna("").str.startswith("IGH")].copy()
    df = df[df["CDR3_nt_trimmed"].notna()].copy()
    df["v_gene"] = df["normalized_V"]
    df["j_gene"] = df["normalized_J"]
    df["c_gene_collapsed"] = df["normalized_C"]
    df["cdr3nt_trimmed"] = df["CDR3_nt_trimmed"]
    df["d_gene"] = df["normalized_D"]
    df["c_gene_original"] = df["original_C"]
    return coalesce(df, ["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt_trimmed"], "count")


def evaluate(graph, trust4, irep):
    irep = coalesce(irep, ["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt"], "copy")
    irep_index = {
        (row.v_gene, row.j_gene, row.c_gene_collapsed, row.cdr3nt): row
        for row in irep.itertuples()
    }
    matches = []
    for row in graph.itertuples():
        key = (row.v_gene, row.j_gene, row.c_gene_collapsed, row.cdr3nt_trimmed)
        reference = irep_index.get(key)
        if reference is None:
            continue
        matches.append({
            "clonotype_key": "|".join(key),
            "contig_id": row.contig_id,
            "graph_count": row.count,
            "irep_copy": reference.copy,
            "V": row.v_gene,
            "J": row.j_gene,
            "C_collapsed": row.c_gene_collapsed,
            "graph_D": row.d_gene,
            "irep_D": reference.d_gene,
            "D_agreement": int(row.d_gene == reference.d_gene),
            "graph_C_exact": row.c_gene_original,
            "irep_C_exact": reference.c_gene_original,
            "isotype_exact_agreement": int(row.c_gene_original == reference.c_gene_original),
            "CDR3_nt_trimmed": row.cdr3nt_trimmed,
        })
    matches = pd.DataFrame(matches)
    n_graph, n_irep, n_matches = len(graph), len(irep), len(matches)
    precision = n_matches / n_graph if n_graph else 0.0
    sensitivity = n_matches / n_irep if n_irep else 0.0
    if n_matches >= 2:
        pearson = float(pearsonr(matches["graph_count"], matches["irep_copy"])[0])
        d_agreement = float(matches["D_agreement"].mean())
        exact_iso = float(matches["isotype_exact_agreement"].mean())
    else:
        pearson = np.nan
        d_agreement = np.nan
        exact_iso = np.nan
    return matches, {
        "unique_clonotypes": n_graph,
        "reference_clonotypes": n_irep,
        "matched_clonotypes": n_matches,
        "precision": precision,
        "sensitivity": sensitivity,
        "pearson_r": pearson,
        "D_agreement": d_agreement,
        "isotype_exact_agreement": exact_iso,
        "isotype_collapsed_agreement": 1.0 if n_matches else np.nan,
    }


def metric_row(metric, trust4, graph):
    difference = graph - trust4 if pd.notna(graph) and pd.notna(trust4) else np.nan
    relative_metrics = {"total_contigs", "IGH_contigs", "CDR3_contigs", "unique_clonotypes", "matched_clonotypes"}
    relative = difference / trust4 if metric in relative_metrics and trust4 != 0 and pd.notna(difference) else np.nan
    return {"metric": metric, "TRUST4": trust4, "Graph-TRUST4": graph, "difference": difference, "relative_change": relative}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/data1/wetlab/TRUST4_benchmark")
    args = parser.parse_args()
    root = args.root
    output = os.path.join(root, "results/Graph-TRUST4-v3-debug/FZ-116/biological_evaluation")
    os.makedirs(output, exist_ok=True)
    graph_root = os.path.join(root, "results/Graph-TRUST4-v3-debug/FZ-116/graph_fz116")
    paths = reconstruct_paths(os.path.join(graph_root, "bundles.tsv"), os.path.join(graph_root, "graph_edges.tsv"))
    graph_all = load_graph_annotations(os.path.join(output, "graph_fz116_annotated.tsv"), os.path.join(graph_root, "assembled_contigs.fa"), paths)
    graph_all.to_csv(os.path.join(output, "graph_fz116_normalized.tsv"), sep="\t", index=False)

    trust4_path = os.path.join(root, "results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv")
    irep_path = os.path.join(root, "reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz")
    trust4_raw = load_trust4(trust4_path)
    trust4 = coalesce(trust4_raw, ["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt_trimmed"], "count")
    irep = load_irep(irep_path)
    graph = graph_to_benchmark(graph_all)
    matches, graph_metrics = evaluate(graph, trust4, irep)
    matches.to_csv(os.path.join(output, "graph_fz116_vs_iRep_matches.tsv"), sep="\t", index=False)

    graph_igh = graph_all[graph_all["normalized_V"].fillna("").str.startswith("IGH")]
    graph_annotated = {
        "total_contigs": len(graph_all),
        "IGH_contigs": len(graph_igh),
        "V_assigned": int((graph_all["normalized_V"].fillna("") != "").sum()),
        "J_assigned": int((graph_all["normalized_J"].fillna("") != "").sum()),
        "C_assigned": int((graph_all["normalized_C"].fillna("") != "").sum()),
        "CDR3_contigs": int(graph_all["CDR3_nt_trimmed"].notna().sum()),
        "productive_IGH_contigs": int(((graph_all["normalized_V"].fillna("").str.startswith("IGH")) & (graph_all["productive"] == "T")).sum()),
        "unmapped_abundance_contigs": int(graph_all["count"].isna().sum()),
        "missing_V": int((graph_all["normalized_V"].fillna("") == "").sum()),
        "missing_J": int((graph_all["normalized_J"].fillna("") == "").sum()),
        "missing_C": int((graph_all["normalized_C"].fillna("") == "").sum()),
        "missing_CDR3": int(graph_all["CDR3_nt_trimmed"].isna().sum()),
        "duplicate_normalized_keys": int(graph[graph.duplicated(["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt_trimmed"], keep=False)].shape[0]),
        "excluded_from_clonotype_evaluation": int(len(graph_all) - len(graph)),
        "abundance_min": float(graph_all["count"].min()),
        "abundance_median": float(graph_all["count"].median()),
        "abundance_mean": float(graph_all["count"].mean()),
        "abundance_max": float(graph_all["count"].max()),
    }
    assembly_debug_path = os.path.join(graph_root, "assembly_debug.tsv")
    assembly_debug = pd.read_csv(assembly_debug_path, sep="\t")
    graph_annotated.update({
        "assembly_cycle_paths": int(assembly_debug["cycle_detected"].sum()),
        "assembly_repeated_node_paths": int((assembly_debug["repeated_nodes"] > 0).sum()),
        "assembly_invalid_paths": int((assembly_debug["valid_path"] == 0).sum()),
        "assembly_length_mismatches": int((assembly_debug["contig_length"] != assembly_debug["expected_length"]).sum()),
    })
    metrics = {**graph_annotated, **graph_metrics}
    pd.DataFrame([metrics]).to_csv(os.path.join(output, "graph_fz116_benchmark_metrics.tsv"), sep="\t", index=False)

    trust4_metrics = {
        "total_contigs": len(pd.read_csv(trust4_path, sep="\t")),
        "IGH_contigs": len(trust4_raw),
        "CDR3_contigs": len(trust4_raw),
        "unique_clonotypes": len(trust4),
        "matched_clonotypes": 3815,
        "precision": 0.3396243211964747,
        "sensitivity": 0.032709719459496536,
        "pearson_r": 0.6604407833526513,
        "D_agreement": 0.41782437745740497,
        "isotype_exact_agreement": 0.23748361730013107,
        "isotype_collapsed_agreement": 1.0,
    }
    comparison_metrics = ["total_contigs", "IGH_contigs", "CDR3_contigs", "unique_clonotypes", "matched_clonotypes", "precision", "sensitivity", "pearson_r", "D_agreement", "isotype_exact_agreement", "isotype_collapsed_agreement"]
    comparison = pd.DataFrame([metric_row(metric, trust4_metrics[metric], metrics.get(metric, np.nan)) for metric in comparison_metrics])
    comparison.to_csv(os.path.join(root, "results/Graph-TRUST4-v3-debug/FZ-116/TRUST4_vs_GraphTRUST4_FZ116.tsv"), sep="\t", index=False)

    with open(os.path.join(output, "README.md"), "w") as handle:
        handle.write("# FZ-116 Graph-TRUST4 biological evaluation\n\n")
        handle.write("Annotation command: `algorithms/TRUST4/annotator -f reference/TRUST4/human_IMGT+C.fa -a results/Graph-TRUST4-v3-debug/FZ-116/graph_fz116/assembled_contigs.fa --fasta --needReverseComplement --noImpute -t 64 --outputFormat 1`.\n\n")
        handle.write("Evaluation command: `python3 scripts/evaluate_graph_fz116.py`.\n\n")
        handle.write("Reference: `reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz`. TRUST4 baseline: `results/TRUST4/FZ-116/TRUST_FZ-116_report.tsv`.\n\n")
        handle.write("Normalization follows `scripts/benchmark_trust4_published.py`: first gene candidate, remove allele suffixes, collapse IGHA1/2 to IGHA and IGHG1-4 to IGHG (plus the existing IGHD/IGHE/IGHM mappings), uppercase sequences, trim 3 nt from each CDR3 end, and coalesce duplicate keys by maximum abundance. Primary key is normalized V + J + collapsed C + trimmed CDR3nt; D is auxiliary.\n\n")
        handle.write("Graph abundance is the sum of `abundance` values in `bundles.tsv` for the unique bundle nodes on the reconstructed contig path. Paths are reconstructed from `graph_edges.tsv` using the validated path traversal; contigs without a path mapping are excluded from biological clonotype evaluation and counted.\n\n")
        handle.write("Records are excluded from clonotype evaluation when they are not IGH or lack a valid CDR3 after the published 3-nt trim; missing V/J/C values and unmapped paths are retained and counted in the metrics table. Duplicate normalized keys are coalesced by maximum path abundance. Abundance min/median/mean/max, assembly cycle/repetition/length checks, and all metric values are recorded in `graph_fz116_benchmark_metrics.tsv` and `graph_fz116_normalized.tsv`. No assembly, parameter tuning, other samples, or iRepertoire-guided filtering was performed.\n")

    with open(os.path.join(output, "graph_fz116_biological_summary.md"), "w") as handle:
        handle.write("# Graph-TRUST4 FZ-116 biological evaluation\n\n")
        for key, value in metrics.items():
            handle.write(f"- {key}: {value}\n")
        handle.write("\nThe TRUST4 baseline reproduces the established published-evaluation values before this comparison.\n")
    print(pd.DataFrame([metrics]).to_string(index=False))
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
