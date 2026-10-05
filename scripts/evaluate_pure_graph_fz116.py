#!/usr/bin/env python3
"""Evaluate frozen Graph-TRUST4-pure FZ-116 outputs against the published benchmark."""

import argparse
import csv
import hashlib
from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr


TRUST4_BASELINE = {
    "unique_clonotypes": 11233,
    "reference_clonotypes": 116632,
    "matched_clonotypes": 3815,
    "precision": 0.339624321196,
    "sensitivity": 0.032709719459,
    "pearson_r": 0.660440783353,
}

OUTPUT_NAMES = (
    "pure_annotation_yield.tsv",
    "pure_benchmark_metrics.tsv",
    "pure_vs_iRep_matches.tsv",
    "pure_vs_TRUST4.tsv",
    "PURE_GRAPH_FZ116_REPORT.md",
)


def remove_allele(value):
    if pd.isna(value):
        return None
    return str(value).split("*")[0]


def collapse_isotype(value):
    groups = {
        "IGHA1": "IGHA",
        "IGHA2": "IGHA",
        "IGHG1": "IGHG",
        "IGHG2": "IGHG",
        "IGHG3": "IGHG",
        "IGHG4": "IGHG",
        "IGHD1": "IGHD",
        "IGHD2": "IGHD",
        "IGHD3": "IGHD",
        "IGHD4": "IGHD",
        "IGHE1": "IGHE",
        "IGHE2": "IGHE",
        "IGHM1": "IGHM",
        "IGHM2": "IGHM",
    }
    return groups.get(value, value)


def is_assigned(value):
    return str(value).strip() not in ("", ".", "NA", "nan", "None")


def trim_cdr3(value):
    if not is_assigned(value):
        return None
    sequence = str(value).upper()
    return sequence[3:-3] if len(sequence) >= 6 else None


def read_stats(path):
    with path.open(newline="") as handle:
        return {row["statistic"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def fasta_summary(path):
    lengths = []
    sequence_hashes = {}
    current_id = None
    sequence_parts = []

    def save_record():
        if current_id is None:
            return
        sequence = "".join(sequence_parts).upper()
        if current_id in sequence_hashes:
            raise ValueError(f"Duplicate FASTA identifier: {current_id}")
        complement = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
        reverse_complement = sequence.translate(complement)[::-1]
        sequence_hashes[current_id] = (
            len(sequence),
            hashlib.sha256(sequence.encode()).hexdigest(),
            hashlib.sha256(reverse_complement.encode()).hexdigest(),
        )
        lengths.append(len(sequence))

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                save_record()
                current_id = line[1:].split()[0]
                sequence_parts = []
            elif line:
                if current_id is None:
                    raise ValueError(f"Sequence encountered before FASTA header in {path}")
                sequence_parts.append(line)
    save_record()

    total_bases = sum(lengths)
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= total_bases:
            n50 = length
            break
    return sequence_hashes, total_bases, n50


def graph_fasta_summary(path):
    """Read TRUST4's existing header/sequence output without interpreting its coverage rows."""
    lengths = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                sequence = next(handle, "").strip()
                if not sequence or sequence.startswith(">"):
                    raise ValueError(f"Missing TRUST4 sequence after header in {path}")
                lengths.append(len(sequence))
    total_bases = sum(lengths)
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative * 2 >= total_bases:
            n50 = length
            break
    return len(lengths), total_bases, n50


def load_graph_nodes(path):
    abundance = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            node_id = row["node_id"]
            if node_id in abundance:
                raise ValueError(f"Duplicate graph node: {node_id}")
            abundance[node_id] = int(row["abundance"])
    return abundance


def validate_node_read_provenance(path, node_abundance):
    read_rows_by_node = {}
    total_rows = 0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            node_id = row["node_id"]
            if node_id not in node_abundance:
                raise ValueError(f"node_reads.tsv references unknown node {node_id}")
            if not is_assigned(row["raw_read_id"]):
                raise ValueError(f"Empty raw read identifier for graph node {node_id}")
            read_rows_by_node[node_id] = read_rows_by_node.get(node_id, 0) + 1
            total_rows += 1
    if set(read_rows_by_node) != set(node_abundance):
        raise ValueError("node_reads.tsv does not represent every graph node")
    if max(read_rows_by_node.values(), default=0) > 64:
        raise ValueError("node_reads.tsv exceeds the source implementation's 64-ID per-node cap")
    return total_rows, max(read_rows_by_node.values(), default=0)


def load_path_abundances(path, node_abundance, fasta_ids):
    abundance_by_contig = {}
    path_id_by_contig = {}
    node_count_by_contig = {}
    current_contig = None
    current_path = None
    current_nodes = set()
    current_abundance = 0
    current_step = 0
    expected_next_node = None
    terminal_seen = False

    def finish_path():
        if current_contig is None:
            return
        fasta_id = f"contig_{current_contig}"
        if fasta_id not in fasta_ids:
            raise ValueError(f"Path references contig absent from FASTA: {fasta_id}")
        if fasta_id in abundance_by_contig:
            raise ValueError(f"Multiple graph paths emitted for contig {fasta_id}")
        abundance_by_contig[fasta_id] = current_abundance
        path_id_by_contig[fasta_id] = current_path
        node_count_by_contig[fasta_id] = len(current_nodes)

    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig_id = row["contig_id"]
            if contig_id != current_contig:
                finish_path()
                current_contig = contig_id
                current_path = row["path_id"]
                current_nodes = set()
                current_abundance = 0
                current_step = 0
                expected_next_node = None
                terminal_seen = False
                if current_path != contig_id:
                    raise ValueError(f"Path ID does not match contig ID {contig_id}")
            elif row["path_id"] != current_path:
                raise ValueError(f"Interleaved paths for contig {contig_id}")

            if terminal_seen:
                raise ValueError(f"Path {current_path} continues after its terminal step")
            if int(row["step"]) != current_step:
                raise ValueError(f"Invalid step sequence in path {current_path}")

            node_id = row["node_id"]
            if node_id not in node_abundance:
                raise ValueError(f"Path {current_path} references unknown node {node_id}")
            if node_id in current_nodes:
                raise ValueError(f"Path {current_path} repeats node {node_id}")
            if expected_next_node is not None and node_id != expected_next_node:
                raise ValueError(f"Broken node transition in path {current_path}")

            step_abundance = int(row["bundle_abundance"])
            if step_abundance != node_abundance[node_id]:
                raise ValueError(f"Bundle abundance mismatch for node {node_id}")
            current_nodes.add(node_id)
            current_abundance += step_abundance

            edge_id = row["edge_id"]
            next_node_id = row["next_node_id"]
            if edge_id == "NOT AVAILABLE":
                if next_node_id != "NOT AVAILABLE":
                    raise ValueError(f"Terminal step has a next node in path {current_path}")
                terminal_seen = True
                expected_next_node = None
            else:
                if next_node_id == "NOT AVAILABLE":
                    raise ValueError(f"Non-terminal edge lacks next node in path {current_path}")
                expected_next_node = next_node_id
            current_step += 1

    finish_path()
    if set(abundance_by_contig) != fasta_ids:
        missing = len(fasta_ids - set(abundance_by_contig))
        extra = len(set(abundance_by_contig) - fasta_ids)
        raise ValueError(f"Path/FASTA coverage mismatch: missing={missing}, extra={extra}")
    return abundance_by_contig, path_id_by_contig, node_count_by_contig


def load_annotations(path, fasta_hashes, path_abundance):
    annotations = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "sequence_id", "sequence", "v_call", "d_call", "j_call", "c_call",
        "junction", "junction_aa", "productive",
    }
    missing = required - set(annotations.columns)
    if missing:
        raise ValueError(f"AIRR annotation is missing columns: {sorted(missing)}")
    if annotations["sequence_id"].duplicated().any():
        raise ValueError("AIRR annotation contains duplicate sequence IDs")
    if set(annotations["sequence_id"]) != set(fasta_hashes):
        raise ValueError("AIRR sequence IDs do not exactly match assembled_contigs.fa")

    orientation_counts = {"same_orientation": 0, "reverse_complement": 0}
    complement = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
    for row in annotations.itertuples(index=False):
        fasta_length, fasta_hash, reverse_hash = fasta_hashes[row.sequence_id]
        sequence = str(row.sequence).upper()
        if len(sequence) != fasta_length:
            raise ValueError(f"AIRR sequence does not match FASTA record {row.sequence_id}")
        sequence_hash = hashlib.sha256(sequence.encode()).hexdigest()
        if sequence_hash == fasta_hash:
            orientation_counts["same_orientation"] += 1
        elif sequence_hash == reverse_hash:
            orientation_counts["reverse_complement"] += 1
        else:
            raise ValueError(f"AIRR sequence does not match FASTA record {row.sequence_id}")

    annotations["v_gene"] = annotations["v_call"].map(remove_allele)
    annotations["j_gene"] = annotations["j_call"].map(remove_allele)
    annotations["c_gene_original"] = annotations["c_call"].map(remove_allele)
    annotations["c_gene_collapsed"] = annotations["c_gene_original"].map(collapse_isotype)
    annotations["d_gene"] = annotations["d_call"].map(remove_allele)
    annotations["cdr3nt_original"] = annotations["junction"].str.upper()
    annotations["cdr3nt_trimmed"] = annotations["cdr3nt_original"].map(trim_cdr3)
    annotations["abundance"] = annotations["sequence_id"].map(path_abundance)
    if annotations["abundance"].isna().any():
        raise ValueError("Some annotated contigs lack graph-derived path abundance")
    return annotations, orientation_counts


def coalesce(df, fields, abundance):
    usable = df.dropna(subset=fields).copy()
    usable = usable[usable[fields].apply(lambda col: col.astype(str).str.len() > 0).all(axis=1)]
    if usable.empty:
        return usable
    selected = usable.groupby(fields, sort=False, dropna=True)[abundance].idxmax()
    return usable.loc[selected].reset_index(drop=True)


def load_reference(path):
    reference = pd.read_csv(path, compression="gzip")
    required = {"V", "D", "J", "C", "CDR3(nuc)", "copy"}
    missing = required - set(reference.columns)
    if missing:
        raise ValueError(f"iRepertoire reference is missing columns: {sorted(missing)}")
    reference["v_gene"] = reference["V"].map(lambda value: remove_allele(str(value).lstrip("h")))
    reference["j_gene"] = reference["J"].map(lambda value: remove_allele(str(value).lstrip("h")))
    reference["c_gene_original"] = reference["C"].map(lambda value: remove_allele(str(value).lstrip("h")))
    reference["c_gene_collapsed"] = reference["c_gene_original"].map(collapse_isotype)
    reference["cdr3nt"] = reference["CDR3(nuc)"].astype(str).str.upper()
    reference["d_gene"] = reference["D"].map(lambda value: remove_allele(str(value).lstrip("h")))
    reference["copy"] = pd.to_numeric(reference["copy"], errors="raise")
    return reference


def load_trust4_annotations(path):
    trust4 = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"sequence_id", "v_call", "productive"}
    missing = required - set(trust4.columns)
    if missing:
        raise ValueError(f"TRUST4 AIRR output is missing columns: {sorted(missing)}")
    v_gene = trust4["v_call"].map(remove_allele)
    igh = v_gene.fillna("").str.startswith("IGH")
    productive = trust4["productive"].str.upper().isin(("T", "TRUE", "1"))
    return len(trust4), int(igh.sum()), int((igh & productive).sum())


def collect_matched_path_nodes(path, matched_contigs):
    node_ids_by_contig = {contig_id: [] for contig_id in matched_contigs}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig_id = f"contig_{row['contig_id']}"
            if contig_id in node_ids_by_contig:
                node_ids_by_contig[contig_id].append(row["node_id"])
    return node_ids_by_contig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz"),
    )
    parser.add_argument(
        "--trust4-airr",
        type=Path,
        default=Path("results/TRUST4/FZ-116/TRUST_FZ-116_airr.tsv"),
    )
    parser.add_argument(
        "--trust4-final",
        type=Path,
        default=Path("results/TRUST4/FZ-116/TRUST_FZ-116_final.out"),
    )
    parser.add_argument(
        "--trust4-published-metrics",
        type=Path,
        default=Path("results/benchmark/FZ-116/published_eval/metrics.tsv"),
    )
    args = parser.parse_args()
    output = args.output

    existing = [name for name in OUTPUT_NAMES if (output / name).exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing evaluation output(s): " + ", ".join(existing)
        )

    fasta_path = output / "assembled_contigs.fa"
    annotation_path = output / "pure_raw_airr.tsv"
    nodes_path = output / "graph_nodes.tsv"
    node_reads_path = output / "node_reads.tsv"
    paths_path = output / "contig_paths.tsv"
    graph_stats_path = output / "graph_statistics.tsv"
    assembly_stats_path = output / "assembly_statistics.tsv"
    qaos_stats_path = output / "QAOS_statistics.tsv"

    fasta_hashes, fasta_bases, fasta_n50 = fasta_summary(fasta_path)
    assembly_stats = read_stats(assembly_stats_path)
    if int(assembly_stats["contigs"]) != len(fasta_hashes):
        raise ValueError("FASTA record count disagrees with assembly_statistics.tsv")
    if int(assembly_stats["assembled_bases"]) != fasta_bases:
        raise ValueError("FASTA base count disagrees with assembly_statistics.tsv")
    if int(assembly_stats["N50"]) != fasta_n50:
        raise ValueError("FASTA N50 disagrees with assembly_statistics.tsv")

    graph_stats = read_stats(graph_stats_path)
    node_abundance = load_graph_nodes(nodes_path)
    if int(graph_stats["bundles"]) != len(node_abundance):
        raise ValueError("graph_nodes.tsv row count disagrees with graph_statistics.tsv")
    if sum(node_abundance.values()) != int(graph_stats["raw_reads"]):
        raise ValueError("Graph bundle abundance sum disagrees with raw_reads")
    node_read_rows, max_sampled_read_ids = validate_node_read_provenance(
        node_reads_path, node_abundance
    )
    path_abundance, path_ids, path_node_counts = load_path_abundances(
        paths_path, node_abundance, set(fasta_hashes)
    )
    if len(path_abundance) != int(graph_stats["contigs"]):
        raise ValueError("contig_paths.tsv path count disagrees with graph_statistics.tsv")

    annotations, annotation_orientations = load_annotations(
        annotation_path, fasta_hashes, path_abundance
    )
    reference = load_reference(args.reference)
    igh_mask = annotations["v_gene"].fillna("").str.startswith("IGH")
    primary_fields = ["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt_trimmed"]
    reference_fields = ["v_gene", "j_gene", "c_gene_collapsed", "cdr3nt"]
    candidates = coalesce(
        annotations.loc[igh_mask & annotations["cdr3nt_trimmed"].notna()],
        primary_fields,
        "abundance",
    )
    reference_unique = coalesce(reference, reference_fields, "copy")
    if len(reference_unique) != TRUST4_BASELINE["reference_clonotypes"]:
        raise ValueError(
            "Reference clonotype count does not match the frozen published TRUST4 baseline: "
            f"{len(reference_unique)}"
        )

    reference_index = {
        tuple(row[field] for field in reference_fields): row
        for _, row in reference_unique.iterrows()
    }
    matched = []
    for _, row in candidates.iterrows():
        key = tuple(row[field] for field in primary_fields)
        reference_row = reference_index.get(key)
        if reference_row is None:
            continue
        matched.append(
            {
                "contig_id": row["sequence_id"],
                "path_id": path_ids[row["sequence_id"]],
                "path_node_count": path_node_counts[row["sequence_id"]],
                "V": row["v_gene"],
                "D": row["d_gene"],
                "J": row["j_gene"],
                "C": row["c_gene_original"],
                "C_collapsed": row["c_gene_collapsed"],
                "CDR3_nt": row["cdr3nt_original"],
                "CDR3_nt_trimmed": row["cdr3nt_trimmed"],
                "graph_abundance": int(row["abundance"]),
                "iRep_copy": int(reference_row["copy"]),
                "D_agreement": row["d_gene"] == reference_row["d_gene"],
                "exact_isotype_agreement": (
                    row["c_gene_original"] == reference_row["c_gene_original"]
                ),
            }
        )

    matches = pd.DataFrame(
        matched,
        columns=[
            "contig_id", "path_id", "path_node_count", "V", "D", "J", "C",
            "C_collapsed", "CDR3_nt", "CDR3_nt_trimmed", "graph_abundance",
            "iRep_copy", "D_agreement", "exact_isotype_agreement",
        ],
    )
    path_nodes = collect_matched_path_nodes(paths_path, set(matches["contig_id"]))
    if not matches.empty:
        matches["path_node_ids"] = matches["contig_id"].map(
            lambda contig_id: ",".join(path_nodes[contig_id])
        )
        if (matches["path_node_count"] != matches["contig_id"].map(
            lambda contig_id: len(path_nodes[contig_id])
        )).any():
            raise ValueError("Matched contig path node count does not match provenance")
        matches = matches[
            [
                "contig_id", "path_id", "path_node_count", "path_node_ids", "V", "D",
                "J", "C", "C_collapsed", "CDR3_nt", "CDR3_nt_trimmed",
                "graph_abundance", "iRep_copy", "D_agreement",
                "exact_isotype_agreement",
            ]
        ]

    match_count = len(matches)
    pure_unique_count = len(candidates)
    reference_count = len(reference_unique)
    precision = match_count / pure_unique_count if pure_unique_count else "NOT AVAILABLE"
    sensitivity = match_count / reference_count if reference_count else "NOT AVAILABLE"
    pearson = "NOT AVAILABLE"
    if match_count >= 2:
        graph_counts = pd.to_numeric(matches["graph_abundance"], errors="raise")
        irep_copies = pd.to_numeric(matches["iRep_copy"], errors="raise")
        if graph_counts.nunique() > 1 and irep_copies.nunique() > 1:
            pearson = float(pearsonr(graph_counts, irep_copies).statistic)

    annotation_metrics = [
        ("total_contigs", pd.Series(True, index=annotations.index)),
        ("IGH_contigs", igh_mask),
        ("V_assigned", annotations["v_call"].map(is_assigned)),
        ("D_assigned", annotations["d_call"].map(is_assigned)),
        ("J_assigned", annotations["j_call"].map(is_assigned)),
        ("C_assigned", annotations["c_call"].map(is_assigned)),
        ("CDR3_assigned", annotations["junction"].map(is_assigned)),
        ("productive_IGH", igh_mask & annotations["productive"].str.upper().eq("T")),
    ]
    yield_rows = [
        {
            "metric": metric,
            "count": int(mask.sum()),
            "percent_all": 100.0 * float(mask.mean()),
        }
        for metric, mask in annotation_metrics
    ]

    benchmark_metrics = [
        ("total_contigs", len(annotations)),
        ("IGH_contigs", int(igh_mask.sum())),
        ("V_assigned", int(annotations["v_call"].map(is_assigned).sum())),
        ("D_assigned", int(annotations["d_call"].map(is_assigned).sum())),
        ("J_assigned", int(annotations["j_call"].map(is_assigned).sum())),
        ("C_assigned", int(annotations["c_call"].map(is_assigned).sum())),
        ("CDR3_assigned", int(annotations["junction"].map(is_assigned).sum())),
        (
            "productive_IGH",
            int((igh_mask & annotations["productive"].str.upper().eq("T")).sum()),
        ),
        ("pure_unique_clonotypes", pure_unique_count),
        ("reference_clonotypes", reference_count),
        ("matched_clonotypes", match_count),
        ("precision", precision),
        ("sensitivity", sensitivity),
        ("pearson_r", pearson),
        ("abundance_method", "sum of graph bundle abundances along each path"),
    ]

    trust4_contigs, trust4_bases, trust4_n50 = graph_fasta_summary(args.trust4_final)
    trust4_airr_records, trust4_igh, trust4_productive_igh = load_trust4_annotations(
        args.trust4_airr
    )
    with args.trust4_published_metrics.open(newline="") as handle:
        frozen_metric_row = next(csv.DictReader(handle, delimiter="\t"))
    for key, expected in TRUST4_BASELINE.items():
        if key in ("unique_clonotypes", "reference_clonotypes", "matched_clonotypes"):
            observed = int(frozen_metric_row[{
                "unique_clonotypes": "n_trust4",
                "reference_clonotypes": "n_irep",
                "matched_clonotypes": "n_matches",
            }[key]])
            if observed != expected:
                raise ValueError(f"Frozen TRUST4 baseline mismatch for {key}: {observed}")
        elif key in ("precision", "sensitivity") and abs(
            float(frozen_metric_row[key]) - expected
        ) > 1e-12:
            raise ValueError(f"Frozen TRUST4 baseline mismatch for {key}")

    pure_values = {
        "contigs": len(fasta_hashes),
        "N50": fasta_n50,
        "IGH_contigs": int(igh_mask.sum()),
        "productive_IGH": int(
            (igh_mask & annotations["productive"].str.upper().eq("T")).sum()
        ),
        "unique_clonotypes": pure_unique_count,
        "matches": match_count,
        "precision": precision,
        "sensitivity": sensitivity,
        "Pearson r": pearson,
    }
    trust4_values = {
        "contigs": trust4_contigs,
        "N50": trust4_n50,
        "IGH_contigs": trust4_igh,
        "productive_IGH": trust4_productive_igh,
        "unique_clonotypes": TRUST4_BASELINE["unique_clonotypes"],
        "matches": TRUST4_BASELINE["matched_clonotypes"],
        "precision": TRUST4_BASELINE["precision"],
        "sensitivity": TRUST4_BASELINE["sensitivity"],
        "Pearson r": TRUST4_BASELINE["pearson_r"],
    }
    comparison_rows = []
    for metric, trust4_value in trust4_values.items():
        pure_value = pure_values[metric]
        difference = (
            pure_value - trust4_value
            if isinstance(pure_value, (int, float))
            and isinstance(trust4_value, (int, float))
            else "NOT AVAILABLE"
        )
        comparison_rows.append(
            {
                "metric": metric,
                "TRUST4": trust4_value,
                "Graph-TRUST4-pure": pure_value,
                "difference": difference,
            }
        )

    metrics_path = output / "pure_benchmark_metrics.tsv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(benchmark_metrics)
    pd.DataFrame(yield_rows).to_csv(
        output / "pure_annotation_yield.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    matches.to_csv(
        output / "pure_vs_iRep_matches.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    pd.DataFrame(comparison_rows).to_csv(
        output / "pure_vs_TRUST4.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )

    qaos_stats = read_stats(qaos_stats_path)
    metadata = {}
    with (output / "run_metadata.txt").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) == 2:
                metadata[fields[0]] = fields[1]

    report = f"""# PURE Graph FZ-116 Biological Evaluation

## Experiment status

- Assembly reused from the frozen output; Graph-TRUST4-pure was not rerun.
- Annotation reused from `pure_raw_airr.tsv`; no annotation rerun was needed.
- All {len(fasta_hashes):,} FASTA IDs have exactly one AIRR record. Every AIRR sequence matches its FASTA sequence in the same orientation or as its reverse complement ({annotation_orientations["same_orientation"]:,} same-orientation; {annotation_orientations["reverse_complement"]:,} reverse-complement).
- Matching uses the established published TRUST4 key: trimmed CDR3 nucleotide + normalized V + normalized J + normalized/collapsed C.
- CDR3 nucleotide normalization removes the first and last 3 nt from the TRUST4/AIRR junction. V/J allele suffixes are removed; constant-region subclasses are collapsed as in the published benchmark.
- This does not use the earlier CDR3 amino-acid comparison that produced zero matches.

## COMPUTATIONAL ASSEMBLY RESULT

| Metric | Graph-TRUST4-pure |
|---|---:|
| Contigs | {len(fasta_hashes):,} |
| N50 | {fasta_n50:,} bp |
| Assembled bases | {fasta_bases:,} |
| Graph nodes / bundles | {len(node_abundance):,} |
| Accepted edges | {int(graph_stats["accepted_edges"]):,} |
| Paths | {int(graph_stats["contigs"]):,} |
| Branches | {int(graph_stats["branches"]):,} |
| Components | {int(graph_stats["components"]):,} |
| Cycles | {int(graph_stats["cycles"]):,} |
| Raw read records | {int(graph_stats["raw_reads"]):,} |
| Mean overlap | {float(qaos_stats["mean_overlap"]):.4f} bp |
| Mean identity | {float(qaos_stats["mean_identity"]):.6f} |
| Mean QAOS | {float(qaos_stats["mean_QAOS"]):.6f} |

Frozen parameters include k={metadata.get("k")}, minimum overlap={metadata.get("minimum_overlap")}, minimum identity={metadata.get("minimum_identity")}, minimum QAOS={metadata.get("minimum_QAOS")}, Q_HIGH={metadata.get("Q_HIGH")}, candidate cap={metadata.get("candidate_cap")}, and {metadata.get("threads")} threads. Assembly wall time was {float(metadata["wall_time_seconds"]):.1f} seconds; peak RSS was {metadata.get("peak_RSS_MB")} MB.

## Annotation yield and completeness

`pure_raw_airr.tsv` has {len(annotations):,} unique annotation records for {len(fasta_hashes):,} contigs. Every sequence matches its FASTA record either directly or by reverse complement.

| Metric | Count | Percent of contigs |
|---|---:|---:|
"""
    report += "".join(
        f"| {row['metric']} | {row['count']:,} | {row['percent_all']:.3f}% |\n"
        for row in yield_rows
    )
    report += f"""
## BIOLOGICAL RECONSTRUCTION RESULT

| Metric | Graph-TRUST4-pure |
|---|---:|
| Pure unique clonotypes | {pure_unique_count:,} |
| iRepertoire reference clonotypes | {reference_count:,} |
| Matched clonotypes | {match_count:,} |
| Precision | {precision:.6f} |
| Sensitivity (primary recall measure) | {sensitivity:.6f} |
| Pearson abundance correlation | {pearson if isinstance(pearson, str) else f"{pearson:.6f}"} |

The iRepertoire reference has {reference_count:,} unique clonotypes under the frozen CDR3nt + V + J + C benchmark key. No extra CDR3aa-stop filter was applied; this preserves the denominator used by the published TRUST4 comparison.

## Direct comparison with frozen TRUST4 baseline

| Metric | TRUST4 | Graph-TRUST4-pure | Difference (pure - TRUST4) |
|---|---:|---:|---:|
"""
    for row in comparison_rows:
        trust4_value = row["TRUST4"]
        pure_value = row["Graph-TRUST4-pure"]
        difference = row["difference"]
        if isinstance(trust4_value, float):
            trust4_text = f"{trust4_value:.6f}"
            pure_text = (
                pure_value if isinstance(pure_value, str) else f"{pure_value:.6f}"
            )
            difference_text = (
                difference if isinstance(difference, str) else f"{difference:.6f}"
            )
        else:
            trust4_text = f"{trust4_value:,}"
            pure_text = (
                pure_value if isinstance(pure_value, str) else f"{pure_value:,}"
            )
            difference_text = (
                difference if isinstance(difference, str) else f"{difference:,}"
            )
        report += f"| {row['metric']} | {trust4_text} | {pure_text} | {difference_text} |\n"

    report += f"""
The TRUST4 assembly figures come from its existing `TRUST_FZ-116_final.out` and AIRR outputs. Its published benchmark metrics were reused without rerunning TRUST4: {TRUST4_BASELINE["unique_clonotypes"]:,} unique clonotypes, {TRUST4_BASELINE["matched_clonotypes"]:,} matches, precision {TRUST4_BASELINE["precision"]:.6f}, sensitivity {TRUST4_BASELINE["sensitivity"]:.6f}, and Pearson r {TRUST4_BASELINE["pearson_r"]:.6f}.

## Abundance and graph provenance

Per-contig graph abundance is the sum of `bundle_abundance` over nodes in that contig's path. Every path-step abundance matches `graph_nodes.tsv`; the sum of graph-node abundances ({sum(node_abundance.values()):,}) matches the graph's {int(graph_stats["raw_reads"]):,} raw read records. `node_reads.tsv` has {node_read_rows:,} provenance rows, with at most {max_sampled_read_ids} stored raw read IDs per node. Therefore bundle counts are complete, while the listed raw read IDs are capped samples.

Matched clonotypes link to contigs and path-node IDs in `pure_vs_iRep_matches.tsv`; follow those node IDs in `contig_paths.tsv` and `node_reads.tsv`, with retained alternatives documented in `branch_decisions.tsv`. Alternative paths may share graph nodes and thus read support; graph abundance should be interpreted as path support, not as independent molecule counts.

## Limitations and interpretation

- Only {int(igh_mask.sum()):,} of {len(annotations):,} assembled contigs have an IGH V assignment; the biological benchmark evaluates the clonotypes that have the required key fields.
- The graph-derived abundance is a sum of constituent bundle read counts. Alternative paths can share nodes, so support across different contigs is not independent.
- Paired-read links are not represented as paired support in this experimental graph output.
- Contig count, N50, and assembled bases are computational assembly measures; they do not establish biological recovery.
- The primary recall comparison is sensitivity: Graph-TRUST4-pure recovers {match_count:,} reference clonotypes ({sensitivity:.6f}), versus TRUST4's frozen {TRUST4_BASELINE["matched_clonotypes"]:,} matches ({TRUST4_BASELINE["sensitivity"]:.6f} sensitivity). The graph result therefore does not show improved reference-clonotype recall, despite its larger contig count and N50.

The annotation-yield counts and full benchmark metrics are in `pure_annotation_yield.tsv` and `pure_benchmark_metrics.tsv`; side-by-side values are in `pure_vs_TRUST4.tsv`.
"""
    (output / "PURE_GRAPH_FZ116_REPORT.md").write_text(report)

    print(f"ASSEMBLY_REUSED=true")
    print(f"ASSEMBLY_RERUN=false")
    print()
    print(f"CONTIGS={len(fasta_hashes)}")
    print(f"N50={fasta_n50}")
    print(f"ASSEMBLED_BASES={fasta_bases}")
    print()
    print(f"IGH_CONTIGS={int(igh_mask.sum())}")
    print(f"V_ASSIGNED={int(annotations['v_call'].map(is_assigned).sum())}")
    print(f"J_ASSIGNED={int(annotations['j_call'].map(is_assigned).sum())}")
    print(f"C_ASSIGNED={int(annotations['c_call'].map(is_assigned).sum())}")
    print(f"CDR3_ASSIGNED={int(annotations['junction'].map(is_assigned).sum())}")
    print(
        "PRODUCTIVE_IGH="
        f"{int((igh_mask & annotations['productive'].str.upper().eq('T')).sum())}"
    )
    print()
    print(f"PURE_UNIQUE_CLONOTYPES={pure_unique_count}")
    print(f"REFERENCE_CLONOTYPES={reference_count}")
    print(f"MATCHES={match_count}")
    print(f"PRECISION={precision}")
    print(f"SENSITIVITY={sensitivity}")
    print(f"PEARSON_R={pearson}")
    print()
    print(f"TRUST4_UNIQUE={TRUST4_BASELINE['unique_clonotypes']}")
    print(f"TRUST4_MATCHES={TRUST4_BASELINE['matched_clonotypes']}")
    print(f"TRUST4_PRECISION={TRUST4_BASELINE['precision']}")
    print(f"TRUST4_SENSITIVITY={TRUST4_BASELINE['sensitivity']}")
    print(f"TRUST4_PEARSON_R={TRUST4_BASELINE['pearson_r']}")
    print()
    print("BIOLOGICAL_EVALUATION_COMPLETE=true")


if __name__ == "__main__":
    main()
