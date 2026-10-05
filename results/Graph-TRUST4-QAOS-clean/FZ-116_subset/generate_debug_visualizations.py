#!/usr/bin/env python3
"""Generate interactive raw-graph debugging views from an existing clean run."""

import csv
import itertools
import math
import random
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EDGE_FILE = ROOT / "graph_edges.tsv"
NODE_FILE = ROOT / "graph_nodes.tsv"
STATS_FILE = ROOT / "graph_statistics.tsv"
PATHS_FILE = ROOT / "graph_paths.bed"
PATH_STATS_FILE = ROOT / "path_statistics.tsv"
FASTA_FILE = ROOT / "assembled_contigs.fa"
BEDPE_FILE = ROOT / "graph_edges.bed"
PNG_DPI = 120
EDGE_COLORS = {
    ("0", "0"): "#277da1",
    ("0", "1"): "#43aa8b",
    ("1", "0"): "#f8961e",
    ("1", "1"): "#d1495b",
}


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header in {path}")
        return reader.fieldnames, list(reader)


class DisjointSet:
    def __init__(self, size):
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, node):
        parent = self.parent
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]


def read_nodes():
    names = {}
    with NODE_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"node_id", "read_id"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected graph_nodes.tsv schema: {reader.fieldnames}")
        for row in reader:
            node_id = int(row["node_id"])
            if node_id in names:
                raise ValueError(f"Duplicate node_id {node_id}")
            names[node_id] = row["read_id"]
    if sorted(names) != list(range(len(names))):
        raise ValueError("Node IDs are not contiguous zero-based integers")
    return names


def edge_headers():
    with EDGE_FILE.open(newline="") as handle:
        headers = csv.DictReader(handle, delimiter="\t").fieldnames
    required = {
        "edge_id", "source_read", "target_read", "source_orientation",
        "target_orientation", "overlap_length", "identity", "QAOS",
    }
    if not headers or not required.issubset(headers):
        raise ValueError(f"Unexpected graph_edges.tsv schema: {headers}")
    return headers


def edge_reader(names=None):
    with EDGE_FILE.open(newline="") as edge_handle, BEDPE_FILE.open(newline="") as bed_handle:
        reader = csv.DictReader(edge_handle, delimiter="\t")
        bed_reader = csv.reader(bed_handle, delimiter="\t")
        bed_header = next(bed_reader, None)
        if not bed_header or not bed_header[0].startswith("#"):
            raise ValueError("Unexpected graph_edges.bed header")
        for row, bed in itertools.zip_longest(reader, bed_reader):
            if row is None or bed is None:
                raise ValueError("graph_edges.tsv and graph_edges.bed row counts differ")
            if len(bed) != 10:
                raise ValueError(f"Unexpected graph_edges.bed row: {bed}")
            source_node = int(bed[0].removeprefix("read_"))
            target_node = int(bed[3].removeprefix("read_"))
            if (
                row["source_orientation"] != bed[7]
                or row["target_orientation"] != bed[8]
                or not math.isclose(float(row["identity"]), float(bed[6]), abs_tol=1e-5)
                or not math.isclose(float(row["QAOS"]), float(bed[9]), abs_tol=1e-5)
            ):
                raise ValueError(f"graph_edges.tsv and BED edge row disagree at edge {row['edge_id']}")
            if names is not None and (
                names[source_node] != row["source_read"]
                or names[target_node] != row["target_read"]
            ):
                raise ValueError(f"Endpoint read names disagree at edge {row['edge_id']}")
            row["source_node_id"] = source_node
            row["target_node_id"] = target_node
            yield row


def get_stats():
    headers, reader = read_tsv(STATS_FILE)
    if headers != ["statistic", "value"]:
        raise ValueError(f"Unexpected graph_statistics.tsv schema: {headers}")
    return {row["statistic"]: row["value"] for row in reader}


def summarize_graph(names):
    dsu = DisjointSet(len(names))
    degree = [0] * len(names)
    count = 0
    identity_sum = 0.0
    qaos_sum = 0.0
    overlap_sum = 0
    qaos_min = math.inf
    qaos_max = -math.inf
    identity_min = math.inf
    identity_max = -math.inf
    overlap_min = math.inf
    overlap_max = -math.inf
    sums = {"identity_sq": 0.0, "qaos_sq": 0.0, "identity_qaos": 0.0}
    headers = edge_headers()
    orientation_counts = Counter()
    for row in edge_reader(names):
        source = row["source_node_id"]
        target = row["target_node_id"]
        dsu.union(source, target)
        degree[source] += 1
        degree[target] += 1
        identity = float(row["identity"])
        qaos = float(row["QAOS"])
        overlap = int(row["overlap_length"])
        count += 1
        identity_sum += identity
        qaos_sum += qaos
        overlap_sum += overlap
        sums["identity_sq"] += identity * identity
        sums["qaos_sq"] += qaos * qaos
        sums["identity_qaos"] += identity * qaos
        identity_min = min(identity_min, identity)
        identity_max = max(identity_max, identity)
        qaos_min = min(qaos_min, qaos)
        qaos_max = max(qaos_max, qaos)
        overlap_min = min(overlap_min, overlap)
        overlap_max = max(overlap_max, overlap)
        orientation_counts[(row["source_orientation"], row["target_orientation"])] += 1
    roots = [dsu.find(node) for node in range(len(names))]
    members = defaultdict(list)
    for node_id, root in enumerate(roots):
        members[root].append(node_id)
    components = {
        root: {
            "nodes": node_ids,
            "component_id": min(node_ids),
            "edges": 0,
            "raw_branches": sum(degree[node] > 2 for node in node_ids),
            "raw_degree_sum": sum(degree[node] for node in node_ids),
        }
        for root, node_ids in members.items()
    }
    for row in edge_reader(names):
        root = roots[row["source_node_id"]]
        components[root]["edges"] += 1

    degrees_sorted = sorted(degree)
    mean_degree = 2 * count / len(names) if names else 0.0
    median_degree = statistics.median(degrees_sorted) if degrees_sorted else 0
    mean_identity = identity_sum / count if count else 0.0
    mean_qaos = qaos_sum / count if count else 0.0
    mean_overlap = overlap_sum / count if count else 0.0
    covariance = sums["identity_qaos"] / count - mean_identity * mean_qaos if count else 0.0
    identity_variance = sums["identity_sq"] / count - mean_identity**2 if count else 0.0
    qaos_variance = sums["qaos_sq"] / count - mean_qaos**2 if count else 0.0
    correlation = (
        covariance / math.sqrt(identity_variance * qaos_variance)
        if identity_variance > 0 and qaos_variance > 0
        else float("nan")
    )
    return {
        "names": names,
        "roots": roots,
        "degree": degree,
        "components": components,
        "edge_count": count,
        "mean_degree": mean_degree,
        "median_degree": median_degree,
        "max_degree": max(degree, default=0),
        "isolated_nodes": sum(value == 0 for value in degree),
        "branch_nodes": sum(value > 2 for value in degree),
        "mean_identity": mean_identity,
        "mean_qaos": mean_qaos,
        "mean_overlap": mean_overlap,
        "identity_min": identity_min,
        "identity_max": identity_max,
        "qaos_min": qaos_min,
        "qaos_max": qaos_max,
        "overlap_min": overlap_min,
        "overlap_max": overlap_max,
        "identity_qaos_correlation": correlation,
        "orientation_counts": orientation_counts,
        "headers": headers,
    }


def select_small_components(graph):
    components = graph["components"]
    degree = graph["degree"]
    small = [
        (root, value) for root, value in components.items()
        if 10 <= len(value["nodes"]) <= 30
    ]
    branched = [
        (root, value) for root, value in components.items()
        if 10 <= len(value["nodes"]) <= 50
    ]
    if not small:
        raise RuntimeError("No 10-30 node connected component is available")
    if not branched:
        raise RuntimeError("No 10-50 node connected component is available")

    candidate_roots = {root for root, _ in branched}
    neighbors = {root: defaultdict(set) for root in candidate_roots}
    node_root = graph["roots"]
    for row in edge_reader(graph["names"]):
        source = row["source_node_id"]
        target = row["target_node_id"]
        root = node_root[source]
        if root in candidate_roots:
            neighbors[root][source].add(target)
            neighbors[root][target].add(source)

    def distinct_branch_count(root):
        comp = components[root]
        return sum(len(neighbors[root].get(node, ())) > 2 for node in comp["nodes"])

    def distinct_max_degree(root):
        return max(
            (len(neighbors[root].get(node, ())) for node in components[root]["nodes"]),
            default=0,
        )

    linear_root = min(
        (root for root, _ in small),
        key=lambda root: (
            distinct_max_degree(root),
            distinct_branch_count(root),
            components[root]["edges"],
            abs(len(components[root]["nodes"]) - 20),
            components[root]["component_id"],
        ),
    )
    branched_candidates = [
        root for root, comp in branched
        if distinct_branch_count(root) > 0 and root != linear_root
    ]
    if not branched_candidates:
        raise RuntimeError("No distinct-neighbor branch in a 10-50 node component")
    readable = [root for root in branched_candidates if components[root]["edges"] <= 180]
    if readable:
        branched_candidates = readable
    branched_root = min(
        branched_candidates,
        key=lambda root: (
            -distinct_branch_count(root),
            -distinct_max_degree(root),
            components[root]["edges"],
            abs(len(components[root]["nodes"]) - 30),
            components[root]["component_id"],
        ),
    )
    for root, comp in ((linear_root, components[linear_root]), (branched_root, components[branched_root])):
        comp["distinct_branch_nodes"] = distinct_branch_count(root)
        comp["distinct_max_degree"] = distinct_max_degree(root)
    return linear_root, branched_root, neighbors


def select_high_degree_neighborhood(graph):
    degree = graph["degree"]
    center = max(range(len(degree)), key=degree.__getitem__)
    neighbors = Counter()
    for row in edge_reader(graph["names"]):
        if row["source_node_id"] == center:
            neighbors[row["target_node_id"]] += 1
        elif row["target_node_id"] == center:
            neighbors[row["source_node_id"]] += 1
    chosen = [center] + sorted(
        neighbors, key=lambda node: (-neighbors[node], -degree[node], node)
    )[:99]
    if len(chosen) < 20:
        frontier = set(chosen)
        expanded = set(chosen)
        for row in edge_reader(graph["names"]):
            source = row["source_node_id"]
            target = row["target_node_id"]
            if source in frontier and target not in expanded:
                expanded.add(target)
            elif target in frontier and source not in expanded:
                expanded.add(source)
        chosen = [center] + sorted(
            expanded - {center}, key=lambda node: (-degree[node], node)
        )[:99]
    return center, set(chosen), degree[center], len(neighbors)


def sample_qaos_edges(sample_size=100):
    rng = random.Random(116)
    reservoir = []
    seen = 0
    min_qaos_edge = None
    max_qaos_edge = None
    for row in edge_reader():
        seen += 1
        item = (float(row["identity"]), float(row["QAOS"]), int(row["overlap_length"]))
        if min_qaos_edge is None or item[1] < min_qaos_edge[1]:
            min_qaos_edge = item
        if max_qaos_edge is None or item[1] > max_qaos_edge[1]:
            max_qaos_edge = item
        if len(reservoir) < 10000:
            reservoir.append(item)
        else:
            slot = rng.randrange(seen)
            if slot < len(reservoir):
                reservoir[slot] = item
    reservoir.sort(key=lambda item: (item[1], item[0], item[2]))
    if len(reservoir) <= sample_size:
        return reservoir
    indices = [round(i * (len(reservoir) - 1) / (sample_size - 1)) for i in range(sample_size)]
    sample = [reservoir[index] for index in indices]
    sample[0] = min_qaos_edge
    sample[-1] = max_qaos_edge
    sample.sort(key=lambda item: (item[1], item[0], item[2]))
    return sample


def read_path_metadata():
    _, reader = read_tsv(PATH_STATS_FILE)
    path_stats = {int(row["contig_id"]): row for row in reader}
    groups = defaultdict(list)
    with PATHS_FILE.open(newline="") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields:
                continue
            if len(fields) != 6:
                raise ValueError(f"Unexpected graph_paths.bed row: {fields}")
            contig = int(fields[0].removeprefix("contig_"))
            node_id = int(fields[3].removeprefix("read_"))
            strand = fields[5]
            if strand not in ("+", "-"):
                raise ValueError(f"Unexpected path strand: {strand}")
            groups[contig].append({
                "node_id": node_id,
                "start": int(fields[1]),
                "end": int(fields[2]),
                "qaos": float(fields[4]),
                "orientation": "0" if strand == "+" else "1",
            })
    eligible = [
        (contig, rows) for contig, rows in groups.items()
        if 4 <= len(rows) <= 12 and int(path_stats[contig]["edges"]) >= 3
    ]
    if not eligible:
        eligible = [
            (contig, rows) for contig, rows in groups.items()
            if len(rows) >= 2 and int(path_stats[contig]["edges"]) >= 1
        ]
    contig, rows = min(eligible, key=lambda item: (len(item[1]), item[0]))
    rows.sort(key=lambda row: (row["start"], row["end"]))
    return contig, rows, path_stats[contig]


def edge_records_for(nodes):
    selected = []
    for row in edge_reader():
        if row["source_node_id"] in nodes and row["target_node_id"] in nodes:
            selected.append(row)
    return selected


def quote(value):
    return str(value).replace('"', '\\"')


def render_dot(dot_path, stem, layout="neato"):
    for format_name in ("svg", "png"):
        command = [layout, f"-T{format_name}", str(dot_path), "-o", str(ROOT / f"{stem}.{format_name}")]
        if format_name == "png":
            command.insert(1, f"-Gdpi={PNG_DPI}")
        subprocess.run(command, check=True)


def render_component(stem, nodes, edges, names, degree, component_id, description, layout="neato"):
    selected_edges = [
        row for row in edges
        if row["source_node_id"] in nodes and row["target_node_id"] in nodes
    ]
    dot_path = ROOT / f"{stem}.dot"
    with dot_path.open("w") as handle:
        handle.write(f'digraph "{stem}" {{\n  graph [overlap=false, splines=true, label="{quote(description)}", labelloc=t];\n')
        handle.write('  node [shape=box, fontsize=9, style="rounded,filled", fillcolor="#edf2f4"];\n')
        handle.write('  edge [fontsize=7, arrowsize=0.55];\n')
        for node in sorted(nodes):
            label = f"{names[node]}\\nnode={node} degree={degree[node]}"
            handle.write(f'  n{node} [label="{quote(label)}"];\n')
        for row in selected_edges:
            source = row["source_node_id"]
            target = row["target_node_id"]
            source_orientation = row["source_orientation"]
            target_orientation = row["target_orientation"]
            label = (
                f"ov={row['overlap_length']} id={row['identity']} "
                f"QAOS={row['QAOS']} ori={source_orientation}>{target_orientation}"
            )
            color = EDGE_COLORS[(source_orientation, target_orientation)]
            handle.write(
                f'  n{source} -> n{target} [label="{quote(label)}", color="{color}"];\n'
            )
        handle.write("}\n")
    render_dot(dot_path, stem, layout)
    return len(nodes), len(selected_edges)


def render_path(stem, contig, rows, path_stats, names, edges, degree):
    chosen_edges = []
    for previous, current in zip(rows, rows[1:]):
        candidates = [
            row for row in edges
            if row["source_node_id"] == previous["node_id"]
            and row["target_node_id"] == current["node_id"]
            if row["source_orientation"] == previous["orientation"]
            and row["target_orientation"] == current["orientation"]
        ]
        if not candidates:
            candidates = [
                row for row in edges
                if row["source_node_id"] == current["node_id"]
                and row["target_node_id"] == previous["node_id"]
                if row["source_orientation"] == current["orientation"]
                and row["target_orientation"] == previous["orientation"]
            ]
        if not candidates:
            raise ValueError(
                f"No graph_edges.tsv link for emitted path transition "
                f"{previous['node_id']} -> {current['node_id']}"
            )
        candidates.sort(
            key=lambda row: (
                abs(float(row["QAOS"]) - previous["qaos"]),
                -int(row["overlap_length"]),
            )
        )
        chosen_edges.append(candidates[0])
    dot_path = ROOT / f"{stem}.dot"
    with dot_path.open("w") as handle:
        handle.write(
            f'digraph "{stem}" {{\n  graph [rankdir=LR, label="contig_{contig} '
            f'(component {path_stats["component_id"]})", labelloc=t];\n'
        )
        handle.write('  node [shape=box, fontsize=9, style="rounded,filled", fillcolor="#edf2f4"];\n')
        for index, row in enumerate(rows):
            node_id = row["node_id"]
            orientation = "forward" if row["orientation"] == "0" else "reverse"
            label = f"{names[node_id]}\\n{orientation} degree={degree[node_id]}"
            handle.write(f'  p{index} [label="{quote(label)}"];\n')
        for index, edge in enumerate(chosen_edges):
            label = (
                f"ov={edge['overlap_length']} id={edge['identity']} "
                f"QAOS={edge['QAOS']} ori={edge['source_orientation']}>{edge['target_orientation']}"
            )
            color = EDGE_COLORS[(edge["source_orientation"], edge["target_orientation"])]
            handle.write(f'  p{index} -> p{index + 1} [label="{quote(label)}", color="{color}"];\n')
        handle.write("}\n")
    render_dot(dot_path, stem, "dot")
    return len(rows), len(chosen_edges)


def render_scatter(sample):
    data_path = ROOT / "debug_identity_vs_qaos_sample.tsv"
    with data_path.open("w") as handle:
        handle.write("identity\tQAOS\toverlap_length\n")
        for identity, qaos, overlap in sample:
            handle.write(f"{identity}\t{qaos}\t{overlap}\n")
    plot_path = ROOT / "debug_identity_vs_qaos.gnuplot"
    plot_path.write_text(
        "set terminal pngcairo size 1100,800 enhanced font 'Sans,11'\n"
        f"set output '{ROOT / 'debug_identity_vs_qaos.png'}'\n"
        "set title 'Accepted overlap edges: identity vs QAOS (quantile sample)'\n"
        "set xlabel 'Sequence identity'\n"
        "set ylabel 'QAOS'\n"
        "set cblabel 'Overlap length (bp)'\n"
        "set xrange [0.89:1.005]\n"
        "set yrange [0.89:1.005]\n"
        "set grid\n"
        f"plot '{data_path}' using 1:2:3 with points pt 7 ps 1.2 palette notitle\n"
    )
    subprocess.run(["gnuplot", str(plot_path)], check=True)


def write_summary(graph, stats, selected, path_info, sample):
    degree = graph["degree"]
    output = ROOT / "debug_graph_summary.txt"
    lines = [
        f"total nodes: {len(degree)}",
        f"total edges: {graph['edge_count']}",
        f"connected components: {len(graph['components'])}",
        f"largest component size: {max(len(comp['nodes']) for comp in graph['components'].values())}",
        f"mean degree: {graph['mean_degree']:.6f}",
        f"median degree: {graph['median_degree']}",
        f"maximum degree: {graph['max_degree']}",
        f"number of isolated nodes: {graph['isolated_nodes']}",
        f"number of branch nodes (raw edge degree > 2): {graph['branch_nodes']}",
        f"number of cycles if available: {stats.get('cycles', 'NOT AVAILABLE')}",
        f"mean overlap: {stats.get('mean_overlap', graph['mean_overlap'])}",
        f"median overlap: {stats.get('median_overlap', 'NOT AVAILABLE')}",
        f"mean identity: {stats.get('mean_identity', graph['mean_identity'])}",
        f"median identity: {stats.get('median_identity', 'NOT AVAILABLE')}",
        f"mean QAOS: {stats.get('mean_QAOS', graph['mean_qaos'])}",
        f"median QAOS: {stats.get('median_QAOS', 'NOT AVAILABLE')}",
        f"linear component ID (minimum node_id): {selected['linear']['component_id']}",
        f"linear component nodes / edge records / distinct branch nodes / max distinct degree: {len(selected['linear']['nodes'])} / {selected['linear']['edges']} / {selected['linear']['distinct_branch_nodes']} / {selected['linear']['distinct_max_degree']}",
        f"branched component ID (minimum node_id): {selected['branched']['component_id']}",
        f"branched component nodes / edge records / distinct branch nodes / max distinct degree: {len(selected['branched']['nodes'])} / {selected['branched']['edges']} / {selected['branched']['distinct_branch_nodes']} / {selected['branched']['distinct_max_degree']}",
        f"high-degree center node_id: {selected['high']['center']}",
        f"high-degree center raw edge-record degree: {selected['high']['degree']}",
        f"high-degree center distinct neighboring nodes: {selected['high']['distinct_neighbors']}",
        f"high-degree neighborhood component ID: {selected['high']['component_id']}",
        f"example path contig_id: {path_info['contig_id']}",
        f"example path component ID: {path_info['component_id']}",
        f"accepted-edge QAOS sample size: {len(sample)}",
        f"QAOS scatter-sample range: {min(item[1] for item in sample):.6f} .. {max(item[1] for item in sample):.6f}",
        f"accepted-edge identity range: {graph['identity_min']:.6f} .. {graph['identity_max']:.6f}",
        f"accepted-edge QAOS range: {graph['qaos_min']:.6f} .. {graph['qaos_max']:.6f}",
        f"accepted-edge overlap range: {graph['overlap_min']} .. {graph['overlap_max']}",
        f"Pearson correlation, all accepted-edge identity vs QAOS: {graph['identity_qaos_correlation']:.6f}",
        "",
        f"duplicate read IDs in graph_nodes.tsv: {selected['duplicate_read_ids']} names occur twice; edge endpoints are disambiguated by row-wise graph_edges.bed read_N IDs (coordinates unused)",
        "",
        "Orientation counts in graph_edges.tsv (source_orientation>target_orientation):",
    ]
    for orientation, count in sorted(graph["orientation_counts"].items()):
        lines.append(f"  {orientation[0]}>{orientation[1]}: {count}")
    lines.extend([
        "",
        'BED orientation requires correction. graph_edges.bed carries orientation '
        "fields, but its source/target intervals are emitted as suffix/prefix "
        "coordinates without reflecting reverse-oriented endpoints. Graphviz "
        "uses edge labels/orientations from graph_edges.tsv and the row-aligned "
        "BED endpoint IDs only; BED coordinates are ignored.",
    ])
    output.write_text("\n".join(lines) + "\n")


def write_report(graph, selected, path_info, counts, sample):
    linear = selected["linear"]
    branched = selected["branched"]
    high = selected["high"]
    high_edges = counts["high"][1]
    orientation_counts = ", ".join(
        f"{left}>{right}: {count}"
        for (left, right), count in sorted(graph["orientation_counts"].items())
    )
    qaos_diff = graph["qaos_max"] - graph["qaos_min"]
    report = f"""# Raw Graph Debug Visualization Report

These figures are for manual inspection of graph topology only. They do not
establish that graph construction is correct and use no V/J/CDR3 annotations.

## Selected views

| View | Component / path | Selection reason | Nodes shown | Edges shown |
|---|---|---|---:|---:|
| Simple linear | Component {linear['component_id']} (minimum node_id) | 10–30-node component selected for lowest maximum distinct-neighbor degree, then fewest distinct-neighbor branch nodes and raw overlap edges | {counts['linear'][0]} | {counts['linear'][1]} |
| Branched | Component {branched['component_id']} (minimum node_id) | 10–50-node component selected for distinct-neighbor branches, preferring at most 180 edge records for legibility | {counts['branched'][0]} | {counts['branched'][1]} |
| High connectivity | Center node {high['center']}; local component {high['component_id']} | The node with maximum incident edge-record degree ({high['degree']}); shown with up to 99 neighbors prioritized by parallel edge count | {counts['high'][0]} | {high_edges} |
| Example assembled path | contig_{path_info['contig_id']}, component {path_info['component_id']} | Shortest emitted path in the selected 4–12-read range; ordered from graph_paths.bed | {counts['path'][0]} | {counts['path'][1]} |

Component IDs shown for raw graph components are the minimum numeric node_id
in that connected component. The path component ID is the assembler's
component_id from path_statistics.tsv.

SVG and PNG files are generated for the linear component, branched component,
high-degree neighborhood, and ordered example path. Their edge labels include
overlap length, identity, QAOS, and source-to-target orientation bits. Edge
colors distinguish orientation combinations: {orientation_counts}.

## QAOS inspection

`debug_identity_vs_qaos.png` plots {len(sample)} accepted edges selected at
evenly spaced QAOS quantiles from a deterministic reservoir sample of up to
10,000 edges. The full accepted-edge ranges are:

- identity: {graph['identity_min']:.6f} to {graph['identity_max']:.6f}
- QAOS: {graph['qaos_min']:.6f} to {graph['qaos_max']:.6f}
- overlap length: {graph['overlap_min']} to {graph['overlap_max']} bp
- Pearson correlation for all accepted edges, identity vs QAOS:
  {graph['identity_qaos_correlation']:.6f}

QAOS varies across the accepted edges (range width {qaos_diff:.6f}); it is not
numerically identical to identity. The plot includes the exact minimum- and
maximum-QAOS edges. Its actual sampled QAOS range is
{min(item[1] for item in sample):.6f} to {max(item[1] for item in sample):.6f}.
The scatter plot and its sampled data are
available as `debug_identity_vs_qaos.png` and
`debug_identity_vs_qaos_sample.tsv`.

## Orientation and BED coordinates

`graph_edges.tsv` explicitly represents `source_orientation` and
`target_orientation` as 0/1 values. `graph_paths.bed` represents path
orientation as `+`/`-`. The edge BED has orientation columns but its
coordinates are always written as source suffix and target prefix, with no
coordinate reflection for reverse-oriented endpoints. **BED orientation requires correction.**
Also, `graph_edges.tsv` endpoint names are not unique:
all 100,000 paired read IDs occur twice in `graph_nodes.tsv`. The generator
uses edge measurements/orientations from `graph_edges.tsv` and pairs each
edge row with the same-order BED record solely to recover its unique `read_N`
node IDs; BED coordinates are never used. Endpoint names are cross-checked.

## Patterns to inspect

- The high-degree center has {high['degree']} incident edge records and
  {high['distinct_neighbors']} distinct neighboring reads. The displayed
  100-node local induced view contains {high_edges} edge records; neighbors
  were prioritized by number of parallel links to the center. Inspect for
  dense parallel links.
- The branched example has {counts['branched'][1]} edge records over
  {counts['branched'][0]} nodes. Edge-record degree and distinct-neighbor
  branching are both visible from repeated labels and node degree annotations.
- The simple-component selection minimizes maximum distinct-neighbor degree
  among eligible components, but it still retains parallel overlap edges.
- These are visual patterns to review manually, not biological conclusions
  or proof that overlaps or paths are correct.
"""
    (ROOT / "DEBUG_VISUALIZATION_REPORT.md").write_text(report)


def main():
    names = read_nodes()
    graph = summarize_graph(names)
    stats = get_stats()
    if int(stats["graph_nodes"]) != len(names):
        raise ValueError("graph_statistics.tsv and graph_nodes.tsv node counts disagree")
    if int(stats["graph_edges"]) != graph["edge_count"]:
        raise ValueError("graph_statistics.tsv and graph_edges.tsv edge counts disagree")

    linear_root, branched_root, _ = select_small_components(graph)
    linear_component = graph["components"][linear_root]
    branched_component = graph["components"][branched_root]

    center, high_nodes, center_degree, distinct_neighbors = select_high_degree_neighborhood(graph)
    high_component_root = graph["roots"][center]
    high_component = graph["components"][high_component_root]

    contig, path_rows, path_stats = read_path_metadata()
    path_stats["contig_id"] = contig
    with FASTA_FILE.open() as handle:
        fasta_records = sum(line.startswith(">") for line in handle)
    path_stats_ids = {
        int(row["contig_id"]) for row in read_tsv(PATH_STATS_FILE)[1]
    }
    if fasta_records != len(path_stats_ids):
        raise ValueError("FASTA and path_statistics.tsv record counts disagree")

    selected_nodes = {
        "linear": set(linear_component["nodes"]),
        "branched": set(branched_component["nodes"]),
        "high": high_nodes,
        "path": {row["node_id"] for row in path_rows},
    }
    component_edges = {key: [] for key in ("linear", "branched", "high", "path")}
    for row in edge_reader(names):
        for key, node_ids in selected_nodes.items():
            if row["source_node_id"] in node_ids and row["target_node_id"] in node_ids:
                component_edges[key].append(row)

    counts = {}
    counts["linear"] = render_component(
        "debug_linear_component",
        selected_nodes["linear"],
        component_edges["linear"],
        names,
        graph["degree"],
        linear_component["component_id"],
        f"Simple component {linear_component['component_id']}",
    )
    counts["branched"] = render_component(
        "debug_branched_component",
        selected_nodes["branched"],
        component_edges["branched"],
        names,
        graph["degree"],
        branched_component["component_id"],
        f"Branched component {branched_component['component_id']}",
    )
    counts["high"] = render_component(
        "debug_high_degree",
        selected_nodes["high"],
        component_edges["high"],
        names,
        graph["degree"],
        high_component["component_id"],
        f"High-degree neighborhood centered at node {center}",
        "sfdp",
    )
    counts["path"] = render_path(
        "debug_example_path",
        contig,
        path_rows,
        path_stats,
        names,
        component_edges["path"],
        graph["degree"],
    )
    sample = sample_qaos_edges()
    render_scatter(sample)

    selected = {
        "linear": linear_component,
        "branched": branched_component,
        "high": {
            "center": center,
            "component_id": high_component["component_id"],
            "nodes": selected_nodes["high"],
            "degree": center_degree,
            "distinct_neighbors": distinct_neighbors,
        },
    }
    name_counts = Counter(names.values())
    selected["duplicate_read_ids"] = sum(count > 1 for count in name_counts.values())
    write_summary(graph, stats, selected, path_stats, sample)
    write_report(graph, selected, path_stats, counts, sample)
    print(f"linear component={linear_component['component_id']} nodes/edges={counts['linear']}")
    print(f"branched component={branched_component['component_id']} nodes/edges={counts['branched']}")
    print(f"high-degree center={center} nodes/edges={counts['high']}")
    print(f"example path=contig_{contig} nodes/edges={counts['path']}")
    print(f"QAOS plot sample={len(sample)} range={graph['qaos_min']}..{graph['qaos_max']}")


if __name__ == "__main__":
    main()
