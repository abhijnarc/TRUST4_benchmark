#!/usr/bin/env python3
"""Write small diagnostic DOT subgraphs selected from clean graph outputs."""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root = args.output
    node_count = 0
    names = {}
    lengths = {}
    with (root / "graph_nodes.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            node = int(row["node_id"])
            node_count = max(node_count, node + 1)
            names[node] = row["read_id"]
            lengths[node] = int(row["length"])

    parent = list(range(node_count))
    rank = [0] * node_count
    degree = [0] * node_count

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root == right_root:
            return
        if rank[left_root] < rank[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        if rank[left_root] == rank[right_root]:
            rank[left_root] += 1

    edge_path = root / "graph_edges.bed"
    with edge_path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            source, target = int(row[0][5:]), int(row[3][5:])
            union(source, target)
            degree[source] += 1
            degree[target] += 1

    roots = [find(node) for node in range(node_count)]
    sizes = Counter(roots)
    max_degree = Counter()
    for node, component in enumerate(roots):
        max_degree[component] = max(max_degree[component], degree[node])
    component_nodes = defaultdict(list)
    for node, component in enumerate(roots):
        if len(component_nodes[component]) < 100:
            component_nodes[component].append(node)

    bcr_path = None
    with (root / "annotated_contigs.tsv").open() as annotations:
        for row in csv.DictReader(annotations, delimiter="\t"):
            if row.get("v_call") and row.get("j_call") and row.get("junction"):
                bcr_path = int(row["sequence_id"].rsplit("_", 1)[1])
                break
    path_counts = Counter()
    path_roots = {}
    bcr_path_nodes = []
    bcr_component = None
    with (root / "graph_paths.bed").open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            path_id = int(row[0].split("_", 1)[1])
            node = int(row[3].split("_", 1)[1])
            component = roots[node]
            if row[1] == "0":
                path_roots[path_id] = component
                path_counts[component] += 1
            if bcr_path is not None and path_id == bcr_path and len(bcr_path_nodes) < 100:
                bcr_path_nodes.append(node)
    if bcr_path is not None:
        bcr_component = path_roots.get(bcr_path)

    representatives = {}
    for component, count in sizes.items():
        if "linear" not in representatives and 1 < count <= 100 and max_degree[component] <= 2:
            representatives["linear"] = component
        if "branched" not in representatives and max_degree[component] > 2:
            representatives["branched"] = component
    multi = next((component for component, count in path_counts.items() if count > 1), None)
    if multi is not None:
        representatives["multiple_paths"] = multi
    if bcr_component is not None:
        representatives["annotated_bcr_like"] = bcr_component
    if not representatives and sizes:
        representatives["largest_component"] = max(sizes, key=sizes.get)

    out_dir = root / "visualization"
    out_dir.mkdir(exist_ok=True)
    for category, component in representatives.items():
        selected_nodes = set(bcr_path_nodes) if category == "annotated_bcr_like" else set(component_nodes[component])
        if category == "branched":
            selected_nodes.update(sorted(
                (node for node, comp in enumerate(roots) if comp == component),
                key=lambda node: degree[node], reverse=True
            )[:100])
        path = out_dir / f"{category}_component_{component}.dot"
        emitted = 0
        with path.open("w") as dot:
            dot.write(f"digraph {category} {{\nrankdir=LR;\n")
            for node in sorted(selected_nodes):
                dot.write(f'r{node} [label="{node}\\n{lengths[node]}bp\\n{names[node]}"];\n')
            with edge_path.open() as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if not row or row[0].startswith("#"):
                        continue
                    source, target = int(row[0][5:]), int(row[3][5:])
                    if source in selected_nodes and target in selected_nodes:
                        dot.write(f'r{source} -> r{target} [label="ov={row[5]} id={row[6]} QAOS={row[9]}"];\n')
                        emitted += 1
                        if emitted >= 300:
                            break
            dot.write("}\n")
    with (out_dir / "README.txt").open("w") as handle:
        handle.write("Selected small overlap-graph subgraphs; coordinates are read-relative, not genomic.\n")
        handle.write("Categories identify linear, branched, multiple-path, or TRUST4-annotated BCR-like examples where available.\n")
if __name__ == "__main__":
    main()
