# Raw Graph Debug Visualization Report

These figures are for manual inspection of graph topology only. They do not
establish that graph construction is correct and use no V/J/CDR3 annotations.

## Selected views

| View | Component / path | Selection reason | Nodes shown | Edges shown |
|---|---|---|---:|---:|
| Simple linear | Component 155529 (minimum node_id) | 10–30-node component selected for lowest maximum distinct-neighbor degree, then fewest distinct-neighbor branch nodes and raw overlap edges | 10 | 12 |
| Branched | Component 198262 (minimum node_id) | 10–50-node component selected for distinct-neighbor branches, preferring at most 180 edge records for legibility | 48 | 178 |
| High connectivity | Center node 8926; local component 7706 | The node with maximum incident edge-record degree (3417); shown with up to 99 neighbors prioritized by parallel edge count | 100 | 146 |
| Example assembled path | contig_1460, component 360 | Shortest emitted path in the selected 4–12-read range; ordered from graph_paths.bed | 4 | 3 |

Component IDs shown for raw graph components are the minimum numeric node_id
in that connected component. The path component ID is the assembler's
component_id from path_statistics.tsv.

SVG and PNG files are generated for the linear component, branched component,
high-degree neighborhood, and ordered example path. Their edge labels include
overlap length, identity, QAOS, and source-to-target orientation bits. Edge
colors distinguish orientation combinations: 0>0: 2365816, 0>1: 4292538, 1>0: 35, 1>1: 196.

## QAOS inspection

`debug_identity_vs_qaos.png` plots 100 accepted edges selected at
evenly spaced QAOS quantiles from a deterministic reservoir sample of up to
10,000 edges. The full accepted-edge ranges are:

- identity: 0.935484 to 1.000000
- QAOS: 0.900000 to 0.999841
- overlap length: 31 to 150 bp
- Pearson correlation for all accepted edges, identity vs QAOS:
  0.831756

QAOS varies across the accepted edges (range width 0.099841); it is not
numerically identical to identity. The plot includes the exact minimum- and
maximum-QAOS edges. Its actual sampled QAOS range is
0.900000 to 0.999841.
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

- The high-degree center has 3417 incident edge records and
  3370 distinct neighboring reads. The displayed
  100-node local induced view contains 146 edge records; neighbors
  were prioritized by number of parallel links to the center. Inspect for
  dense parallel links.
- The branched example has 178 edge records over
  48 nodes. Edge-record degree and distinct-neighbor
  branching are both visible from repeated labels and node degree annotations.
- The simple-component selection minimizes maximum distinct-neighbor degree
  among eligible components, but it still retains parallel overlap edges.
- These are visual patterns to review manually, not biological conclusions
  or proof that overlaps or paths are correct.
