#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/r1.fq" <<'EOF'
@r1/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2/1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
cat > "$tmp/r2.fq" <<'EOF'
@r1/2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@r2/2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
printf 'metric\tvalue\nraw_read_pairs\t2\n' > "$tmp/candidates.tsv"
make -C "$root/algorithms/Graph-TRUST4-QAOS-clean" >/dev/null
"$root/algorithms/Graph-TRUST4-QAOS-clean/graph-trust4-qaos-clean" \
  -1 "$tmp/r1.fq" -2 "$tmp/r2.fq" -o "$tmp/out" -s "$tmp/candidates.tsv" \
  -k 9 -m 31 -t 2 --max-candidates 8 --max-paths 10 >/dev/null
test -s "$tmp/out/assembled_contigs.fa"
test -s "$tmp/out/graph_edges.tsv"
test -s "$tmp/out/path_statistics.tsv"
test -s "$tmp/out/consensus_support.tsv"
test -s "$tmp/out/graph_nodes.bed"
test -s "$tmp/out/graph_edges.bed"
test -s "$tmp/out/graph_paths.bed"
test -s "$tmp/out/visualization/README.txt"
grep -q 'graph_nodes' "$tmp/out/graph_statistics.tsv"
echo "graph_qaos_clean_smoke=PASS"
