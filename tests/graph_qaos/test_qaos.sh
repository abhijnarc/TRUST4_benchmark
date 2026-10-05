#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
bin="$root/algorithms/Graph-TRUST4-QAOS/graph-trust4-qaos"
tmp="$root/results/Graph-TRUST4-QAOS/test_unit"
mkdir -p "$tmp"

make -C "$root/algorithms/Graph-TRUST4-QAOS" graph-trust4-qaos >/dev/null
"$bin" --self-test >/dev/null

cat > "$tmp/r1.fq" <<'EOF'
@pair1/1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/1
ACGTACGTACGTACGTACGTACGTACGTACGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
cat > "$tmp/r2.fq" <<'EOF'
@pair1/2
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@pair2/2
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

"$bin" -1 "$tmp/r1.fq" -2 "$tmp/r2.fq" -o "$tmp/out" \
  -k 9 -m 31 -t 2 -M 32 >/dev/null
test -s "$tmp/out/assembled_contigs.fa"
test -s "$tmp/out/graph_statistics.tsv"
test -s "$tmp/out/qaos_statistics.tsv"
test -s "$tmp/out/path_statistics.tsv"
test -s "$tmp/out/runtime.tsv"
grep -q $'paired_links\t' "$tmp/out/graph_statistics.tsv"
grep -q $'paired_support' "$tmp/out/path_statistics.tsv"
echo "graph_qaos_tests=PASS"
