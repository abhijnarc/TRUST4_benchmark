#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
r1="$1"
r2="$2"
out="$3"
mode="${4:-qaos}"
threads="${5:-64}"
mkdir -p "$out"

ref="$root/reference/TRUST4/human_IMGT+C.fa"
prefix="$out/trust4_candidates"
"$root/algorithms/TRUST4/fastq-extractor" -t "$threads" -f "$ref" -o "$prefix" -1 "$r1" -2 "$r2" \
  >"$out/candidate_extraction.stdout" 2>"$out/candidate_extraction.stderr"

candidate_r1="${prefix}_1.fq"
candidate_r2="${prefix}_2.fq"
raw_reads=$(( $(wc -l < "$r1") / 4 + $(wc -l < "$r2") / 4 ))
candidate_r1_reads=$(( $(wc -l < "$candidate_r1") / 4 ))
candidate_r2_reads=$(( $(wc -l < "$candidate_r2") / 4 ))
candidate_reads=$(( candidate_r1_reads + candidate_r2_reads ))
candidate_pairs=$(( candidate_r1_reads < candidate_r2_reads ? candidate_r1_reads : candidate_r2_reads ))
python - "$out/candidate_statistics.tsv" "$raw_reads" "$candidate_reads" "$candidate_pairs" <<'PY'
import sys
out, raw, reads, pairs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
with open(out, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"raw_reads\t{raw}\n")
    f.write(f"candidate_reads\t{reads}\n")
    f.write(f"candidate_pairs\t{pairs}\n")
    f.write(f"candidate_fraction\t{reads / raw if raw else 0:.12g}\n")
PY

args=(-1 "$candidate_r1" -2 "$candidate_r2" -o "$out" -k 9 -m 31 -t "$threads" -M 32)
if [[ "$mode" == "identity-only" ]]; then
  args+=(--identity-only)
fi
if [[ "$mode" == "v3" ]]; then
  args+=(--v3-resolution)
fi
exec "$root/algorithms/Graph-TRUST4-QAOS/graph-trust4-qaos" "${args[@]}"
