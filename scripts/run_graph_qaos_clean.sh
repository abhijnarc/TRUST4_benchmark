#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
r1="${1:?raw R1 FASTQ required}"
r2="${2:?raw R2 FASTQ required}"
sample="${3:-FZ-116}"
out="${4:-$root/results/Graph-TRUST4-QAOS-clean/$sample}"
threads="${5:-64}"
extract_prefix="$out/candidates/candidate"
mkdir -p "$out/candidates" "$root/logs/Graph-TRUST4-QAOS-clean/$sample"

ref="$root/reference/TRUST4/human_IMGT+C.fa"
extract_log="$root/logs/Graph-TRUST4-QAOS-clean/$sample/candidate_extraction.log"
"$root/algorithms/TRUST4/fastq-extractor" -t "$threads" -f "$ref" -o "$extract_prefix" \
  -1 "$r1" -2 "$r2" >"$extract_log" 2>&1

candidate_r1="${extract_prefix}_1.fq"
candidate_r2="${extract_prefix}_2.fq"
raw_pairs=$(( $(wc -l < "$r1") / 4 ))
candidate_pairs=$(( $(wc -l < "$candidate_r1") / 4 ))
python3 - "$out/candidate_statistics.tsv" "$raw_pairs" "$candidate_pairs" <<'PY'
import sys
path, raw_pairs, candidate_pairs = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with open(path, "w") as handle:
    handle.write("metric\tvalue\n")
    handle.write(f"raw_read_pairs\t{raw_pairs}\n")
    handle.write(f"candidate_read_pairs\t{candidate_pairs}\n")
    handle.write(f"candidate_fraction\t{candidate_pairs / raw_pairs if raw_pairs else 0:.12g}\n")
    handle.write(f"raw_reads\t{raw_pairs * 2}\n")
    handle.write(f"candidate_reads\t{candidate_pairs * 2}\n")
PY

make -C "$root/algorithms/Graph-TRUST4-QAOS-clean"
/usr/bin/time -v "$root/algorithms/Graph-TRUST4-QAOS-clean/graph-trust4-qaos-clean" \
  -1 "$candidate_r1" -2 "$candidate_r2" -o "$out" -s "$out/candidate_statistics.tsv" \
  -k 9 -m 31 -t "$threads" >"$root/logs/Graph-TRUST4-QAOS-clean/$sample/assembly.stdout" \
  2>"$root/logs/Graph-TRUST4-QAOS-clean/$sample/assembly.time"

python3 "$root/scripts/annotate_fasta_batched.py" \
  --annotator "$root/algorithms/TRUST4/annotator" --reference "$ref" \
  --fasta "$out/assembled_contigs.fa" --output "$out/annotated_contigs.tsv" \
  --log-dir "$root/logs/Graph-TRUST4-QAOS-clean/$sample/annotation" --threads "$threads"
python3 "$root/scripts/visualize_graph_qaos_clean.py" --output "$out"
