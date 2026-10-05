#!/bin/bash

IREP="/data1/wetlab/TRUST4_benchmark/reference/iRep/ng-bcr-validate/iRep"
OUT="/data1/wetlab/TRUST4_benchmark/reference/iRep/iRep_reference_summary.tsv"

echo -e "Sample\tTotal_records\tUnique_CDR3_nuc\tUnique_CDR3_pep\tCopy1_records\tMax_copy\tUnique_CDR3_V_J_C" > "$OUT"

for f in "$IREP"/FZ-*.csv.gz
do
    echo "Found file: $f"

    sample=$(basename "$f" .csv.gz)

    echo "Processing $sample..."

    tmp="/tmp/${sample}_irep_$$"
    mkdir -p "$tmp"

    zcat "$f" | tail -n +2 > "$tmp/data"

    total=$(wc -l < "$tmp/data")

    unique_nuc=$(cut -d',' -f6 "$tmp/data" | tr '[:lower:]' '[:upper:]' | sort -u | wc -l)

    unique_pep=$(cut -d',' -f1 "$tmp/data" | sort -u | wc -l)

    copy1=$(awk -F',' '$7 == 1 {n++} END {print n+0}' "$tmp/data")

    maxcopy=$(cut -d',' -f7 "$tmp/data" | sort -n | tail -1)

    combo=$(awk -F',' '
    {
        nuc=toupper($6)

        v=$2
        j=$4
        c=$5

        sub(/\*.*/, "", v)
        sub(/\*.*/, "", j)
        sub(/\*.*/, "", c)

        print nuc "|" v "|" j "|" c
    }' "$tmp/data" | sort -u | wc -l)

    echo -e "$sample\t$total\t$unique_nuc\t$unique_pep\t$copy1\t$maxcopy\t$combo" >> "$OUT"

    echo "Completed $sample"
    echo

    rm -rf "$tmp"
done

echo "======================================"
echo "ALL SAMPLES COMPLETE"
echo "======================================"
echo
cat "$OUT"
