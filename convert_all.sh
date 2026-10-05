#!/bin/bash

BASE="/data1/wetlab/TRUST4_benchmark"
OUT="$BASE/data/fastq"
TMP="$BASE/data/tmp"

SAMPLES=(
    SRR7882935
    SRR7882936
    SRR7882937
    SRR7882938
    SRR7882939
    SRR7882940
)

for SRR in "${SAMPLES[@]}"; do

    echo "=========================================="
    echo "STARTING $SRR - $(date)"
    echo "=========================================="

    SRA=$(find "$BASE/data/raw_sra/$SRR" -type f -name "$SRR.sra" | head -1)

    if [ -z "$SRA" ]; then
        echo "ERROR: SRA file not found: $SRR"
        continue
    fi

    echo "Using: $SRA"

    # Skip conversion if both compressed FASTQs already exist
    if [ -f "$OUT/${SRR}_1.fastq.gz" ] && [ -f "$OUT/${SRR}_2.fastq.gz" ]; then
        echo "$SRR already converted. Skipping."
        continue
    fi

    fasterq-dump "$SRA" \
        --split-files \
        --threads 12 \
        --temp "$TMP" \
        --outdir "$OUT"

    if [ $? -ne 0 ]; then
        echo "ERROR: fasterq-dump failed for $SRR"
        continue
    fi

    echo "Compressing $SRR..."

    pigz -p 12 "$OUT/${SRR}_1.fastq"
    pigz -p 12 "$OUT/${SRR}_2.fastq"

    echo "FINISHED $SRR - $(date)"
done

echo "=========================================="
echo "ALL SAMPLES PROCESSED - $(date)"
echo "=========================================="