#!/bin/bash

BASE="/data1/wetlab/TRUST4_benchmark"
RAW="$BASE/data/raw_sra"
OUT="$BASE/data/fastq"
TMP="$BASE/data/tmp"
LOG="$BASE/logs/convert_remaining.log"

mkdir -p "$OUT" "$TMP" "$BASE/logs"

SAMPLES=(
    SRR7882936
    SRR7882937
    SRR7882938
    SRR7882939
    SRR7882940
)

echo "Queue started: $(date)" >> "$LOG"

# Wait for the currently running SRR7882935 fasterq-dump
while pgrep -x fasterq-dump > /dev/null; do
    echo "Waiting for current fasterq-dump: $(date)" >> "$LOG"
    sleep 60
done

echo "Current fasterq-dump finished: $(date)" >> "$LOG"

# Compress SRR7882935 if its FASTQs were successfully produced
if [[ -f "$OUT/SRR7882935_1.fastq" && -f "$OUT/SRR7882935_2.fastq" ]]; then
    echo "Compressing SRR7882935: $(date)" >> "$LOG"
    pigz -p 12 "$OUT/SRR7882935_1.fastq" "$OUT/SRR7882935_2.fastq"
else
    echo "WARNING: SRR7882935 FASTQs not found; skipping compression." >> "$LOG"
fi

# Process remaining five sequentially
for SRR in "${SAMPLES[@]}"; do

    echo "========================================" >> "$LOG"
    echo "Starting $SRR: $(date)" >> "$LOG"

    SRA=$(find "$RAW/$SRR" -type f -name "${SRR}.sra" | head -1)

    if [[ -z "$SRA" ]]; then
        echo "ERROR: $SRR .sra not found. Skipping." >> "$LOG"
        continue
    fi

    echo "Using: $SRA" >> "$LOG"

    fasterq-dump "$SRA" \
        --split-files \
        --threads 12 \
        --temp "$TMP" \
        --outdir "$OUT" >> "$LOG" 2>&1

    STATUS=$?

    if [[ $STATUS -ne 0 ]]; then
        echo "ERROR: $SRR fasterq-dump failed with exit code $STATUS" >> "$LOG"
        echo "Moving to next sample." >> "$LOG"
        continue
    fi

    if [[ -f "$OUT/${SRR}_1.fastq" && -f "$OUT/${SRR}_2.fastq" ]]; then
        echo "Dump successful. Compressing $SRR: $(date)" >> "$LOG"

        pigz -p 12 \
            "$OUT/${SRR}_1.fastq" \
            "$OUT/${SRR}_2.fastq" >> "$LOG" 2>&1

        echo "Finished $SRR: $(date)" >> "$LOG"
    else
        echo "ERROR: expected FASTQ files for $SRR were not produced." >> "$LOG"
    fi

done

echo "========================================" >> "$LOG"
echo "Queue finished: $(date)" >> "$LOG"