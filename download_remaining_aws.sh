#!/bin/bash

BASE="/data1/wetlab/TRUST4_benchmark/data/raw_sra"

download_sample() {
    SRR="$1"

    echo "========================================"
    echo "Starting $SRR: $(date)"
    echo "========================================"

    mkdir -p "$BASE/$SRR"

    aws s3 cp \
        "s3://sra-pub-run-odp/sra/$SRR/$SRR" \
        "$BASE/$SRR/$SRR.sra" \
        --no-sign-request

    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "===== FINISHED $SRR: $(date) ====="
    else
        echo "===== FAILED $SRR (exit code $STATUS): $(date) ====="
    fi
}

export -f download_sample
export BASE

printf "%s\n" \
    SRR7882940 \
    SRR7882939 \
    SRR7882938 \
    SRR7882937 \
    | xargs -n 1 -P 2 bash -c 'download_sample "$1"' _
