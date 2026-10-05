#!/bin/bash
set -e

# TRUST4 test run for FZ-116 (SRR7882936)
# Date: 2026-09-18

PROJECT_ROOT="/data1/wetlab/TRUST4_benchmark"
TRUST4_DIR="${PROJECT_ROOT}/algorithms/TRUST4"
FASTQ_DIR="${PROJECT_ROOT}/data/fastq"
REFERENCE_DIR="${PROJECT_ROOT}/reference/TRUST4"
RESULTS_DIR="${PROJECT_ROOT}/results/TRUST4"
LOG_DIR="${PROJECT_ROOT}/logs"
THREADS=8
SAMPLE="FZ-116"
SRR_ID="SRR7882936"

# Create directories
mkdir -p "${RESULTS_DIR}/${SAMPLE}"
mkdir -p "${LOG_DIR}"

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh || source ~/miniforge3/etc/profile.d/conda.sh
conda activate trust4_benchmark

echo "================================================================"
echo "TRUST4 Test Run: ${SAMPLE}"
echo "================================================================"
echo "Sample: ${SAMPLE} (${SRR_ID})"
echo "Threads: ${THREADS}"
echo "Input FASTQ:"
echo "  Read 1: ${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"
echo "  Read 2: ${FASTQ_DIR}/${SRR_ID}_2.fastq.gz"
echo "Reference: ${REFERENCE_DIR}/human_IMGT+C.fa"
echo "Output directory: ${RESULTS_DIR}/${SAMPLE}"
echo "Log file: ${LOG_DIR}/TRUST4_${SAMPLE}.log"
echo "================================================================"

# Run TRUST4
cd "${RESULTS_DIR}/${SAMPLE}"
${TRUST4_DIR}/run-trust4 \
    -f "${REFERENCE_DIR}/human_IMGT+C.fa" \
    --ref "${REFERENCE_DIR}/human_IMGT+C.fa" \
    -1 "${FASTQ_DIR}/${SRR_ID}_1.fastq.gz" \
    -2 "${FASTQ_DIR}/${SRR_ID}_2.fastq.gz" \
    -t ${THREADS} \
    -o "TRUST_${SAMPLE}" \
    2>&1 | tee "${LOG_DIR}/TRUST4_${SAMPLE}.log"

echo ""
echo "================================================================"
echo "TRUST4 Run Completed"
echo "================================================================"
echo "Output files:"
ls -lh "${RESULTS_DIR}/${SAMPLE}"/TRUST_${SAMPLE}* 2>/dev/null || echo "No output files found"
echo ""
