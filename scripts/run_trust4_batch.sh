#!/bin/bash

###############################################################################
# TRUST4 Batch Runner - Remaining 5 Samples
#
# Purpose: Run TRUST4 sequentially on FZ-20, FZ-83, FZ-94, FZ-97, FZ-122
# Uses exact same parameters as FZ-116 baseline run
# Safe for nohup execution
###############################################################################

set -u  # Exit on undefined variable

# Configuration
PROJECT_ROOT="/data1/wetlab/TRUST4_benchmark"
TRUST4_DIR="${PROJECT_ROOT}/algorithms/TRUST4"
FASTQ_DIR="${PROJECT_ROOT}/data/fastq"
REFERENCE_DIR="${PROJECT_ROOT}/reference/TRUST4"
RESULTS_DIR="${PROJECT_ROOT}/results/TRUST4"
LOG_DIR="${PROJECT_ROOT}/logs/TRUST4"
THREADS=8

# Sample mapping (FZ sample -> SRR accession)
declare -A SAMPLES=(
    [FZ-20]="SRR7882940"
    [FZ-83]="SRR7882939"
    [FZ-94]="SRR7882938"
    [FZ-97]="SRR7882937"
    [FZ-122]="SRR7882935"
)

# Create directories
mkdir -p "${LOG_DIR}"
BATCH_LOG="${LOG_DIR}/batch_runner.log"
SUMMARY_FILE="${LOG_DIR}/batch_summary.tsv"

# Start batch logging
{
    echo "================================================================================"
    echo "TRUST4 Batch Runner - 5 Remaining Samples"
    echo "================================================================================"
    echo "Start Time: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo "Project Root: ${PROJECT_ROOT}"
    echo "TRUST4 Binary: ${TRUST4_DIR}/run-trust4"
    echo "Reference: ${REFERENCE_DIR}/human_IMGT+C.fa"
    echo "Threads per sample: ${THREADS}"
    echo "Samples to process: ${!SAMPLES[@]}"
    echo "================================================================================"
} | tee "${BATCH_LOG}"

# Activate conda environment
{
    echo ""
    echo "Activating conda environment..."
    source /opt/conda/etc/profile.d/conda.sh 2>/dev/null || source ~/miniforge3/etc/profile.d/conda.sh 2>/dev/null
    conda activate trust4_benchmark 2>&1
} >> "${BATCH_LOG}" 2>&1

# Initialize summary file
cat > "${SUMMARY_FILE}" << 'EOF'
sample	sra_run	start_time	end_time	exit_status	report_exists	report_size_bytes	status
EOF

# Process each sample
for SAMPLE in "${!SAMPLES[@]}"; do
    SRR_ID="${SAMPLES[$SAMPLE]}"
    SAMPLE_RESULT_DIR="${RESULTS_DIR}/${SAMPLE}"
    SAMPLE_LOG="${LOG_DIR}/TRUST4_${SAMPLE}.log"

    REPORT_FILE="${SAMPLE_RESULT_DIR}/TRUST_${SAMPLE}_report.tsv"

    START_TIME=$(date -u '+%Y-%m-%d %H:%M:%S UTC')
    START_EPOCH=$(date +%s)

    {
        echo ""
        echo "================================================================================"
        echo "Processing Sample: ${SAMPLE} (${SRR_ID})"
        echo "================================================================================"
        echo "Start: ${START_TIME}"
        echo "Output directory: ${SAMPLE_RESULT_DIR}"
        echo "Log file: ${SAMPLE_LOG}"
        echo "Expected report: ${REPORT_FILE}"
        echo "================================================================================"
        echo ""
    } | tee -a "${BATCH_LOG}"

    # Verify input FASTQ files exist
    R1="${FASTQ_DIR}/${SRR_ID}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR_ID}_2.fastq.gz"

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        END_TIME=$(date -u '+%Y-%m-%d %H:%M:%S UTC')
        {
            echo "ERROR: Input FASTQ files not found:"
            echo "  R1: $R1 ($([ -f "$R1" ] && echo 'EXISTS' || echo 'MISSING'))"
            echo "  R2: $R2 ($([ -f "$R2" ] && echo 'EXISTS' || echo 'MISSING'))"
            echo ""
        } | tee -a "${BATCH_LOG}"

        echo "${SAMPLE}	${SRR_ID}	${START_TIME}	${END_TIME}	ERROR_INPUT	N/A	0	INPUT_FILES_MISSING" >> "${SUMMARY_FILE}"
        continue
    fi

    # Create output directory
    mkdir -p "${SAMPLE_RESULT_DIR}"

    # Run TRUST4
    EXIT_STATUS=0
    {
        cd "${SAMPLE_RESULT_DIR}"

        echo "Executing TRUST4..."
        echo "Command: ${TRUST4_DIR}/run-trust4 \\"
        echo "  -f ${REFERENCE_DIR}/human_IMGT+C.fa \\"
        echo "  --ref ${REFERENCE_DIR}/human_IMGT+C.fa \\"
        echo "  -1 ${R1} \\"
        echo "  -2 ${R2} \\"
        echo "  -t ${THREADS} \\"
        echo "  -o TRUST_${SAMPLE}"
        echo ""

        "${TRUST4_DIR}/run-trust4" \
            -f "${REFERENCE_DIR}/human_IMGT+C.fa" \
            --ref "${REFERENCE_DIR}/human_IMGT+C.fa" \
            -1 "${R1}" \
            -2 "${R2}" \
            -t ${THREADS} \
            -o "TRUST_${SAMPLE}" \
            2>&1 | tee "${SAMPLE_LOG}"

    } || EXIT_STATUS=$?

    # Check for successful completion
    END_TIME=$(date -u '+%Y-%m-%d %H:%M:%S UTC')
    END_EPOCH=$(date +%s)
    ELAPSED=$((END_EPOCH - START_EPOCH))

    {
        echo ""
        echo "Elapsed time: ${ELAPSED} seconds"
        echo "Exit status: ${EXIT_STATUS}"
    } | tee -a "${BATCH_LOG}"

    # Verify output report exists
    REPORT_SIZE=0
    REPORT_EXISTS="NO"
    STATUS="FAILED"

    if [ -f "${REPORT_FILE}" ]; then
        REPORT_SIZE=$(stat -f%z "${REPORT_FILE}" 2>/dev/null || stat -c%s "${REPORT_FILE}" 2>/dev/null || echo 0)
        if [ "${REPORT_SIZE}" -gt 100 ]; then
            REPORT_EXISTS="YES"
            if [ ${EXIT_STATUS} -eq 0 ]; then
                STATUS="SUCCESS"
            else
                STATUS="COMPLETED_WITH_ERROR"
            fi
        else
            STATUS="REPORT_EMPTY"
        fi
    else
        STATUS="REPORT_MISSING"
    fi

    {
        echo "Report file: ${REPORT_FILE}"
        echo "Report exists: ${REPORT_EXISTS}"
        echo "Report size: ${REPORT_SIZE} bytes"
        echo "Status: ${STATUS}"
        echo ""
    } | tee -a "${BATCH_LOG}"

    # Append to summary
    echo "${SAMPLE}	${SRR_ID}	${START_TIME}	${END_TIME}	${EXIT_STATUS}	${REPORT_EXISTS}	${REPORT_SIZE}	${STATUS}" >> "${SUMMARY_FILE}"
done

# Final summary
{
    echo ""
    echo "================================================================================"
    echo "BATCH PROCESSING COMPLETE"
    echo "================================================================================"
    echo "End Time: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo ""
    echo "Summary saved to: ${SUMMARY_FILE}"
    echo ""
    echo "Results:"
} | tee -a "${BATCH_LOG}"

# Print summary table
cat "${SUMMARY_FILE}" | tee -a "${BATCH_LOG}"

# Print status summary
{
    echo ""
    echo "================================================================================"
    TOTAL=$(tail -n +2 "${SUMMARY_FILE}" | wc -l)
    SUCCESS=$(tail -n +2 "${SUMMARY_FILE}" | grep "SUCCESS" | wc -l)
    FAILED=$(tail -n +2 "${SUMMARY_FILE}" | grep -v "SUCCESS" | wc -l)

    echo "Total samples: ${TOTAL}"
    echo "Successful: ${SUCCESS}"
    echo "Failed: ${FAILED}"
    echo "================================================================================"
    echo ""
} | tee -a "${BATCH_LOG}"

# List all output directories
{
    echo "Output directories created:"
    ls -d "${RESULTS_DIR}"/FZ-* 2>/dev/null | sort
    echo ""
    echo "Batch log: ${BATCH_LOG}"
    echo ""
} | tee -a "${BATCH_LOG}"

exit 0
