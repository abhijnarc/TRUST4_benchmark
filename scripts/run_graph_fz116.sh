#!/bin/bash

###############################################################################
# Graph-TRUST4 FZ-116 Pilot Runner
#
# Background-safe execution script for full FZ-116 sample
# Safe to use with nohup and survives SSH disconnection
###############################################################################

set -u

# Configuration
PROJECT_ROOT="/data1/wetlab/TRUST4_benchmark"
GRAPH_BINARY="${PROJECT_ROOT}/algorithms/Graph-TRUST4/graph-trust4"
CANDIDATE_R1="${PROJECT_ROOT}/results/TRUST4/FZ-116/TRUST_FZ-116_toassemble_1.fq"
CANDIDATE_R2="${PROJECT_ROOT}/results/TRUST4/FZ-116/TRUST_FZ-116_toassemble_2.fq"
OUTPUT_DIR="${PROJECT_ROOT}/results/Graph-TRUST4/FZ-116"
LOG_DIR="${PROJECT_ROOT}/logs/Graph-TRUST4"
OUTPUT_PREFIX="${OUTPUT_DIR}/graph_fz116"

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# Timestamp
START_TIME=$(date -u '+%Y-%m-%d %H:%M:%S UTC')
START_EPOCH=$(date +%s)

{
    echo "================================================================================"
    echo "Graph-TRUST4 FZ-116 Pilot Execution"
    echo "================================================================================"
    echo "Start Time: ${START_TIME}"
    echo "Project Root: ${PROJECT_ROOT}"
    echo "Binary: ${GRAPH_BINARY}"
    echo "Candidate R1: ${CANDIDATE_R1}"
    echo "Candidate R2: ${CANDIDATE_R2}"
    echo "Output Directory: ${OUTPUT_DIR}"
    echo "Output Prefix: ${OUTPUT_PREFIX}"
    echo "================================================================================"
    echo ""

    # Verify inputs
    if [ ! -f "${GRAPH_BINARY}" ]; then
        echo "ERROR: Graph binary not found: ${GRAPH_BINARY}"
        exit 1
    fi

    if [ ! -f "${CANDIDATE_R1}" ] || [ ! -f "${CANDIDATE_R2}" ]; then
        echo "ERROR: Candidate read files not found"
        exit 1
    fi

    echo "Input verification: OK"
    echo "Candidate R1: $(wc -l < ${CANDIDATE_R1}) lines"
    echo "Candidate R2: $(wc -l < ${CANDIDATE_R2}) lines"
    echo ""
    echo "Starting assembly..."
    echo "================================================================================"
    echo ""

    # Run Graph-TRUST4 with resource monitoring
    /usr/bin/time -v \
        "${GRAPH_BINARY}" \
            -1 "${CANDIDATE_R1}" \
            -2 "${CANDIDATE_R2}" \
            -o "${OUTPUT_PREFIX}" \
            -k 9 \
            -m 20 \
            -i 0.90 2>&1

    EXIT_STATUS=$?

    END_TIME=$(date -u '+%Y-%m-%d %H:%M:%S UTC')
    END_EPOCH=$(date +%s)
    ELAPSED=$((END_EPOCH - START_EPOCH))

    echo ""
    echo "================================================================================"
    echo "Execution Complete"
    echo "================================================================================"
    echo "End Time: ${END_TIME}"
    echo "Elapsed: ${ELAPSED} seconds"
    echo "Exit Status: ${EXIT_STATUS}"
    echo ""

    # Verify outputs
    if [ ${EXIT_STATUS} -eq 0 ]; then
        echo "Output verification:"
        if [ -f "${OUTPUT_PREFIX}_contigs.fa" ]; then
            CONTIG_COUNT=$(grep -c "^>" "${OUTPUT_PREFIX}_contigs.fa")
            echo "  Contigs FASTA: $(wc -c < ${OUTPUT_PREFIX}_contigs.fa) bytes (${CONTIG_COUNT} contigs)"
        fi
        if [ -f "${OUTPUT_PREFIX}_stats.tsv" ]; then
            echo "  Statistics TSV: $(wc -c < ${OUTPUT_PREFIX}_stats.tsv) bytes"
        fi
    else
        echo "ERROR: Assembly failed with exit status ${EXIT_STATUS}"
    fi

    echo ""
    echo "================================================================================"

} | tee "${LOG_DIR}/graph_fz116_execution.log"

# Record timing
{
    echo "sample	status	start_time	end_time	elapsed_seconds	exit_status"
    echo "FZ-116	$([ $EXIT_STATUS -eq 0 ] && echo 'SUCCESS' || echo 'FAILED')	${START_TIME}	${END_TIME}	${ELAPSED}	${EXIT_STATUS}"
} > "${LOG_DIR}/graph_fz116_timing.tsv"

# Create config record
{
    echo "GRAPH_CONFIG"
    echo "sample	FZ-116"
    echo "kmer_size	9"
    echo "min_overlap_length	20"
    echo "min_identity_threshold	0.90"
    echo "candidate_input_r1	${CANDIDATE_R1}"
    echo "candidate_input_r2	${CANDIDATE_R2}"
    echo "output_prefix	${OUTPUT_PREFIX}"
} > "${OUTPUT_DIR}/graph_config.txt"

exit ${EXIT_STATUS}
