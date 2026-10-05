#!/bin/bash

cd /data1/wetlab/TRUST4_benchmark

mkdir -p results/Graph-TRUST4-v2-memoryfix/benchmarks

LOG_FILE="results/Graph-TRUST4-v2-memoryfix/benchmark_all.log"
RESULTS_FILE="results/Graph-TRUST4-v2-memoryfix/benchmark_threads.tsv"

{
    echo "threads,runtime_seconds,peak_rss_mb,candidate_pairs,accepted_edges" > "$RESULTS_FILE"

    for threads in 1 8 16 32 64; do
        echo "====== BENCHMARK: $threads threads ======"
        echo ""

        OUTDIR="results/Graph-TRUST4-v2-memoryfix/benchmark_${threads}t"
        mkdir -p "$OUTDIR"

        /usr/bin/time -v \
          algorithms/Graph-TRUST4-v2-memoryfix/graph-trust4-v2-memoryfix \
          -1 test_data/FZ116_test_1.fq \
          -2 test_data/FZ116_test_2.fq \
          -o "$OUTDIR/test" \
          -k 9 -m 31 -i 0.90 -q 0.90 -H 20 \
          -c 500 -M 32 -t $threads \
          > "$OUTDIR/run.log" 2>&1

        # Extract metrics from timing output
        if grep -q "Elapsed (wall clock)" "$OUTDIR/run.log"; then
            RUNTIME=$(grep "Elapsed (wall clock)" "$OUTDIR/run.log" | awk '{print $6}')
            # Convert m:ss to seconds
            if [[ $RUNTIME == *:* ]]; then
                RUNTIME=$(echo "$RUNTIME" | awk -F: '{print int($1)*60+$2}')
            fi
        else
            RUNTIME="0"
        fi

        PEAK_RSS=$(grep "Maximum resident set size" "$OUTDIR/run.log" | awk '{print int($6/1024)}')
        PAIRS=$(grep "Generated.*pairs" "$OUTDIR/run.log" | head -1 | awk '{print $2}' | sed 's/,//g')

        # Try to get edges from graph stats
        EDGES=0
        if [ -f "$OUTDIR/test_graph_stats.tsv" ]; then
            EDGES=$(grep "total_edges" "$OUTDIR/test_graph_stats.tsv" | awk '{print $2}')
        fi

        echo "$threads,$RUNTIME,$PEAK_RSS,$PAIRS,$EDGES" >> "$RESULTS_FILE"

        echo "Completed: $threads threads"
        echo "  Runtime: $RUNTIME sec"
        echo "  Peak RSS: $PEAK_RSS MB"
        echo "  Pairs: $PAIRS"
        echo "  Edges: $EDGES"
        echo ""
    done

    echo "====== BENCHMARK RESULTS ======"
    cat "$RESULTS_FILE"

} | tee "$LOG_FILE"
