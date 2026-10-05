# Graph-TRUST4-QAOS usage

Build and run the focused tests:

```bash
make -C algorithms/Graph-TRUST4-QAOS
make -C algorithms/Graph-TRUST4-QAOS self-test
tests/graph_qaos/test_qaos.sh
```

Run on paired candidate FASTQs:

```bash
/usr/bin/time -v algorithms/Graph-TRUST4-QAOS/graph-trust4-qaos \
  -1 candidate_R1.fq -2 candidate_R2.fq \
  -o results/Graph-TRUST4-QAOS/test \
  -k 9 -m 31 -t 64 -M 32
```

The output directory contains `assembled_contigs.fa`, graph/node/edge
provenance, `graph_statistics.tsv`, `qaos_statistics.tsv`,
`path_statistics.tsv`, and `runtime.tsv`. Annotate without changing the
assembly:

```bash
algorithms/TRUST4/annotator \
  -f reference/iRep/ng-bcr-validate/iRep/FZ-116.csv.gz \
  -a results/Graph-TRUST4-QAOS/test/assembled_contigs.fa \
  --fasta --needReverseComplement --noImpute --outputFormat 1 \
  > results/Graph-TRUST4-QAOS/test/annotated_contigs.tsv
```

The annotator command is deliberately separate from graph assembly so that
the existing TRUST4 annotation implementation is reused unchanged.
