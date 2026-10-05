#!/usr/bin/env python3
"""Run unchanged TRUST4 annotation in bounded FASTA batches."""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def records(path):
    name = None
    sequence = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(sequence)
                name, sequence = line, []
            elif name is not None:
                sequence.append(line)
    if name is not None:
        yield name, "".join(sequence)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotator", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--log-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=64)
    parser.add_argument("--batch-records", type=int, default=1000)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="trust4-annotation-") as temp_name:
        temp = Path(temp_name)
        output = args.output.open("w")
        header_written = False
        batch_id = 0
        batch = []

        def annotate_batch(rows):
            nonlocal batch_id, header_written
            if not rows:
                return
            fasta_path = temp / f"batch_{batch_id}.fa"
            raw_out = temp / f"batch_{batch_id}.tsv"
            with fasta_path.open("w") as handle:
                for name, sequence in rows:
                    handle.write(f"{name}\n{sequence}\n")
            command = [
                "/usr/bin/time", "-v", str(args.annotator),
                "-f", str(args.reference), "-a", str(fasta_path),
                "--fasta", "-t", str(args.threads),
                "--needReverseComplement", "--noImpute", "--outputFormat", "1",
            ]
            with raw_out.open("w") as stdout, (args.log_dir / f"batch_{batch_id}.time").open("w") as stderr:
                subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
            with raw_out.open() as annotated:
                for line_number, line in enumerate(annotated):
                    if line_number == 0 and header_written:
                        continue
                    output.write(line)
                    if line_number == 0:
                        header_written = True
            batch_id += 1

        try:
            for row in records(args.fasta):
                batch.append(row)
                if len(batch) >= args.batch_records:
                    annotate_batch(batch)
                    batch.clear()
            annotate_batch(batch)
        finally:
            output.close()

    if not header_written:
        raise RuntimeError(f"Annotator produced no AIRR header for {args.fasta}")
    with (args.log_dir / "summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["annotation_batches", batch_id])
        writer.writerow(["records_per_batch", args.batch_records])
        writer.writerow(["records", sum(1 for _ in records(args.fasta))])


if __name__ == "__main__":
    main()
