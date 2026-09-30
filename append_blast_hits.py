#!/usr/bin/env python3

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

BLAST_FIELDS = [
    "search_database",
    "match_type",
    "segment",
    "gene",
    "accession",
    "identity",
    "alignment_length",
    "query_coverage",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "title",
]

HIT_FIELDS = [
    f"{prefix}_{field}"
    for prefix in ("top", "second")
    for field in BLAST_FIELDS
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ranked")
    parser.add_argument("contigs_dir")
    parser.add_argument("output")
    parser.add_argument(
        "--blast-script",
        default=str(Path(__file__).with_name("blast_and_pick.py"))
    )
    args = parser.parse_args()

    ranked_path = Path(args.ranked)
    contigs_dir = Path(args.contigs_dir)
    output_path = Path(args.output)

    with ranked_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        original_fields = reader.fieldnames
        rows = list(reader)

    if not original_fields:
        raise ValueError("Input ranked table has no header")

    fields = original_fields + HIT_FIELDS

    with output_path.open("w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            contig = row["contig"]
            contig_fasta = contigs_dir / f"{contig}.fa"

            if not contig_fasta.exists():
                raise FileNotFoundError(contig_fasta)

            with tempfile.TemporaryDirectory() as tmp:
                blast_output = Path(tmp) / "blast.tsv"

                command = [
                    sys.executable,
                    args.blast_script,
                    str(contig_fasta),
                    "--output",
                    str(blast_output),
                ]

                subprocess.run(
                    command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                with blast_output.open() as f:
                    hits = list(
                        csv.DictReader(f, delimiter="\t")
                    )

            result = dict(row)

            for prefix, hit in zip(
                ("top", "second"), hits[:2]
            ):
                for field in BLAST_FIELDS:
                    result[f"{prefix}_{field}"] = hit.get(
                        field, ""
                    )

            writer.writerow(result)

            print(
                f"{contig}: retained {min(len(hits), 2)} "
                f"BLAST hit(s)",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
