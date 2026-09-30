
#!/usr/bin/env python3

import argparse
import csv
import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

def normalize_accession(accession):
    """Normalize BLAST identifiers such as gb|JX272491.1."""
    accession = accession.strip()

    if "|" in accession:
        parts = accession.split("|")
        if len(parts) > 1:
            accession = parts[1]

    return accession


def read_fasta(path):
    """Yield normalized accession, full header and sequence."""
    header = None
    sequence = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield (
                        normalize_accession(header.split()[0]),
                        header,
                        "".join(sequence),
                    )

                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        yield (
            normalize_accession(header.split()[0]),
            header,
            "".join(sequence),
        )

def write_fasta(records, path):
    with open(path, "w") as out:
        for accession, sequence in records.items():
            out.write(f">{accession}\n")
            for i in range(0, len(sequence), 80):
                out.write(sequence[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Map reads to the top BLAST reference per gene."
    )
    parser.add_argument("ranked_blast")
    parser.add_argument("reference_fasta")
    parser.add_argument("output_pdf")
    parser.add_argument("output_tsv")
    parser.add_argument("--r1", nargs="+", required=True)
    parser.add_argument("--r2", nargs="+", required=True)
    parser.add_argument("--threads", type=int, default=10)
    args = parser.parse_args()

    if len(args.r1) != len(args.r2):
        parser.error("The number of R1 and R2 files must match")

    ranked = pd.read_csv(args.ranked_blast, sep="\t",
                         dtype=str, keep_default_na=False)

    required = {"genes", "top_accession"}
    missing = required - set(ranked.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    # Preserve input order: first row for each genes value.
    first_rows = ranked.drop_duplicates(
        subset=["genes"], keep="first"
    )

    # Select only rows with a usable top accession.
    first_rows = first_rows[
        (first_rows["genes"] != "") &
        (first_rows["top_accession"] != "")
    ]

    if first_rows.empty:
        raise ValueError("No top BLAST accessions found")

    # An accession can occur for more than one gene.
    # Map each accession only once, while retaining gene labels.
    accession_to_genes = OrderedDict()
    for _, row in first_rows.iterrows():
        accession = normalize_accession(row["top_accession"])
        accession_to_genes.setdefault(accession, []).append(
            row["genes"]
        )

    wanted = set(accession_to_genes)
    reference_sequences = {}

    for accession, header, sequence in read_fasta(
        args.reference_fasta
    ):
        if accession in wanted:
            reference_sequences[accession] = sequence

    missing_accessions = wanted - set(reference_sequences)
    if missing_accessions:
        raise ValueError(
            "Accessions absent from reference FASTA: "
            + ", ".join(sorted(missing_accessions))
        )

    output_pdf = Path(args.output_pdf)
    output_tsv = Path(args.output_tsv)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix="reference_coverage_",
        dir=output_pdf.parent,
    ) as tmp:
        tmp = Path(tmp)
        reference_fasta = tmp / "top_references.fasta"
        bam = tmp / "reads.bam"

        write_fasta(reference_sequences, reference_fasta)

        # Align all read pairs against the selected references.
        command = [
            "minimap2",
            "-ax", "sr",
            "-t", str(args.threads),
            str(reference_fasta),
        ]
        for r1, r2 in zip(args.r1, args.r2):
            command.extend([r1, r2])

        with bam.open("wb") as bam_out:
            minimap = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
            )
            sort = subprocess.run(
                ["samtools", "sort", "-@", str(args.threads),
                 "-o", str(bam), "-"],
                stdin=minimap.stdout,
                check=True,
            )
            minimap.stdout.close()
            return_code = minimap.wait()
            if return_code:
                raise subprocess.CalledProcessError(
                    return_code, command
                )

        subprocess.run(
            ["samtools", "index", str(bam)],
            check=True,
        )

        depth_file = tmp / "depth.tsv"
        with depth_file.open("w") as out:
            subprocess.run(
                ["samtools", "depth", "-aa", str(bam)],
                stdout=out,
                check=True,
            )

        depth = pd.read_csv(
            depth_file,
            sep="\t",
            names=["accession", "position", "depth"],
            dtype={"accession": str, "position": int,
                   "depth": int},
        )

    # Produce a row for every reference position, including zero depth.
    depth_by_accession = {
        accession: group.sort_values("position")
        for accession, group in depth.groupby("accession")
    }

    all_depth_rows = []

    with PdfPages(output_pdf) as pdf:
        for _, row in first_rows.iterrows():
            gene = row["genes"]
            accession = normalize_accession(row["top_accession"])
            sequence_length = len(reference_sequences[accession])

            gene_depth = depth_by_accession.get(
                accession,
                pd.DataFrame(columns=["position", "depth"]),
            )

            # samtools depth -aa should emit all positions, but
            # reindex to ensure a complete 1..length series.
            series = gene_depth.set_index("position")[
                "depth"
            ].reindex(range(1, sequence_length + 1), fill_value=0)

            plot_df = pd.DataFrame({
                "gene": gene,
                "accession": accession,
                "position": range(1, sequence_length + 1),
                "depth": series.to_numpy(),
            })
            all_depth_rows.append(plot_df)

            breadth = (plot_df["depth"] > 0).mean()
            mean_depth = plot_df["depth"].mean()
            median_depth = plot_df["depth"].median()

            fig, ax = plt.subplots(figsize=(12, 4.5))
            ax.plot(
                plot_df["position"],
                plot_df["depth"],
                linewidth=0.6,
            )
            ax.set_title(
                f"{gene} | {accession}\n"
                f"Length={sequence_length:,} bp  "
                f"Breadth={breadth:.1%}  "
                f"Mean depth={mean_depth:.1f}  "
                f"Median depth={median_depth:.1f}"
            )
            ax.set_xlabel("Reference position")
            ax.set_ylabel("Read depth")
            ax.grid(alpha=0.3, linewidth=0.5)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    pd.concat(all_depth_rows, ignore_index=True).to_csv(
        output_tsv, sep="\t", index=False
    )

    print(f"Wrote coverage plot: {output_pdf}")
    print(f"Wrote coverage data: {output_tsv}")


if __name__ == "__main__":
    main()
