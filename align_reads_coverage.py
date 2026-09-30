#!/usr/bin/env python3

import argparse
import csv
import html
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def accession_candidates(value):
    """Generate common accession representations for FASTA lookup."""
    if value is None or pd.isna(value):
        return []

    value = str(value).strip()
    if not value or value.lower() in {"nan", "none", "na"}:
        return []

    # BLAST identifiers may look like gb|KP821098.1| or ref|NC_...|.
    parts = [p for p in value.split("|") if p]
    candidates = [value]

    candidates.extend(parts)
    if parts:
        candidates.append(parts[-1])

    # Also try accession without its version suffix.
    expanded = []
    for candidate in candidates:
        candidate = candidate.strip()
        if candidate:
            expanded.append(candidate)
            expanded.append(re.sub(r"\.\d+$", "", candidate))

    return list(dict.fromkeys(expanded))


def read_fasta(path):
    """Read FASTA records, keyed by their first header token."""
    records = {}
    header = None
    sequence = []

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence)

                header = line[1:].split()[0]
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence)

    return records


def build_accession_index(records):
    """Index FASTA records by header token and normalized accession."""
    index = {}

    for header, sequence in records.items():
        for candidate in accession_candidates(header):
            index.setdefault(candidate, (header, sequence))

    return index


def find_reference(accession, accession_index):
    for candidate in accession_candidates(accession):
        if candidate in accession_index:
            return accession_index[candidate]
    return None


def run_command(command, stdout=None):
    subprocess.run(
        [str(arg) for arg in command],
        check=True,
        stdout=stdout,
    )
import os

def make_depth_plot(depth_tsv, output_png, title, reference_length=1000):
    positions = []
    depths = []

    with open(depth_tsv, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                positions.append(int(fields[1]))
                depths.append(int(fields[2]))
            except ValueError:
                continue

    if not positions:
        print(f"WARNING: No depth values produced for {title}; plotting zero coverage.")
        positions = [1, reference_length]
        depths = [0, 0]

    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.fill_between(positions, depths, step="mid", alpha=0.35)
    ax.plot(positions, depths, linewidth=0.7)
    ax.set_xlim(1, max(positions))
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Reference position (bp)")
    ax.set_ylabel("Read depth")
    ax.set_title(title)
    ax.grid(True, linewidth=0.3, alpha=0.5)

    fig.tight_layout()

    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    fig.savefig(output_png, dpi=150)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description="Align paired reads to BLAST top-accession references "
                    "and generate full-reference depth plots."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--hits", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--reads1", required=True)
    parser.add_argument("--reads2", required=True)
    parser.add_argument("--plot-dir", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    hits_path = Path(args.hits)
    reads1 = Path(args.reads1)
    reads2 = Path(args.reads2)
    plot_dir = Path(args.plot_dir)
    manifest_path = Path(args.manifest)

    for path in (hits_path, reads1, reads2, Path(args.reference_fasta)):
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")

    hits = pd.read_csv(hits_path, sep="\t", dtype=str).fillna("")

    required = {"contig", "top_accession"}
    missing_columns = required - set(hits.columns)
    if missing_columns:
        raise ValueError(
            f"{hits_path}: missing columns {sorted(missing_columns)}"
        )

    if "genes" not in hits.columns:
        hits["genes"] = ""

    # Retain the first row per contig, consistent with the report workflow.
    hits["contig"] = hits["contig"].astype(str).str.strip()
    hits = hits[hits["contig"].ne("")].drop_duplicates(
        subset="contig", keep="first"
    )

    records = read_fasta(args.reference_fasta)
    accession_index = build_accession_index(records)

    plot_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    missing_accessions = []

    # Align once per unique accession, then create a plot for each contig.
    by_accession = {}
    for _, row in hits.iterrows():
        accession = row["top_accession"].strip()
        by_accession.setdefault(accession, []).append(row)

    with tempfile.TemporaryDirectory(prefix="coverage_") as temp_name:
        temp_dir = Path(temp_name)

        for accession, rows in by_accession.items():
            reference = find_reference(accession, accession_index)

            if not accession or reference is None:
                missing_accessions.append(accession or "(empty)")
                for row in rows:
                    manifest_rows.append({
                        "sample": args.sample,
                        "genes": row["genes"],
                        "contig": row["contig"],
                        "top_accession": accession,
                        "coverage_plot": "",
                        "status": "missing_accession",
                    })
                continue

            fasta_header, sequence = reference
            safe_id = "ref"
            reference_fasta = temp_dir / f"{safe_id}.fa"
            sam_file = temp_dir / f"{safe_id}.sam"
            bam_file = temp_dir / f"{safe_id}.bam"
            sorted_bam = temp_dir / f"{safe_id}.sorted.bam"
            depth_file = temp_dir / f"{safe_id}.depth.tsv"
            reference_fasta.write_text(
                f">{safe_id}\n" +
                "\n".join(
                    sequence[i:i + 80]
                    for i in range(0, len(sequence), 80)
                ) + "\n",
                encoding="utf-8",
            )

            # Short-read paired-end alignment.
            with open(sam_file, "w", encoding="utf-8") as sam_out:
                run_command([
                    "minimap2", "-ax", "sr",
                    reference_fasta,
                    reads1,
                    reads2,
                ], stdout=sam_out)

            run_command(["samtools", "view", "-b", "-o", bam_file, sam_file])
            run_command(["samtools", "sort", "-o", sorted_bam, bam_file])

            # -aa emits zero-depth positions; -d 0 disables depth cap.
            with open(depth_file, "w", encoding="utf-8") as depth_out:
                run_command([
                    "samtools", "depth", "-aa", "-d", "0", sorted_bam
                ], stdout=depth_out)

            for row in rows:
                contig = row["contig"]
                # Keep filenames portable and prevent path traversal.
                safe_contig = re.sub(r"[^A-Za-z0-9_.-]", "_", contig)
                relative_plot = (
                    Path("coverage_plots") / args.sample /
                    f"{safe_contig}.png"
                )
                output_png = plot_dir / f"{safe_contig}.png"

                make_depth_plot(
                    depth_file,
                    output_png,
                    f"{args.sample} | {contig} | {accession}",
                )

                manifest_rows.append({
                    "sample": args.sample,
                    "genes": row["genes"],
                    "contig": contig,
                    "top_accession": accession,
                    "coverage_plot": relative_plot.as_posix(),
                    "status": "ok",
                })

    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample", "genes", "contig", "top_accession",
                "coverage_plot", "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    if missing_accessions:
        print(
            "Warning: reference accessions not found in FASTA: "
            + ", ".join(sorted(set(missing_accessions)))
        )


if __name__ == "__main__":
    main()
