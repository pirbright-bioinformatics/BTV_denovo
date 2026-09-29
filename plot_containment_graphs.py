#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import SeqIO
import subprocess

# -----------------------------
# Inputs
# -----------------------------

FASTA_DIR = Path("gene_fastas")
OUTDIR = Path("containment_maps")
OUTDIR.mkdir(exist_ok=True)

TMPDIR = Path("tmp_paf")
TMPDIR.mkdir(exist_ok=True)

# -----------------------------
# Helper: run minimap2
# -----------------------------

def run_minimap(ref, query, out_paf):

    cmd = [
        "minimap2",
        "-x", "asm5",
        ref,
        query
    ]

    with open(out_paf, "w") as f:
        subprocess.run(cmd, stdout=f, check=True)

# -----------------------------
# Parse PAF
# -----------------------------

def parse_paf(paf_file):

    cols = [
        "qname","qlen","qstart","qend","strand",
        "tname","tlen","tstart","tend",
        "matches","alen","mapq"
    ]

    return pd.read_csv(paf_file, sep="\t", header=None, names=cols)

# -----------------------------
# Main loop
# -----------------------------

for fa in FASTA_DIR.glob("*.fa"):

    gene = fa.stem

    records = list(SeqIO.parse(fa, "fasta"))

    if len(records) < 2:
        continue

    # -----------------------------
    # Pick longest contig as reference
    # -----------------------------

    lengths = {r.id: len(r.seq) for r in records}
    ref = max(lengths, key=lengths.get)

    ref_fa = TMPDIR / f"{gene}_ref.fa"
    qry_fa = TMPDIR / f"{gene}_qry.fa"
    paf_file = TMPDIR / f"{gene}.paf"

    SeqIO.write(
        [r for r in records if r.id == ref],
        ref_fa,
        "fasta"
    )

    SeqIO.write(records, qry_fa, "fasta")

    # -----------------------------
    # Align
    # -----------------------------

    run_minimap(str(ref_fa), str(qry_fa), paf_file)

    paf = parse_paf(paf_file)

    # remove self alignment
    paf = paf[paf.qname != paf.tname]

    # -----------------------------
    # Plot
    # -----------------------------

    fig, ax = plt.subplots(figsize=(12, 1.2 * len(records)))

    y_positions = {}
    y = len(records)

    for r in sorted(records, key=lambda x: -len(x.seq)):

        y_positions[r.id] = y
        y -= 1

    ref_len = lengths[ref]

    # draw reference
    ax.plot([0, ref_len], [y_positions[ref], y_positions[ref]],
            color="black", lw=4)
    ax.text(-0.02 * ref_len, y_positions[ref], ref,
            ha="right", va="center", fontsize=10)

    # draw others
    for r in records:

        if r.id == ref:
            continue

        y = y_positions[r.id]

        # full contig
        ax.plot([0, len(r.seq)], [y, y],
                color="lightgray", lw=4)

        hits = paf[paf.qname == r.id]

        for _, h in hits.iterrows():

            start = h.tstart
            end = h.tend

            color = "steelblue" if h.strand == "+" else "darkorange"

            ax.plot([start, end], [y, y],
                    color=color, lw=8, solid_capstyle="butt")

        ax.text(-0.02 * ref_len, y, r.id,
                ha="right", va="center", fontsize=9)

    ax.set_xlim(0, ref_len * 1.02)
    ax.set_ylim(0.5, len(records) + 0.5)

    ax.set_yticks([])
    ax.set_xlabel("Position on longest contig (bp)")
    ax.set_title(f"{gene} containment map")

    plt.tight_layout()

    out_pdf = OUTDIR / f"{gene}.pdf"
    plt.savefig(out_pdf)
    plt.close()

    print(f"Saved {out_pdf}")

print("Done.")
