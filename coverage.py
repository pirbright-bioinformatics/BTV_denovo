#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path


DEPTH_DIR = sys.argv[1] 
RANKED_TSV = sys.argv[2]
OUTPUT_PDF = sys.argv[3]


def load_depth(depth_file):

    return pd.read_csv(
        depth_file,
        sep="\t",
        names=["contig", "position", "depth"]
    )


def coverage_stats(df):

    length = len(df)

    breadth = (df.depth > 0).sum() / length

    mean_depth = df.depth.mean()

    median_depth = df.depth.median()

    edge = max(10, int(length * 0.05))

    five_prime = df.head(edge).depth.mean()

    three_prime = df.tail(edge).depth.mean()

    return {
        "length": length,
        "breadth": breadth,
        "mean_depth": mean_depth,
        "median_depth": median_depth,
        "five_prime_depth": five_prime,
        "three_prime_depth": three_prime,
    }


ranked = pd.read_csv(
    RANKED_TSV,
    sep=r"\s+"

)

all_stats = []

for _, row in ranked.iterrows():

    contig = row["contig"]

    depth_file = Path(
        DEPTH_DIR,
        f"{contig}.fa.depth.txt"
    )

    if not depth_file.exists():
        print("No Depth{}".format(DEPTH_DIR))
        continue

    depth = load_depth(depth_file)

    stats = coverage_stats(depth)

    stats["contig"] = contig

    all_stats.append(stats)

stats_df = pd.DataFrame(all_stats)

ranked = ranked.merge(
    stats_df,
    on="contig",
    how="left"
)

with PdfPages(OUTPUT_PDF) as pdf:

    #
    # One section per gene
    #
    for gene in sorted(ranked["genes"].unique()):
        gene_df = (
            ranked[ranked["genes"] == gene]
            .sort_values(
                ["confidence", "score"],
                ascending=False
            )
        )
        print(gene_df)

        #
        # Summary page
        #
        fig, ax = plt.subplots(
            figsize=(11, 8.5)
        )

        ax.axis("off")

        ax.set_title(
            f"{gene}\nCandidate contigs",
            fontsize=16,
            pad=20
        )

        table_df = gene_df[
            [
                "contig",
                "segment",
                "score",
                "confidence",
                "length",
                "breadth",
                "median_depth",
                "five_prime_depth",
                "three_prime_depth"
            ]
        ].copy()

        table_df = table_df.round(3)

        table = ax.table(
            cellText=table_df.values,
            colLabels=table_df.columns,
            loc="center"
        )

        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.5)

        pdf.savefig(fig)
        plt.close()

        #
        # Coverage pages
        #
        for _, row in gene_df.iterrows():

            contig = row["contig"]

            depth_file = Path(
                DEPTH_DIR,
                f"{contig}.fa.depth.txt"
            )

            if not depth_file.exists():
                continue

            depth = load_depth(depth_file)

            fig, ax = plt.subplots(
                figsize=(11, 4)
            )

            ax.plot(
                depth.position,
                depth.depth,
                linewidth=0.6
            )

            ax.set_title(
                f"{gene} | {contig}\n"
                f"score={row.score:.3f}  "
                f"confidence={row.confidence:.3f}  "
                f"breadth={row.breadth:.3f}  "
                f"median={row.median_depth:.1f}"
            )

            ax.set_xlabel("Position")
            ax.set_ylabel("Depth")

            ax.grid(
                alpha=0.3,
                linewidth=0.5
            )

            pdf.savefig(fig)
            plt.close()

print(f"Wrote {OUTPUT_PDF}")
