#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

INPUT = sys.argv[1] 
OUTPUT = sys.argv[2]

# Remove obvious length outliers
ABS_Z_DISCARD = 3.5

# Maximum number of contigs to retain per gene
MAX_CONTIGS = 3

# Read whitespace-delimited table
df = pd.read_csv(INPUT, sep="\t")
df = df[ df["top_title"].str.contains( "Bluetongue", case=False, na=False) ]

selected = []

for gene, group in df.groupby("genes", sort=False):

    group = group.copy()

    # -----------------------------
    # Stage 1: Remove extreme length outliers
    # -----------------------------
    filtered = group[group["z"].abs() <= ABS_Z_DISCARD].copy()

    # If everything was removed, keep the least extreme contig
    if filtered.empty:
        filtered = group.loc[[group["z"].abs().idxmin()]]

    # -----------------------------
    # Stage 2: Keep all if <= MAX_CONTIGS
    # -----------------------------
    if len(filtered) <= MAX_CONTIGS:
        selected.append(filtered)
        continue

    # -----------------------------
    # Stage 3: Rank remaining contigs
    # -----------------------------
    filtered["rank_score"] = (
        filtered["confidence"] /
        (1 + filtered["z"].abs())
    )

    filtered = filtered.sort_values(
        by=["rank_score", "confidence"],
        ascending=False
    )

    selected.append(filtered.head(MAX_CONTIGS))

result = pd.concat(selected)

# Remove helper column if present
if "rank_score" in result.columns:
    result = result.drop(columns="rank_score")

result.to_csv(OUTPUT, sep="\t", index=False)

print(f"Selected {len(result)} contigs.")
print(f"Results written to {OUTPUT}")
