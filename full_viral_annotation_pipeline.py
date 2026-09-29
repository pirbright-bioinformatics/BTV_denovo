#!/usr/bin/env python3

"""
FULL VIRAL PROTEIN ANNOTATION PIPELINE
=====================================

Features:
- Ontology-aware protein normalization
- Stopword learning + filtering
- Genome-aware validation
- Graph-based clustering with fuzzy matching
* Todo: UniProt ontology generation (optional), at the moment we use a curated ontology curated_ontology.json
* QC metric calculation
* QC dashboard generation
input: tab-separated file with two columns: entry and counts (dictionary of genome:protein -> count)
e.g. WEF49311.1|Bluetongue_virus_1|2|VP2_protein {'2:VP2_protein': 1, '2:major_virus_neutralization_protein': 7, '2:VP2': 13, '2:OC1': 2, '2:outer_capsid_protein_VP2': 3, 'unknown_segment:VP2': 1, 'unknown_segment:VP2_protein': 1}
in analysis/cleaned.tsv, we output the canonical protein name and its count for each entry, along with the list of canonical proteins
Usage:
    python full_pipeline.py input.txt ontology
"""

import re
import ast
import math
import json
from collections import Counter
from rapidfuzz import fuzz
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# CONFIG
# -----------------------------
SIM_THRESHOLD = 85
PARTIAL_THRESHOLD = 80
RARE_THRESHOLD = 2

BASE_STOPWORDS = {
    "protein", "enzyme", "gene", "putative",
    "hypothetical", "like", "family"
}

# -----------------------------
# NORMALIZATION
# -----------------------------
def normalize(text):
    text = text.lower()
    text = text.replace("-", "_")
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^a-z0-9_]", "", text)
    return text


def split_entry(entry):
    genome, protein = entry.split(":", 1)
    return genome, normalize(protein)

# -----------------------------
# STOPWORDS
# -----------------------------
def tokenize(p):
    return p.split("_")


def build_stopwords(proteins):
    counts = Counter()
    for p in proteins:
        counts.update(tokenize(p))

    total = len(proteins)
    dynamic = {t for t, c in counts.items() if c / total > 0.6}
    return BASE_STOPWORDS.union(dynamic)


def remove_stopwords(p, stopwords):
    tokens = [t for t in tokenize(p) if t not in stopwords]
    return "_".join(tokens) if tokens else p

# -----------------------------
# STRUCTURED RULES
# -----------------------------
def extract_structured(p):
    m = re.search(r"(vp|ns|orf)(\d+[a-z]?)", p)
    return m.group() if m else None


def number_conflict(a, b):
    ta = extract_structured(a)
    tb = extract_structured(b)
    if ta and tb and ta != tb and ta[:2] == tb[:2]:
        return True
    return False

# -----------------------------
# SIMILARITY
# -----------------------------
def similar(a, b, stopwords):
    if number_conflict(a, b):
        return False

    a_clean = remove_stopwords(a, stopwords)
    b_clean = remove_stopwords(b, stopwords)

    return (
        fuzz.token_set_ratio(a_clean, b_clean) >= SIM_THRESHOLD or
        fuzz.partial_ratio(a_clean, b_clean) >= PARTIAL_THRESHOLD
    )

# -----------------------------
# GRAPH CLUSTERING
# -----------------------------
def build_graph(names, stopwords):
    graph = {n: set() for n in names}
    names = list(names)

    for i, a in enumerate(names):
        for b in names[i+1:]:
            if similar(a, b, stopwords):
                graph[a].add(b)
                graph[b].add(a)

    return graph


def connected_components(graph):
    visited = set()
    clusters = []

    for node in graph:
        if node in visited:
            continue

        stack = [node]
        comp = []

        while stack:
            n = stack.pop()
            if n in visited:
                continue
            visited.add(n)
            comp.append(n)
            stack.extend(graph[n])

        clusters.append(comp)

    return clusters

# -----------------------------
# ONTOLOGY (simple local)
# -----------------------------
def load_ontology(file="curated_ontology.json"):
    try:
        return json.load(open(file))
    except:
        return {}


def build_lookup(ontology):
    lookup = {}
    for k, v in ontology.items():
        segment = v.get("segment", "unknown_segment")
        k = f"{segment}:{k}"
        lookup[k] = k
        for a in v.get("aliases", []):
            lookup[normalize(a)] = k
    return lookup


def ontology_match(p, lookup):
    if p in lookup:
        return lookup[p], "exact"

    best, best_score = None, 0
    for a in lookup:
        score = fuzz.partial_ratio(p, a)
        if score > best_score:
            best, best_score = a, score

    if best_score > 85:
        return lookup[best], "fuzzy"

    return p, "unknown"

# -----------------------------
# QC
# -----------------------------
def compute_qc(genome_counts, protein_counts, pair_counts, dominant, outliers):
    total = sum(pair_counts.values())

    entropy = 0
    for c in protein_counts.values():
        p = c / total
        entropy -= p * math.log2(p)

    return {
        "genome_frac": genome_counts[dominant] / total,
        "protein_frac": max(protein_counts.values()) / total,
        "entropy": entropy,
        "conflict": (total - genome_counts[dominant]) / total,
        "unknown": genome_counts.get("unknown_segment", 0) / total,
        "outlier": len(outliers) / max(len(pair_counts), 1)
    }

# -----------------------------
# ENTRY PROCESSING
# -----------------------------
def process_entry(entry_counts, stopwords, lookup):

    genome_counts = Counter()
    protein_counts = Counter()
    pair_counts = Counter()

    for raw, count in entry_counts.items():
        g, p = split_entry(raw)

        genome_counts[g] += count
        protein_counts[p] += count
        pair_counts[(g, p)] += count

    dominant = genome_counts.most_common(1)[0][0]

    graph = build_graph(protein_counts.keys(), stopwords)
    clusters = connected_components(graph)
    #print(f"clusters: {clusters}")

    canonical = []
    outliers = []
    canonical_counts = Counter()

    for cluster in clusters:
        #print(f"\tcluster: {cluster}")
        for p in cluster:
            canon, src = ontology_match(p, lookup)
            canonical_counts[canon] += protein_counts[p] 
            #print(f"\t\t{p} -> {canon} ({src})")
            canonical.append(f"{canon}")

            if src == "unknown":
                outliers.append((dominant, p, "unknown"))

    qc = compute_qc(genome_counts, protein_counts, pair_counts, dominant, outliers)

    return canonical, outliers, qc, canonical_counts

# -----------------------------
# DASHBOARD
# -----------------------------
def make_dashboard(qc_file="qc_metrics.tsv"):
    df = pd.read_csv(qc_file, sep="\t")

    for col in df.columns[1:]:
        plt.figure()
        df[col].hist(bins=30)
        plt.title(col)
        plt.savefig(f"qc_{col}.png")
        plt.close()

    plt.scatter(df["genome_frac"], df["protein_frac"], alpha=0.5)
    plt.savefig("qc_scatter.png")

# -----------------------------
# MAIN
# -----------------------------
def run(input_file,onto_file):

    # gather proteins
    proteins = []
    lines = []

    for line in open(input_file):
        if not line.strip():
            continue
        entry, d = line.strip().split("\t", 1)
        #d = "{" + d
        counts = ast.literal_eval(d)

        for k in counts:
            _, p = split_entry(k)
            proteins.append(p)

        lines.append((entry, counts))

    stopwords = build_stopwords(proteins)
    #print("STOPWORDS:", stopwords)

    ontology = load_ontology(onto_file)
    lookup = build_lookup(ontology)
    print(f"lookup: {lookup}")

    qc_rows = []

    with open("cleaned.tsv", "w") as out, \
         open("outliers.tsv", "w") as o, \
         open("qc_metrics.tsv", "w") as q:

        q.write("entry\tgenome_frac\tprotein_frac\tentropy\tconflict\tunknown\toutlier\n")

        for entry, counts in lines:
            canon, outliers, qc, canonical_counts= process_entry(counts, stopwords, lookup)
            canonical_protein = "";canonical_protein_count = 0

            for c in canonical_counts:
                if canonical_counts[c] > canonical_protein_count:
                    canonical_protein = c
                    canonical_protein_count = canonical_counts[c]

            out.write(f"{entry}\t{canonical_protein}\t{canonical_protein_count}\t{canon}\n")

            for g, p, r in outliers:
                o.write(f"{entry}\t{g}:{p}\t{r}\n")

            q.write(f"{entry}\t{qc['genome_frac']}\t{qc['protein_frac']}\t{qc['entropy']}\t{qc['conflict']}\t{qc['unknown']}\t{qc['outlier']}\n")

    make_dashboard()


if __name__ == "__main__":
    import sys
    run(sys.argv[1],sys.argv[2])
