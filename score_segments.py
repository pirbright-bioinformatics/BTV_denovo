#!/usr/bin/env python3

import csv
from collections import defaultdict

diamond_file = "analysis/diamond.tsv"
metadata_file = "analysis/cluster_metadata.tsv"

def clean_id(x):
    return x.split("|")[0]

###############################################################
# Read cluster metadata
###############################################################

metadata = {}

with open(metadata_file) as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata[row["protein_id"]] = {
            "cluster": row["cluster_id"],
            "segment": row["segment"],
            "gene": row["gene"],
            "protein_length": float(row["protein_length"])
        }

###############################################################
# Keep best hit per (contig,segment,cluster)
###############################################################

best_hits = {}

for row in csv.reader(open(diamond_file), delimiter="\t"):
    (
        qseqid,
        sseqid,
        stitle,
        pident,
        aln_length,
        evalue,
        bitscore,
        qstart,
        qend,
        sstart,
        send,
        qlen,
        slen

    ) = row

    sseqid = clean_id(sseqid)
    if sseqid not in metadata:
        continue

    pident = float(pident)
    bitscore = float(bitscore)
    aln_length = float(aln_length)

    cluster = metadata[sseqid]["cluster"]
    gene = metadata[sseqid]["gene"]
    segment = gene.split(":")[0] 
    protein_length = metadata[sseqid]["protein_length"]

    coverage = aln_length / protein_length

    if coverage > 1.0:
        coverage = 1.0

    # higher weight for higher bitscore, pident as a percentage, and coverage
    score = bitscore * (pident / 100.0) * coverage/protein_length
    key = (qseqid, segment, cluster)

    if key not in best_hits:
        best_hits[key] = {
            "score": score,
            "gene": gene,
            "cluster": cluster
        }

    elif score > best_hits[key]["score"]:
        best_hits[key] = {
            "score": score,
            "gene": gene,
            "cluster": cluster
        }

###############################################################
# Aggregate by segment
###############################################################

segment_scores = defaultdict(lambda: defaultdict(float))
segment_genes = defaultdict(lambda: defaultdict(set))
segment_cluster = defaultdict(lambda: defaultdict(str))

for (contig, segment, cluster), info in best_hits.items():
    segment_scores[contig][segment] += info["score"]
    segment_genes[contig][segment].add(info["gene"])
    segment_cluster[contig][segment] = cluster
    

###############################################################
# Produce ranked candidate table
###############################################################

with open("ranked_candidates.tsv", "w") as out:

    out.write("genes\tcontig\tsegment\tscore\tconfidence\tcluster\n")
    
    best_contigs = defaultdict(list)
    best_contigs_stat = defaultdict()
    for contig in sorted(segment_scores):
        total = sum(segment_scores[contig].values())

        ranked = sorted(
            segment_scores[contig].items(),
            key=lambda x: x[1],
            reverse=True
        )

        for segment, score in ranked:
            cluster = segment_cluster[contig][segment]
            confidence = score / total if total > 0 else 0
            genes = ",".join(
                sorted(segment_genes[contig][segment])
            )
            if contig not in best_contigs[genes]:
               best_contigs[genes].append(contig)
               best_contigs_stat[contig] =  (
                                                "{}\t{}\t{}\t{:.3f}\t{:.3f}\t{}".format(
                                                    genes,
                                                    contig,
                                                    segment,
                                                    score,
                                                    confidence,
                                                    cluster
                                                )
                                            )


    for gene,contigs in best_contigs.items():
        for contig in contigs:
            out.write(f"{best_contigs_stat[contig]}\n") 
