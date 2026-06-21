#!/usr/bin/env python3

import csv

def clean_id(x):
    return x.split("|")[0]

# protein metadata
metadata = {}

with open("ref/protein_metadata.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata[row["protein_id"]] = row

canonical_map = {}
with open("analysis/cleaned.tsv") as f:
    reader = csv.reader(f, delimiter="\t")
    for id,canonical,count,cluster in reader:
        canonical_map[clean_id(id)] = canonical
    
with open("analysis/cluster_metadata.tsv", "w") as out:

    out.write("protein_id\tcluster_id\tsegment\tgene\tprotein_length\n")

    with open("analysis/clusters.tsv") as f:
        reader = csv.reader(f, delimiter="\t")

        for rep, member in reader:

            member_id = clean_id(member)
            rep_id = clean_id(rep)

            if member_id not in metadata:
                continue

            m = metadata[member_id]

            out.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    member_id,
                    rep_id,
                    m["segment"],
                    canonical_map[rep_id],
                    m["protein_length"]
                )
            )
