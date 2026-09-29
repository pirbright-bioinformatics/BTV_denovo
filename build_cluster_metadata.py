#!/usr/bin/env python3

import csv
import sys

def clean_id(x):
    return x.split("|")[0]

# protein metadata
metadata = {}

with open(sys.argv[2]) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata[row["protein_id"]] = row

canonical_map = {}
with open(sys.argv[3]) as f:
    reader = csv.reader(f, delimiter="\t")
    for id,canonical,count,cluster in reader:
        canonical_map[clean_id(id)] = canonical
    
with open(sys.argv[4], "w") as out:

    out.write("protein_id\tcluster_id\tsegment\tgene\tprotein_length\n")

    with open(sys.argv[1]) as f:
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
