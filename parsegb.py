#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
import csv
import sys


gb_dir = Path(sys.argv[1])

nt_out = open(sys.argv[2], "w")
aa_out = open(sys.argv[3], "w")
meta = open(sys.argv[4], "w", newline="")
writer = csv.writer(meta, delimiter="\t")

writer.writerow([
    "protein_id",
    "virus",
    "segment",
    "gene",
    "protein_length",
    "gb_accession"
])

for gb_file in gb_dir.glob("*.gb*"):
    for record in SeqIO.parse(gb_file, "genbank"):

        organism = record.annotations.get("organism", "unknown")
        accession = record.id

        # Extract segment information if available
        segment = "unknown_segment"
        for feature in record.features:
            if feature.type == "source":
                segment = feature.qualifiers.get("segment", [segment])[0]
                break

        # Write nucleotide sequence
        nt_header = f"{accession}|{organism}|{segment}"
        nt_out.write(f">{nt_header}\n{record.seq}\n")

        segment = "unknown_segment"

        for feature in record.features:
            if feature.type == "source":
                segment = feature.qualifiers.get("segment", ["unknown_segment"])[0]
                break

        # Extract CDS proteins
        for i, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers

            product = qualifiers.get("product", ["unknown_protein"])[0]

            protein_id = qualifiers.get(
                "protein_id",
                [f"{accession}_cds{i}"]
            )[0]

            translation = qualifiers.get("translation")[0]

            if translation:
                aa_header = (
                    f"{protein_id}|{organism}|"
                    f"{segment}|{product}"
                ).replace(" ","_")
                aa_out.write(
                    f">{aa_header}\n{translation}\n"
                )
                product = qualifiers.get("product", ["unknown"])[0]

                writer.writerow([
                    protein_id,
                    record.annotations.get("organism", "Unknown"),
                    segment,
                    product,
                    len(translation),
                    record.id
                ])

nt_out.close()
aa_out.close()
