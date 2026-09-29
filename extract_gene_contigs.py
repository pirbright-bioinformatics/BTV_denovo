import pandas as pd
from Bio import SeqIO
from pathlib import Path
import sys

CAND = sys.argv[1] #"selected_candidates.tsv"
FASTA = sys.argv[2] #"all_contigs.fa"
OUTDIR = Path(sys.argv[3])
OUTDIR.mkdir(exist_ok=True)

df = pd.read_csv(CAND, sep="\t")

records = SeqIO.to_dict(SeqIO.parse(FASTA, "fasta"))

for gene, group in df.groupby("genes"):

    out_path = OUTDIR / f"{gene}.fa"

    with open(out_path, "w") as out:

        for contig in group["contig"]:
            SeqIO.write(records[contig], out, "fasta")

print("FASTA files created per gene.")
