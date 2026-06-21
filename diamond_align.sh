mkdir diamonddb
diamond makedb \
    --in ref/references.aa.fasta \
    -d diamonddb/viral_proteins

mkdir analysis -p
diamond blastx \
    --query contig.fa \
    --db diamonddb/viral_proteins \
    --out analysis/diamond.tsv \
    --outfmt 6 \
    qseqid sseqid stitle pident length evalue bitscore \
    qstart qend sstart send \
    qlen slen \
    --evalue 1e-5 \
    --max-target-seqs 100 \
#    --top 10
