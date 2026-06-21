mkdir mmseqdb -p
mmseqs createdb ref/references.aa.fasta mmseqdb/refsDB

mkdir tmp
mmseqs cluster \
    mmseqdb/refsDB \
    tmp/clusterDB \
    tmp \
    --min-seq-id 0.7 \
    -c 0.8

mmseqs createtsv \
    mmseqdb/refsDB \
    mmseqdb/refsDB \
    tmp/clusterDB \
    analysis/clusters.tsv

