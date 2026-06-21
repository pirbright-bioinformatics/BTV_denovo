# =============================================================================
# Snakemake Pipeline: Viral Genome Segment Matching
# =============================================================================
#
# Usage:
#   snakemake --configfile config.yaml --cores <N>
#
# For a dry-run:
#   snakemake --configfile config.yaml --cores 1 -n
#
# One-time index building only:
#   snakemake --configfile config.yaml --cores <N> build_index
#
# Directory layout produced
# ─────────────────────────
# ref/                          ← parsegb.py outputs (one-time)
#   references.nt.fasta
#   references.aa.fasta
#   protein_metadata.tsv
#   clusters.tsv                ← mmseqs cluster output (shared)
#
# diamonddb/                    ← diamond makedb output (one-time)
#   viral_proteins.dmnd
#
# mmseqs_tmp/                   ← mmseqs working files (one-time)
#
# {sample}/                     ← per-sample outputs
#   contig.fa                   ← INPUT (must exist)
#   analysis/
#     diamond.tsv
#     seg_miss.txt
#     protein_miss.txt
#     cluster_metadata.tsv
#     cleaned.tsv
#     outliers.tsv
#     qc_metrics.tsv
#     qc_*.png
#   ranked_candidates.tsv
# =============================================================================

import os

configfile: "config.yaml"

SAMPLES        = config["samples"]
GB_DIR         = config["paths"]["gb_dir"]
REF_DIR        = config["paths"]["ref_dir"]
DIAMOND_DB_DIR = config["paths"]["diamond_db_dir"]
MMSEQS_DB_DIR  = config["paths"]["mmseqs_db_dir"]
MMSEQS_TMP     = config["paths"]["mmseqs_tmp_dir"]
CLUSTERS_TSV   = config["paths"]["clusters_tsv"]
SAMPLEDIR      = config["paths"]["sample_dir"]
ONTOLOGY       = config["ontology"]

THREADS_DIAMOND = config["tools"]["diamond_threads"]
THREADS_MMSEQS  = config["tools"]["mmseqs_threads"]
DIAMOND_EVALUE  = config["tools"]["diamond_evalue"]
DIAMOND_MAXHITS = config["tools"]["diamond_max_target_seqs"]
MMSEQS_MINID    = config["tools"]["mmseqs_min_seq_id"]
MMSEQS_COV      = config["tools"]["mmseqs_coverage"]

# Scripts live next to the Snakefile
SCRIPT_DIR = workflow.basedir
LOGDIR = workflow.basedir + "/logs"
# =============================================================================
# Helper: resolve the ontology path relative to a sample dir
# =============================================================================
def ontology_for_sample(sample):
    """Return the ontology path as seen from {sample}/ working directory."""
    p = ONTOLOGY
    if not os.path.isabs(p):
        # Make it absolute so it survives chdir into sample dir
        p = os.path.normpath(os.path.join(workflow.basedir, sample, p))
    return p

def log_path(*parts):
    return os.path.join(LOGDIR, *parts)

# =============================================================================
# Target rules
# =============================================================================

rule all:
    """Build the index and run the full per-sample pipeline for every sample."""
    input:
        # Index sentinels
        DIAMOND_DB_DIR + "/viral_proteins.dmnd",
        CLUSTERS_TSV,
        # Per-sample final outputs
        expand(SAMPLEDIR + "/{sample}/ranked_candidates.tsv", sample=SAMPLES),


rule build_index:
    """Convenience target: build all index files without running per-sample steps."""
    input:
        DIAMOND_DB_DIR + "/viral_proteins.dmnd",
        CLUSTERS_TSV,


# =============================================================================
# ONE-TIME INDEX RULES
# =============================================================================

rule parsegb:
    """
    Parse GenBank files → ref FASTAs + protein_metadata.tsv.
    Triggered once; sentinel file marks completion.
    """
    input:
        gb_dir = GB_DIR,
    output:
        nt_fasta  = REF_DIR + "/references.nt.fasta",
        aa_fasta  = REF_DIR + "/references.aa.fasta",
        metadata  = REF_DIR + "/protein_metadata.tsv",
    log:
        log_path("parsegb.log")
    shell:
        """
        mkdir -p {REF_DIR}
        python {SCRIPT_DIR}/parsegb.py {input.gb_dir} {output.nt_fasta} {output.aa_fasta} {output.metadata} > {log} 2>&1
        """


rule diamond_makedb:
    """Build Diamond protein database from reference AA FASTA (one-time)."""
    input:
        aa_fasta = REF_DIR + "/references.aa.fasta",
    output:
        db = DIAMOND_DB_DIR + "/viral_proteins.dmnd",
    log:
        log_path("diamond_makedb.log")
    threads: THREADS_DIAMOND
    shell:
        """
        mkdir -p {DIAMOND_DB_DIR}
        diamond makedb \
            --in {input.aa_fasta} \
            -d {DIAMOND_DB_DIR}/viral_proteins \
            --threads {threads} \
            > {log} 2>&1
        """


rule mmseqs_cluster:
    """
    Build MMseqs2 reference DB and cluster it (one-time).
    Writes global clusters.tsv shared by all samples.
    """
    input:
        aa_fasta = REF_DIR + "/references.aa.fasta",
    output:
        clusters = CLUSTERS_TSV,
    log:
        log_path("mmseqs_cluster.log")
    threads: THREADS_MMSEQS
    shell:
        """
        mkdir -p {MMSEQS_DB_DIR} {MMSEQS_TMP}
        mmseqs createdb {input.aa_fasta} {MMSEQS_DB_DIR}/refsDB \
            >> {log} 2>&1

        mmseqs cluster \
            {MMSEQS_DB_DIR}/refsDB \
            {MMSEQS_TMP}/clusterDB \
            {MMSEQS_TMP} \
            --min-seq-id {MMSEQS_MINID} \
            -c {MMSEQS_COV} \
            --threads {threads} \
            >> {log} 2>&1

        mmseqs createtsv \
            {MMSEQS_DB_DIR}/refsDB \
            {MMSEQS_DB_DIR}/refsDB \
            {MMSEQS_TMP}/clusterDB \
            {output.clusters} \
            >> {log} 2>&1
        rm -r {MMSEQS_TMP} {MMSEQS_DB_DIR}
        """


# =============================================================================
# PER-SAMPLE RULES
# Each rule changes into {sample}/ so that scripts using relative paths
# (analysis/, cleaned.tsv, etc.) work unchanged.
# =============================================================================

rule diamond_blastx:
    """
    Align sample contigs against the viral protein database with Diamond BLASTx.
    """
    input:
        contigs = SAMPLEDIR + "/{sample}/assembly.fasta",
        db      = DIAMOND_DB_DIR + "/viral_proteins.dmnd",
    output:
        hits =  SAMPLEDIR + "/{sample}/analysis/diamond.tsv",
    log:
        log_path("{sample}/diamond_blastx.log")
    threads: THREADS_DIAMOND
    shell:
        """
        cd {SAMPLEDIR}/{wildcards.sample}
        mkdir -p analysis
        diamond blastx \
            --query {input.contigs} \
            --db {input.db} \
            --out {output.hits} \
            --outfmt 6 qseqid sseqid stitle pident length evalue bitscore \
                       qstart qend sstart send qlen slen \
            --evalue {DIAMOND_EVALUE} \
            --max-target-seqs {DIAMOND_MAXHITS} \
            --threads {threads} \
            > {log} 2>&1
        """


rule qc_cluster:
    """
    Run QC_cluster.py inside the sample directory.
    Reads analysis/clusters.tsv (symlinked from the global one).
    Writes analysis/seg_miss.txt and analysis/protein_miss.txt.
    """
    input:
        clusters = CLUSTERS_TSV,
        # diamond output must exist so the analysis dir is present
        diamond  = SAMPLEDIR + "/{sample}/analysis/diamond.tsv",
    output:
        seg_miss     = SAMPLEDIR +"/{sample}/analysis/seg_miss.txt",
        protein_miss = SAMPLEDIR +"/{sample}/analysis/protein_miss.txt",
    log:
        log_path("{sample}/qc_cluster.log")
    shell:
        """
        # Make a local symlink so QC_cluster.py can open analysis/clusters.tsv

        cd {SAMPLEDIR}/{wildcards.sample}
        ln -sf $(realpath {input.clusters}) analysis/clusters.tsv
        python {SCRIPT_DIR}/QC_cluster.py > {log} 2>&1
        """


rule full_viral_annotation:
    """
    Run the full viral annotation pipeline on seg_miss.txt.
    Writes cleaned.tsv, outliers.tsv, qc_metrics.tsv, and QC PNGs.
    """
    input:
        seg_miss = SAMPLEDIR + "/{sample}/analysis/seg_miss.txt",
    output:
        cleaned  = SAMPLEDIR + "/{sample}/analysis/cleaned.tsv",
        outliers = SAMPLEDIR + "/{sample}/analysis/outliers.tsv",
        qc_tsv   = SAMPLEDIR + "/{sample}/analysis/qc_metrics.tsv",
    log:
        log_path("{sample}/full_viral_annotation.log")
    shell:
        """
        cd {SAMPLEDIR}/{wildcards.sample}
        python {SCRIPT_DIR}/full_viral_annotation_pipeline.py \
            analysis/seg_miss.txt {ONTOLOGY} \
            > {log} 2>&1
        """


rule build_cluster_metadata:
    """
    Build cluster_metadata.tsv from protein metadata + cleaned.tsv + clusters.tsv.
    """
    input:
        protein_meta = REF_DIR + "/protein_metadata.tsv",
        cleaned      = SAMPLEDIR + "/{sample}/analysis/cleaned.tsv",
        clusters     = CLUSTERS_TSV,
    output:
        cluster_meta = SAMPLEDIR + "/{sample}/analysis/cluster_metadata.tsv",
    log:
        log_path("{sample}/build_cluster_metadata.log")
    shell:
        """
        # Ensure symlinks are in place (clusters already linked in qc_cluster)
        # protein_metadata.tsv is read as ref/protein_metadata.tsv from sample dir
        cd {SAMPLEDIR}/{wildcards.sample}
        mkdir -p ref
        ln -sf $(realpath {input.protein_meta}) \
            ref/protein_metadata.tsv 2>/dev/null || true

        python {SCRIPT_DIR}/build_cluster_metadata.py \
            > {log} 2>&1
        """


rule score_segments:
    """
    Score and rank contig-to-segment assignments.
    Writes ranked_candidates.tsv inside the sample directory.
    """
    input:
        diamond      = SAMPLEDIR + "/{sample}/analysis/diamond.tsv",
        cluster_meta = SAMPLEDIR + "/{sample}/analysis/cluster_metadata.tsv",
    output:
        ranked = SAMPLEDIR + "/{sample}/ranked_candidates.tsv",
    log:
        log_path("{sample}/score_segments.log")
    shell:
        """
        cd {SAMPLEDIR}/{wildcards.sample}
        python {SCRIPT_DIR}/score_segments.py \
            > {log} 2>&1
        """
