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
from glob import glob


configfile: "config.yaml"

SAMPLES        = config["samples"]
GB_DIR         = config["paths"]["gb_dir"]
REF_DIR        = config["paths"]["ref_dir"]
DIAMOND_DB_DIR = config["paths"]["diamond_db_dir"]
BLAST_DB_DIR   = config["paths"]["blast_db_dir"]
MMSEQS_DB_DIR  = config["paths"]["mmseqs_db_dir"]
MMSEQS_TMP     = config["paths"]["mmseqs_tmp_dir"]
CLUSTER_DIR   = config["paths"]["cluster_dir"]
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
CLUSTERS_TSV = os.path.join(CLUSTER_DIR, "clusters.tsv")
BTV_FASTA = os.path.expanduser("~/BTV/BTV.fasta")

def log_path(*parts):
    return os.path.join(LOGDIR, *parts)

import glob
import os

from pathlib import Path
def trimmed_reads1(wc):
    return [r1 for r1, r2 in trimmed_read_pair(wc)]

def trimmed_reads2(wc):
    return [r2 for r1, r2 in trimmed_read_pair(wc)]

def trimmed_read_pair(wc):
    sample = wc.sample
    trimmed_dir = Path(SAMPLEDIR) / sample / "trimmed"

    r1_files = sorted(trimmed_dir.glob("*_R1_*_val_1.fq"))
    if not r1_files:
        raise ValueError(f"No R1 files found for {sample} in {trimmed_dir}")

    pairs = []
    for r1 in r1_files:
        r2 = trimmed_dir / r1.name.replace("_R1_", "_R2_").replace(
            "_val_1.fq", "_val_2.fq"
        )
        if not r2.exists():
            raise ValueError(f"Paired read file not found for {sample}: {r2}")

        pairs.append((str(r1), str(r2)))

    return pairs
# =============================================================================
# Target rules
# =============================================================================

rule all:
    """Build the index and run the full per-sample pipeline for every sample."""
    input:
        # Index sentinels
        DIAMOND_DB_DIR + "/viral_proteins.dmnd",
        CLUSTER_DIR + "/cluster_metadata.tsv",

        # Per-sample final outputs
        expand(SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted", sample=SAMPLES),
        #expand(SAMPLEDIR + "/{sample}/selected_contigs/coverage_report.pdf", sample=SAMPLES),
        expand( SAMPLEDIR + "/{sample}/analysis/reference_coverage/coverage.pdf", sample=SAMPLES),
        expand( SAMPLEDIR + "/{sample}/analysis/reference_coverage/coverage.tsv", sample=SAMPLES),
        SAMPLEDIR + "/blast_summary.tsv",
        SAMPLEDIR + "/blast_summary.html"

rule build_index:
    """Convenience target: build all index files without running per-sample steps."""
    input:
        DIAMOND_DB_DIR + "/viral_proteins.dmnd",
        CLUSTER_DIR + "/cluster_metadata.tsv",
        REF_DIR + "/gene_stats.txt",


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
        genome_lengths = REF_DIR + "/genome_lengths.txt",
    log:
        log_path("parsegb.log")
    shell:
        """
        mkdir -p {REF_DIR}
        python {SCRIPT_DIR}/parsegb.py {input.gb_dir} {output.nt_fasta} {output.aa_fasta} {output.metadata} > {log} 2>&1
        awk '{{if(NR%2){{A=substr($1,2);sub(/\|.*/, "", A)}}else{{print A,length($0)}}}}' {output.nt_fasta} > {output.genome_lengths}
        """


rule diamond_makedb:
    """Build Diamond protein database from reference AA FASTA (one-time)."""
    input:
        aa_fasta = REF_DIR + "/references.aa.fasta",
        nt_fasta = REF_DIR + "/references.nt.fasta",
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
        mkdir -p {BLAST_DB_DIR}
        makeblastdb -in {input.nt_fasta} -dbtype nucl -out {BLAST_DB_DIR}/blastdb.ref
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
        mkdir -p {CLUSTER_DIR} {MMSEQS_DB_DIR} {MMSEQS_TMP}
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
        rm -rf {MMSEQS_TMP} {MMSEQS_DB_DIR}
        """


# =============================================================================
# PER-SAMPLE RULES
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
                       qstart qend sstart send qlen slen qframe \
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
    output:
        cluster_proteins = CLUSTER_DIR +"/clusters.txt",
    log:
        log_path("qc_cluster.log")
    shell:
        """
        mkdir -p {CLUSTER_DIR}
        python {SCRIPT_DIR}/QC_cluster.py {output.cluster_proteins} {input.clusters} > {log} 2>&1
        """


rule full_viral_annotation:
    """
    Run the full viral annotation pipeline on seg_miss.txt.
    Writes cleaned.tsv, outliers.tsv, qc_metrics.tsv, and QC PNGs.
    """
    input:
        cluster_proteins = CLUSTER_DIR +"/clusters.txt",
    output:
        cleaned  = CLUSTER_DIR + "/cleaned.tsv",
        outliers = CLUSTER_DIR + "/outliers.tsv",
        qc_tsv   = CLUSTER_DIR + "/qc_metrics.tsv",
    log:
        log_path("full_viral_annotation.log")
    shell:
        """
        cd {CLUSTER_DIR}
        python {SCRIPT_DIR}/full_viral_annotation_pipeline.py \
            {input.cluster_proteins} {ONTOLOGY} \
            > {log} 2>&1
        """


rule build_cluster_metadata:
    """
    Build cluster_metadata.tsv from protein metadata + cleaned.tsv + clusters.tsv.
    """
    input:
        protein_metadata = REF_DIR + "/protein_metadata.tsv",
        cleaned      = CLUSTER_DIR + "/cleaned.tsv",
        clusters     = CLUSTERS_TSV,
    output:
        cluster_meta = CLUSTER_DIR + "/cluster_metadata.tsv",
    log:
        log_path("build_cluster_metadata.log")
    shell:
        """
        python {SCRIPT_DIR}/build_cluster_metadata.py {input.clusters} {input.protein_metadata} {input.cleaned} {output.cluster_meta} \
            > {log} 2>&1
        """
rule build_length_stats:
    """
    Build std and length of general segments 
    """
    input:
        cluster_meta = CLUSTER_DIR + "/cluster_metadata.tsv",
        genome_lengths = REF_DIR + "/genome_lengths.txt",
        protein_metadata = REF_DIR + "/protein_metadata.tsv",
    output:
        gene_stats = REF_DIR + "/gene_stats.txt"
    shell:
        """
        Rscript generate_length_stats.R {input.cluster_meta} {input.genome_lengths} {input.protein_metadata} > {output.gene_stats} 
        """

rule score_segments:
    """
    Score and rank contig-to-segment assignments.

    BLAST validation is downstream because selected contig FASTAs do not exist
    until the extract_contigs checkpoint has completed.
    """
    input:
        diamond      = SAMPLEDIR + "/{sample}/analysis/diamond.tsv",
        cluster_meta = CLUSTER_DIR + "/cluster_metadata.tsv",
        gene_stats   = REF_DIR + "/gene_stats.txt"
    output:
        ranked = SAMPLEDIR + "/{sample}/ranked_candidates.tsv",
    log:
        log_path("{sample}/score_segments.log")
    shell:
        """
        cd {SAMPLEDIR}/{wildcards.sample}
        python {SCRIPT_DIR}/score_segments.py {input.cluster_meta} \
            > {log} 2>&1

        Rscript {SCRIPT_DIR}/filter_by_length.R \
            {output.ranked} {input.gene_stats} {output.ranked}.filtered
        """


checkpoint extract_contigs:
    input:
        ranked = SAMPLEDIR + "/{sample}/ranked_candidates.tsv.filtered",
    output:
        contigs_dir = directory(SAMPLEDIR + "/{sample}/analysis/contigs")
    shell:
        """
        mkdir -p {output.contigs_dir}
        awk '(NR>1) {{print $2}}' {input.ranked} | while read -r F; do
            grep -A1 "$F" {SAMPLEDIR}/{wildcards.sample}/assembly.fasta \
                > {output.contigs_dir}/"$F".fa
        done
        """

def get_selected_contigs_dir(wc):
    return checkpoints.extract_contigs.get(
        sample=wc.sample
    ).output.contigs_dir


rule blast_selected_contigs:
    """
    BLAST selected contigs and retain the top two hits
    with their details in separate columns.
    """
    input:
        contigs = get_selected_contigs_dir,
        ranked = SAMPLEDIR + "/{sample}/ranked_candidates.tsv.filtered"
    output:
        ranked_blast = SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted"
    log:
        log_path("{sample}/blast_selected_contigs.log")
    shell:
        r"""
        set -euo pipefail

        cd {SAMPLEDIR}/{wildcards.sample}

        python {SCRIPT_DIR}/append_blast_hits.py \
            {input.ranked} \
            {input.contigs} \
            {output.ranked_blast} \
            >> {log} 2>&1
        """


rule reference_coverage:
    """
    Map reads to the top BLAST reference for the first
    candidate of each gene and generate coverage plots.
    """
    input:
        ranked_blast = SAMPLEDIR +
            "/{sample}/ranked_candidates.tsv.blasted",
        reference = BTV_FASTA,
        R1 = lambda wc: sorted(glob.glob(
            f"{SAMPLEDIR}/{wc.sample}/trimmed/*_val_1.fq"
        )),
        R2 = lambda wc: sorted(glob.glob(
            f"{SAMPLEDIR}/{wc.sample}/trimmed/*_val_2.fq"
        )),
    output:
        coverage_pdf = (
            SAMPLEDIR +
            "/{sample}/analysis/reference_coverage/coverage.pdf"
        ),
        coverage_tsv = (
            SAMPLEDIR +
            "/{sample}/analysis/reference_coverage/coverage.tsv"
        )
    log:
        log_path("{sample}/reference_coverage.log")
    threads: 10
    params:
        script = SCRIPT_DIR + "/reference_coverage.py"
    shell:
        """
        set -euo pipefail

        python {params.script} \
            {input.ranked_blast} \
            {input.reference} \
            {output.coverage_pdf} \
            {output.coverage_tsv} \
            --r1 {input.R1} \
            --r2 {input.R2} \
            --threads {threads} \
            > {log} 2>&1
        """

rule retain_selected_tophits:
    """
    Select the final contigs from the BLAST-augmented candidate table.
    """
    input:
        ranked_blast = SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted",
    output:
        selected_hits = SAMPLEDIR + "/{sample}/selected_hits.tsv",
    log:
        log_path("{sample}/retain_tophits.log")
    shell:
        """
        python {SCRIPT_DIR}/retain_tophits.py \
            {input.ranked_blast} {output.selected_hits} \
            > {log} 2>&1
        """


rule align_contig:
    input:
        contig = SAMPLEDIR + "/{sample}/analysis/contigs/{contig}.fa",
        R1 = lambda wc: glob(f"{SAMPLEDIR}/{wc.sample}/trimmed/*_val_1.fq"),
        R2 = lambda wc: glob(f"{SAMPLEDIR}/{wc.sample}/trimmed/*_val_2.fq"),
    output:
        bam = SAMPLEDIR + "/{sample}/analysis/bam/{contig}.fa.bam",
        bai = SAMPLEDIR + "/{sample}/analysis/bam/{contig}.fa.bam.bai",
        depth = SAMPLEDIR + "/{sample}/analysis/{contig}.fa.depth.txt"
    threads: 10
    shell:
        """
        mkdir -p {SAMPLEDIR}/{wildcards.sample}/analysis/bam
        minimap2 -ax sr -t {threads} {input.contig} {input.R1} {input.R2} | samtools sort -o {output.bam}
        samtools index {output.bam}
        samtools depth -aa {output.bam} > {output.depth}
        """

def get_contig_depths(wc):
    contigs_dir = checkpoints.extract_contigs.get(sample=wc.sample).output.contigs_dir
    contigs = glob_wildcards(os.path.join(contigs_dir, "{contig}.fa")).contig
    return expand(
        SAMPLEDIR + "/{sample}/analysis/{contig}.fa.depth.txt",
        sample=wc.sample,
        contig=contigs
    )

rule align_to_contigs:
    input:
        selected = SAMPLEDIR + "/{sample}/selected_hits.tsv",
        depths = get_contig_depths,
        nt_fasta = REF_DIR + "/references.nt.fasta",
    output:
        report = SAMPLEDIR + "/{sample}/analysis/coverage_report.pdf",
        containment = directory(SAMPLEDIR + "/{sample}/analysis/containment")

    shell:
        """
        python {SCRIPT_DIR}/coverage.py {SAMPLEDIR}/{wildcards.sample}/analysis {input.selected} {output.report}
        viral-containment   --hits {input.selected} --contigs {SAMPLEDIR}/{wildcards.sample}/assembly.fasta --refs {input.nt_fasta} --outdir {output.containment}
        """

rule extract_selected_contigs:
    input:
        selected = SAMPLEDIR + "/{sample}/selected_hits.tsv",
        contigs = SAMPLEDIR + "/{sample}/assembly.fasta",
        coverage_report = SAMPLEDIR + "/{sample}/analysis/coverage_report.pdf",
        containment = directory(SAMPLEDIR + "/{sample}/analysis/containment"),
    output:
        selected_contigs = directory(SAMPLEDIR + "/{sample}/selected_contigs"),
        coverage_report = SAMPLEDIR + "/{sample}/selected_contigs/coverage_report.pdf"
    shell:
        """
        python {SCRIPT_DIR}/extract_gene_contigs.py {input.selected} {input.contigs} {output.selected_contigs}
        pdfunite {SAMPLEDIR}/{wildcards.sample}/analysis/containment/plots/*.pdf {output.selected_contigs}/containment_plots.pdf 
        cp {SAMPLEDIR}/{wildcards.sample}/analysis/coverage_report.pdf {output.selected_contigs}
        """


#rule summarise_blast:
#    input:
#        protein=expand(
#            SAMPLEDIR + "/{sample}/selected_hits.tsv",
#            sample=SAMPLES
#        ),
#        nucleotide=expand(
#            SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted",
#            sample=SAMPLES
#        )
#    output:
#        tsv=SAMPLEDIR + "/blast_summary.tsv",
#        html=SAMPLEDIR + "/blast_summary.html"
#    params:
#        sample_dir=SAMPLEDIR,
#        samples=" ".join(SAMPLES),
#        script=SCRIPT_DIR + "/summarise_blast.py",
#    log:
#        log_path("summarise_blast.log")
#    shell:
#        """
#        python {params.script:q} \
#            --sample-dir {params.sample_dir:q} \
#            --samples {params.samples} \
#            --output-tsv {output.tsv:q} \
#            --output-html {output.html:q} \
#            > {log:q} 2>&1
#        """
rule align_reads_coverage:
    input:
        hits=SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted",
        reads1=trimmed_reads1,
        reads2=trimmed_reads2
    output:
        manifest=SAMPLEDIR + "/{sample}/coverage_manifest.tsv",
        plots=directory(SAMPLEDIR + "/coverage_plots/{sample}")
    params:
        reference=BTV_FASTA,
        script=SCRIPT_DIR + "/align_reads_coverage.py"
    log:
        log_path("align_reads_coverage", "{sample}.log")
    threads: 4
    #conda:
    #    "envs/coverage.yaml"
    shell:
        r"""
        python {params.script:q} \
            --sample {wildcards.sample:q} \
            --hits {input.hits:q} \
            --reference-fasta {params.reference:q} \
            --reads1 {input.reads1:q} \
            --reads2 {input.reads2:q} \
            --plot-dir {output.plots:q} \
            --manifest {output.manifest:q} \
            > {log:q} 2>&1
        """

rule summarise_blast:
    input:
        diamond=expand(
            SAMPLEDIR + "/{sample}/analysis/diamond.tsv",
            sample=SAMPLES
        ),
        nucleotide=expand(
            SAMPLEDIR + "/{sample}/ranked_candidates.tsv.blasted",
            sample=SAMPLES
        ),
        coverage=expand(
            SAMPLEDIR + "/{sample}/coverage_manifest.tsv",
            sample=SAMPLES
        ),
        metadata=CLUSTER_DIR + "/cluster_metadata.tsv"
    output:
        tsv=SAMPLEDIR + "/blast_summary.tsv",
        html=SAMPLEDIR + "/blast_summary.html"
    params:
        sample_dir=SAMPLEDIR,
        samples=" ".join(SAMPLES),
        script=SCRIPT_DIR + "/summarise_blast.py"
    log:
        log_path("summarise_blast.log")
    shell:
        """
        python {params.script:q} \
            --sample-dir {params.sample_dir:q} \
            --samples {params.samples} \
            --metadata {input.metadata:q} \
            --output-tsv {output.tsv:q} \
            --output-html {output.html:q} \
            > {log:q} 2>&1
        """
