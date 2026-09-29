#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path


# ----------------------------------------------------------------------
# Defaults
# ----------------------------------------------------------------------

DEFAULT_WORD_SIZE = 7
DEFAULT_EVALUE = 100
DEFAULT_MAX_TARGET_SEQS = 100
TOP_HITS = 2


# ----------------------------------------------------------------------
# Header parsing
# ----------------------------------------------------------------------

SEGMENT_RE = re.compile(
    r"\bsegment\s+(10|[1-9])\b",
    re.IGNORECASE,
)

GENE_RE = re.compile(
    r"\b("
    r"VP1|VP2|VP3|VP4|VP5|VP6|VP7|"
    r"NS1|NS2|NS3|NS3A"
    r")\s+gene\b",
    re.IGNORECASE,
)


def get_segment(title):
    match = SEGMENT_RE.search(title)

    if match:
        return int(match.group(1))

    return ""


def get_gene(title):
    match = GENE_RE.search(title)

    if match:
        return match.group(1).upper()

    return ""


# ----------------------------------------------------------------------
# Read FASTA
# ----------------------------------------------------------------------

def read_fasta(path):
    """
    Yield:

        (sequence_id, header, sequence)

    for each FASTA record.
    """

    header = None
    sequence = []

    with open(path) as f:

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if header is not None:

                    sequence_id = header.split()[0]

                    yield (
                        sequence_id,
                        header,
                        "".join(sequence),
                    )

                header = line[1:].strip()
                sequence = []

            else:

                sequence.append(line.strip())

    if header is not None:

        sequence_id = header.split()[0]

        yield (
            sequence_id,
            header,
            "".join(sequence),
        )


# ----------------------------------------------------------------------
# Write temporary single-query FASTA
# ----------------------------------------------------------------------

def write_query_fasta(query_id, header, sequence, directory):

    path = Path(directory) / f"{query_id}.fasta"

    with path.open("w") as f:

        f.write(f">{header}\n")

        # Wrap sequence at 80 characters
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i + 80] + "\n")

    return path


# ----------------------------------------------------------------------
# Run BLAST
# ----------------------------------------------------------------------

def run_blast(
    query_fasta,
    database,
    word_size,
    evalue,
    max_target_seqs,
    verbose=False,
):
    """
    Run BLASTN against one query sequence.

    Returns a list of parsed hits.
    """

    outfmt = (
        "6 "
        "qseqid "
        "sseqid "
        "pident "
        "length "
        "qstart "
        "qend "
        "sstart "
        "send "
        "evalue "
        "bitscore "
        "qcovs "
        "stitle"
    )

    command = [
        "blastn",
        "-task", "blastn",
        "-query", str(query_fasta),
        "-db", str(database),
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-reward", "2",
        "-penalty", "-3",
        "-gapopen", "5",
        "-gapextend", "2",
        "-dust", "yes",
        "-outfmt", outfmt,
    ]

    if verbose:

        print(
            "\nBLAST command:",
            " ".join(command),
            file=sys.stderr,
        )

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:

        raise RuntimeError(
            "BLASTN failed:\n\n"
            + result.stderr
        )

    if verbose and result.stderr.strip():

        print(
            "BLAST stderr:",
            result.stderr,
            file=sys.stderr,
        )

    hits = []

    for line in result.stdout.splitlines():

        if not line.strip():
            continue

        fields = line.split("\t")

        if len(fields) != 12:

            if verbose:

                print(
                    f"WARNING: unexpected BLAST output "
                    f"({len(fields)} fields): {line}",
                    file=sys.stderr,
                )

            continue

        (
            qseqid,
            sseqid,
            pident,
            length,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
            qcovs,
            stitle,
        ) = fields

        hits.append(
            {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": float(pident),
                "length": int(length),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "qcovs": float(qcovs),
                "stitle": stitle,
            }
        )

    # Highest bitscore first
    hits.sort(
        key=lambda h: (
            -h["bitscore"],
            h["evalue"],
        )
    )

    return hits


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Search BTV nucleotide databases hierarchically: "
            "complete segments first, then genes/CDS."
        )
    )

    parser.add_argument(
        "query",
        help="Query FASTA",
    )

    parser.add_argument(
        "--complete-db",
        default="/mnt/lustre/RDS-live/tennakoon/BTV/BTV_complete_segments",
        help="Complete-segment BLAST database",
    )

    parser.add_argument(
        "--gene-db",
        default="/mnt/lustre/RDS-live/tennakoon/BTV/BTV_genes",
        help="Gene/CDS BLAST database",
    )

    parser.add_argument(
        "--word-size",
        type=int,
        default=DEFAULT_WORD_SIZE,
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=DEFAULT_EVALUE,
    )

    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=DEFAULT_MAX_TARGET_SEQS,
    )

    parser.add_argument(
        "--output",
        default="BTV_blast_results.tsv",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
    )

    args = parser.parse_args()

    query = Path(args.query)

    if not query.exists():

        sys.exit(
            f"ERROR: query file does not exist: {query}"
        )

    queries = list(read_fasta(query))

    if not queries:

        sys.exit(
            "ERROR: no FASTA sequences found."
        )

    print(
        f"Found {len(queries)} query sequence(s).",
        file=sys.stderr,
    )


    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    fields = [
        "query",
        "search_database",
        "match_type",
        "segment",
        "gene",
        "accession",
        "identity",
        "alignment_length",
        "query_coverage",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "title",
    ]


    with tempfile.TemporaryDirectory() as tmpdir:

        with open(args.output, "w", newline="") as output:

            writer = csv.DictWriter(
                output,
                fieldnames=fields,
                delimiter="\t",
            )

            writer.writeheader()


            # ----------------------------------------------------------
            # Process each query independently
            # ----------------------------------------------------------

            for query_id, header, sequence in queries:

                print(
                    f"\nSearching {query_id} "
                    f"against complete BTV segments...",
                    file=sys.stderr,
                )


                # ------------------------------------------------------
                # Make a one-sequence FASTA
                # ------------------------------------------------------

                query_fasta = write_query_fasta(
                    query_id,
                    header,
                    sequence,
                    tmpdir,
                )


                # ------------------------------------------------------
                # Search complete segments
                # ------------------------------------------------------

                hits = run_blast(
                    query_fasta,
                    args.complete_db,
                    args.word_size,
                    args.evalue,
                    args.max_target_seqs,
                    args.verbose,
                )


                print(
                    f"  BLAST returned {len(hits)} "
                    f"complete-segment hit(s).",
                    file=sys.stderr,
                )


                # ------------------------------------------------------
                # Complete segment hits
                # ------------------------------------------------------

                if hits:

                    top_hits = hits[:TOP_HITS]

                    for hit in top_hits:

                        writer.writerow(
                            {
                                "query": query_id,
                                "search_database":
                                    "complete_segments",
                                "match_type":
                                    "complete_segment",
                                "segment":
                                    get_segment(hit["stitle"]),
                                "gene":
                                    "",
                                "accession":
                                    hit["sseqid"],
                                "identity":
                                    hit["pident"],
                                "alignment_length":
                                    hit["length"],
                                "query_coverage":
                                    hit["qcovs"],
                                "qstart":
                                    hit["qstart"],
                                "qend":
                                    hit["qend"],
                                "sstart":
                                    hit["sstart"],
                                "send":
                                    hit["send"],
                                "evalue":
                                    hit["evalue"],
                                "bitscore":
                                    hit["bitscore"],
                                "title":
                                    hit["stitle"],
                            }
                        )

                        print(
                            f"    segment "
                            f"{get_segment(hit['stitle'])}: "
                            f"{hit['sseqid']} "
                            f"identity={hit['pident']:.2f}% "
                            f"length={hit['length']} "
                            f"coverage={hit['qcovs']:.1f}% "
                            f"bitscore={hit['bitscore']:.1f}",
                            file=sys.stderr,
                        )

                    # --------------------------------------------------
                    # IMPORTANT:
                    #
                    # We found complete-segment hits, so do not search
                    # genes/CDS.
                    # --------------------------------------------------

                    continue


                # ------------------------------------------------------
                # No complete segment hit
                # ------------------------------------------------------

                print(
                    "  No complete-segment hit; "
                    "searching genes/CDS...",
                    file=sys.stderr,
                )


                gene_hits = run_blast(
                    query_fasta,
                    args.gene_db,
                    args.word_size,
                    args.evalue,
                    args.max_target_seqs,
                    args.verbose,
                )


                print(
                    f"  BLAST returned {len(gene_hits)} "
                    f"gene/CDS hit(s).",
                    file=sys.stderr,
                )


                if gene_hits:

                    top_hits = gene_hits[:TOP_HITS]

                    for hit in top_hits:

                        writer.writerow(
                            {
                                "query": query_id,
                                "search_database": "genes",
                                "match_type": "gene",
                                "segment": "",
                                "gene":
                                    get_gene(hit["stitle"]),
                                "accession":
                                    hit["sseqid"],
                                "identity":
                                    hit["pident"],
                                "alignment_length":
                                    hit["length"],
                                "query_coverage":
                                    hit["qcovs"],
                                "qstart":
                                    hit["qstart"],
                                "qend":
                                    hit["qend"],
                                "sstart":
                                    hit["sstart"],
                                "send":
                                    hit["send"],
                                "evalue":
                                    hit["evalue"],
                                "bitscore":
                                    hit["bitscore"],
                                "title":
                                    hit["stitle"],
                            }
                        )

                        print(
                            f"    gene "
                            f"{get_gene(hit['stitle'])}: "
                            f"{hit['sseqid']} "
                            f"identity={hit['pident']:.2f}% "
                            f"length={hit['length']} "
                            f"coverage={hit['qcovs']:.1f}% "
                            f"bitscore={hit['bitscore']:.1f}",
                            file=sys.stderr,
                        )

                else:

                    print(
                        "  No BTV nucleotide hit found.",
                        file=sys.stderr,
                    )


    print(
        f"\nResults written to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
