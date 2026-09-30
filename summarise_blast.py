#!/usr/bin/env python3

import argparse
from pathlib import Path
import html

import pandas as pd


HIT_FIELDS = [
    "top_accession",
    "top_identity",
    "top_alignment_length",
    "top_query_coverage",
    "top_bitscore",
    "second_accession",
    "second_identity",
    "second_alignment_length",
    "second_query_coverage",
    "second_bitscore",
]


def read_blast_table(path, sample, prefix):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing BLAST table: {path}")

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    required = {"contig"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns: {sorted(missing)}")

    # A contig should identify one result row within a sample. Avoid an
    # accidental many-to-many merge if an input contains duplicate contigs.
    #duplicated = df["contig"].duplicated(keep=False)
    #if duplicated.any():
    #    examples = df.loc[duplicated, "contig"].unique()[:5].tolist()
    #    raise ValueError(
    #        f"{path}: duplicate contig rows detected, examples: {examples}"
    #    )

    # Keep the first row for each contig and discard subsequent duplicates.
    duplicate_count = df["contig"].duplicated(keep="first").sum()
    if duplicate_count:
        print(
            f"Warning: {path}: discarded {duplicate_count} duplicate "
            "contig row(s); keeping the first occurrence."
        )

    df = df.drop_duplicates(subset="contig", keep="first").copy()

    keep = ["contig"]
    if "genes" in df.columns:
        keep.append("genes")
    keep.extend(field for field in HIT_FIELDS if field in df.columns)
    df = df[keep].copy()

    rename = {
        field: f"{prefix}_{field}"
        for field in df.columns
        if field not in ("contig", "genes")
    }
    df = df.rename(columns=rename)

    if "genes" in df.columns:
        df = df.rename(columns={"genes": f"{prefix}_genes"})

    return df


def combine_sample(sample, sample_dir):
    sample_path = Path(sample_dir) / sample

    protein = read_blast_table(
        sample_path / "selected_hits.tsv",
        sample,
        "protein",
    )
    nucleotide = read_blast_table(
        sample_path / "ranked_candidates.tsv.blasted",
        sample,
        "nt",
    )

    combined = protein.merge(
        nucleotide,
        on="contig",
        how="outer",
        validate="one_to_one",
    )

    protein_genes = combined.get("protein_genes", pd.Series("", index=combined.index))
    nt_genes = combined.get("nt_genes", pd.Series("", index=combined.index))
    combined.insert(
        0,
        "genes",
        protein_genes.where(protein_genes.ne(""), nt_genes),
    )
    combined.insert(0, "sample", sample)
    # Remove rows where gene assignment is missing or empty.
    combined["genes"] = combined["genes"].fillna("").astype(str).str.strip()
    combined = combined[combined["genes"].ne("")].copy()

    return combined


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-html", required=True)
    args = parser.parse_args()

    tables = [
        combine_sample(sample, args.sample_dir)
        for sample in args.samples
    ]
    result = pd.concat(tables, ignore_index=True, sort=False)

    # Add per-sample coverage paths from the manifests produced by
    # align_reads_coverage.py. Paths are relative to the report directory.
    coverage_tables = []
    for sample in args.samples:
        manifest_path = Path(args.sample_dir) / sample / "coverage_manifest.tsv"
        if not manifest_path.is_file():
            raise FileNotFoundError(
                f"Missing coverage manifest for {sample}: {manifest_path}"
            )
        coverage = pd.read_csv(
            manifest_path, sep="\\t", dtype=str
        ).fillna("")
        required_manifest_columns = {
            "sample", "contig", "coverage_plot", "status"
        }
        missing_manifest_columns = required_manifest_columns - set(coverage.columns)
        if missing_manifest_columns:
            raise ValueError(
                f"{manifest_path}: missing columns "
                f"{sorted(missing_manifest_columns)}"
            )
        coverage_tables.append(coverage[
            ["sample", "contig", "coverage_plot", "status"]
        ])

    coverage_result = pd.concat(coverage_tables, ignore_index=True)
    if coverage_result.duplicated(["sample", "contig"]).any():
        duplicates = coverage_result.loc[
            coverage_result.duplicated(["sample", "contig"], keep=False),
            ["sample", "contig"]
        ].drop_duplicates().head(5).to_dict("records")
        raise ValueError(
            f"Duplicate sample/contig entries in coverage manifests: {duplicates}"
        )

    result = result.merge(
        coverage_result,
        on=["sample", "contig"],
        how="left",
        validate="one_to_one",
    )
    result["coverage_plot"] = result["coverage_plot"].fillna("")
    result["status"] = result["status"].fillna("no_coverage_manifest_entry")

    result = result.sort_values(
        ["genes", "sample", "contig"],
        na_position="last",
    )

    output_tsv = Path(args.output_tsv)
    output_html = Path(args.output_html)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    result.to_csv(output_tsv, sep="\t", index=False)

    # Keep the TSV machine-readable; create links only in the HTML copy.
    report = result.copy()
###########################
    column_groups = {
        "General": [
            "sample",
            "genes",
            "contig",
            "coverage_plot",
        ],
        "Protein — top hit": [
            "protein_top_accession",
            "protein_top_identity",
            "protein_top_alignment_length",
            "protein_top_query_coverage",
            "protein_top_bitscore",
        ],
        "Protein — second hit": [
            "protein_second_accession",
            "protein_second_identity",
            "protein_second_alignment_length",
            "protein_second_query_coverage",
            "protein_second_bitscore",
        ],
        "Nucleotide — top hit": [
            "nt_top_accession",
            "nt_top_identity",
            "nt_top_alignment_length",
            "nt_top_query_coverage",
            "nt_top_bitscore",
        ],
        "Nucleotide — second hit": [
            "nt_second_accession",
            "nt_second_identity",
            "nt_second_alignment_length",
            "nt_second_query_coverage",
            "nt_second_bitscore",
        ],
    }

    # Retain only columns present in this report.
    column_groups = {
        group: [column for column in columns if column in report.columns]
        for group, columns in column_groups.items()
    }
    column_groups = {
        group: columns for group, columns in column_groups.items() if columns
    }

    # Include any report columns not explicitly assigned above.
    assigned_columns = {
        column for columns in column_groups.values() for column in columns
    }
    other_columns = [
        column for column in report.columns if column not in assigned_columns
    ]
    if other_columns:
        column_groups["Other fields"] = other_columns

    column_indexes = {
        column: index for index, column in enumerate(report.columns)
    }

    group_controls = []

    for group, columns in column_groups.items():
        group_id = f"group-{len(group_controls)}"
        indexes = [column_indexes[column] for column in columns]

        individual_controls = "\n".join(
            f"""
            <label class="column-option">
                <input type="checkbox"
                       class="column-checkbox"
                       data-column="{column_indexes[column]}"
                       checked>
                {html.escape(column)}
            </label>
            """
            for column in columns
        )

        group_controls.append(
            f"""
            <details class="column-group" open>
                <summary>
                    <input type="checkbox"
                           class="group-checkbox"
                           data-group="{group_id}"
                           checked
                           aria-label="Toggle all columns in {html.escape(group)}">
                    <strong>{html.escape(group)}</strong>
                </summary>
                <div class="group-columns"
                     data-group-columns="{group_id}">
                    {individual_controls}
                </div>
            </details>
            """
        )

    column_controls = "\n".join(group_controls)

    # Escape all cell content, then render only the trusted, generated
    # relative coverage path as an anchor. This avoids interpreting BLAST
    # titles/accessions as HTML.
    for column in report.columns:
        if column != "coverage_plot":
            report[column] = report[column].map(
                lambda value: html.escape(str(value))
            )

    def coverage_link(path):
        if not path or pd.isna(path):
            return "Not available"
        # Manifest paths are generated by the pipeline, but reject absolute
        # paths and parent traversal so links remain within the report bundle.
        relative = Path(str(path))
        if relative.is_absolute() or ".." in relative.parts:
            return "Invalid coverage path"
        href = html.escape(relative.as_posix(), quote=True)
        return f'<a href="{href}" target="_blank">View PNG</a>'

    report["coverage_plot"] = result["coverage_plot"].map(coverage_link)

    table_html = report.to_html(
        index=False,
        escape=False,
        classes="blast-table",
        table_id="blastTable",
    )

    page = f"""<!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Cross-sample BLAST summary</title>
    <style>
    body {{
        font-family: sans-serif;
        margin: 2rem;
    }}
    .note {{
        color: #555;
        margin-bottom: 1rem;
    }}
    input[type="search"] {{
        padding: 0.5rem;
        width: min(30rem, 90%);
        margin: 0.5rem 0 1rem;
    }}
    .actions {{
        display: flex;
        flex-wrap: wrap;
        gap: 0.5rem;
        margin: 0.5rem 0;
    }}
    button {{
        padding: 0.4rem 0.7rem;
        cursor: pointer;
    }}
    .column-groups {{
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
        gap: 0.7rem;
        margin: 0.7rem 0 1.2rem;
    }}
    .column-group {{
        border: 1px solid #ddd;
        border-radius: 6px;
        padding: 0.6rem 0.75rem;
    }}
    .column-group summary {{
        cursor: pointer;
        display: flex;
        align-items: center;
        gap: 0.5rem;
    }}
    .group-columns {{
        display: flex;
        flex-direction: column;
        gap: 0.35rem;
        padding: 0.7rem 0.2rem 0.2rem 1.5rem;
    }}
    .column-option {{
        font-size: 0.88rem;
        overflow-wrap: anywhere;
    }}
    .table-wrap {{
        overflow: auto;
        max-height: 75vh;
    }}
    table {{
        border-collapse: collapse;
        font-size: 0.85rem;
    }}
    th, td {{
        border: 1px solid #ddd;
        padding: 0.4rem 0.55rem;
        text-align: left;
        white-space: nowrap;
    }}
    th {{
        position: sticky;
        top: 0;
        background: #f2f2f2;
        z-index: 1;
    }}
    tr:nth-child(even) {{
        background: #fafafa;
    }}
    </style>
    </head>
    <body>

    <h1>Cross-sample BLAST summary</h1>
    <p class="note">
    One row per sample and contig. Protein and nucleotide results are joined
    by contig within each sample. Expand a column group to select individual
    fields, or use the group checkbox to toggle all fields in that group.
    </p>

    <label for="filter">Filter rows:</label><br>
    <input id="filter" type="search"
           placeholder="Search gene, sample, contig, accession...">

    <h3>Visible columns</h3>
    <div class="actions">
        <button type="button" id="showAll">Show all</button>
        <button type="button" id="hideAll">Hide all</button>
    </div>

    <div class="column-groups">
        {column_controls}
    </div>

    <div class="table-wrap">
        {table_html}
    </div>

    <script>
    const table = document.getElementById("blastTable");
    const filter = document.getElementById("filter");
    const columnCheckboxes = [
        ...document.querySelectorAll(".column-checkbox")
    ];
    const groupCheckboxes = [
        ...document.querySelectorAll(".group-checkbox")
    ];

    function setColumnVisibility(index, visible) {{
        for (const row of table.rows) {{
            const cell = row.cells[index];
            if (cell) {{
                cell.style.display = visible ? "" : "none";
            }}
        }}
    }}

    function updateGroupState(groupId) {{
        const group = document.querySelector(
            `.group-checkbox[data-group="${{groupId}}"]`
        );
        const columns = document.querySelector(
            `.group-columns[data-group-columns="${{groupId}}"]`
        );
        const checkboxes = [
            ...columns.querySelectorAll(".column-checkbox")
        ];

        const checkedCount = checkboxes.filter(box => box.checked).length;
        group.checked = checkedCount === checkboxes.length;
        group.indeterminate =
            checkedCount > 0 && checkedCount < checkboxes.length;
    }}

    columnCheckboxes.forEach(checkbox => {{
        checkbox.addEventListener("change", () => {{
            setColumnVisibility(
                Number(checkbox.dataset.column),
                checkbox.checked
            );

            const groupId = checkbox.closest(".group-columns")
                .dataset.groupColumns;
            updateGroupState(groupId);
        }});
    }});

    groupCheckboxes.forEach(groupCheckbox => {{
        groupCheckbox.addEventListener("change", () => {{
            const groupId = groupCheckbox.dataset.group;
            const columns = document.querySelector(
                `.group-columns[data-group-columns="${{groupId}}"]`
            );
            const checkboxes = columns.querySelectorAll(".column-checkbox");

            checkboxes.forEach(checkbox => {{
                checkbox.checked = groupCheckbox.checked;
                setColumnVisibility(
                    Number(checkbox.dataset.column),
                    checkbox.checked
                );
            }});

            groupCheckbox.indeterminate = false;
        }});
    }});

    document.getElementById("showAll").addEventListener("click", () => {{
        columnCheckboxes.forEach(checkbox => {{
            checkbox.checked = true;
            setColumnVisibility(Number(checkbox.dataset.column), true);
        }});

        groupCheckboxes.forEach(checkbox => {{
            checkbox.checked = true;
            checkbox.indeterminate = false;
        }});
    }});

    document.getElementById("hideAll").addEventListener("click", () => {{
        columnCheckboxes.forEach(checkbox => {{
            checkbox.checked = false;
            setColumnVisibility(Number(checkbox.dataset.column), false);
        }});

        groupCheckboxes.forEach(checkbox => {{
            checkbox.checked = false;
            checkbox.indeterminate = false;
        }});
    }});

    filter.addEventListener("input", () => {{
        const query = filter.value.toLowerCase();

        for (const row of table.tBodies[0].rows) {{
            row.style.display =
                row.textContent.toLowerCase().includes(query) ? "" : "none";
        }}
    }});
    </script>

    </body>
    </html>
    """

    output_html.write_text(page, encoding="utf-8")


if __name__ == "__main__":
    main()
