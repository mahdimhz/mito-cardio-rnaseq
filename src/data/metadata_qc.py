from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path


METADATA_COLUMNS = [
    "sample_id",
    "sample_title",
    "expression_column",
    "group",
    "disease_label",
    "mutation",
    "paired_patient_id",
    "cell_line",
    "genotype",
    "cell_type",
    "organism",
    "replicate",
    "source_name",
    "instrument_model",
    "library_name",
    "geo_accession",
    "notes",
]

SUMMARY_COLUMNS = [
    "group",
    "sample_count",
    "patient_count",
    "cell_lines",
    "description",
]


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value.replace('""', '"').strip()


def parse_key_value(value: str) -> tuple[str, str] | None:
    if ":" not in value:
        return None
    key, raw = value.split(":", 1)
    key = key.strip().lower().replace(" ", "_")
    return key, raw.strip()


def parse_expression_column(descriptions: list[str]) -> str:
    for description in descriptions:
        match = re.search(r"Column name .*?:\s*([A-Za-z0-9_.-]+)", description)
        if match:
            return match.group(1)
    return ""


def parse_library_name(descriptions: list[str]) -> str:
    for description in descriptions:
        match = re.search(r"Library name:\s*(.+)", description)
        if match:
            return match.group(1).strip()
    return ""


def parse_series_matrix(path: Path) -> list[dict[str, str]]:
    """Parse sample-level metadata from a GEO series matrix file."""
    sample_rows: list[tuple[str, list[str]]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].replace("!Sample_", "").lower()
            values = [strip_quotes(part) for part in parts[1:]]
            sample_rows.append((key, values))

    if not sample_rows:
        return []

    sample_count = max(len(values) for _, values in sample_rows)
    records: list[dict[str, str | list[str]]] = [{"descriptions": []} for _ in range(sample_count)]
    for key, values in sample_rows:
        for index, value in enumerate(values):
            if key == "characteristics_ch1":
                parsed = parse_key_value(value)
                if parsed:
                    records[index][parsed[0]] = parsed[1]
            elif key == "description":
                records[index].setdefault("descriptions", []).append(value)
            else:
                records[index][f"sample_{key}"] = value

    parsed_records: list[dict[str, str]] = []
    for record in records:
        descriptions = [str(item) for item in record.get("descriptions", [])]
        parsed_records.append(
            {
                "sample_geo_accession": str(record.get("sample_geo_accession", "")),
                "sample_title": str(record.get("sample_title", "")),
                "source_name": str(record.get("sample_source_name_ch1", "")),
                "organism": str(record.get("sample_organism_ch1", "")),
                "patient": str(record.get("patient", "")),
                "cell_line": str(record.get("cell_line", "")),
                "cell_type": str(record.get("cell_type", "")),
                "genotype": str(record.get("genotype", "")),
                "instrument_model": str(record.get("sample_instrument_model", "")),
                "expression_column": parse_expression_column(descriptions),
                "library_name": parse_library_name(descriptions),
            }
        )
    return parsed_records


def infer_group(sample_title: str, genotype: str, cell_line: str) -> str:
    text = f"{sample_title} {genotype} {cell_line}".lower()
    if "healthy" in text or "corrected" in text or cell_line.lower().endswith("ic"):
        return "control"
    if "disease" in text or genotype.lower() == "wt":
        return "disease"
    return "unknown"


def infer_replicate(sample_title: str) -> str:
    parts = [part.strip() for part in sample_title.split(",")]
    return parts[-1] if parts and parts[-1].isdigit() else ""


def clean_sample_metadata(records: list[dict[str, str]]) -> list[dict[str, str]]:
    clean = []
    for record in records:
        group = infer_group(record["sample_title"], record["genotype"], record["cell_line"])
        row = {
            "sample_id": record["sample_geo_accession"],
            "sample_title": record["sample_title"],
            "expression_column": record["expression_column"],
            "group": group,
            "disease_label": "Friedreich ataxia" if group == "disease" else "isogenic corrected control",
            "mutation": "FXN GAA repeat expansion" if group == "disease" else "CRISPR-corrected FXN",
            "paired_patient_id": record["patient"],
            "cell_line": record["cell_line"],
            "genotype": record["genotype"],
            "cell_type": record["cell_type"],
            "organism": record["organism"],
            "replicate": infer_replicate(record["sample_title"]),
            "source_name": record["source_name"],
            "instrument_model": record["instrument_model"],
            "library_name": record["library_name"],
            "geo_accession": "GSE305638",
            "notes": "Disease samples are FRDA patient iPSC-cardiomyocytes; controls are paired isogenic CRISPR-corrected lines.",
        }
        clean.append(row)
    return clean


def summarize_samples(clean: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = []
    groups = sorted({row["group"] for row in clean})
    for group in groups:
        group_rows = [row for row in clean if row["group"] == group]
        patient_count = len({row["paired_patient_id"] for row in group_rows if row["paired_patient_id"]})
        cell_lines = "; ".join(sorted({row["cell_line"] for row in group_rows if row["cell_line"]}))
        description = (
            "FRDA patient iPSC-cardiomyocytes"
            if group == "disease"
            else "Isogenic corrected iPSC-cardiomyocyte controls"
        )
        rows.append(
            {
                "group": group,
                "sample_count": str(len(group_rows)),
                "patient_count": str(patient_count),
                "cell_lines": cell_lines,
                "description": description,
            }
        )
    return rows


def write_csv(rows: list[dict[str, str]], path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_metadata_outputs(clean: list[dict[str, str]], metadata_path: Path, summary_path: Path) -> None:
    write_csv(clean, metadata_path, METADATA_COLUMNS)
    write_csv(summarize_samples(clean), summary_path, SUMMARY_COLUMNS)


def make_notebook(notebook_path: Path) -> None:
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Step 5 - Download and Metadata QC\n",
                "\n",
                "Clean sample metadata for `GSE305638`, the selected FRDA iPSC-cardiomyocyte RNA-seq dataset.\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "from pathlib import Path\n",
                "import sys\n",
                "\n",
                "PROJECT_ROOT = Path.cwd().parent if Path.cwd().name == 'notebooks' else Path.cwd()\n",
                "sys.path.insert(0, str(PROJECT_ROOT / 'src'))\n",
                "\n",
                "from data.metadata_qc import parse_series_matrix, clean_sample_metadata, write_metadata_outputs\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "series_matrix = PROJECT_ROOT / 'data/raw/GSE305638_series_matrix.txt.gz'\n",
                "metadata_path = PROJECT_ROOT / 'data/processed/metadata_clean.csv'\n",
                "summary_path = PROJECT_ROOT / 'reports/tables/sample_summary.csv'\n",
                "\n",
                "records = parse_series_matrix(series_matrix)\n",
                "metadata = clean_sample_metadata(records)\n",
                "write_metadata_outputs(metadata, metadata_path, summary_path)\n",
                "len(metadata)\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "import pandas as pd\n",
                "\n",
                "metadata_df = pd.read_csv(metadata_path)\n",
                "summary_df = pd.read_csv(summary_path)\n",
                "summary_df\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "Metadata notes: disease samples are FRDA patient iPSC-cardiomyocytes with FXN GAA repeat expansion; controls are paired isogenic CRISPR-corrected iPSC-cardiomyocytes. This is MERRF-relevant mitochondrial cardiomyopathy biology, not a MERRF-specific dataset.\n",
            ],
        },
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    notebook_path.parent.mkdir(parents=True, exist_ok=True)
    notebook_path.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean GSE305638 sample metadata.")
    parser.add_argument("--series-matrix", type=Path, default=Path("data/raw/GSE305638_series_matrix.txt.gz"))
    parser.add_argument("--metadata-output", type=Path, default=Path("data/processed/metadata_clean.csv"))
    parser.add_argument("--summary-output", type=Path, default=Path("reports/tables/sample_summary.csv"))
    parser.add_argument("--notebook-output", type=Path, default=Path("notebooks/02_download_and_metadata_qc.ipynb"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records = parse_series_matrix(args.series_matrix)
    clean = clean_sample_metadata(records)
    write_metadata_outputs(clean, args.metadata_output, args.summary_output)
    make_notebook(args.notebook_output)
    counts = Counter(row["group"] for row in clean)
    print(f"wrote {len(clean)} samples to {args.metadata_output}")
    print(f"wrote sample summary to {args.summary_output}")
    print(f"wrote notebook to {args.notebook_output}")
    print("group_counts=" + ", ".join(f"{group}:{count}" for group, count in sorted(counts.items())))


if __name__ == "__main__":
    main()
