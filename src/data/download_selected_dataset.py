from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path
from typing import Any
from urllib.request import urlopen

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from data.search_ncbi_geo import entrez_json


SOURCE_REGISTER_COLUMNS = [
    "source_name",
    "source_type",
    "accession",
    "url",
    "local_path",
    "date_accessed",
    "notes",
]


def gse_prefix(accession: str) -> str:
    if not accession.startswith("GSE") or not accession[3:].isdigit():
        raise ValueError(f"Expected a GEO Series accession like GSE305638, got {accession!r}")
    digits = accession[3:]
    return f"GSE{digits[:-3]}nnn"


def build_download_plan(accession: str) -> list[dict[str, str]]:
    """Build the Step 4 download plan for the selected GEO accession."""
    if accession != "GSE305638":
        raise ValueError("Step 4 is locked to selected dataset GSE305638.")

    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_prefix(accession)}/{accession}"
    return [
        {
            "source_name": f"{accession} series matrix",
            "source_type": "geo_series_matrix",
            "accession": accession,
            "url": f"{base}/matrix/{accession}_series_matrix.txt.gz",
            "local_path": f"data/raw/{accession}_series_matrix.txt.gz",
            "notes": "GEO series matrix containing sample metadata and expression matrix annotations when available.",
        },
        {
            "source_name": f"{accession} processed gene counts",
            "source_type": "geo_processed_counts",
            "accession": accession,
            "url": f"{base}/suppl/{accession}_Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv.gz",
            "local_path": f"data/raw/{accession}_Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv.gz",
            "notes": "Processed supplementary RNA-seq gene-count table from GEO; not raw FASTQ.",
        },
        {
            "source_name": f"{accession} MINiML metadata",
            "source_type": "geo_miniml_metadata",
            "accession": accession,
            "url": f"{base}/miniml/{accession}_family.xml.tgz",
            "local_path": f"data/interim/{accession}_family.xml.tgz",
            "notes": "GEO MINiML archive for sample and platform metadata.",
        },
        {
            "source_name": f"{accession} NCBI GDS summary",
            "source_type": "ncbi_gds_summary_json",
            "accession": accession,
            "url": (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                f"?db=gds&id=200305638&retmode=json"
            ),
            "local_path": f"data/interim/{accession}_gds_summary.json",
            "notes": "NCBI Entrez GDS summary metadata, including linked publication fields when available.",
        },
    ]


def download_binary(url: str, output_path: Path, timeout: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=timeout) as response:
        output_path.write_bytes(response.read())


def download_gds_summary(output_path: Path, timeout: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = entrez_json("esummary", {"db": "gds", "id": "200305638"}, timeout)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def source_register_rows(plan: list[dict[str, str]], date_accessed: str) -> list[dict[str, str]]:
    rows = []
    for item in plan:
        rows.append(
            {
                "source_name": item["source_name"],
                "source_type": item["source_type"],
                "accession": item["accession"],
                "url": item["url"],
                "local_path": item["local_path"],
                "date_accessed": date_accessed,
                "notes": item["notes"],
            }
        )
    return rows


def write_source_register(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_REGISTER_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in SOURCE_REGISTER_COLUMNS})


def execute_download(plan: list[dict[str, str]], timeout: int) -> None:
    for item in plan:
        output_path = Path(item["local_path"])
        if item["source_type"] == "ncbi_gds_summary_json":
            download_gds_summary(output_path, timeout)
        else:
            download_binary(item["url"], output_path, timeout)
        print(f"downloaded {item['url']} -> {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download lightweight files for GSE305638.")
    parser.add_argument("--accession", default="GSE305638")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--source-register", type=Path, default=Path("references/source_register.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plan = build_download_plan(args.accession)
    execute_download(plan, args.timeout)
    rows = source_register_rows(plan, date.today().isoformat())
    write_source_register(rows, args.source_register)
    print(f"wrote source register -> {args.source_register}")


if __name__ == "__main__":
    main()
