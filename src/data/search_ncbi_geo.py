from __future__ import annotations

import argparse
import csv
import json
import os
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import urlopen


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

SEARCH_QUERIES = [
    "MERRF iPSC RNA-seq",
    "MERRF induced pluripotent stem cell transcriptome",
    "MERRF m.8344A>G RNA-seq",
    "MERRF A8344G expression",
    "MT-TK mitochondrial disease iPSC",
    "mitochondrial disease iPSC RNA-seq",
    "mitochondrial disorder iPSC transcriptome",
    "MELAS iPSC RNA-seq",
    "mtDNA mutation iPSC RNA-seq",
    "mitochondrial cardiomyopathy iPSC cardiomyocyte RNA-seq",
    "Barth syndrome hiPSC cardiomyocyte RNA-seq",
    "OXPHOS hiPSC cardiomyocyte RNA-seq",
    "mitochondrial dysfunction cardiomyocyte RNA-seq",
    "human cardiomyopathy cardiomyocyte RNA-seq",
    "dilated cardiomyopathy single nucleus RNA-seq human heart",
    "heart failure cardiomyocyte single nucleus RNA-seq",
]

REGISTRY_COLUMNS = [
    "accession",
    "source_database",
    "title",
    "organism",
    "disease",
    "mutation",
    "cell_type",
    "assay_type",
    "platform",
    "sample_count",
    "control_count",
    "case_count",
    "processed_matrix_available",
    "raw_fastq_available",
    "metadata_available",
    "publication",
    "publication_year",
    "download_url",
    "keep_status",
    "priority_score",
    "reason_for_keep_or_exclude",
    "notes",
]


def entrez_json(endpoint: str, params: dict[str, Any], timeout: int) -> dict[str, Any]:
    query = dict(params)
    query.setdefault("retmode", "json")
    url = f"{EUTILS_BASE}/{endpoint}.fcgi?{urlencode(query)}"
    with urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def search_gds(
    query: str,
    retmax: int,
    timeout: int,
    email: str | None = None,
    api_key: str | None = None,
) -> list[str]:
    params: dict[str, Any] = {
        "db": "gds",
        "term": f"GSE[ETYP] AND ({query})",
        "retmax": retmax,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    payload = entrez_json("esearch", params, timeout)
    return payload.get("esearchresult", {}).get("idlist", [])


def fetch_gds_summaries(
    uids: list[str],
    timeout: int,
    email: str | None = None,
    api_key: str | None = None,
) -> list[dict[str, Any]]:
    if not uids:
        return []

    params: dict[str, Any] = {"db": "gds", "id": ",".join(uids)}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    payload = entrez_json("esummary", params, timeout)
    result = payload.get("result", {})
    return [result[uid] for uid in result.get("uids", []) if uid in result]


def has_processed_supplement(suppfile: str) -> str:
    """Return a conservative draft flag for processed table availability."""
    text = (suppfile or "").upper()
    if not text:
        return "unknown"
    if text in {"NONE", "NA", "N/A"}:
        return "no"
    processed_extensions = ("CSV", "TSV", "TXT", "XLS", "XLSX", "SOFT", "SERIES_MATRIX")
    if any(extension in text for extension in processed_extensions):
        return "yes"
    return "unknown"


def normalize_platform(gpl: Any) -> str:
    value = str(gpl or "").strip()
    if not value:
        return ""
    if value.upper().startswith("GPL"):
        return value.upper()
    if value.isdigit():
        return f"GPL{value}"
    return value


def publication_from_summary(record: dict[str, Any]) -> str:
    pubmedids = record.get("pubmedids") or []
    if isinstance(pubmedids, str):
        pubmedids = [pubmedids]
    return "; ".join(f"PMID:{pmid}" for pmid in pubmedids if str(pmid).strip())


def year_from_public_date(record: dict[str, Any]) -> str:
    pdat = str(record.get("pdat") or "")
    return pdat[:4] if len(pdat) >= 4 and pdat[:4].isdigit() else ""


def row_from_summary(record: dict[str, Any], matched_queries: list[str]) -> dict[str, str]:
    accession = str(record.get("accession") or "").strip()
    samples = record.get("samples") or []
    ftplink = str(record.get("ftplink") or "").strip()
    download_url = ftplink or f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"

    notes = [
        f"Entrez UID: {record.get('uid', '')}",
        f"matched_queries: {' | '.join(matched_queries)}",
    ]
    if record.get("suppfile"):
        notes.append(f"suppfile: {record.get('suppfile')}")
    if record.get("bioproject"):
        notes.append(f"bioproject: {record.get('bioproject')}")

    row = {column: "" for column in REGISTRY_COLUMNS}
    row.update(
        {
            "accession": accession,
            "source_database": "NCBI GEO/GDS",
            "title": str(record.get("title") or "").strip(),
            "organism": str(record.get("taxon") or "").strip(),
            "assay_type": str(record.get("gdstype") or "").strip(),
            "platform": normalize_platform(record.get("gpl")),
            "sample_count": str(record.get("n_samples") or ""),
            "processed_matrix_available": has_processed_supplement(str(record.get("suppfile") or "")),
            "raw_fastq_available": "unknown",
            "metadata_available": "yes" if samples else "unknown",
            "publication": publication_from_summary(record),
            "publication_year": year_from_public_date(record),
            "download_url": download_url,
            "keep_status": "candidate_unscreened",
            "reason_for_keep_or_exclude": "Found by NCBI GEO/GDS Entrez search; requires Step 3 screening.",
            "notes": "; ".join(part for part in notes if part),
        }
    )
    return row


def collect_candidate_rows(
    queries: list[str],
    retmax_per_query: int,
    timeout: int,
    sleep_seconds: float,
    email: str | None = None,
    api_key: str | None = None,
) -> list[dict[str, str]]:
    uid_to_queries: dict[str, list[str]] = {}
    for query in queries:
        for uid in search_gds(query, retmax_per_query, timeout, email=email, api_key=api_key):
            uid_to_queries.setdefault(uid, []).append(query)
        time.sleep(sleep_seconds)

    rows: list[dict[str, str]] = []
    uids = sorted(uid_to_queries)
    for start in range(0, len(uids), 100):
        batch = uids[start : start + 100]
        summaries = fetch_gds_summaries(batch, timeout, email=email, api_key=api_key)
        for summary in summaries:
            uid = str(summary.get("uid") or "")
            rows.append(row_from_summary(summary, uid_to_queries.get(uid, [])))
        time.sleep(sleep_seconds)

    return sorted(rows, key=lambda row: (row["accession"], row["title"]))


def write_registry(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGISTRY_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in REGISTRY_COLUMNS})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Search NCBI GEO/GDS with mitochondrial-disease and cardiomyocyte queries and draft dataset_registry.csv."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("references/dataset_registry.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--retmax-per-query",
        type=int,
        default=20,
        help="Maximum GEO/GDS records to retrieve per search query.",
    )
    parser.add_argument("--timeout", type=int, default=30, help="HTTP timeout in seconds.")
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.34,
        help="Delay between Entrez requests to respect NCBI rate limits.",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL"),
        help="Optional NCBI contact email. Defaults to NCBI_EMAIL.",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="Optional NCBI API key. Defaults to NCBI_API_KEY.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = collect_candidate_rows(
        SEARCH_QUERIES,
        retmax_per_query=args.retmax_per_query,
        timeout=args.timeout,
        sleep_seconds=args.sleep_seconds,
        email=args.email,
        api_key=args.api_key,
    )
    write_registry(rows, args.output)
    print(f"Wrote {len(rows)} candidate records to {args.output}")


if __name__ == "__main__":
    main()
