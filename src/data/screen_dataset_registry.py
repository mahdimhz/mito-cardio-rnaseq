from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from data.search_ncbi_geo import REGISTRY_COLUMNS, entrez_json, fetch_gds_summaries, write_registry


MERRF_TERMS = ("merrf", "m.8344a>g", "a8344g", "mt-tk", "mttk")
IPSC_TERMS = ("ipsc", "hipsc", "induced pluripotent")
CARDIAC_TERMS = (
    "cardiomyocyte",
    "cardiomyocytes",
    "cardiomyopathy",
    "heart",
    "cardiac",
    "myocard",
    "ventricle",
    "dilated cardiomyopathy",
)
EXPRESSION_TERMS = ("expression profiling", "rna-seq", "transcriptome", "single cell", "single nucleus")


def text_blob(row: dict[str, str], summary: dict[str, Any] | None = None) -> str:
    summary = summary or {}
    sample_titles = " ".join(str(sample.get("title", "")) for sample in summary.get("samples", []) or [])
    parts = [
        row.get("title", ""),
        row.get("organism", ""),
        row.get("assay_type", ""),
        row.get("notes", ""),
        str(summary.get("summary", "")),
        sample_titles,
    ]
    return " ".join(parts).lower()


def contains_any(text: str, terms: tuple[str, ...]) -> bool:
    return any(term in text for term in terms)


def infer_disease(text: str) -> str:
    if "merrf" in text or "m.8344a>g" in text or "a8344g" in text or "mt-tk" in text:
        return "MERRF syndrome"
    if "melas" in text or "m.3243a>g" in text or "a3243g" in text:
        return "MELAS / m.3243A>G mitochondrial disease"
    if "barth" in text:
        return "Barth syndrome"
    if "friedreich" in text or "frataxin" in text or "gaa repeat" in text:
        return "Friedreich ataxia / frataxin-related mitochondrial disease"
    if "mitochondrial cardiomyopathy" in text:
        return "mitochondrial cardiomyopathy"
    if "heart failure" in text:
        return "heart failure"
    if "dilated cardiomyopathy" in text:
        return "dilated cardiomyopathy"
    if "cardiomyopathy" in text:
        return "cardiomyopathy"
    if "mitochondrial disease" in text or "mitochondrial disorder" in text:
        return "mitochondrial disease"
    if "mitochondrial dysfunction" in text or "oxphos" in text:
        return "mitochondrial dysfunction / OXPHOS"
    return ""


def infer_mutation(text: str) -> str:
    mutations = []
    mutation_patterns = [
        ("m.8344A>G", ("m.8344a>g", "8344a>g")),
        ("A8344G", ("a8344g",)),
        ("MT-TK", ("mt-tk", "mttk")),
        ("m.3243A>G", ("m.3243a>g", "3243a>g")),
        ("A3243G", ("a3243g",)),
        ("MT-TL1", ("mt-tl1", "mttl1")),
        ("TAZ", (" taz ", "taz mutation")),
        ("FXN / frataxin", (" frataxin ", " fxn ", "gaa repeat")),
    ]
    padded = f" {text} "
    for label, terms in mutation_patterns:
        if any(term in padded for term in terms):
            mutations.append(label)
    return "; ".join(dict.fromkeys(mutations))


def infer_cell_type(text: str, organism: str) -> str:
    human = "homo sapiens" in organism.lower() or "human" in text
    prefix = "human " if human else ""
    has_ipsc = contains_any(text, IPSC_TERMS)
    has_cardiomyocyte = "cardiomyocyte" in text or "cardiomyocytes" in text

    if has_ipsc and has_cardiomyocyte:
        return f"{prefix}iPSC-derived cardiomyocytes".strip()
    if has_ipsc and ("neuron" in text or "neuronal" in text or "cortical" in text):
        return f"{prefix}iPSC-derived neurons / neural cells".strip()
    if has_ipsc and "organoid" in text:
        return f"{prefix}iPSC-derived organoid model".strip()
    if has_ipsc:
        return f"{prefix}iPSC-derived cells".strip()
    if has_cardiomyocyte:
        return f"{prefix}cardiomyocytes".strip()
    if "single nucleus" in text and ("heart" in text or "cardiac" in text):
        return f"{prefix}heart tissue, single-nucleus".strip()
    if "single cell" in text and ("heart" in text or "cardiac" in text):
        return f"{prefix}heart tissue, single-cell".strip()
    if "heart" in text or "cardiac" in text:
        return f"{prefix}heart tissue".strip()
    if "fibroblast" in text:
        return f"{prefix}fibroblasts".strip()
    if "neuron" in text or "neuronal" in text:
        return f"{prefix}neuronal cells".strip()
    return ""


def count_case_control(summary: dict[str, Any] | None) -> tuple[str, str, bool]:
    """Estimate case/control counts from GSM titles when labels are obvious."""
    summary = summary or {}
    control_terms = ("control", "ctrl", "healthy", "wild type", "wildtype", "vehicle", "untreated")
    case_terms = (
        "patient",
        "merrf",
        "melas",
        "barth",
        "mutant",
        "mutation",
        "disease",
        "treated",
        "failing",
        "heart failure",
        "dcm",
        "cardiomyopathy",
    )
    control_count = 0
    case_count = 0
    for sample in summary.get("samples", []) or []:
        title = str(sample.get("title", "")).lower()
        is_control = any(term in title for term in control_terms)
        is_case = any(term in title for term in case_terms)
        if is_control and not is_case:
            control_count += 1
        elif is_case and not is_control:
            case_count += 1

    if control_count and case_count:
        return str(control_count), str(case_count), True
    return "", "", False


def infer_keep_status(score: int, row: dict[str, str]) -> str:
    if score >= 7:
        return "keep_high_priority"
    if score >= 4:
        return "keep_candidate"
    if row.get("processed_matrix_available") == "yes" and row.get("assay_type"):
        return "maybe_secondary"
    return "exclude_low_priority"


def score_candidate(row: dict[str, str]) -> int:
    text = text_blob(row)
    score = 0
    if contains_any(text, MERRF_TERMS):
        score += 3
    if contains_any(text, IPSC_TERMS) and "homo sapiens" in row.get("organism", "").lower():
        score += 2
    if contains_any(text, CARDIAC_TERMS):
        score += 2
    if contains_any(text, EXPRESSION_TERMS):
        score += 1
    if row.get("processed_matrix_available") == "yes":
        score += 1
    if row.get("control_count") and row.get("case_count") and row.get("metadata_available") == "yes":
        score += 1
    if row.get("raw_fastq_available") == "yes" and row.get("processed_matrix_available") != "yes":
        score -= 2
    if not row.get("control_count") or not row.get("case_count"):
        score -= 2
    organism = row.get("organism", "").lower()
    if "mus musculus" in organism and "homo sapiens" not in organism:
        score -= 2
    if not contains_any(text, EXPRESSION_TERMS):
        score -= 3
    return score


def reason_for_row(row: dict[str, str]) -> str:
    reasons = []
    text = text_blob(row)
    if contains_any(text, MERRF_TERMS):
        reasons.append("contains MERRF/m.8344A>G/MT-TK terms")
    if contains_any(text, IPSC_TERMS):
        reasons.append("uses iPSC-related model")
    if contains_any(text, CARDIAC_TERMS):
        reasons.append("cardiac/cardiomyopathy relevance")
    if row.get("processed_matrix_available") == "yes":
        reasons.append("processed supplementary table likely available")
    if row.get("control_count") and row.get("case_count"):
        reasons.append("case/control labels detected in sample titles")
    if not reasons:
        reasons.append("limited relevance based on automated screening")
    return "; ".join(reasons)


def infer_screening_fields(
    row: dict[str, str],
    summary: dict[str, Any] | None,
    raw_fastq_available: str,
) -> dict[str, str]:
    screened = dict(row)
    summary = summary or {}
    text = text_blob(screened, summary)
    organism = screened.get("organism") or str(summary.get("taxon") or "")

    control_count, case_count, counted = count_case_control(summary)
    if control_count and case_count:
        screened["control_count"] = control_count
        screened["case_count"] = case_count

    screened["organism"] = organism
    screened["disease"] = infer_disease(text)
    screened["mutation"] = infer_mutation(text)
    screened["cell_type"] = infer_cell_type(text, organism)
    screened["raw_fastq_available"] = raw_fastq_available
    screened["metadata_available"] = "yes" if summary.get("samples") else screened.get("metadata_available", "unknown")
    if summary.get("gdstype"):
        screened["assay_type"] = str(summary.get("gdstype"))
    if summary.get("n_samples"):
        screened["sample_count"] = str(summary.get("n_samples"))
    if summary.get("pubmedids") and not screened.get("publication"):
        pubmedids = summary.get("pubmedids") or []
        screened["publication"] = "; ".join(f"PMID:{pmid}" for pmid in pubmedids)

    score = score_candidate(screened)
    screened["priority_score"] = str(score)
    screened["keep_status"] = infer_keep_status(score, screened)
    screened["reason_for_keep_or_exclude"] = reason_for_row(screened)

    notes = [note.strip() for note in screened.get("notes", "").split(";") if note.strip()]
    if counted:
        notes.append("case/control counts are heuristic from GEO sample titles")
    if summary.get("summary"):
        notes.append("screened_with_gds_summary")
    screened["notes"] = "; ".join(dict.fromkeys(note for note in notes if note))
    return {column: screened.get(column, "") for column in REGISTRY_COLUMNS}


def read_registry(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def uid_from_notes(notes: str) -> str:
    match = re.search(r"Entrez UID:\s*(\d+)", notes or "")
    return match.group(1) if match else ""


def sra_available_for_accession(accession: str, timeout: int) -> str:
    payload = entrez_json(
        "esearch",
        {"db": "sra", "term": f"{accession}[All Fields]", "retmax": 0},
        timeout,
    )
    count = int(payload.get("esearchresult", {}).get("count", "0"))
    return "yes" if count > 0 else "no"


def screen_registry(input_path: Path, timeout: int, sleep_seconds: float) -> list[dict[str, str]]:
    rows = read_registry(input_path)
    uid_to_row = {uid_from_notes(row.get("notes", "")): row for row in rows if uid_from_notes(row.get("notes", ""))}
    summaries = {}
    uids = list(uid_to_row)
    for start in range(0, len(uids), 100):
        for summary in fetch_gds_summaries(uids[start : start + 100], timeout):
            summaries[str(summary.get("uid", ""))] = summary
        time.sleep(sleep_seconds)

    screened_rows = []
    for row in rows:
        uid = uid_from_notes(row.get("notes", ""))
        summary = summaries.get(uid, {})
        raw_fastq_available = sra_available_for_accession(row.get("accession", ""), timeout)
        screened_rows.append(infer_screening_fields(row, summary, raw_fastq_available))
        time.sleep(sleep_seconds)

    return sorted(
        screened_rows,
        key=lambda item: (int(item.get("priority_score") or 0), item.get("accession", "")),
        reverse=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Screen candidate GEO datasets in dataset_registry.csv.")
    parser.add_argument("--input", type=Path, default=Path("references/dataset_registry.csv"))
    parser.add_argument("--output", type=Path, default=Path("references/dataset_registry.csv"))
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--sleep-seconds", type=float, default=0.34)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    screened_rows = screen_registry(args.input, args.timeout, args.sleep_seconds)
    write_registry(screened_rows, args.output)
    print(f"Wrote {len(screened_rows)} screened candidate records to {args.output}")


if __name__ == "__main__":
    main()
