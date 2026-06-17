from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def metric_lookup(path: Path) -> dict[str, str]:
    df = pd.read_csv(path)
    return dict(zip(df["metric"].astype(str), df["value"].astype(str), strict=False))


def load_summary_inputs(project_root: Path) -> dict[str, Any]:
    registry = pd.read_csv(project_root / "references/dataset_registry.csv")
    selected = registry.loc[registry["accession"] == "GSE305638"].iloc[0].fillna("")
    metadata = pd.read_csv(project_root / "data/processed/metadata_clean.csv")
    qc = metric_lookup(project_root / "reports/tables/expression_qc_summary.csv")
    de = metric_lookup(project_root / "reports/tables/differential_expression_summary.csv")
    mito_results = pd.read_csv(project_root / "reports/tables/mitochondrial_gene_results.csv")
    mito_summary = pd.read_csv(project_root / "reports/tables/mitochondrial_pathway_summary.csv")
    top_mito = mito_results[mito_results["result_status"] == "tested"].sort_values("p_value").head(5)

    return {
        "selected_dataset": selected.to_dict(),
        "metadata_counts": {
            "disease": int((metadata["group"] == "disease").sum()),
            "control": int((metadata["group"] == "control").sum()),
            "patients": int(metadata["paired_patient_id"].nunique()),
        },
        "qc": qc,
        "de": de,
        "mito": {
            "gene_count": int(mito_results["gene_symbol"].nunique()),
            "tested": int((mito_results["result_status"] == "tested").sum()),
            "min_padj": float(mito_results["padj"].min()),
            "top_direction": "mtDNA-encoded OXPHOS genes trend lower in disease",
            "top_genes": ", ".join(top_mito["gene_symbol"].tolist()),
            "pathways": mito_summary,
        },
    }


def build_summary_sections(inputs: dict[str, Any]) -> list[dict[str, str]]:
    selected = inputs["selected_dataset"]
    counts = inputs["metadata_counts"]
    qc = inputs["qc"]
    de = inputs["de"]
    mito = inputs["mito"]
    top_mito_genes = mito.get("top_genes", "no FDR-significant mitochondrial genes")

    sections = [
        {
            "section": "Dataset choice",
            "body": (
                f"Selected dataset: {selected['accession']} - {selected['title']}. "
                f"This is a {selected['disease']} model using {selected['cell_type']}. "
                "It combines human iPSC-cardiomyocytes, frataxin-deficiency disease biology, cardiac dysfunction, "
                "processed RNA-seq counts, and paired disease/control metadata."
            ),
        },
        {
            "section": "Methods",
            "body": (
                f"Metadata curation produced {counts['disease']} disease and {counts['control']} control samples "
                f"across {counts['patients']} patient/isogenic-control pairs. Expression QC aligned "
                f"{qc.get('sample_count', '')} samples and {qc.get('gene_count', '')} gene rows with "
                f"{qc.get('missing_value_count', '')} missing count values. Exploratory paired expression analysis used log2 CPM "
                "normalization, replicate averaging within patient/group, paired disease-control tests across "
                "three patient pairs, and Benjamini-Hochberg FDR correction."
            ),
        },
        {
            "section": "Main expression findings",
            "body": (
                f"The exploratory comparison tested {de.get('genes_tested', '')} expressed genes. "
                f"{de.get('nominal_p_lt_0.05', '')} genes had nominal p < 0.05, but "
                f"{de.get('fdr_lt_0.10', '')} genes passed FDR < 0.10 and {de.get('fdr_lt_0.05', '')} "
                "passed FDR < 0.05. This supports reporting ranked exploratory signals rather than claiming "
                "genome-wide findings after FDR correction."
            ),
        },
        {
            "section": "Mitochondrial gene findings",
            "body": (
                f"The mitochondrial reference contains {mito['gene_count']} genes, with {mito['tested']} tested "
                f"in the exploratory expression table. The minimum mitochondrial adjusted p-value was "
                f"{mito['min_padj']:.3f}, so no mitochondrial gene is FDR-significant. Descriptively, "
                f"{mito['top_direction']}; the strongest nominal mitochondrial genes include {top_mito_genes}."
            ),
        },
        {
            "section": "Limitations",
            "body": (
                "The exploratory paired expression analysis is limited by the effective paired sample size of "
                "three patient pairs and the Python implementation is not a full DESeq2 negative-binomial model. "
                "Mitochondrial pathway summaries are descriptive selected-gene summaries, not formal enrichment tests."
            ),
        },
        {
            "section": "Project scope",
            "body": (
                "The repository focuses on public dataset screening, biomedical metadata curation, RNA-seq count QC, "
                "paired expression analysis, mitochondrial gene-set interpretation, and reproducible reporting for "
                "a Friedreich ataxia iPSC-cardiomyocyte model."
            ),
        },
    ]
    return sections


def build_artifact_inventory(project_root: Path) -> pd.DataFrame:
    rows = []
    for artifact_type, folder in [("table", "reports/tables"), ("figure", "reports/figures")]:
        for path in sorted((project_root / folder).glob("*")):
            if path.is_file() and path.name != ".gitkeep":
                rows.append(
                    {
                        "artifact_type": artifact_type,
                        "relative_path": path.relative_to(project_root).as_posix(),
                        "bytes": path.stat().st_size,
                    }
                )
    return pd.DataFrame(rows, columns=["artifact_type", "relative_path", "bytes"])


def write_summary_tables(
    sections: list[dict[str, str]],
    inventory: pd.DataFrame,
    table_dir: Path,
) -> dict[str, str]:
    table_dir.mkdir(parents=True, exist_ok=True)
    summary_path = table_dir / "project_summary_table.csv"
    inventory_path = table_dir / "final_artifact_inventory.csv"
    pd.DataFrame(sections).to_csv(summary_path, index=False)
    inventory.to_csv(inventory_path, index=False)
    return {
        "project_summary_table": str(summary_path),
        "artifact_inventory": str(inventory_path),
    }


def make_notebook(notebook_path: Path, sections: list[dict[str, str]]) -> None:
    markdown = ["# Step 9 - Final Results Summary\n\n"]
    for section in sections:
        markdown.append(f"## {section['section']}\n\n{section['body']}\n\n")
    cells = [
        {"cell_type": "markdown", "metadata": {}, "source": markdown},
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "from pathlib import Path\n",
                "import pandas as pd\n",
                "\n",
                "PROJECT_ROOT = Path.cwd().parent if Path.cwd().name == 'notebooks' else Path.cwd()\n",
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/project_summary_table.csv')\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/final_artifact_inventory.csv')\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "README wording is kept separate from these generated summary tables.\n",
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


def run_results_summary(project_root: Path) -> dict[str, str]:
    inputs = load_summary_inputs(project_root)
    sections = build_summary_sections(inputs)
    inventory = build_artifact_inventory(project_root)
    outputs = write_summary_tables(sections, inventory, project_root / "reports/tables")
    notebook_path = project_root / "notebooks/06_results_summary.ipynb"
    make_notebook(notebook_path, sections)
    outputs["notebook"] = str(notebook_path)
    return outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Step 9 final summary artifacts.")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_results_summary(args.project_root.resolve())
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
