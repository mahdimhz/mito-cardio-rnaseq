from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GENE_SETS = [
    ("MT-TK", "mitochondrial translation / mtDNA maintenance", "Mitochondrial tRNA lysine gene."),
    ("MT-ND1", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND2", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND3", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND4", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND4L", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND5", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-ND6", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex I subunit."),
    ("MT-CO1", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex IV subunit."),
    ("MT-CO2", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex IV subunit."),
    ("MT-CO3", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex IV subunit."),
    ("MT-CYB", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded complex III subunit."),
    ("MT-ATP6", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded ATP synthase subunit."),
    ("MT-ATP8", "mtDNA-encoded OXPHOS genes", "mtDNA-encoded ATP synthase subunit."),
    ("NDUFA1", "respiratory chain complex I", "Nuclear-encoded complex I subunit."),
    ("NDUFB8", "respiratory chain complex I", "Nuclear-encoded complex I subunit."),
    ("SDHA", "respiratory chain complex II", "Nuclear-encoded complex II subunit."),
    ("UQCRC2", "respiratory chain complex III", "Nuclear-encoded complex III subunit."),
    ("COX4I1", "respiratory chain complex IV", "Nuclear-encoded complex IV subunit."),
    ("ATP5F1A", "ATP synthase", "Nuclear-encoded ATP synthase subunit."),
    ("TFAM", "mitochondrial translation / mtDNA maintenance", "mtDNA packaging and transcription factor."),
    ("POLG", "mitochondrial translation / mtDNA maintenance", "mtDNA polymerase catalytic subunit."),
    ("POLG2", "mitochondrial translation / mtDNA maintenance", "mtDNA polymerase accessory subunit."),
    ("TWNK", "mitochondrial translation / mtDNA maintenance", "mtDNA helicase."),
    ("MFN1", "fusion/fission", "Mitochondrial fusion."),
    ("MFN2", "fusion/fission", "Mitochondrial fusion."),
    ("OPA1", "fusion/fission", "Mitochondrial inner-membrane fusion/cristae structure."),
    ("DNM1L", "fusion/fission", "Mitochondrial fission."),
    ("FIS1", "fusion/fission", "Mitochondrial fission adaptor."),
    ("MFF", "fusion/fission", "Mitochondrial fission factor."),
    ("PINK1", "mitophagy", "Mitophagy initiation."),
    ("PRKN", "mitophagy", "Parkin-mediated mitophagy."),
    ("BNIP3", "mitophagy", "Hypoxia-associated mitophagy."),
    ("BNIP3L", "mitophagy", "NIX-mediated mitophagy."),
    ("SOD2", "oxidative stress / ROS", "Mitochondrial superoxide dismutase."),
    ("GPX1", "oxidative stress / ROS", "Glutathione peroxidase."),
    ("CAT", "oxidative stress / ROS", "Catalase."),
    ("PRDX3", "oxidative stress / ROS", "Mitochondrial peroxiredoxin."),
    ("PPARGC1A", "mitochondrial biogenesis", "PGC-1 alpha biogenesis regulator."),
    ("NRF1", "mitochondrial biogenesis", "Nuclear respiratory factor 1."),
    ("ESRRA", "mitochondrial biogenesis", "ERR alpha mitochondrial transcriptional regulator."),
]


def build_default_gene_sets() -> pd.DataFrame:
    return pd.DataFrame(GENE_SETS, columns=["gene_symbol", "category", "notes"])


def write_gene_sets(path: Path) -> pd.DataFrame:
    gene_sets = build_default_gene_sets()
    path.parent.mkdir(parents=True, exist_ok=True)
    gene_sets.to_csv(path, index=False)
    return gene_sets


def filter_mitochondrial_de_results(de_results: pd.DataFrame, gene_sets: pd.DataFrame) -> pd.DataFrame:
    """Annotate requested mitochondrial genes with DE results when tested."""
    de_columns = [
        "gene_id",
        "gene_symbol",
        "mean_log2_cpm_disease",
        "mean_log2_cpm_control",
        "mean_log2_cpm",
        "log2_fold_change",
        "p_value",
        "padj",
        "n_pairs",
        "method",
        "comparison",
    ]
    de_subset = de_results.copy()
    for column in de_columns:
        if column not in de_subset.columns:
            de_subset[column] = np.nan
    de_subset = de_subset[de_columns].copy()
    merged = gene_sets.merge(de_subset, on="gene_symbol", how="left")
    merged["result_status"] = np.where(merged["gene_id"].notna(), "tested", "not_tested_in_de_filter")
    merged["direction"] = np.select(
        [merged["log2_fold_change"] > 0, merged["log2_fold_change"] < 0],
        ["up_in_disease", "down_in_disease"],
        default="not_tested_or_no_change",
    )
    merged["rank_by_padj"] = merged["padj"].rank(method="min", na_option="bottom").astype(int)
    return merged.sort_values(["result_status", "padj", "p_value", "gene_symbol"], na_position="last").reset_index(drop=True)


def summarize_pathways(mito_results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for category, group in mito_results.groupby("category", sort=True):
        tested = group[group["result_status"] == "tested"].copy()
        if tested.empty:
            rows.append(
                {
                    "category": category,
                    "gene_count": str(len(group)),
                    "tested_gene_count": "0",
                    "mean_log2_fold_change": "",
                    "median_log2_fold_change": "",
                    "min_p_value": "",
                    "min_padj": "",
                    "direction": "not_tested",
                    "notes": "No genes from this category passed the Step 7 expression filter.",
                }
            )
            continue
        mean_lfc = float(tested["log2_fold_change"].mean())
        median_lfc = float(tested["log2_fold_change"].median())
        if mean_lfc > 0.05:
            direction = "up_in_disease"
        elif mean_lfc < -0.05:
            direction = "down_in_disease"
        else:
            direction = "mixed_or_near_zero"
        rows.append(
            {
                "category": category,
                "gene_count": str(len(group)),
                "tested_gene_count": str(len(tested)),
                "mean_log2_fold_change": f"{mean_lfc:.6f}",
                "median_log2_fold_change": f"{median_lfc:.6f}",
                "min_p_value": f"{float(tested['p_value'].min()):.6g}",
                "min_padj": f"{float(tested['padj'].min()):.6g}",
                "direction": direction,
                "notes": "Pathway direction is descriptive; no enrichment test is implied.",
            }
        )
    return pd.DataFrame(rows)


def prepare_heatmap_matrix(expression_log2: pd.DataFrame, gene_sets: pd.DataFrame) -> pd.DataFrame:
    """Return row-z-scored log2 expression for mitochondrial genes."""
    sample_cols = [col for col in expression_log2.columns if col not in {"gene_id", "gene_symbol"}]
    requested = gene_sets["gene_symbol"].drop_duplicates().tolist()
    subset = expression_log2[expression_log2["gene_symbol"].isin(requested)].copy()
    if subset.empty:
        return pd.DataFrame()
    collapsed = subset.groupby("gene_symbol", sort=False)[sample_cols].mean()
    collapsed = collapsed.reindex([gene for gene in requested if gene in collapsed.index])
    means = collapsed.mean(axis=1)
    stds = collapsed.std(axis=1).replace(0, 1)
    return collapsed.sub(means, axis=0).div(stds, axis=0).fillna(0)


def plot_heatmap(heatmap: pd.DataFrame, metadata: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if heatmap.empty:
        return
    sample_order = metadata["sample_id"].tolist()
    heatmap = heatmap[[sample for sample in sample_order if sample in heatmap.columns]]
    groups = metadata.set_index("sample_id").loc[heatmap.columns, "group"]
    fig_height = max(8, len(heatmap) * 0.22)
    fig, ax = plt.subplots(figsize=(11, fig_height))
    im = ax.imshow(heatmap.values, aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(heatmap.index)))
    ax.set_yticklabels(heatmap.index, fontsize=7)
    ax.set_xticks(range(len(heatmap.columns)))
    x_labels = [f"{sample}\n{groups.loc[sample]}" for sample in heatmap.columns]
    ax.set_xticklabels(x_labels, rotation=90, fontsize=7)
    ax.set_title("Mitochondrial Gene log2 Expression Z-scores")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Gene")
    fig.colorbar(im, ax=ax, label="row z-score")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def make_notebook(notebook_path: Path) -> None:
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Step 8 - Mitochondrial Gene-Set Analysis\n",
                "\n",
                "Mitochondrial gene-level summaries and sample-level pathway scores from the exploratory paired expression results.\n",
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
                "from analysis.mitochondrial_gene_set_analysis import run_mitochondrial_gene_set_analysis\n",
                "from analysis.mitochondrial_pathway_scores import run_mitochondrial_pathway_scores\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "outputs = run_mitochondrial_gene_set_analysis(\n",
                "    de_results_path=PROJECT_ROOT / 'reports/tables/differential_expression_results.csv',\n",
                "    expression_log2_path=PROJECT_ROOT / 'data/processed/expression_log2.csv',\n",
                "    metadata_path=PROJECT_ROOT / 'data/processed/metadata_clean.csv',\n",
                "    gene_sets_path=PROJECT_ROOT / 'references/mitochondrial_gene_sets.csv',\n",
                "    table_dir=PROJECT_ROOT / 'reports/tables',\n",
                "    figure_dir=PROJECT_ROOT / 'reports/figures',\n",
                "    notebook_path=PROJECT_ROOT / 'notebooks/05_mitochondrial_gene_set_analysis.ipynb',\n",
                ")\n",
                "outputs\n",
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
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/mitochondrial_gene_results.csv').head(20)\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/mitochondrial_pathway_summary.csv')\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "score_outputs = run_mitochondrial_pathway_scores(\n",
                "    expression_log2_path=PROJECT_ROOT / 'data/processed/expression_log2.csv',\n",
                "    metadata_path=PROJECT_ROOT / 'data/processed/metadata_clean.csv',\n",
                "    gene_sets_path=PROJECT_ROOT / 'references/mitochondrial_gene_sets.csv',\n",
                "    table_dir=PROJECT_ROOT / 'reports/tables',\n",
                "    figure_dir=PROJECT_ROOT / 'reports/figures',\n",
                ")\n",
                "score_outputs\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/mitochondrial_pathway_score_tests.csv')\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "Interpretation note: pathway scores are mean z-scored expression summaries for selected genes. They are not formal pathway enrichment results and should be read alongside the limitation that no genes passed FDR correction.\n",
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


def run_mitochondrial_gene_set_analysis(
    de_results_path: Path,
    expression_log2_path: Path,
    metadata_path: Path,
    gene_sets_path: Path,
    table_dir: Path,
    figure_dir: Path,
    notebook_path: Path,
) -> dict[str, str]:
    gene_sets = write_gene_sets(gene_sets_path)
    de_results = pd.read_csv(de_results_path)
    expression_log2 = pd.read_csv(expression_log2_path)
    metadata = pd.read_csv(metadata_path)
    mito_results = filter_mitochondrial_de_results(de_results, gene_sets)
    pathway_summary = summarize_pathways(mito_results)
    heatmap = prepare_heatmap_matrix(expression_log2, gene_sets)

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    results_path = table_dir / "mitochondrial_gene_results.csv"
    summary_path = table_dir / "mitochondrial_pathway_summary.csv"
    heatmap_matrix_path = table_dir / "mitochondrial_gene_heatmap_matrix.csv"
    heatmap_path = figure_dir / "mitochondrial_gene_heatmap.png"

    mito_results.to_csv(results_path, index=False)
    pathway_summary.to_csv(summary_path, index=False)
    heatmap.to_csv(heatmap_matrix_path)
    plot_heatmap(heatmap, metadata, heatmap_path)
    make_notebook(notebook_path)
    return {
        "gene_sets": str(gene_sets_path),
        "mitochondrial_gene_results": str(results_path),
        "mitochondrial_pathway_summary": str(summary_path),
        "heatmap_matrix": str(heatmap_matrix_path),
        "heatmap": str(heatmap_path),
        "notebook": str(notebook_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run mitochondrial gene-set analysis.")
    parser.add_argument("--de-results", type=Path, default=Path("reports/tables/differential_expression_results.csv"))
    parser.add_argument("--expression-log2", type=Path, default=Path("data/processed/expression_log2.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/metadata_clean.csv"))
    parser.add_argument("--gene-sets", type=Path, default=Path("references/mitochondrial_gene_sets.csv"))
    parser.add_argument("--table-dir", type=Path, default=Path("reports/tables"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--notebook", type=Path, default=Path("notebooks/05_mitochondrial_gene_set_analysis.ipynb"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_mitochondrial_gene_set_analysis(
        de_results_path=args.de_results,
        expression_log2_path=args.expression_log2,
        metadata_path=args.metadata,
        gene_sets_path=args.gene_sets,
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
        notebook_path=args.notebook,
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
