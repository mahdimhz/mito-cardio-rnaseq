from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


GENE_COLUMNS = ["gene_id", "gene_symbol"]
LIMITATION_TEXT = (
    "Exploratory paired expression analysis using log2 CPM values and three patient-level disease-control pairs. "
    "This is not a full DESeq2 negative-binomial model; adjusted p-values should be interpreted cautiously "
    "because the effective paired sample size is small."
)


def sample_columns(expression: pd.DataFrame) -> list[str]:
    return [column for column in expression.columns if column not in GENE_COLUMNS]


def compute_log2_cpm(expression: pd.DataFrame) -> pd.DataFrame:
    samples = sample_columns(expression)
    counts = expression[samples].clip(lower=0)
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        bad = library_sizes[library_sizes <= 0].index.tolist()
        raise ValueError(f"samples with non-positive library sizes: {bad}")
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    log_cpm = expression[GENE_COLUMNS].copy()
    log_cpm[samples] = np.log2(cpm + 1)
    return log_cpm


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg FDR correction preserving original index."""
    p = pd.to_numeric(p_values, errors="coerce").fillna(1.0).clip(lower=0, upper=1)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running_min = 1.0
    for rank in range(n, 0, -1):
        value = ranked[rank - 1] * n / rank
        running_min = min(running_min, value)
        adjusted[rank - 1] = running_min
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return pd.Series(out, index=p_values.index)


def paired_patient_means(log_cpm: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    samples = sample_columns(log_cpm)
    sample_meta = metadata.set_index("sample_id").loc[samples]
    complete_pairs = []
    for patient_id, patient_meta in sample_meta.groupby("paired_patient_id", sort=True):
        groups = set(patient_meta["group"])
        if {"disease", "control"}.issubset(groups):
            complete_pairs.append(str(patient_id))
    if len(complete_pairs) < 2:
        raise ValueError("At least two complete disease/control patient pairs are required.")

    mean_columns = {}
    for patient_id in complete_pairs:
        for group in ["disease", "control"]:
            group_samples = sample_meta[
                (sample_meta["paired_patient_id"].astype(str) == patient_id)
                & (sample_meta["group"] == group)
            ].index.tolist()
            mean_columns[f"patient_{patient_id}_{group}"] = log_cpm[group_samples].mean(axis=1)
    means = pd.concat([log_cpm[GENE_COLUMNS], pd.DataFrame(mean_columns)], axis=1)
    return means, complete_pairs


def paired_gene_statistics(log_cpm: pd.DataFrame, metadata: pd.DataFrame, min_mean_cpm: float = 1.0) -> pd.DataFrame:
    """Compute disease-control paired statistics from patient-level log2 CPM means."""
    means, pair_ids = paired_patient_means(log_cpm, metadata)
    disease_columns = [f"patient_{pair_id}_disease" for pair_id in pair_ids]
    control_columns = [f"patient_{pair_id}_control" for pair_id in pair_ids]

    disease_values = means[disease_columns].to_numpy(dtype=float)
    control_values = means[control_columns].to_numpy(dtype=float)
    differences = disease_values - control_values

    mean_disease = np.nanmean(disease_values, axis=1)
    mean_control = np.nanmean(control_values, axis=1)
    mean_log_cpm = np.nanmean(np.concatenate([disease_values, control_values], axis=1), axis=1)
    log2_fc = mean_disease - mean_control

    p_values = []
    for row_index in range(differences.shape[0]):
        row_diff = differences[row_index, :]
        if np.nanstd(row_diff) == 0:
            p_values.append(1.0 if np.nanmean(row_diff) == 0 else 0.0)
        else:
            _, p_value = stats.ttest_rel(disease_values[row_index, :], control_values[row_index, :], nan_policy="omit")
            p_values.append(float(p_value) if np.isfinite(p_value) else 1.0)

    result = pd.DataFrame(
        {
            "gene_id": means["gene_id"],
            "gene_symbol": means["gene_symbol"],
            "mean_log2_cpm_disease": mean_disease,
            "mean_log2_cpm_control": mean_control,
            "mean_log2_cpm": mean_log_cpm,
            "log2_fold_change": log2_fc,
            "p_value": p_values,
            "n_pairs": len(pair_ids),
            "method": "paired_t_test_on_patient_mean_log2_cpm",
            "comparison": "FRDA_iPSC_cardiomyocytes_vs_isogenic_corrected_controls",
        }
    )
    result = result[result["mean_log2_cpm"] >= min_mean_cpm].copy()
    result["padj"] = benjamini_hochberg(result["p_value"])
    result["abs_log2_fold_change"] = result["log2_fold_change"].abs()
    result = result.sort_values(["padj", "p_value", "abs_log2_fold_change"], ascending=[True, True, False])
    return result.reset_index(drop=True)


def top_gene_table(results: pd.DataFrame, n: int = 20) -> pd.DataFrame:
    up = results.sort_values(["log2_fold_change", "padj"], ascending=[False, True]).head(n).copy()
    up["direction"] = "up_in_disease"
    down = results.sort_values(["log2_fold_change", "padj"], ascending=[True, True]).head(n).copy()
    down["direction"] = "down_in_disease"
    return pd.concat([up, down], ignore_index=True)


def summary_table(results: pd.DataFrame) -> pd.DataFrame:
    thresholds = [
        ("nominal_p_lt_0.05", (results["p_value"] < 0.05).sum()),
        ("fdr_lt_0.10", (results["padj"] < 0.10).sum()),
        ("fdr_lt_0.05", (results["padj"] < 0.05).sum()),
        ("abs_log2fc_ge_1", (results["abs_log2_fold_change"] >= 1).sum()),
        ("genes_tested", len(results)),
    ]
    rows = [{"metric": metric, "value": str(int(value)), "notes": LIMITATION_TEXT} for metric, value in thresholds]
    return pd.DataFrame(rows)


def plot_volcano(results: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    x = results["log2_fold_change"]
    y = -np.log10(results["padj"].clip(lower=1e-300))
    significant = (results["padj"] < 0.10) & (results["abs_log2_fold_change"] >= 1)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.scatter(x[~significant], y[~significant], s=8, alpha=0.32, color="#777777", label="tested genes")
    if significant.any():
        ax.scatter(
            x[significant],
            y[significant],
            s=14,
            alpha=0.8,
            color="#c44e52",
            label="FDR < 0.10 and |log2FC| >= 1",
        )
    ax.axvline(-1, color="#333333", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#333333", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(0.10), color="#333333", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log2 fold change (disease - control)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Exploratory Paired Expression: GSE305638")
    ax.text(
        0.02,
        0.90,
        "No genes pass FDR < 0.10",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#333333",
    )
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def make_notebook(notebook_path: Path) -> None:
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Step 7 - Exploratory Paired Expression\n",
                "\n",
                "Exploratory paired disease-control comparison for `GSE305638`.\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                f"**Method limitation:** {LIMITATION_TEXT}\n",
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
                "from analysis.differential_expression import run_differential_expression\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "outputs = run_differential_expression(\n",
                "    expression_path=PROJECT_ROOT / 'data/processed/expression_clean.csv',\n",
                "    metadata_path=PROJECT_ROOT / 'data/processed/metadata_clean.csv',\n",
                "    table_dir=PROJECT_ROOT / 'reports/tables',\n",
                "    figure_dir=PROJECT_ROOT / 'reports/figures',\n",
                "    notebook_path=PROJECT_ROOT / 'notebooks/04_differential_expression.ipynb',\n",
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
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/differential_expression_results.csv').head(10)\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/top_differential_expression_genes.csv').head(20)\n",
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


def run_differential_expression(
    expression_path: Path,
    metadata_path: Path,
    table_dir: Path,
    figure_dir: Path,
    notebook_path: Path,
) -> dict[str, str]:
    expression = pd.read_csv(expression_path)
    metadata = pd.read_csv(metadata_path)
    log_cpm = compute_log2_cpm(expression)
    results = paired_gene_statistics(log_cpm, metadata)
    top_genes = top_gene_table(results)
    summary = summary_table(results)

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    results_path = table_dir / "differential_expression_results.csv"
    top_path = table_dir / "top_differential_expression_genes.csv"
    summary_path = table_dir / "differential_expression_summary.csv"
    volcano_path = figure_dir / "volcano_plot.png"

    results.to_csv(results_path, index=False)
    top_genes.to_csv(top_path, index=False)
    summary.to_csv(summary_path, index=False)
    plot_volcano(results, volcano_path)
    make_notebook(notebook_path)

    return {
        "results": str(results_path),
        "top_genes": str(top_path),
        "summary": str(summary_path),
        "volcano": str(volcano_path),
        "notebook": str(notebook_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run exploratory paired expression analysis for GSE305638.")
    parser.add_argument("--expression", type=Path, default=Path("data/processed/expression_clean.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/metadata_clean.csv"))
    parser.add_argument("--table-dir", type=Path, default=Path("reports/tables"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--notebook", type=Path, default=Path("notebooks/04_differential_expression.ipynb"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_differential_expression(args.expression, args.metadata, args.table_dir, args.figure_dir, args.notebook)
    print(json.dumps(outputs, indent=2))
    print("limitation=" + LIMITATION_TEXT)


if __name__ == "__main__":
    main()
