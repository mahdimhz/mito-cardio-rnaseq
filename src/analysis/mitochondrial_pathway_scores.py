from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


GENE_COLUMNS = {"gene_id", "gene_symbol"}


def sample_columns(expression: pd.DataFrame) -> list[str]:
    return [column for column in expression.columns if column not in GENE_COLUMNS]


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce").fillna(1.0).clip(lower=0, upper=1)
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running_min = 1.0
    for rank in range(n, 0, -1):
        running_min = min(running_min, ranked[rank - 1] * n / rank)
        adjusted[rank - 1] = running_min
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return pd.Series(out, index=p_values.index)


def gene_z_scores(expression_log2: pd.DataFrame) -> pd.DataFrame:
    samples = sample_columns(expression_log2)
    collapsed = expression_log2.groupby("gene_symbol", sort=False)[samples].mean()
    means = collapsed.mean(axis=1)
    stds = collapsed.std(axis=1).replace(0, 1)
    return collapsed.sub(means, axis=0).div(stds, axis=0).fillna(0)


def compute_pathway_scores(
    expression_log2: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_sets: pd.DataFrame,
) -> pd.DataFrame:
    z_scores = gene_z_scores(expression_log2)
    metadata_by_sample = metadata.set_index("sample_id", drop=False)
    rows = []
    for category, category_genes in gene_sets.groupby("category", sort=True):
        requested = category_genes["gene_symbol"].drop_duplicates().tolist()
        available = [gene for gene in requested if gene in z_scores.index]
        if not available:
            continue
        category_scores = z_scores.loc[available].mean(axis=0)
        for sample_id, score in category_scores.items():
            if sample_id not in metadata_by_sample.index:
                continue
            sample_meta = metadata_by_sample.loc[sample_id]
            rows.append(
                {
                    "category": category,
                    "sample_id": sample_id,
                    "group": sample_meta.get("group", ""),
                    "paired_patient_id": str(sample_meta.get("paired_patient_id", "")),
                    "cell_line": sample_meta.get("cell_line", ""),
                    "replicate": sample_meta.get("replicate", ""),
                    "pathway_score": float(score),
                    "available_gene_count": len(available),
                    "requested_gene_count": len(requested),
                    "available_genes": ";".join(available),
                }
            )
    return pd.DataFrame(rows).sort_values(["category", "paired_patient_id", "group", "sample_id"]).reset_index(drop=True)


def paired_pathway_tests(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for category, category_scores in scores.groupby("category", sort=True):
        patient_means = (
            category_scores.groupby(["paired_patient_id", "group"], as_index=False)["pathway_score"].mean()
        )
        pivot = patient_means.pivot(index="paired_patient_id", columns="group", values="pathway_score")
        if not {"disease", "control"}.issubset(pivot.columns):
            continue
        paired = pivot[["disease", "control"]].dropna()
        n_pairs = len(paired)
        if n_pairs < 2:
            p_value = np.nan
            statistic = np.nan
        else:
            differences = paired["disease"] - paired["control"]
            if float(differences.std(ddof=1)) == 0.0:
                statistic = np.nan
                p_value = 1.0 if float(differences.mean()) == 0.0 else np.nan
            else:
                statistic, p_value = stats.ttest_rel(paired["disease"], paired["control"])
        mean_difference = float((paired["disease"] - paired["control"]).mean()) if n_pairs else np.nan
        if mean_difference > 0.05:
            direction = "higher_in_disease"
        elif mean_difference < -0.05:
            direction = "lower_in_disease"
        else:
            direction = "near_zero"
        rows.append(
            {
                "category": category,
                "n_pairs": n_pairs,
                "mean_score_disease": float(paired["disease"].mean()) if n_pairs else np.nan,
                "mean_score_control": float(paired["control"].mean()) if n_pairs else np.nan,
                "mean_difference_disease_minus_control": mean_difference,
                "statistic": float(statistic) if np.isfinite(statistic) else np.nan,
                "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                "padj": np.nan,
                "direction": direction,
                "method": "paired_t_test_on_patient_mean_pathway_scores",
                "notes": "Exploratory pathway-score summary; scores are mean z-scored expression of available category genes.",
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["padj"] = benjamini_hochberg(result["p_value"])
    return result


def plot_paired_scores(scores: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if scores.empty:
        return
    categories = sorted(scores["category"].unique())
    fig_width = max(10, len(categories) * 0.82)
    fig, ax = plt.subplots(figsize=(fig_width, 5.8))
    colors = {"control": "#4c72b0", "disease": "#c44e52"}
    offsets = {"control": -0.14, "disease": 0.14}
    patient_scores = scores.groupby(["category", "paired_patient_id", "group"], as_index=False)["pathway_score"].mean()
    for x_index, category in enumerate(categories):
        category_scores = patient_scores[patient_scores["category"] == category]
        for patient_id, pair in category_scores.groupby("paired_patient_id"):
            pair = pair.set_index("group")
            if {"control", "disease"}.issubset(pair.index):
                ax.plot(
                    [x_index + offsets["control"], x_index + offsets["disease"]],
                    [pair.loc["control", "pathway_score"], pair.loc["disease", "pathway_score"]],
                    color="#888888",
                    linewidth=0.9,
                    alpha=0.55,
                )
        for group in ["control", "disease"]:
            group_scores = category_scores[category_scores["group"] == group]["pathway_score"]
            ax.scatter(
                np.full(len(group_scores), x_index + offsets[group]),
                group_scores,
                color=colors[group],
                s=32,
                alpha=0.9,
                label=group if x_index == 0 else None,
            )
    ax.axhline(0, color="#333333", linewidth=0.8, linestyle=":")
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Mean gene z-score")
    ax.set_title("Mitochondrial Pathway Scores by Patient Pair")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.7)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_mitochondrial_pathway_scores(
    expression_log2_path: Path,
    metadata_path: Path,
    gene_sets_path: Path,
    table_dir: Path,
    figure_dir: Path,
) -> dict[str, str]:
    expression_log2 = pd.read_csv(expression_log2_path)
    metadata = pd.read_csv(metadata_path)
    gene_sets = pd.read_csv(gene_sets_path)
    scores = compute_pathway_scores(expression_log2, metadata, gene_sets)
    tests = paired_pathway_tests(scores)

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    scores_path = table_dir / "mitochondrial_pathway_scores.csv"
    tests_path = table_dir / "mitochondrial_pathway_score_tests.csv"
    plot_path = figure_dir / "mitochondrial_pathway_scores_paired.png"

    scores.to_csv(scores_path, index=False)
    tests.to_csv(tests_path, index=False)
    plot_paired_scores(scores, plot_path)
    return {
        "pathway_scores": str(scores_path),
        "pathway_score_tests": str(tests_path),
        "paired_plot": str(plot_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run mitochondrial pathway score analysis.")
    parser.add_argument("--expression-log2", type=Path, default=Path("data/processed/expression_log2.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/metadata_clean.csv"))
    parser.add_argument("--gene-sets", type=Path, default=Path("references/mitochondrial_gene_sets.csv"))
    parser.add_argument("--table-dir", type=Path, default=Path("reports/tables"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_mitochondrial_pathway_scores(
        expression_log2_path=args.expression_log2,
        metadata_path=args.metadata,
        gene_sets_path=args.gene_sets,
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
