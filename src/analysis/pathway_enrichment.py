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


def gene_level_de_table(de_results: pd.DataFrame) -> pd.DataFrame:
    required = {"gene_symbol", "p_value", "log2_fold_change"}
    missing = sorted(required - set(de_results.columns))
    if missing:
        raise ValueError(f"expression result table is missing required columns: {missing}")
    de = de_results.copy()
    de["gene_symbol"] = de["gene_symbol"].astype(str)
    de["p_value"] = pd.to_numeric(de["p_value"], errors="coerce")
    de["abs_log2_fold_change"] = pd.to_numeric(de["log2_fold_change"], errors="coerce").abs()
    de = de[de["gene_symbol"].str.len() > 0].dropna(subset=["p_value"])
    de = de.sort_values(["p_value", "abs_log2_fold_change"], ascending=[True, False])
    return de.drop_duplicates("gene_symbol", keep="first")


def test_local_gene_set_enrichment(
    de_results: pd.DataFrame,
    gene_sets: pd.DataFrame,
    nominal_p_threshold: float = 0.05,
) -> pd.DataFrame:
    de = gene_level_de_table(de_results)
    background = set(de["gene_symbol"])
    nominal = set(de.loc[de["p_value"] < nominal_p_threshold, "gene_symbol"])
    rows = []
    for category, category_genes in gene_sets.groupby("category", sort=True):
        genes = set(category_genes["gene_symbol"].astype(str)) & background
        if not genes:
            continue
        overlap = genes & nominal
        table = [
            [len(overlap), len(genes - nominal)],
            [len(nominal - genes), len(background - genes - nominal)],
        ]
        odds_ratio, p_value = stats.fisher_exact(table, alternative="greater")
        category_de = de[de["gene_symbol"].isin(genes)]
        mean_lfc = float(category_de["log2_fold_change"].mean())
        rows.append(
            {
                "category": category,
                "tested_gene_count": len(genes),
                "nominal_gene_count": len(overlap),
                "background_gene_count": len(background),
                "background_nominal_gene_count": len(nominal),
                "overlap_genes": ";".join(sorted(overlap)),
                "odds_ratio": float(odds_ratio) if np.isfinite(odds_ratio) else np.inf,
                "p_value": float(p_value),
                "padj": np.nan,
                "mean_log2_fold_change": mean_lfc,
                "direction": "higher_in_disease" if mean_lfc > 0.05 else "lower_in_disease" if mean_lfc < -0.05 else "near_zero",
                "method": "fisher_exact_over_representation_nominal_p_lt_0.05",
                "notes": "Exploratory local mitochondrial-category over-representation; not broad pathway discovery.",
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["padj"] = benjamini_hochberg(result["p_value"])
        result = result.sort_values(["padj", "p_value", "category"]).reset_index(drop=True)
    return result


def plot_enrichment_dotplot(enrichment: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if enrichment.empty:
        return
    plot_df = enrichment.sort_values(["padj", "p_value"]).head(15).copy()
    plot_df = plot_df.sort_values("p_value", ascending=False)
    y = np.arange(len(plot_df))
    x = -np.log10(plot_df["p_value"].clip(lower=1e-300))
    sizes = 35 + 18 * pd.to_numeric(plot_df["nominal_gene_count"], errors="coerce").fillna(0)
    fig_height = max(4.5, len(plot_df) * 0.45)
    fig, ax = plt.subplots(figsize=(8.4, fig_height))
    scatter = ax.scatter(x, y, s=sizes, c=plot_df["padj"], cmap="viridis_r", edgecolor="#333333", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["category"], fontsize=8)
    ax.set_xlabel("-log10 enrichment p-value")
    ax.set_title("Exploratory Mitochondrial Category Enrichment")
    ax.text(
        0.98,
        0.04,
        "No category passes FDR correction",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#333333",
    )
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("adjusted p-value")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", linewidth=0.6, alpha=0.7)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def run_pathway_enrichment(
    de_results_path: Path,
    gene_sets_path: Path,
    table_dir: Path,
    figure_dir: Path,
) -> dict[str, str]:
    de_results = pd.read_csv(de_results_path)
    gene_sets = pd.read_csv(gene_sets_path)
    enrichment = test_local_gene_set_enrichment(de_results, gene_sets)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    results_path = table_dir / "pathway_enrichment_results.csv"
    dotplot_path = figure_dir / "pathway_enrichment_dotplot.png"
    enrichment.to_csv(results_path, index=False)
    plot_enrichment_dotplot(enrichment, dotplot_path)
    return {"enrichment_results": str(results_path), "dotplot": str(dotplot_path)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run exploratory local mitochondrial-category enrichment.")
    parser.add_argument("--de-results", type=Path, default=Path("reports/tables/differential_expression_results.csv"))
    parser.add_argument("--gene-sets", type=Path, default=Path("references/mitochondrial_gene_sets.csv"))
    parser.add_argument("--table-dir", type=Path, default=Path("reports/tables"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_pathway_enrichment(
        de_results_path=args.de_results,
        gene_sets_path=args.gene_sets,
        table_dir=args.table_dir,
        figure_dir=args.figure_dir,
    )
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
