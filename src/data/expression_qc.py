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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


GENE_COLUMNS = ["gene_id", "gene_symbol"]


def split_gene_identifier(value: str) -> tuple[str, str]:
    """Split GEO gene row names like ENSG00000000003__TSPAN6."""
    text = str(value)
    if "__" in text:
        gene_id, gene_symbol = text.split("__", 1)
        return gene_id, gene_symbol
    return text, text


def load_count_matrix(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", compression="gzip")
    first_column = counts.columns[0]
    gene_parts = counts[first_column].map(split_gene_identifier)
    counts.insert(0, "gene_symbol", [symbol for _, symbol in gene_parts])
    counts.insert(0, "gene_id", [gene_id for gene_id, _ in gene_parts])
    counts = counts.drop(columns=[first_column])
    numeric_columns = [column for column in counts.columns if column not in GENE_COLUMNS]
    counts[numeric_columns] = counts[numeric_columns].apply(pd.to_numeric, errors="coerce")
    return counts


def load_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path)
    required = {"sample_id", "expression_column", "group", "paired_patient_id", "cell_line"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"metadata is missing required columns: {missing}")
    return metadata


def align_expression_to_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    missing = [column for column in metadata["expression_column"].tolist() if column not in counts.columns]
    if missing:
        raise ValueError(f"expression columns missing from count matrix: {missing}")

    aligned = counts[GENE_COLUMNS + metadata["expression_column"].tolist()].copy()
    rename_map = dict(zip(metadata["expression_column"], metadata["sample_id"], strict=True))
    aligned = aligned.rename(columns=rename_map)
    return aligned


def sample_columns(expression: pd.DataFrame) -> list[str]:
    return [column for column in expression.columns if column not in GENE_COLUMNS]


def metadata_annotations(metadata: pd.DataFrame) -> pd.DataFrame:
    columns = ["sample_id", "group", "paired_patient_id", "cell_line", "replicate"]
    annotations = metadata.copy()
    for column in columns:
        if column not in annotations.columns:
            annotations[column] = ""
    return annotations[columns]


def compute_qc_tables(expression: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = sample_columns(expression)
    values = expression[samples]
    missing_values = int(values.isna().sum().sum())
    duplicated_gene_ids = int(expression["gene_id"].duplicated().sum())
    duplicated_gene_symbols = int(expression["gene_symbol"].duplicated().sum())
    negative_values = int((values < 0).sum().sum())

    summary = pd.DataFrame(
        [
            {"metric": "gene_count", "value": str(len(expression)), "notes": "Rows in aligned expression matrix."},
            {"metric": "sample_count", "value": str(len(samples)), "notes": "Samples aligned with cleaned metadata."},
            {"metric": "missing_value_count", "value": str(missing_values), "notes": "Missing values in sample count columns."},
            {"metric": "duplicate_gene_id_count", "value": str(duplicated_gene_ids), "notes": "Duplicate Ensembl gene IDs."},
            {"metric": "duplicate_gene_symbol_count", "value": str(duplicated_gene_symbols), "notes": "Duplicate gene symbols; gene_id remains the primary key."},
            {"metric": "negative_value_count", "value": str(negative_values), "notes": "Negative count values should be zero for count data."},
            {"metric": "matrix_type", "value": "raw_or_estimated_counts", "notes": "GEO describes the supplementary TSV as raw gene counts; decimal values indicate imported/estimated counts."},
        ]
    )

    library_sizes = values.sum(axis=0)
    detected_genes = (values > 0).sum(axis=0)
    sample_metrics = pd.DataFrame({"sample_id": samples})
    sample_metrics["library_size"] = sample_metrics["sample_id"].map(library_sizes).astype(float)
    sample_metrics["detected_genes"] = sample_metrics["sample_id"].map(detected_genes).astype(int)
    sample_metrics["missing_values"] = sample_metrics["sample_id"].map(values.isna().sum(axis=0)).astype(int)
    sample_metrics = sample_metrics.merge(metadata_annotations(metadata), on="sample_id", how="left")
    return summary, sample_metrics


def log2_transform(expression: pd.DataFrame) -> pd.DataFrame:
    transformed = expression.copy()
    samples = sample_columns(transformed)
    transformed[samples] = np.log2(transformed[samples].clip(lower=0) + 1)
    return transformed


def run_pca(log_expression: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, list[float]]:
    """Run PCA on log2(count + 1) expression, returning sample coordinates."""
    samples = sample_columns(log_expression)
    matrix = log_expression[samples].T
    variable_genes = matrix.loc[:, matrix.var(axis=0) > 0]
    if variable_genes.shape[0] < 2 or variable_genes.shape[1] < 2:
        return pd.DataFrame(), []

    scaled = StandardScaler().fit_transform(variable_genes)
    pca = PCA(n_components=2, random_state=0)
    coordinates = pca.fit_transform(scaled)
    pca_df = pd.DataFrame(
        {
            "sample_id": samples,
            "PC1": coordinates[:, 0],
            "PC2": coordinates[:, 1],
        }
    )
    pca_df = pca_df.merge(metadata_annotations(metadata), on="sample_id", how="left")
    return pca_df, [float(value) for value in pca.explained_variance_ratio_]


def save_tables(
    expression: pd.DataFrame,
    log_expression: pd.DataFrame,
    summary: pd.DataFrame,
    sample_metrics: pd.DataFrame,
    pca_df: pd.DataFrame,
    output_dir: Path,
    table_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    expression.to_csv(output_dir / "expression_clean.csv", index=False)
    log_expression.to_csv(output_dir / "expression_log2.csv", index=False)
    summary.to_csv(table_dir / "expression_qc_summary.csv", index=False)
    sample_metrics.to_csv(table_dir / "sample_library_sizes.csv", index=False)
    pca_df.to_csv(table_dir / "pca_coordinates.csv", index=False)


def plot_library_sizes(sample_metrics: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    colors = sample_metrics["group"].map({"disease": "#c44e52", "control": "#4c72b0"}).fillna("#777777")
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.bar(sample_metrics["sample_id"], sample_metrics["library_size"] / 1_000_000, color=colors, width=0.78)
    ax.set_ylabel("Library size (million counts)")
    ax.set_xlabel("Sample")
    ax.set_title("GSE305638 RNA-seq Library Sizes")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_expression_distributions(log_expression: pd.DataFrame, metadata: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    samples = sample_columns(log_expression)
    labels = metadata.set_index("sample_id").loc[samples, "group"].tolist()
    fig, ax = plt.subplots(figsize=(10.5, 4.2))
    ax.boxplot(
        [log_expression[sample].dropna().values for sample in samples],
        tick_labels=samples,
        showfliers=False,
        widths=0.55,
        medianprops={"color": "#333333", "linewidth": 1.1},
    )
    for tick, group in zip(ax.get_xticklabels(), labels, strict=True):
        tick.set_color("#c44e52" if group == "disease" else "#4c72b0")
    ax.set_ylabel("log2(count + 1)")
    ax.set_xlabel("Sample")
    ax.set_title("GSE305638 Expression Distributions")
    ax.tick_params(axis="x", rotation=90, labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, alpha=0.8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_pca(pca_df: pd.DataFrame, explained: list[float], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for group, color in [("disease", "#c44e52"), ("control", "#4c72b0")]:
        subset = pca_df[pca_df["group"] == group]
        ax.scatter(subset["PC1"], subset["PC2"], label=group, color=color, s=58, alpha=0.9)
        for _, row in subset.iterrows():
            ax.annotate(
                str(row["cell_line"]),
                (row["PC1"], row["PC2"]),
                xytext=(4, 3),
                textcoords="offset points",
                fontsize=7,
                alpha=0.75,
            )
    x_label = f"PC1 ({explained[0] * 100:.1f}% variance)" if explained else "PC1"
    y_label = f"PC2 ({explained[1] * 100:.1f}% variance)" if len(explained) > 1 else "PC2"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title("GSE305638 PCA of log2(count + 1) Expression")
    ax.axhline(0, color="#dddddd", linewidth=0.7, zorder=0)
    ax.axvline(0, color="#dddddd", linewidth=0.7, zorder=0)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_figures(
    log_expression: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_metrics: pd.DataFrame,
    pca_df: pd.DataFrame,
    explained: list[float],
    figure_dir: Path,
) -> None:
    plot_library_sizes(sample_metrics, figure_dir / "library_size_by_sample.png")
    plot_expression_distributions(log_expression, metadata, figure_dir / "expression_distribution_boxplot.png")
    if not pca_df.empty:
        plot_pca(pca_df, explained, figure_dir / "pca_plot.png")


def make_notebook(notebook_path: Path) -> None:
    cells = [
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "# Step 6 - Expression Matrix QC\n",
                "\n",
                "Load the selected `GSE305638` processed gene-count matrix, align samples to cleaned metadata, run count-matrix QC, log-transform counts, and generate PCA and QC figures.\n",
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
                "from data.expression_qc import run_expression_qc\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "outputs = run_expression_qc(\n",
                "    counts_path=PROJECT_ROOT / 'data/raw/GSE305638_Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv.gz',\n",
                "    metadata_path=PROJECT_ROOT / 'data/processed/metadata_clean.csv',\n",
                "    processed_dir=PROJECT_ROOT / 'data/processed',\n",
                "    table_dir=PROJECT_ROOT / 'reports/tables',\n",
                "    figure_dir=PROJECT_ROOT / 'reports/figures',\n",
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
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/expression_qc_summary.csv')\n",
            ],
        },
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "pd.read_csv(PROJECT_ROOT / 'reports/tables/pca_coordinates.csv').head()\n",
            ],
        },
        {
            "cell_type": "markdown",
            "metadata": {},
            "source": [
                "The matrix is retained at gene-ID level for downstream analysis. Duplicate gene symbols are reported rather than collapsed at this stage to avoid silently merging Ensembl IDs.\n",
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


def run_expression_qc(
    counts_path: Path,
    metadata_path: Path,
    processed_dir: Path,
    table_dir: Path,
    figure_dir: Path,
) -> dict[str, str]:
    metadata = load_metadata(metadata_path)
    counts = load_count_matrix(counts_path)
    expression = align_expression_to_metadata(counts, metadata)
    summary, sample_metrics = compute_qc_tables(expression, metadata)
    log_expression = log2_transform(expression)
    pca_df, explained = run_pca(log_expression, metadata)
    if not pca_df.empty:
        summary = pd.concat(
            [
                summary,
                pd.DataFrame(
                    [
                        {"metric": "pca_pc1_variance_ratio", "value": f"{explained[0]:.6f}", "notes": "PCA on log2(count + 1), variable genes only."},
                        {"metric": "pca_pc2_variance_ratio", "value": f"{explained[1]:.6f}", "notes": "PCA on log2(count + 1), variable genes only."},
                    ]
                ),
            ],
            ignore_index=True,
        )
    save_tables(expression, log_expression, summary, sample_metrics, pca_df, processed_dir, table_dir)
    save_figures(log_expression, metadata, sample_metrics, pca_df, explained, figure_dir)
    return {
        "expression_clean": str(processed_dir / "expression_clean.csv"),
        "expression_log2": str(processed_dir / "expression_log2.csv"),
        "qc_summary": str(table_dir / "expression_qc_summary.csv"),
        "sample_metrics": str(table_dir / "sample_library_sizes.csv"),
        "pca_coordinates": str(table_dir / "pca_coordinates.csv"),
        "library_size_figure": str(figure_dir / "library_size_by_sample.png"),
        "distribution_figure": str(figure_dir / "expression_distribution_boxplot.png"),
        "pca_figure": str(figure_dir / "pca_plot.png"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Step 6 expression matrix QC.")
    parser.add_argument("--counts", type=Path, default=Path("data/raw/GSE305638_Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv.gz"))
    parser.add_argument("--metadata", type=Path, default=Path("data/processed/metadata_clean.csv"))
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--table-dir", type=Path, default=Path("reports/tables"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--notebook-output", type=Path, default=Path("notebooks/03_expression_qc.ipynb"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_expression_qc(args.counts, args.metadata, args.processed_dir, args.table_dir, args.figure_dir)
    make_notebook(args.notebook_output)
    print(json.dumps(outputs, indent=2))
    print(f"wrote notebook -> {args.notebook_output}")


if __name__ == "__main__":
    main()
