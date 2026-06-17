import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from data.expression_qc import (
    align_expression_to_metadata,
    compute_qc_tables,
    load_count_matrix,
    run_pca,
)


COUNTS_FIXTURE = """\traw_a\traw_b\traw_c
ENSG000001__GENE1\t10\t20\t30
ENSG000002__GENE2\t0\t5\t15
ENSG000003__GENE2\t7\t8\t9
"""


class ExpressionQcTests(unittest.TestCase):
    def test_load_count_matrix_splits_gene_id_and_symbol(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "counts.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(COUNTS_FIXTURE)

            counts = load_count_matrix(path)

        self.assertEqual(list(counts.columns), ["gene_id", "gene_symbol", "raw_a", "raw_b", "raw_c"])
        self.assertEqual(counts.loc[0, "gene_id"], "ENSG000001")
        self.assertEqual(counts.loc[0, "gene_symbol"], "GENE1")
        self.assertEqual(float(counts.loc[1, "raw_b"]), 5.0)

    def test_align_expression_to_metadata_renames_columns_to_sample_ids(self):
        counts = pd.DataFrame(
            {
                "gene_id": ["ENSG1"],
                "gene_symbol": ["GENE1"],
                "raw_b": [2],
                "raw_a": [1],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["GSM_A", "GSM_B"],
                "expression_column": ["raw_a", "raw_b"],
                "group": ["disease", "control"],
            }
        )

        aligned = align_expression_to_metadata(counts, metadata)

        self.assertEqual(list(aligned.columns), ["gene_id", "gene_symbol", "GSM_A", "GSM_B"])
        self.assertEqual(float(aligned.loc[0, "GSM_A"]), 1.0)
        self.assertEqual(float(aligned.loc[0, "GSM_B"]), 2.0)

    def test_compute_qc_tables_reports_library_size_and_duplicate_symbols(self):
        expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_symbol": ["GENE1", "GENE2", "GENE2"],
                "GSM_A": [1, 2, 3],
                "GSM_B": [4, 5, 6],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["GSM_A", "GSM_B"],
                "group": ["disease", "control"],
                "paired_patient_id": ["1", "1"],
                "cell_line": ["FA1", "FA1ic"],
            }
        )

        summary, sample_metrics = compute_qc_tables(expression, metadata)

        self.assertEqual(summary.loc[summary["metric"] == "gene_count", "value"].item(), "3")
        self.assertEqual(summary.loc[summary["metric"] == "duplicate_gene_symbol_count", "value"].item(), "1")
        self.assertEqual(sample_metrics.loc[sample_metrics["sample_id"] == "GSM_A", "library_size"].item(), 6.0)
        self.assertEqual(sample_metrics.loc[sample_metrics["sample_id"] == "GSM_B", "library_size"].item(), 15.0)

    def test_run_pca_returns_sample_coordinates(self):
        log_expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_symbol": ["GENE1", "GENE2", "GENE3"],
                "GSM_A": [1.0, 2.0, 3.0],
                "GSM_B": [2.0, 3.0, 4.0],
                "GSM_C": [4.0, 3.0, 2.0],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["GSM_A", "GSM_B", "GSM_C"],
                "group": ["disease", "control", "disease"],
                "paired_patient_id": ["1", "1", "2"],
                "cell_line": ["FA1", "FA1ic", "FA2"],
            }
        )

        pca_df, explained = run_pca(log_expression, metadata)

        self.assertEqual(set(pca_df["sample_id"]), {"GSM_A", "GSM_B", "GSM_C"})
        self.assertEqual({"PC1", "PC2"}.issubset(pca_df.columns), True)
        self.assertEqual(len(explained), 2)


if __name__ == "__main__":
    unittest.main()

