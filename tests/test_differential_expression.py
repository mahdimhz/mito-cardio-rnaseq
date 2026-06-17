import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.differential_expression import (
    benjamini_hochberg,
    compute_log2_cpm,
    paired_gene_statistics,
    run_differential_expression,
)


class DifferentialExpressionTests(unittest.TestCase):
    def test_compute_log2_cpm_returns_finite_normalized_values(self):
        expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2"],
                "gene_symbol": ["GENE1", "GENE2"],
                "S1": [100, 900],
                "S2": [200, 1800],
            }
        )

        log_cpm = compute_log2_cpm(expression)

        self.assertEqual(list(log_cpm.columns), ["gene_id", "gene_symbol", "S1", "S2"])
        self.assertTrue((log_cpm[["S1", "S2"]] > 0).all().all())
        self.assertAlmostEqual(float(log_cpm.loc[0, "S1"]), float(log_cpm.loc[0, "S2"]), places=6)

    def test_benjamini_hochberg_is_monotonic_and_bounded(self):
        adjusted = benjamini_hochberg(pd.Series([0.001, 0.02, 0.5, 1.0]))

        self.assertEqual(len(adjusted), 4)
        self.assertTrue(((adjusted >= 0) & (adjusted <= 1)).all())
        self.assertTrue(adjusted.iloc[0] <= adjusted.iloc[1] <= adjusted.iloc[2] <= adjusted.iloc[3])

    def test_paired_gene_statistics_reports_disease_minus_control_effect(self):
        log_cpm = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2"],
                "gene_symbol": ["UP_IN_DISEASE", "FLAT"],
                "D1": [10.0, 5.0],
                "C1": [1.0, 5.0],
                "D2": [11.0, 5.0],
                "C2": [2.0, 5.0],
                "D3": [12.0, 5.0],
                "C3": [3.0, 5.0],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["D1", "C1", "D2", "C2", "D3", "C3"],
                "group": ["disease", "control", "disease", "control", "disease", "control"],
                "paired_patient_id": ["1", "1", "2", "2", "3", "3"],
            }
        )

        result = paired_gene_statistics(log_cpm, metadata, min_mean_cpm=0)

        top = result.iloc[0]
        self.assertEqual(top["gene_symbol"], "UP_IN_DISEASE")
        self.assertGreater(float(top["log2_fold_change"]), 0)
        self.assertEqual(int(top["n_pairs"]), 3)
        self.assertIn("padj", result.columns)

    def test_run_differential_expression_writes_expected_outputs(self):
        expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2"],
                "gene_symbol": ["UP_IN_DISEASE", "FLAT"],
                "D1": [1000, 100],
                "C1": [100, 100],
                "D2": [1200, 100],
                "C2": [100, 100],
                "D3": [1400, 100],
                "C3": [100, 100],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["D1", "C1", "D2", "C2", "D3", "C3"],
                "group": ["disease", "control", "disease", "control", "disease", "control"],
                "paired_patient_id": ["1", "1", "2", "2", "3", "3"],
            }
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            expression_path = tmp / "expression_clean.csv"
            metadata_path = tmp / "metadata_clean.csv"
            expression.to_csv(expression_path, index=False)
            metadata.to_csv(metadata_path, index=False)

            outputs = run_differential_expression(
                expression_path=expression_path,
                metadata_path=metadata_path,
                table_dir=tmp,
                figure_dir=tmp,
                notebook_path=tmp / "notebook.ipynb",
            )

            self.assertTrue(Path(outputs["results"]).exists())
            self.assertTrue(Path(outputs["top_genes"]).exists())
            self.assertTrue(Path(outputs["volcano"]).exists())
            self.assertTrue(Path(outputs["summary"]).exists())
            self.assertTrue((tmp / "notebook.ipynb").exists())


if __name__ == "__main__":
    unittest.main()
