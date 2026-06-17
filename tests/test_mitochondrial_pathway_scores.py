import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.mitochondrial_pathway_scores import (
    compute_pathway_scores,
    paired_pathway_tests,
    run_mitochondrial_pathway_scores,
)


class MitochondrialPathwayScoresTests(unittest.TestCase):
    def test_compute_pathway_scores_uses_available_gene_z_scores_per_sample(self):
        expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_symbol": ["MT-ND1", "MT-ND2", "MFN1"],
                "S1": [1.0, 2.0, 8.0],
                "S2": [3.0, 4.0, 10.0],
            }
        )
        metadata = pd.DataFrame(
            {
                "sample_id": ["S1", "S2"],
                "group": ["disease", "control"],
                "paired_patient_id": ["1", "1"],
            }
        )
        gene_sets = pd.DataFrame(
            {
                "gene_symbol": ["MT-ND1", "MT-ND2", "MFN1", "NOT_IN_MATRIX"],
                "category": ["mtDNA-encoded OXPHOS genes", "mtDNA-encoded OXPHOS genes", "fusion/fission", "fusion/fission"],
                "notes": ["", "", "", ""],
            }
        )

        scores = compute_pathway_scores(expression, metadata, gene_sets)

        self.assertEqual(set(scores["category"]), {"mtDNA-encoded OXPHOS genes", "fusion/fission"})
        self.assertEqual(set(scores["sample_id"]), {"S1", "S2"})
        oxphos = scores[(scores["category"] == "mtDNA-encoded OXPHOS genes") & (scores["sample_id"] == "S1")].iloc[0]
        self.assertEqual(int(oxphos["available_gene_count"]), 2)
        self.assertLess(float(oxphos["pathway_score"]), 0)

    def test_paired_pathway_tests_compare_patient_level_disease_control_scores(self):
        scores = pd.DataFrame(
            {
                "category": ["fusion/fission"] * 4,
                "sample_id": ["D1", "C1", "D2", "C2"],
                "group": ["disease", "control", "disease", "control"],
                "paired_patient_id": ["1", "1", "2", "2"],
                "pathway_score": [2.0, 1.0, 4.0, 2.0],
                "available_gene_count": [2, 2, 2, 2],
            }
        )

        tests = paired_pathway_tests(scores)

        self.assertEqual(len(tests), 1)
        row = tests.iloc[0]
        self.assertEqual(row["category"], "fusion/fission")
        self.assertEqual(int(row["n_pairs"]), 2)
        self.assertAlmostEqual(float(row["mean_difference_disease_minus_control"]), 1.5)
        self.assertIn("paired_t_test", row["method"])

    def test_run_mitochondrial_pathway_scores_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            expression_path = root / "expression_log2.csv"
            metadata_path = root / "metadata.csv"
            gene_sets_path = root / "gene_sets.csv"
            table_dir = root / "tables"
            figure_dir = root / "figures"

            pd.DataFrame(
                {
                    "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                    "gene_symbol": ["MT-ND1", "MT-ND2", "MFN1"],
                    "D1": [1.0, 2.0, 8.0],
                    "C1": [3.0, 4.0, 10.0],
                    "D2": [2.0, 3.0, 9.0],
                    "C2": [4.0, 5.0, 11.0],
                }
            ).to_csv(expression_path, index=False)
            pd.DataFrame(
                {
                    "sample_id": ["D1", "C1", "D2", "C2"],
                    "group": ["disease", "control", "disease", "control"],
                    "paired_patient_id": ["1", "1", "2", "2"],
                }
            ).to_csv(metadata_path, index=False)
            pd.DataFrame(
                {
                    "gene_symbol": ["MT-ND1", "MT-ND2", "MFN1"],
                    "category": ["mtDNA-encoded OXPHOS genes", "mtDNA-encoded OXPHOS genes", "fusion/fission"],
                    "notes": ["", "", ""],
                }
            ).to_csv(gene_sets_path, index=False)

            outputs = run_mitochondrial_pathway_scores(
                expression_log2_path=expression_path,
                metadata_path=metadata_path,
                gene_sets_path=gene_sets_path,
                table_dir=table_dir,
                figure_dir=figure_dir,
            )

            self.assertTrue(Path(outputs["pathway_scores"]).exists())
            self.assertTrue(Path(outputs["pathway_score_tests"]).exists())
            self.assertTrue(Path(outputs["paired_plot"]).exists())
            written_scores = pd.read_csv(outputs["pathway_scores"])
            self.assertIn("pathway_score", written_scores.columns)


if __name__ == "__main__":
    unittest.main()
