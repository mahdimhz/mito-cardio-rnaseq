import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.pathway_enrichment import run_pathway_enrichment, test_local_gene_set_enrichment


class PathwayEnrichmentTests(unittest.TestCase):
    def test_local_gene_set_enrichment_uses_tested_background_and_nominal_gene_list(self):
        de = pd.DataFrame(
            {
                "gene_symbol": ["A", "B", "C", "D", "E", "F"],
                "p_value": [0.01, 0.02, 0.50, 0.80, 0.03, 0.70],
                "log2_fold_change": [1.0, 0.5, -0.2, 0.1, -1.0, 0.2],
            }
        )
        gene_sets = pd.DataFrame(
            {
                "gene_symbol": ["A", "B", "C", "F"],
                "category": ["category one", "category one", "category one", "category two"],
                "notes": ["", "", "", ""],
            }
        )

        enrichment = test_local_gene_set_enrichment(de, gene_sets, nominal_p_threshold=0.05)

        first = enrichment.loc[enrichment["category"] == "category one"].iloc[0]
        self.assertEqual(int(first["tested_gene_count"]), 3)
        self.assertEqual(int(first["nominal_gene_count"]), 2)
        self.assertEqual(int(first["background_gene_count"]), 6)
        self.assertEqual(int(first["background_nominal_gene_count"]), 3)
        self.assertGreaterEqual(float(first["p_value"]), 0)
        self.assertLessEqual(float(first["p_value"]), 1)

    def test_run_pathway_enrichment_writes_table_and_dotplot(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            de_path = root / "de.csv"
            gene_sets_path = root / "gene_sets.csv"
            table_dir = root / "tables"
            figure_dir = root / "figures"
            pd.DataFrame(
                {
                    "gene_symbol": ["A", "B", "C", "D", "E", "F"],
                    "p_value": [0.01, 0.02, 0.50, 0.80, 0.03, 0.70],
                    "log2_fold_change": [1.0, 0.5, -0.2, 0.1, -1.0, 0.2],
                }
            ).to_csv(de_path, index=False)
            pd.DataFrame(
                {
                    "gene_symbol": ["A", "B", "C", "F"],
                    "category": ["category one", "category one", "category one", "category two"],
                    "notes": ["", "", "", ""],
                }
            ).to_csv(gene_sets_path, index=False)

            outputs = run_pathway_enrichment(
                de_results_path=de_path,
                gene_sets_path=gene_sets_path,
                table_dir=table_dir,
                figure_dir=figure_dir,
            )

            self.assertTrue(Path(outputs["enrichment_results"]).exists())
            self.assertTrue(Path(outputs["dotplot"]).exists())
            written = pd.read_csv(outputs["enrichment_results"])
            self.assertIn("padj", written.columns)
            self.assertIn("odds_ratio", written.columns)


if __name__ == "__main__":
    unittest.main()
