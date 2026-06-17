import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.mitochondrial_gene_set_analysis import (
    build_default_gene_sets,
    filter_mitochondrial_de_results,
    prepare_heatmap_matrix,
    summarize_pathways,
    write_gene_sets,
)


class MitochondrialGeneSetAnalysisTests(unittest.TestCase):
    def test_build_default_gene_sets_contains_required_genes_and_categories(self):
        gene_sets = build_default_gene_sets()

        self.assertIn("MT-TK", set(gene_sets["gene_symbol"]))
        self.assertIn("ATP5F1A", set(gene_sets["gene_symbol"]))
        self.assertIn("fusion/fission", set(gene_sets["category"]))
        self.assertIn("mitochondrial biogenesis", set(gene_sets["category"]))
        self.assertEqual(len(set(gene_sets["gene_symbol"])), 41)

    def test_write_gene_sets_uses_expected_columns(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "mitochondrial_gene_sets.csv"
            write_gene_sets(path)

            written = pd.read_csv(path)

        self.assertEqual(list(written.columns), ["gene_symbol", "category", "notes"])
        self.assertEqual(len(set(written["gene_symbol"])), 41)

    def test_filter_mitochondrial_de_results_keeps_missing_de_genes_with_status(self):
        gene_sets = pd.DataFrame(
            {
                "gene_symbol": ["MT-TK", "MFN1"],
                "category": ["mitochondrial translation / mtDNA maintenance", "fusion/fission"],
                "notes": ["", ""],
            }
        )
        de = pd.DataFrame(
            {
                "gene_id": ["ENSG1"],
                "gene_symbol": ["MFN1"],
                "log2_fold_change": [0.5],
                "p_value": [0.01],
                "padj": [0.2],
                "mean_log2_cpm": [5.0],
            }
        )

        filtered = filter_mitochondrial_de_results(de, gene_sets)

        self.assertEqual(set(filtered["gene_symbol"]), {"MT-TK", "MFN1"})
        self.assertEqual(filtered.loc[filtered["gene_symbol"] == "MFN1", "result_status"].item(), "tested")
        self.assertEqual(filtered.loc[filtered["gene_symbol"] == "MT-TK", "result_status"].item(), "not_tested_in_de_filter")

    def test_summarize_pathways_reports_direction(self):
        mito = pd.DataFrame(
            {
                "category": ["fusion/fission", "fusion/fission", "mitophagy"],
                "gene_symbol": ["MFN1", "MFN2", "PINK1"],
                "log2_fold_change": [1.0, 0.5, -0.2],
                "p_value": [0.01, 0.02, 0.5],
                "padj": [0.2, 0.2, 0.8],
                "result_status": ["tested", "tested", "tested"],
            }
        )

        summary = summarize_pathways(mito)

        fusion = summary.loc[summary["category"] == "fusion/fission"].iloc[0]
        self.assertEqual(int(fusion["tested_gene_count"]), 2)
        self.assertEqual(fusion["direction"], "up_in_disease")

    def test_prepare_heatmap_matrix_collapses_duplicate_symbols_and_orders_genes(self):
        expression = pd.DataFrame(
            {
                "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_symbol": ["MFN1", "MFN1", "PINK1"],
                "S1": [1.0, 3.0, 2.0],
                "S2": [5.0, 7.0, 2.0],
            }
        )
        gene_sets = pd.DataFrame(
            {
                "gene_symbol": ["PINK1", "MFN1"],
                "category": ["mitophagy", "fusion/fission"],
                "notes": ["", ""],
            }
        )

        heatmap = prepare_heatmap_matrix(expression, gene_sets)

        self.assertEqual(list(heatmap.index), ["PINK1", "MFN1"])
        self.assertEqual(list(heatmap.columns), ["S1", "S2"])
        self.assertFalse(heatmap.isna().any().any())


if __name__ == "__main__":
    unittest.main()

