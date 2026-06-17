import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.results_summary import (
    build_artifact_inventory,
    build_summary_sections,
    make_notebook,
    write_summary_tables,
)


class ResultsSummaryTests(unittest.TestCase):
    def test_build_summary_sections_mentions_dataset_methods_limitations_and_relevance(self):
        inputs = {
            "selected_dataset": {
                "accession": "GSE305638",
                "title": "Frataxin deficiency drives cardiac dysfunction",
                "disease": "Friedreich ataxia / frataxin-related mitochondrial disease",
                "cell_type": "human iPSC-derived cardiomyocytes",
            },
            "metadata_counts": {"disease": 9, "control": 9, "patients": 3},
            "qc": {"gene_count": "59743", "missing_value_count": "0"},
            "de": {"genes_tested": "15100", "fdr_lt_0.10": "0", "nominal_p_lt_0.05": "859"},
            "mito": {
                "gene_count": 41,
                "tested": 40,
                "min_padj": 0.67,
                "top_direction": "mtDNA-encoded OXPHOS genes trend lower in disease",
            },
        }

        sections = build_summary_sections(inputs)
        joined = "\n".join(section["body"] for section in sections)

        self.assertIn("GSE305638", joined)
        self.assertIn("log2 CPM", joined)
        self.assertIn("not MERRF-specific", joined)
        self.assertIn("related mitochondrial cardiomyopathy model", joined)
        self.assertIn("clear disease-model framing", joined)
        self.assertIn("FDR", joined)

    def test_build_artifact_inventory_lists_tables_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            (root / "reports/tables").mkdir(parents=True)
            (root / "reports/figures").mkdir(parents=True)
            (root / "reports/tables/example.csv").write_text("a,b\n1,2\n", encoding="utf-8")
            (root / "reports/figures/example.png").write_bytes(b"png")

            inventory = build_artifact_inventory(root)

        self.assertEqual(set(inventory["artifact_type"]), {"table", "figure"})
        self.assertIn("reports/tables/example.csv", set(inventory["relative_path"]))
        self.assertIn("reports/figures/example.png", set(inventory["relative_path"]))

    def test_write_summary_tables_creates_expected_csvs(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            table_dir = root / "reports/tables"
            table_dir.mkdir(parents=True)
            sections = [{"section": "Dataset", "body": "Selected dataset summary."}]
            inventory = pd.DataFrame(
                [{"artifact_type": "table", "relative_path": "reports/tables/example.csv", "bytes": 10}]
            )

            outputs = write_summary_tables(sections, inventory, table_dir)

            self.assertTrue(Path(outputs["project_summary_table"]).exists())
            self.assertTrue(Path(outputs["artifact_inventory"]).exists())

    def test_make_notebook_writes_valid_notebook(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            notebook_path = Path(tmp_dir) / "06_results_summary.ipynb"
            sections = [{"section": "Dataset", "body": "Selected dataset summary."}]
            make_notebook(notebook_path, sections)

            text = notebook_path.read_text(encoding="utf-8")

        self.assertIn("Step 9 - Final Results Summary", text)
        self.assertIn("Selected dataset summary.", text)


if __name__ == "__main__":
    unittest.main()
