import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from data.screen_dataset_registry import infer_screening_fields, score_candidate


class ScreenDatasetRegistryTests(unittest.TestCase):
    def test_score_candidate_prioritizes_true_merrf_ipsc_cardiomyocyte_expression(self):
        row = {
            "title": "MERRF m.8344A>G MT-TK patient iPSC-derived cardiomyocyte RNA-seq",
            "organism": "Homo sapiens",
            "assay_type": "Expression profiling by high throughput sequencing",
            "processed_matrix_available": "yes",
            "raw_fastq_available": "yes",
            "metadata_available": "yes",
            "control_count": "3",
            "case_count": "3",
        }

        self.assertEqual(score_candidate(row), 10)

    def test_score_candidate_penalizes_animal_only_unclear_design_without_omics(self):
        row = {
            "title": "Mouse imaging assay",
            "organism": "Mus musculus",
            "assay_type": "Genome variation profiling",
            "processed_matrix_available": "no",
            "raw_fastq_available": "yes",
            "metadata_available": "unknown",
            "control_count": "",
            "case_count": "",
        }

        self.assertEqual(score_candidate(row), -9)

    def test_infer_screening_fields_extracts_mitochondrial_disease_context(self):
        base_row = {
            "title": "Barth syndrome hiPSC cardiomyocyte RNA-seq",
            "organism": "Homo sapiens",
            "assay_type": "Expression profiling by high throughput sequencing",
            "processed_matrix_available": "yes",
            "metadata_available": "yes",
        }
        summary = {
            "summary": "Patient-derived induced pluripotent stem cell cardiomyocytes with TAZ mutation.",
            "samples": [
                {"title": "control cardiomyocyte replicate 1"},
                {"title": "Barth patient cardiomyocyte replicate 1"},
            ],
        }

        inferred = infer_screening_fields(base_row, summary, raw_fastq_available="yes")

        self.assertEqual(inferred["disease"], "Barth syndrome")
        self.assertEqual(inferred["mutation"], "TAZ")
        self.assertEqual(inferred["cell_type"], "human iPSC-derived cardiomyocytes")
        self.assertEqual(inferred["case_count"], "1")
        self.assertEqual(inferred["control_count"], "1")
        self.assertIn("case/control counts are heuristic", inferred["notes"])


if __name__ == "__main__":
    unittest.main()
