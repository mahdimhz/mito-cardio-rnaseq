import csv
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from data.search_ncbi_geo import REGISTRY_COLUMNS, row_from_summary, write_registry


class SearchNcbiGeoTests(unittest.TestCase):
    def test_row_from_summary_maps_ncbi_gds_fields(self):
        record = {
            "uid": "200333536",
            "accession": "GSE333536",
            "title": "Mitochondrial dysfunction in cardiomyocytes",
            "taxon": "Homo sapiens",
            "gdstype": "Expression profiling by high throughput sequencing",
            "gpl": "GPL24676",
            "n_samples": 6,
            "pdat": "2026/05/29",
            "suppfile": "CSV; XLSX",
            "ftplink": "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333536/",
            "bioproject": "PRJNA123456",
            "pubmedids": ["12345678"],
            "samples": [{"accession": "GSM1", "title": "control"}],
        }

        row = row_from_summary(record, ["mitochondrial dysfunction cardiomyocyte RNA-seq"])

        self.assertEqual(row["accession"], "GSE333536")
        self.assertEqual(row["source_database"], "NCBI GEO/GDS")
        self.assertEqual(row["organism"], "Homo sapiens")
        self.assertEqual(row["assay_type"], "Expression profiling by high throughput sequencing")
        self.assertEqual(row["platform"], "GPL24676")
        self.assertEqual(row["sample_count"], "6")
        self.assertEqual(row["processed_matrix_available"], "yes")
        self.assertEqual(row["metadata_available"], "yes")
        self.assertEqual(row["publication"], "PMID:12345678")
        self.assertEqual(row["publication_year"], "2026")
        self.assertTrue(row["download_url"].endswith("/GSE333536/"))
        self.assertEqual(row["keep_status"], "candidate_unscreened")
        self.assertIn("Entrez UID: 200333536", row["notes"])
        self.assertIn("matched_queries: mitochondrial dysfunction cardiomyocyte RNA-seq", row["notes"])

    def test_write_registry_uses_required_column_order(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "dataset_registry.csv"
            row = {column: "" for column in REGISTRY_COLUMNS}
            row["accession"] = "GSE1"

            write_registry([row], output)

            with output.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(reader.fieldnames, REGISTRY_COLUMNS)
                rows = list(reader)

        self.assertEqual(rows, [row])


if __name__ == "__main__":
    unittest.main()
