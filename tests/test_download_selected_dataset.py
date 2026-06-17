import csv
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from data.download_selected_dataset import (
    SOURCE_REGISTER_COLUMNS,
    build_download_plan,
    source_register_rows,
    write_source_register,
)


class DownloadSelectedDatasetTests(unittest.TestCase):
    def test_build_download_plan_targets_only_geo_lightweight_files(self):
        plan = build_download_plan("GSE305638")
        urls = [item["url"] for item in plan]
        local_paths = [item["local_path"] for item in plan]

        self.assertEqual(len(plan), 4)
        self.assertIn("data/raw/GSE305638_series_matrix.txt.gz", local_paths)
        self.assertIn("data/raw/GSE305638_Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv.gz", local_paths)
        self.assertIn("data/interim/GSE305638_family.xml.tgz", local_paths)
        self.assertIn("data/interim/GSE305638_gds_summary.json", local_paths)
        self.assertTrue(all("ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305638" in url or "eutils.ncbi.nlm.nih.gov" in url for url in urls))
        self.assertFalse(any("sra" in url.lower() or "fastq" in url.lower() for url in urls))

    def test_source_register_rows_include_accession_and_date(self):
        rows = source_register_rows(build_download_plan("GSE305638"), "2026-06-16")

        self.assertEqual({row["accession"] for row in rows}, {"GSE305638"})
        self.assertEqual({row["date_accessed"] for row in rows}, {"2026-06-16"})
        self.assertTrue(all(row["source_name"] for row in rows))
        self.assertTrue(all(row["source_type"] for row in rows))
        self.assertTrue(all(row["url"] for row in rows))
        self.assertTrue(all(row["local_path"] for row in rows))

    def test_write_source_register_uses_existing_columns(self):
        rows = source_register_rows(build_download_plan("GSE305638"), "2026-06-16")
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "source_register.csv"
            write_source_register(rows, output)
            with output.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(reader.fieldnames, SOURCE_REGISTER_COLUMNS)
                written_rows = list(reader)

        self.assertEqual(len(written_rows), 4)


if __name__ == "__main__":
    unittest.main()

