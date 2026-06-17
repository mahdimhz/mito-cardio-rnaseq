import gzip
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from data.metadata_qc import (
    clean_sample_metadata,
    parse_series_matrix,
    write_metadata_outputs,
)


SERIES_FIXTURE = """!Series_title\t\"Example\"
!Sample_title\t\"Cardiomyocytes, FA1, Disease, 1\"\t\"Cardiomyocytes, FA1ic, Healthy, 1\"
!Sample_geo_accession\t\"GSM1\"\t\"GSM2\"
!Sample_source_name_ch1\t\"FA1\"\t\"FA1ic\"
!Sample_organism_ch1\t\"Homo sapiens\"\t\"Homo sapiens\"
!Sample_characteristics_ch1\t\"patient: 1\"\t\"patient: 1\"
!Sample_characteristics_ch1\t\"cell line: FA1\"\t\"cell line: FA1ic\"
!Sample_characteristics_ch1\t\"cell type: iPSC-cardiomyocytes\"\t\"cell type: iPSC-cardiomyocytes\"
!Sample_characteristics_ch1\t\"genotype: WT\"\t\"genotype: corrected\"
!Sample_description\t\"Library name: Sample 1\"\t\"Library name: Sample 2\"
!Sample_description\t\"Column name in \"\"Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv\"\": JL_poolC_b_07_S19\"\t\"Column name in \"\"Lees_FRDA_CM_RNAseq_GeneCounts_AllSamples.tsv\"\": JL_poolC_b_10_S22\"
!Sample_instrument_model\t\"NextSeq 2000\"\t\"NextSeq 2000\"
"""


class MetadataQcTests(unittest.TestCase):
    def test_parse_series_matrix_collects_repeated_sample_fields(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "series.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(SERIES_FIXTURE)

            records = parse_series_matrix(path)

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0]["sample_geo_accession"], "GSM1")
        self.assertEqual(records[0]["sample_title"], "Cardiomyocytes, FA1, Disease, 1")
        self.assertEqual(records[0]["patient"], "1")
        self.assertEqual(records[0]["cell_line"], "FA1")
        self.assertEqual(records[0]["genotype"], "WT")
        self.assertEqual(records[0]["expression_column"], "JL_poolC_b_07_S19")

    def test_clean_sample_metadata_assigns_groups_and_ids(self):
        records = parse_series_matrix_from_fixture()

        clean = clean_sample_metadata(records)

        self.assertEqual(clean[0]["sample_id"], "GSM1")
        self.assertEqual(clean[0]["group"], "disease")
        self.assertEqual(clean[0]["disease_label"], "Friedreich ataxia")
        self.assertEqual(clean[0]["mutation"], "FXN GAA repeat expansion")
        self.assertEqual(clean[0]["replicate"], "1")
        self.assertEqual(clean[1]["group"], "control")
        self.assertEqual(clean[1]["disease_label"], "isogenic corrected control")
        self.assertEqual(clean[1]["mutation"], "CRISPR-corrected FXN")
        self.assertEqual(clean[1]["paired_patient_id"], "1")

    def test_write_metadata_outputs_creates_clean_metadata_and_summary(self):
        clean = clean_sample_metadata(parse_series_matrix_from_fixture())
        with tempfile.TemporaryDirectory() as tmp_dir:
            processed = Path(tmp_dir) / "metadata_clean.csv"
            summary = Path(tmp_dir) / "sample_summary.csv"
            write_metadata_outputs(clean, processed, summary)

            metadata_text = processed.read_text(encoding="utf-8")
            summary_text = summary.read_text(encoding="utf-8")

        self.assertIn("sample_id,sample_title,expression_column", metadata_text.splitlines()[0])
        self.assertIn("disease,1,1", summary_text)
        self.assertIn("control,1,1", summary_text)


def parse_series_matrix_from_fixture():
    with tempfile.TemporaryDirectory() as tmp_dir:
        path = Path(tmp_dir) / "series.txt.gz"
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(SERIES_FIXTURE)
        return parse_series_matrix(path)


if __name__ == "__main__":
    unittest.main()

