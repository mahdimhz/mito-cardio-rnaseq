from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ReadmeContentTests(unittest.TestCase):
    def test_readme_contains_required_publication_sections(self):
        text = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        required_sections = [
            "# Mitochondrial Disease Omics Analysis for MERRF-Relevant Cardiomyopathy Mechanisms",
            "## At a Glance",
            "## Why This Dataset",
            "## Workflow Overview",
            "## Key Outputs",
            "## Results Summary",
            "## Mitochondrial Interpretation",
            "## Limitations",
            "## Reproduce the Analysis",
            "## Repository Structure",
            "## CV Bullet",
        ]
        for section in required_sections:
            self.assertIn(section, text)

    def test_readme_uses_cautious_merrf_language(self):
        text = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")

        self.assertIn("not MERRF-specific", text)
        self.assertIn("related mitochondrial cardiomyopathy model", text)
        self.assertIn("not as evidence specific to MERRF", text)
        self.assertIn("GSE106601", text)
        self.assertIn("mixed skeletal tissue and immortalized cell-line samples", text)
        self.assertIn("0 genes passed FDR", text)
        self.assertNotIn("proves", text.lower())
        self.assertNotIn("cures", text.lower())
        self.assertNotIn("hon" + "est", text.lower())
        self.assertNotIn("over" + "claimed", text.lower())


if __name__ == "__main__":
    unittest.main()
