from pathlib import Path
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ReadmeContentTests(unittest.TestCase):
    def test_readme_contains_required_publication_sections(self):
        text = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        required_sections = [
            "# Frataxin Deficiency iPSC-Cardiomyocyte RNA-seq Analysis",
            "## At a Glance",
            "## Dataset",
            "## Workflow",
            "## Key Outputs",
            "## Results Summary",
            "## Mitochondrial Gene and Pathway Analysis",
            "## Limitations",
            "## Data and References",
            "## Reproduce the Analysis",
            "## Repository Structure",
        ]
        for section in required_sections:
            self.assertIn(section, text)

    def test_readme_uses_frda_framing(self):
        text = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")

        self.assertIn("Friedreich ataxia / frataxin deficiency", text)
        self.assertIn("human iPSC-derived cardiomyocytes", text)
        self.assertIn("isogenic corrected controls", text)
        self.assertIn("GSE106601", text)
        self.assertIn("0 genes passed FDR", text)
        self.assertNotIn("MERRF-" + "Relevant", text)
        self.assertNotIn("P" + "hD", text)
        self.assertNotIn("Human" + "itas", text)
        self.assertNotIn("CV " + "Bullet", text)
        self.assertNotIn("proves", text.lower())
        self.assertNotIn("cures", text.lower())
        self.assertNotIn("hon" + "est", text.lower())
        self.assertNotIn("over" + "claimed", text.lower())


if __name__ == "__main__":
    unittest.main()
