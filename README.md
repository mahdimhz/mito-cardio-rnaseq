# Mitochondrial Disease Omics Analysis for MERRF-Relevant Cardiomyopathy Mechanisms

A focused RNA-seq workflow for a related mitochondrial cardiomyopathy model, built from public human iPSC-cardiomyocyte data and reported with cautious MERRF framing.

## At a Glance

| Item | Value |
|---|---|
| Disease model | Friedreich ataxia / frataxin deficiency |
| Cell type | Human iPSC-derived cardiomyocytes |
| Dataset | `GSE305638` |
| Comparison | FRDA iPSC-cardiomyocytes vs isogenic corrected controls |
| Samples | 18 total: 9 disease, 9 control |
| Design | 3 paired patient/control lines with 3 replicates per condition |
| Data type | Processed RNA-seq count matrix from GEO |
| Analysis | Dataset screening, metadata curation, QC, exploratory paired expression analysis, mitochondrial pathway scoring |
| Scope | Related mitochondrial cardiomyopathy model, not direct MERRF data |

## Why This Dataset

MERRF is commonly associated with the mitochondrial `MT-TK` m.8344A>G variant and impaired mitochondrial translation, oxidative phosphorylation, ATP production, and tissue energy balance. The PhD topic is centered on MERRF cardiomyopathy, but public transcriptomic datasets that are direct MERRF, iPSC-cardiomyocyte based, processed, and clearly paired are limited.

`GSE305638` was selected because it provides a usable cardiac mitochondrial disease model:

| Criterion | `GSE305638` |
|---|---|
| Human cells | Yes |
| iPSC-cardiomyocytes | Yes |
| Mitochondrial disease biology | Yes, FRDA/frataxin deficiency |
| Paired disease/control structure | Yes, patient lines with isogenic corrected controls |
| Processed RNA-seq counts | Yes |
| Direct MERRF / `MT-TK` m.8344A>G | No |

The dataset is not MERRF-specific. It is used as a related mitochondrial cardiomyopathy model, not as evidence specific to MERRF.

The dataset registry was re-checked for direct MERRF alternatives. `GSE106601` is a direct MERRF/m.8344A>G dataset with processed supplementary files, but GEO metadata indicate mixed skeletal tissue and immortalized cell-line samples rather than iPSC-cardiomyocytes. `GSE142745` and related subseries include m.8344A>G context, but they are single-cell mitochondrial genotyping/clonal-variation datasets rather than cardiac iPSC disease models. These datasets are retained as secondary registry context.

## Workflow Overview

| Stage | Output |
|---|---|
| Dataset search and screening | `references/dataset_registry.csv` |
| Source tracking | `references/source_register.csv` |
| Metadata curation | `data/processed/metadata_clean.csv` |
| Expression QC | `data/processed/expression_clean.csv`, QC tables and figures |
| Exploratory paired expression analysis | Ranked gene table, top-gene table, volcano plot |
| Mitochondrial gene review | Mitochondrial gene table and heatmap |
| Pathway scoring | Per-sample mitochondrial category scores and paired score plot |
| Local enrichment | Exploratory mitochondrial-category over-representation table |

## Key Outputs

| Type | Path |
|---|---|
| Dataset registry | `references/dataset_registry.csv` |
| Clean metadata | `data/processed/metadata_clean.csv` |
| Clean expression matrix | `data/processed/expression_clean.csv` |
| Expression QC summary | `reports/tables/expression_qc_summary.csv` |
| Ranked expression results | `reports/tables/differential_expression_results.csv` |
| Mitochondrial gene results | `reports/tables/mitochondrial_gene_results.csv` |
| Mitochondrial pathway scores | `reports/tables/mitochondrial_pathway_scores.csv` |
| Pathway score tests | `reports/tables/mitochondrial_pathway_score_tests.csv` |
| Local enrichment results | `reports/tables/pathway_enrichment_results.csv` |

![Library size by sample](reports/figures/library_size_by_sample.png)

![PCA plot](reports/figures/pca_plot.png)

## Results Summary

Metadata curation recovered 18 samples: 9 FRDA disease samples and 9 isogenic corrected controls across 3 paired patient/control lines.

Expression QC found:

| Metric | Value |
|---|---:|
| Gene rows | 59,743 |
| Aligned samples | 18 |
| Missing count values | 0 |
| Duplicate Ensembl gene IDs | 0 |
| Duplicate gene symbols | 3,715 |
| PCA variance, PC1 | 35.3% |
| PCA variance, PC2 | 9.1% |

The exploratory paired expression analysis tested 15,100 expressed genes after filtering. It found 859 genes with nominal `p < 0.05`, but 0 genes passed FDR `< 0.10` and 0 genes passed FDR `< 0.05`.

![Exploratory paired expression volcano plot](reports/figures/volcano_plot.png)

Ranked exploratory signals included `MEG3`, `CBLN2`, `CNTN6`, and `CHCHD2` among disease-up genes, and `TRH`, `DRD1`, `SFRP5`, and `CACNA1G` among disease-down genes. These are ranked signals from a small paired analysis, not genome-wide findings after FDR correction.

## Mitochondrial Interpretation

The mitochondrial reference list contains 41 genes spanning mtDNA-encoded OXPHOS genes, respiratory-chain complexes, ATP synthase, mitochondrial translation and mtDNA maintenance, fusion/fission, mitophagy, oxidative stress / ROS, and mitochondrial biogenesis.

Gene-level review found:

| Metric | Value |
|---|---:|
| Mitochondrial genes in reference | 41 |
| Genes tested in exploratory expression table | 40 |
| `MT-TK` tested after expression filter | No |
| Lowest mitochondrial adjusted p-value | 0.671 |
| FDR-significant mitochondrial genes | 0 |

Descriptively, mtDNA-encoded OXPHOS genes trended lower in disease samples. The strongest nominal mitochondrial genes included `POLG`, `ATP5F1A`, `MT-ND2`, `MT-CO3`, and `PPARGC1A`.

![Mitochondrial gene heatmap](reports/figures/mitochondrial_gene_heatmap.png)

Sample-level pathway scores were computed as mean z-scored expression of available genes in each mitochondrial category. All 11 categories had lower mean disease scores than paired corrected controls. The largest mean differences were seen for mtDNA-encoded OXPHOS genes, mitochondrial biogenesis, and respiratory-chain categories. No mitochondrial pathway-score comparison is interpreted as a formal disease mechanism, and no score-level result is treated as confirmatory given the three-pair design.

![Mitochondrial pathway scores](reports/figures/mitochondrial_pathway_scores_paired.png)

Local mitochondrial-category enrichment used nominal `p < 0.05` genes against the tested-gene background. No category passed FDR correction. The strongest local results were ATP synthase, mtDNA-encoded OXPHOS genes, and mitochondrial translation / mtDNA maintenance, but these remain exploratory.

![Pathway enrichment dotplot](reports/figures/pathway_enrichment_dotplot.png)

## Limitations

`GSE305638` is an FRDA/frataxin-deficiency model, not a MERRF patient dataset and not an `MT-TK` m.8344A>G cardiomyocyte dataset.

The effective paired sample size is 3 patient/control lines. Replicates were averaged within each patient/group before paired testing. This preserves the patient-pair structure, but it does not replace a larger biological cohort.

The paired expression analysis uses log2 CPM values and paired t-tests in Python. It is not a DESeq2, edgeR, or limma/voom model. The results are therefore reported as exploratory.

The mitochondrial pathway scores and local enrichment tables summarize selected gene categories. They do not establish disease mechanism, treatment response, or diagnostic claims specific to MERRF.

## Data and References

| Resource | Role in this project |
|---|---|
| [GSE305638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305638) | Primary iPSC-cardiomyocyte RNA-seq dataset used for analysis |
| [GSE106601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106601) | Direct MERRF-related registry candidate; not selected as the primary dataset |
| [GSE142745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142745) | m.8344A>G-related registry context; not selected as the primary dataset |
| [Frataxin deficiency drives cardiac dysfunction and transcriptional dysregulation in Friedreich ataxia iPSC model](https://www.biorxiv.org/content/10.1101/2025.08.20.671405v1) | Preprint associated with the selected FRDA iPSC-cardiomyocyte dataset |
| [tRNA modification landscape selectively controls mitochondrial translation efficiency in MERRF](https://pubmed.ncbi.nlm.nih.gov/30262910/) | MERRF / m.8344A>G reference associated with the registry review |

## Reproduce the Analysis

Create the environment:

```bash
conda env create -f environment.yml
conda activate merrf-omics
```

Run tests:

```bash
python -m unittest discover -s tests -v
```

Run the workflow:

```bash
python src/data/search_ncbi_geo.py
python src/data/screen_dataset_registry.py
python src/data/download_selected_dataset.py
python src/data/metadata_qc.py
python src/data/expression_qc.py
python src/analysis/differential_expression.py
python src/analysis/mitochondrial_gene_set_analysis.py
python src/analysis/mitochondrial_pathway_scores.py
python src/analysis/pathway_enrichment.py
python src/analysis/results_summary.py
```

The repository uses lightweight GEO files and processed count tables. It does not require downloading raw FASTQ files.

## Repository Structure

~~~text
merrf-mitochondrial-disease-omics-analysis/
|-- data/
|   |-- raw/                 # GEO series matrix and processed count table
|   |-- interim/             # GEO/NCBI metadata files
|   +-- processed/           # cleaned metadata and expression matrices
|-- notebooks/               # notebooks 01-06 for workflow stages
|-- references/              # dataset registry, source register, gene sets
|-- reports/
|   |-- figures/             # QC, PCA, expression, mitochondrial, enrichment figures
|   +-- tables/              # QC, expression, mitochondrial, enrichment, summary tables
|-- src/
|   |-- data/                # dataset search, download, metadata, QC scripts
|   +-- analysis/            # expression, mitochondrial, enrichment, summary scripts
|-- tests/                   # unit tests for workflow logic
|-- environment.yml
|-- requirements.txt
+-- README.md
~~~

## CV Bullet

Built a mitochondrial-disease omics workflow using public human iPSC-cardiomyocyte RNA-seq data, including GEO dataset screening, metadata curation, RNA-seq QC, exploratory paired expression analysis, mitochondrial pathway scoring, local enrichment, and cautious reporting for MERRF-relevant cardiomyopathy mechanisms.
