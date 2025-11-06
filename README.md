# Single-Cell RNA-Seq Analysis: Exploring FDA-Approved Drug Targets in Cancer Cell Lines

## Research Question

How can single-cell RNA sequencing (scRNA-seq) data from diverse cancer cell lines be leveraged to identify potential therapeutic applications of FDA-approved antibody therapies (Trastuzumab and Bevacizumab) in cancer types beyond their current clinical indications?

## Problem Statement

Cancer treatment faces significant challenges due to tumor heterogeneity and the limited applicability of targeted therapies. Traditional bulk RNA sequencing masks cellular diversity within tumors, potentially overlooking subpopulations that could benefit from existing therapies. This analysis employs single-cell RNA sequencing to:

1. **Identify target expression patterns** at the single-cell level across different cancer cell lines
2. **Characterize cellular heterogeneity** within and between cancer cell lines
3. **Assess therapeutic potential** of FDA-approved drugs (Trastuzumab and Bevacizumab) in additional cancer types
4. **Explore cellular trajectories** and state transitions that may influence drug response

## Methodology

### Analysis Pipeline

1. **Data Acquisition**

   - Download scRNA-seq data from public repositories (GEO, SRA)
   - Quality control and filtering of low-quality cells and genes
   - Normalization and batch effect correction
2. **Dimensionality Reduction**

   - Principal Component Analysis (PCA)
   - Uniform Manifold Approximation and Projection (UMAP)
   - t-distributed Stochastic Neighbor Embedding (t-SNE)
3. **Cell Clustering and Annotation**

   - Leiden clustering for cell population identification
   - Marker gene identification for cell type annotation
   - Integration of multiple datasets for comparative analysis
4. **Target Expression Analysis**

   - HER2 (ERBB2) expression profiling for Trastuzumab assessment
   - VEGF-A expression and angiogenesis pathway analysis for Bevacizumab assessment
   - Expression heterogeneity across cell lines and subpopulations
5. **Pseudotime and Trajectory Analysis**

   - Reconstruction of cellular state transitions
   - Identification of developmental trajectories
   - Correlation of trajectory states with drug target expression
6. **Differential Expression and Pathway Enrichment**

   - Comparative analysis between high/low target-expressing populations
   - Gene Ontology and KEGG pathway enrichment
   - Identification of co-expressed therapeutic targets

## Key Features

- **Comprehensive single-cell analysis** using Scanpy workflows with quality control and preprocessing
- **Trajectory and pseudotime analysis** to understand cellular state dynamics and transitions
- **Multi-cancer cell line comparison** for broader therapeutic insights
- **Reproducible pipeline** with documented code and dependencies

## `Project Structure

```
.
├── README.md
├── requirements.txt
├── notebooks/
│   └── scRNA_drug_target_analysis.ipynb  # Main analysis notebook
├── scripts/
│   ├── data_download.sh                    # Data acquisition script
│   ├── preprocessing_utils.py             # Data preprocessing functions
│   ├── advanced_analysis.py               # Unique analytical approaches
│   └── visualization_utils.py             # Visualization functions
└── data/
    ├── raw/                                # Raw sequencing data
    └── processed/                          # Processed analysis-ready data
```

## Setup Instructions

1. **Clone the repository**

   ```bash
   git clone <your-repo-url>
   cd single-cell-RNA_FDA-drugs
   ```
2. **Create virtual environment**

   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. **Install dependencies**

   ```bash
   pip install -r requirements.txt
   ```
4. **Download data**

   ```bash
   chmod +x scripts/data_download.sh
   ./scripts/data_download.sh
   ```
5. **Run analysis**

   - Open `notebooks/scRNA_drug_target_analysis.ipynb` in Jupyter
   - Execute cells sequentially

## Dependencies

- Python 3.8+
- Scanpy 1.9+
- Seurat (via R)
- Standard bioinformatics libraries (pandas, numpy, matplotlib, seaborn)

See `requirements.txt` for complete list.

## Results Overview

Through my analysis of single-cell RNA-seq data from cancer cell lines, I found several interesting patterns that could inform the use of Trastuzumab and Bevacizumab in additional cancer types.

**Target Expression Patterns**

I observed significant heterogeneity in ERBB2 (HER2) and VEGFA expression both within and between different cancer cell lines. Some cell lines showed clear subpopulations with distinct expression levels, while others were more uniform. This variability was particularly interesting because it suggests that bulk RNA-seq might be masking subpopulations that could benefit from targeted therapy.

**Potential New Indications**

My cross-cancer comparison revealed several cancer types with high target expression that aren't currently approved for these drugs. While I haven't validated these findings experimentally, the expression patterns suggest it might be worth investigating. I also noticed that some cell lines had expression levels comparable to or higher than currently approved indications, which was unexpected.

**Combination Therapy Insights**

One of the more interesting findings was the identification of cell populations expressing both targets simultaneously. I developed a synergy scoring approach to quantify this, and found that roughly 15-25% of cells (depending on the cell line) showed high co-expression. This could indicate potential for combination therapy, though this would need functional validation.

**Therapeutic Windows**

I defined therapeutic windows based on expression percentiles to identify optimal expression ranges. Interestingly, very high expression (above the 90th percentile) sometimes correlated with resistance markers, suggesting that extremely high target expression might not always be favorable. The "sweet spot" seemed to be in the moderate expression range for most cell lines I analyzed.

**Regulatory Mechanisms**

I found several transcription factors that correlated with target gene expression, particularly STAT3 and NFKB1 for VEGFA, and MYC and E2F1 for ERBB2. These correlations were statistically significant and might point to upstream regulatory mechanisms that could be exploited therapeutically or might contribute to resistance.

**Treatment Prioritization**

I developed a scoring system that integrates target expression, response markers, and potential resistance indicators. This helped prioritize which cell lines or subpopulations might be most responsive to treatment. However, I should note that these scores are predictive and would need validation with actual drug response data.

## Analytical Approach

### Standard Workflow

1. Quality control and filtering (genes, cells, mitochondrial/ribosomal content)
2. Normalization and log transformation
3. Highly variable gene identification
4. Dimensionality reduction (PCA, UMAP)
5. Clustering (Leiden algorithm)
6. Differential expression analysis
7. Pathway enrichment analysis

### References

### Key Publications

- Kinker, G.S., et al. (2020). Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity. Nature Genetics.
- Haque, A., et al. (2017). A practical guide to single-cell RNA-sequencing for biomedical research. Genome Medicine.
- Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.

### Data Resources

- Cancer Cell Line Encyclopedia (CCLE), Broad Institute
- Gene Expression Omnibus (GEO), NCBI
- Sequence Read Archive (SRA), NCBI
- 10x Genomics Public Datasets

### Drug Information

- Trastuzumab: FDA-approved HER2-targeted therapy for breast and gastric cancers
- Bevacizumab: FDA-approved VEGF-targeted therapy for multiple cancer types

## Computational Methods

### Software Versions

- Python 3.8+
- Scanpy 1.9+
- Key dependencies: pandas, numpy, scipy, matplotlib, seaborn, scikit-learn

### Statistical Methods

- Non-parametric correlation (Spearman) for regulatory network analysis
- Wilcoxon rank-sum test for differential expression
- Z-score normalization for cross-cancer comparisons
- Percentile-based thresholding for therapeutic windows

### Visualization

- UMAP for non-linear dimensionality reduction
- Leiden clustering for cell population identification
- Custom visualization functions for synergy, windows, and comparisons

## License

This project is for educational and research purposes.
