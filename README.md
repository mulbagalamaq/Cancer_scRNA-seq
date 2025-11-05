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

## Unique Contributions

This analysis extends beyond standard single-cell RNA-seq workflows with several novel analytical approaches:

### 1. Drug Synergy Scoring
- Calculates co-expression-based synergy scores for combination therapy
- Identifies cell populations that express both therapeutic targets simultaneously
- Uses geometric mean to emphasize dual-target expression patterns
- Enables identification of optimal candidates for combination treatment strategies

### 2. Therapeutic Window Analysis
- Defines optimal expression ranges for drug targeting
- Distinguishes between responsive and potentially resistant populations
- Identifies expression ranges where drugs are most likely to be effective
- Helps prioritize treatment decisions based on target expression levels

### 3. Cross-Cancer Type Comparison
- Systematic comparison of target expression across diverse cancer cell lines
- Identifies potential new indications beyond current FDA approvals
- Statistical validation of expression differences between cancer types
- Enables evidence-based expansion of therapeutic applications

### 4. Regulatory Network Analysis
- Identifies transcription factors that control drug target expression
- Reveals upstream regulatory mechanisms
- Highlights potential resistance mechanisms through correlation analysis
- Suggests combination therapy targets by identifying co-regulatory factors

### 5. Response Prediction Scoring
- Integrates multiple factors (target expression, response markers, resistance indicators)
- Calculates treatment priority scores for clinical decision support
- Combines expression data with known response-associated markers
- Provides quantitative framework for treatment prioritization

## Project Structure

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

The comprehensive analysis pipeline generates multiple layers of insight:

1. **Target Expression Profiling**
   - Heterogeneous expression of HER2 (ERBB2) and VEGF-A across cancer cell lines
   - Cell-type and subpopulation-specific expression patterns
   - Expression variability within and between cell lines

2. **Therapeutic Potential Assessment**
   - Identification of cancer types with high target expression beyond current indications
   - Therapeutic window definitions for optimal treatment targeting
   - Drug synergy scores for combination therapy candidates

3. **Mechanistic Insights**
   - Regulatory networks controlling target gene expression
   - Potential resistance mechanisms through co-expression analysis
   - Cellular state trajectories correlated with drug target expression

4. **Treatment Prioritization**
   - Response prediction scores integrating multiple factors
   - Statistical validation of potential new indications
   - Priority scoring for clinical decision support

## Analytical Approach

### Standard Workflow
1. Quality control and filtering (genes, cells, mitochondrial/ribosomal content)
2. Normalization and log transformation
3. Highly variable gene identification
4. Dimensionality reduction (PCA, UMAP)
5. Clustering (Leiden algorithm)
6. Differential expression analysis
7. Pathway enrichment analysis

### Unique Extensions
1. **Synergy Analysis**: Geometric mean-based co-expression scoring
2. **Window Analysis**: Percentile-based therapeutic range identification
3. **Cross-Comparison**: Statistical testing across cancer types
4. **Network Analysis**: Transcription factor correlation mapping
5. **Priority Scoring**: Multi-factor integration for treatment decisions

## References

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

## Author

[Your Name]

## License

This project is for educational and research purposes.

