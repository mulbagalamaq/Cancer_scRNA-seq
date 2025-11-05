"""
Utility functions for single-cell RNA-seq data preprocessing
"""

import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

def load_scRNA_data(file_path, file_format='h5'):
    """
    Load single-cell RNA-seq data from various formats
    
    Parameters:
    -----------
    file_path : str
        Path to the data file
    file_format : str
        Format of the file ('h5', 'h5ad', 'csv', 'mtx')
    
    Returns:
    --------
    adata : AnnData object
        Annotated data object containing gene expression matrix
    """
    if file_format == 'h5ad':
        adata = sc.read_h5ad(file_path)
    elif file_format == 'h5':
        adata = sc.read_10x_h5(file_path)
    elif file_format == 'csv':
        adata = sc.read_csv(file_path).T
    elif file_format == 'mtx':
        adata = sc.read_mtx(file_path).T
    else:
        raise ValueError(f"Unsupported file format: {file_format}")
    
    return adata


def quality_control(adata, min_genes=200, min_cells=3, max_genes=5000, 
                   max_mito_pct=20, max_ribo_pct=50):
    """
    Perform quality control filtering on scRNA-seq data
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    min_genes : int
        Minimum number of genes detected per cell
    min_cells : int
        Minimum number of cells expressing a gene
    max_genes : int
        Maximum number of genes per cell (filter doublets)
    max_mito_pct : float
        Maximum percentage of mitochondrial genes
    max_ribo_pct : float
        Maximum percentage of ribosomal genes
    
    Returns:
    --------
    adata : AnnData
        Filtered annotated data object
    """
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    sc.pp.calculate_qc_metrics(adata, percent_top=None, 
                               log1p=False, inplace=True)
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'],
                               percent_top=None, log1p=False, inplace=True)
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter based on QC metrics
    adata = adata[adata.obs['n_genes_by_counts'] < max_genes, :]
    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_mito_pct, :]
    if 'pct_counts_ribo' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_ribo'] < max_ribo_pct, :]
    
    return adata


def normalize_and_log_transform(adata, target_sum=1e4):
    """
    Normalize and log-transform scRNA-seq data
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target_sum : float
        Target sum for normalization (default 10,000)
    
    Returns:
    --------
    adata : AnnData
        Normalized and log-transformed data
    """
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=target_sum)
    
    # Log transform
    sc.pp.log1p(adata)
    
    return adata


def find_highly_variable_genes(adata, n_top_genes=2000, flavor='seurat'):
    """
    Identify highly variable genes for downstream analysis
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    n_top_genes : int
        Number of top variable genes to select
    flavor : str
        Method to use ('seurat', 'cell_ranger', 'seurat_v3')
    
    Returns:
    --------
    adata : AnnData
        Data with highly variable genes marked
    """
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, 
                                flavor=flavor, subset=True)
    return adata


def scale_and_pca(adata, n_comps=50):
    """
    Scale data and perform PCA
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    n_comps : int
        Number of principal components
    
    Returns:
    --------
    adata : AnnData
        Data with PCA results
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver='arpack')
    return adata


def compute_neighbors_and_umap(adata, n_neighbors=15, n_pcs=40, 
                               min_dist=0.5, spread=1.0):
    """
    Compute neighborhood graph and UMAP embedding
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    n_neighbors : int
        Number of neighbors for graph construction
    n_pcs : int
        Number of PCs to use
    min_dist : float
        Minimum distance for UMAP
    spread : float
        Spread parameter for UMAP
    
    Returns:
    --------
    adata : AnnData
        Data with neighbors and UMAP computed
    """
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)
    return adata


def leiden_clustering(adata, resolution=0.5, key_added='leiden'):
    """
    Perform Leiden clustering
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    resolution : float
        Resolution parameter for clustering
    key_added : str
        Key to store cluster labels
    
    Returns:
    --------
    adata : AnnData
        Data with cluster labels
    """
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    return adata


def find_marker_genes(adata, groupby='leiden', n_genes=10):
    """
    Find marker genes for each cluster
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    groupby : str
        Column to group by for marker identification
    n_genes : int
        Number of top marker genes per cluster
    
    Returns:
    --------
    markers_df : DataFrame
        DataFrame with marker genes for each cluster
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, 
                           method='wilcoxon', n_genes=n_genes)
    
    # Convert to DataFrame
    markers_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    return markers_df

