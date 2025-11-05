"""
Visualization utilities for single-cell RNA-seq drug target analysis
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')


def plot_drug_target_comparison(adata, target1='ERBB2', target2='VEGFA', 
                                save_path=None):
    """
    Create comprehensive comparison plots for both drug targets
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target1 : str
        First target gene
    target2 : str
        Second target gene
    save_path : str
        Path to save figure (optional)
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # UMAP plots for each target
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=target1, ax=axes[0, 0], show=False,
                   title=f'{target1} Expression (Trastuzumab Target)')
        sc.pl.umap(adata, color=target2, ax=axes[0, 1], show=False,
                   title=f'{target2} Expression (Bevacizumab Target)')
    
    # Expression distribution
    if target1 in adata.var_names:
        expr1 = adata[:, target1].X.toarray().flatten() if hasattr(
            adata[:, target1].X, 'toarray') else adata[:, target1].X.flatten()
        axes[1, 0].hist(expr1, bins=50, alpha=0.7, edgecolor='black')
        axes[1, 0].set_xlabel('Expression Level')
        axes[1, 0].set_ylabel('Number of Cells')
        axes[1, 0].set_title(f'{target1} Expression Distribution')
        axes[1, 0].axvline(np.percentile(expr1, 75), color='r', 
                          linestyle='--', label='75th percentile')
        axes[1, 0].legend()
    
    if target2 in adata.var_names:
        expr2 = adata[:, target2].X.toarray().flatten() if hasattr(
            adata[:, target2].X, 'toarray') else adata[:, target2].X.flatten()
        axes[1, 1].hist(expr2, bins=50, alpha=0.7, edgecolor='black', color='orange')
        axes[1, 1].set_xlabel('Expression Level')
        axes[1, 1].set_ylabel('Number of Cells')
        axes[1, 1].set_title(f'{target2} Expression Distribution')
        axes[1, 1].axvline(np.percentile(expr2, 75), color='r', 
                          linestyle='--', label='75th percentile')
        axes[1, 1].legend()
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_synergy_analysis(adata, save_path=None):
    """
    Visualize drug synergy scoring results
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with synergy scores
    save_path : str
        Path to save figure (optional)
    """
    if 'synergy_score' not in adata.obs.columns:
        print("Error: Synergy scores not calculated. Run calculate_drug_synergy_score first.")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Synergy score distribution
    axes[0, 0].hist(adata.obs['synergy_score'], bins=50, alpha=0.7, 
                    edgecolor='black', color='purple')
    axes[0, 0].set_xlabel('Synergy Score')
    axes[0, 0].set_ylabel('Number of Cells')
    axes[0, 0].set_title('Drug Synergy Score Distribution')
    axes[0, 0].axvline(adata.obs['synergy_score'].quantile(0.75), 
                      color='r', linestyle='--', label='High synergy threshold')
    axes[0, 0].legend()
    
    # Synergy on UMAP
    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color='synergy_score', ax=axes[0, 1], show=False,
                   title='Synergy Score on UMAP', cmap='viridis')
    
    # Synergy categories
    if 'synergy_category' in adata.obs.columns:
        category_counts = adata.obs['synergy_category'].value_counts()
        axes[1, 0].bar(category_counts.index.astype(str), category_counts.values,
                      color=['blue', 'green', 'orange', 'red'])
        axes[1, 0].set_xlabel('Synergy Category')
        axes[1, 0].set_ylabel('Number of Cells')
        axes[1, 0].set_title('Cell Distribution by Synergy Category')
        axes[1, 0].tick_params(axis='x', rotation=45)
    
    # Synergy by cluster
    if 'leiden' in adata.obs.columns and 'synergy_score' in adata.obs.columns:
        cluster_synergy = adata.obs.groupby('leiden')['synergy_score'].mean()
        axes[1, 1].bar(range(len(cluster_synergy)), cluster_synergy.values,
                      color='teal')
        axes[1, 1].set_xlabel('Cluster')
        axes[1, 1].set_ylabel('Mean Synergy Score')
        axes[1, 1].set_title('Average Synergy Score by Cluster')
        axes[1, 1].set_xticks(range(len(cluster_synergy)))
        axes[1, 1].set_xticklabels(cluster_synergy.index)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_therapeutic_windows(window_stats_list, save_path=None):
    """
    Visualize therapeutic windows for multiple drugs
    
    Parameters:
    -----------
    window_stats_list : list
        List of window statistics dictionaries
    save_path : str
        Path to save figure (optional)
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for i, stats in enumerate(window_stats_list):
        drug = stats['drug']
        window_low = stats['window_low']
        window_high = stats['window_high']
        high_threshold = stats['high_expr_threshold']
        
        # Draw therapeutic window
        rect = Rectangle((i-0.3, window_low), 0.6, window_high - window_low,
                        facecolor=colors[i % len(colors)], alpha=0.3,
                        edgecolor='black', linewidth=2, label=f'{drug} Window')
        ax.add_patch(rect)
        
        # Mark high expression threshold
        ax.axhline(y=high_threshold, xmin=i-0.3, xmax=i+0.3,
                  color='red', linestyle='--', linewidth=2,
                  label=f'{drug} High Expr Threshold' if i == 0 else '')
    
    ax.set_xlabel('Drug')
    ax.set_ylabel('Expression Level')
    ax.set_title('Therapeutic Windows for Drug Targets')
    ax.set_xticks(range(len(window_stats_list)))
    ax.set_xticklabels([s['drug'] for s in window_stats_list])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_cross_cancer_comparison(comparison_df, target_gene, save_path=None):
    """
    Create bar plot comparing target expression across cancer types
    
    Parameters:
    -----------
    comparison_df : DataFrame
        DataFrame from cross-cancer comparison
    target_gene : str
        Target gene name
    save_path : str
        Path to save figure (optional)
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Sort by mean expression
    df_sorted = comparison_df.sort_values('mean_expression', ascending=True)
    
    # Create bar plot
    bars = ax.barh(range(len(df_sorted)), df_sorted['mean_expression'],
                   color='steelblue', alpha=0.7, edgecolor='black')
    
    # Highlight high expression
    high_threshold = df_sorted['mean_expression'].quantile(0.75)
    for i, (idx, row) in enumerate(df_sorted.iterrows()):
        if row['mean_expression'] > high_threshold:
            bars[i].set_color('coral')
    
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted['cell_line'])
    ax.set_xlabel(f'{target_gene} Mean Expression')
    ax.set_title(f'Target Expression Across Cancer Cell Lines: {target_gene}')
    ax.axvline(high_threshold, color='red', linestyle='--', 
              label=f'75th percentile ({high_threshold:.2f})')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


def plot_regulatory_network(regulatory_df, top_n=10, save_path=None):
    """
    Visualize regulatory network analysis results
    
    Parameters:
    -----------
    regulatory_df : DataFrame
        DataFrame from regulatory network analysis
    top_n : int
        Number of top regulators to display
    save_path : str
        Path to save figure (optional)
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Top regulators by correlation
    top_regulators = regulatory_df.nlargest(top_n, 'abs_correlation')
    
    # Bar plot of correlations
    colors = ['green' if x > 0 else 'red' for x in top_regulators['correlation']]
    axes[0].barh(range(len(top_regulators)), top_regulators['correlation'],
                color=colors, alpha=0.7, edgecolor='black')
    axes[0].set_yticks(range(len(top_regulators)))
    axes[0].set_yticklabels(top_regulators['transcription_factor'])
    axes[0].set_xlabel('Correlation Coefficient')
    axes[0].set_title(f'Top {top_n} Regulatory Factors')
    axes[0].axvline(0, color='black', linestyle='-', linewidth=0.5)
    axes[0].grid(True, alpha=0.3, axis='x')
    
    # Significance plot
    significant = regulatory_df[regulatory_df['significant']]
    axes[1].scatter(significant['correlation'], -np.log10(significant['p_value'] + 1e-10),
                   s=100, alpha=0.6, c=significant['abs_correlation'], cmap='viridis')
    axes[1].set_xlabel('Correlation Coefficient')
    axes[1].set_ylabel('-log10(p-value)')
    axes[1].set_title('Significance vs Correlation')
    axes[1].axhline(-np.log10(0.05), color='red', linestyle='--', 
                   label='p = 0.05 threshold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

