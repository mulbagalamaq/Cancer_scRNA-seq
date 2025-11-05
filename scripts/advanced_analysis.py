"""
Advanced analysis functions for single-cell RNA-seq drug target analysis
This module contains unique analytical approaches for therapeutic assessment
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')


def calculate_response_score(adata, target_gene, response_markers=None):
    """
    Calculate a predictive response score based on target expression and
    response-associated markers
    
    Unique Feature: Combines target expression with known response markers
    to predict treatment responsiveness
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target_gene : str
        Target gene for drug therapy
    response_markers : list
        List of genes associated with treatment response (optional)
    
    Returns:
    --------
    adata : AnnData
        Data with response scores added
    """
    if response_markers is None:
        # Default response markers (can be customized)
        response_markers = ['CD274', 'PDCD1', 'CTLA4', 'IFNG', 'GZMB']
    
    # Get target expression
    if target_gene not in adata.var_names:
        matching = [g for g in adata.var_names if g.upper() == target_gene.upper()]
        if matching:
            target_gene = matching[0]
        else:
            return adata
    
    target_expr = adata[:, target_gene].X.toarray().flatten() if hasattr(
        adata[:, target_gene].X, 'toarray') else adata[:, target_gene].X.flatten()
    
    # Calculate response marker scores
    marker_scores = []
    for marker in response_markers:
        if marker in adata.var_names:
            marker_expr = adata[:, marker].X.toarray().flatten() if hasattr(
                adata[:, marker].X, 'toarray') else adata[:, marker].X.flatten()
            marker_scores.append(marker_expr)
    
    if marker_scores:
        marker_array = np.array(marker_scores).T
        marker_score = np.mean(marker_array, axis=1)
    else:
        marker_score = np.zeros(len(target_expr))
    
    # Normalize scores
    target_norm = (target_expr - target_expr.min()) / (target_expr.max() - target_expr.min() + 1e-10)
    marker_norm = (marker_score - marker_score.min()) / (marker_score.max() - marker_score.min() + 1e-10)
    
    # Combined response score (weighted combination)
    response_score = 0.6 * target_norm + 0.4 * marker_norm
    
    adata.obs[f'{target_gene}_response_score'] = response_score
    adata.obs[f'{target_gene}_response_category'] = pd.cut(
        response_score,
        bins=[0, 0.33, 0.66, 1.0],
        labels=['Low Response', 'Moderate Response', 'High Response']
    )
    
    return adata


def identify_resistance_mechanisms(adata, target_gene, resistance_genes=None):
    """
    Identify potential resistance mechanisms based on co-expression patterns
    
    Unique Feature: Identifies genes that might confer resistance when
    co-expressed with drug targets
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target_gene : str
        Target gene for drug therapy
    resistance_genes : list
        List of genes potentially associated with resistance (optional)
    
    Returns:
    --------
    resistance_df : DataFrame
        DataFrame with resistance gene associations
    """
    if resistance_genes is None:
        # Common resistance-associated genes
        resistance_genes = ['ABCB1', 'ABCG2', 'TP53', 'MDM2', 'BCL2', 
                          'BCL2L1', 'MCL1', 'AKT1', 'PIK3CA', 'PTEN']
    
    if target_gene not in adata.var_names:
        matching = [g for g in adata.var_names if g.upper() == target_gene.upper()]
        if matching:
            target_gene = matching[0]
        else:
            return None
    
    target_expr = adata[:, target_gene].X.toarray().flatten() if hasattr(
        adata[:, target_gene].X, 'toarray') else adata[:, target_gene].X.flatten()
    
    # Calculate correlations with resistance genes
    resistance_results = []
    for gene in resistance_genes:
        if gene in adata.var_names:
            gene_expr = adata[:, gene].X.toarray().flatten() if hasattr(
                adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
            
            corr, pval = stats.spearmanr(target_expr, gene_expr)
            
            resistance_results.append({
                'resistance_gene': gene,
                'correlation': corr,
                'p_value': pval,
                'abs_correlation': abs(corr),
                'significant': pval < 0.05
            })
    
    resistance_df = pd.DataFrame(resistance_results)
    resistance_df = resistance_df.sort_values('abs_correlation', ascending=False)
    
    return resistance_df


def calculate_treatment_priority_score(adata, target1='ERBB2', target2='VEGFA'):
    """
    Calculate a priority score for treatment decisions
    
    Unique Feature: Integrates multiple factors to prioritize cells/cell lines
    for treatment:
    - Target expression levels
    - Synergy potential
    - Response markers
    - Resistance indicators
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target1 : str
        First target gene
    target2 : str
        Second target gene
    
    Returns:
    --------
    adata : AnnData
        Data with priority scores added
    """
    # Get target expressions
    targets = [target1, target2]
    target_exprs = []
    
    for target in targets:
        if target in adata.var_names:
            expr = adata[:, target].X.toarray().flatten() if hasattr(
                adata[:, target].X, 'toarray') else adata[:, target].X.flatten()
            expr_norm = (expr - expr.min()) / (expr.max() - expr.min() + 1e-10)
            target_exprs.append(expr_norm)
        else:
            target_exprs.append(np.zeros(adata.n_obs))
    
    # Average target expression
    avg_target = np.mean(target_exprs, axis=0)
    
    # Synergy component (if calculated)
    if 'synergy_score' in adata.obs.columns:
        synergy = adata.obs['synergy_score'].values
    else:
        synergy = np.sqrt(target_exprs[0] * target_exprs[1]) if len(target_exprs) == 2 else avg_target
    
    # Response component (if calculated)
    response_component = np.zeros(adata.n_obs)
    for target in targets:
        col = f'{target}_response_score'
        if col in adata.obs.columns:
            response_component += adata.obs[col].values
    if np.sum(response_component) > 0:
        response_component = response_component / np.max(response_component)
    
    # Combined priority score
    priority_score = (0.4 * avg_target + 
                     0.3 * synergy + 
                     0.3 * response_component)
    
    adata.obs['treatment_priority_score'] = priority_score
    adata.obs['treatment_priority'] = pd.cut(
        priority_score,
        bins=[0, 0.33, 0.66, 1.0],
        labels=['Low Priority', 'Medium Priority', 'High Priority']
    )
    
    return adata


def perform_statistical_validation(comparison_df, current_indications=None):
    """
    Perform statistical validation of potential new indications
    
    Unique Feature: Statistical testing to identify significantly different
    expression patterns between current and potential indications
    
    Parameters:
    -----------
    comparison_df : DataFrame
        DataFrame from cross-cancer comparison
    current_indications : list
        List of cancer types with current FDA approval
    
    Returns:
    --------
    validated_df : DataFrame
        DataFrame with statistical validation results
    """
    if current_indications is None:
        current_indications = []
    
    validated_df = comparison_df.copy()
    
    # Calculate statistics for current indications
    if len(current_indications) > 0:
        current_mask = validated_df['cell_line'].isin(current_indications)
        current_mean = validated_df.loc[current_mask, 'mean_expression'].mean()
        current_std = validated_df.loc[current_mask, 'mean_expression'].std()
    else:
        current_mean = validated_df['mean_expression'].median()
        current_std = validated_df['mean_expression'].std()
    
    # Z-score for each cell line
    validated_df['z_score'] = (validated_df['mean_expression'] - current_mean) / (current_std + 1e-10)
    
    # P-value (two-tailed test)
    validated_df['p_value'] = 2 * (1 - stats.norm.cdf(abs(validated_df['z_score'])))
    
    # Significance
    validated_df['significant'] = validated_df['p_value'] < 0.05
    validated_df['higher_than_current'] = validated_df['mean_expression'] > current_mean
    
    # Potential new indication if significantly higher or similar
    validated_df['potential_indication'] = (
        (validated_df['significant'] & validated_df['higher_than_current']) |
        (validated_df['mean_expression'] >= current_mean * 0.8)
    )
    
    return validated_df

