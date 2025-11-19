"""
Visualization Functions

Plotting functions for SpatioResist results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_resistance_scores(scores, sample_names=None, figsize=(10, 6)):
    """
    Plot resistance risk scores for multiple samples.

    Parameters
    ----------
    scores : array-like
        Resistance risk scores (0-1)
    sample_names : list, optional
        Sample names for labels
    figsize : tuple
        Figure size

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    scores = np.array(scores)

    if sample_names is None:
        sample_names = [f'Sample_{i}' for i in range(len(scores))]

    # Sort by score
    sort_idx = np.argsort(scores)
    scores_sorted = scores[sort_idx]
    names_sorted = [sample_names[i] for i in sort_idx]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Color by risk level
    colors = ['#2ecc71' if s < 0.3 else '#f39c12' if s < 0.7 else '#e74c3c' for s in scores_sorted]

    # Bar plot
    bars = ax.barh(range(len(scores_sorted)), scores_sorted, color=colors)

    # Formatting
    ax.set_yticks(range(len(names_sorted)))
    ax.set_yticklabels(names_sorted)
    ax.set_xlabel('Resistance Risk Score')
    ax.set_title('SpatioResist Risk Scores')
    ax.axvline(0.5, color='black', linestyle='--', alpha=0.5)

    # Add value labels
    for i, (bar, score) in enumerate(zip(bars, scores_sorted)):
        ax.text(score + 0.02, i, f'{score:.2f}', va='center')

    ax.set_xlim(0, 1.1)
    plt.tight_layout()

    return fig


def plot_spatial_risk(adata, risk_scores=None, cmap='RdYlBu_r', figsize=(10, 8)):
    """
    Plot spatial map of resistance risk.

    Parameters
    ----------
    adata : AnnData
        Spatial AnnData object
    risk_scores : array-like, optional
        Risk scores per spot. If None, looks for 'resistance_risk' in obs.
    cmap : str
        Colormap
    figsize : tuple
        Figure size

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Get risk scores
    if risk_scores is not None:
        scores = risk_scores
    elif 'resistance_risk' in adata.obs.columns:
        scores = adata.obs['resistance_risk'].values
    else:
        raise ValueError("No risk scores provided. Pass risk_scores or add 'resistance_risk' to adata.obs")

    # Get spatial coordinates
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")

    # Plot
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=scores,
        cmap=cmap,
        s=20,
        alpha=0.8
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Resistance Risk Score')

    ax.set_xlabel('Spatial X')
    ax.set_ylabel('Spatial Y')
    ax.set_title('Spatial Resistance Risk Map')
    ax.set_aspect('equal')

    plt.tight_layout()

    return fig


def plot_feature_importance(importance, top_n=20, figsize=(10, 8)):
    """
    Plot feature importance from SpatioResist model.

    Parameters
    ----------
    importance : pd.Series
        Feature importance values
    top_n : int
        Number of top features to show
    figsize : tuple
        Figure size

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Get top features
    top_features = importance.head(top_n)

    # Plot
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_features)))
    bars = ax.barh(range(len(top_features)), top_features.values, color=colors)

    # Formatting
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels(top_features.index)
    ax.set_xlabel('Feature Importance')
    ax.set_title(f'Top {top_n} Spatial Features for Resistance Prediction')

    plt.tight_layout()

    return fig


def plot_niche_composition(adata, niche_col='spatial_niche', figsize=(12, 8)):
    """
    Plot cell type composition of spatial niches.

    Parameters
    ----------
    adata : AnnData
        Spatial AnnData with abundance columns
    niche_col : str
        Column name for niche assignments
    figsize : tuple
        Figure size

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    # Get abundance columns
    abundance_cols = [c for c in adata.obs.columns if c.startswith('abundance_')]

    if not abundance_cols:
        raise ValueError("No abundance columns found in adata.obs")

    if niche_col not in adata.obs.columns:
        raise ValueError(f"Niche column '{niche_col}' not found")

    # Calculate mean abundance per niche
    abundance_df = adata.obs[abundance_cols + [niche_col]].copy()
    niche_composition = abundance_df.groupby(niche_col)[abundance_cols].mean()

    # Clean column names
    niche_composition.columns = [c.replace('abundance_', '') for c in niche_composition.columns]

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        niche_composition.T,
        cmap='YlOrRd',
        annot=True,
        fmt='.1f',
        ax=ax
    )

    ax.set_title('Cell Type Composition per Spatial Niche')
    ax.set_xlabel('Niche')
    ax.set_ylabel('Cell Type')

    plt.tight_layout()

    return fig


def plot_survival_km(df, score_col, time_col='OS_time', event_col='OS_event', figsize=(8, 6)):
    """
    Plot Kaplan-Meier survival curves stratified by score.

    Parameters
    ----------
    df : pd.DataFrame
        Clinical data with scores and survival info
    score_col : str
        Column name for the score
    time_col : str
        Column name for survival time
    event_col : str
        Column name for event indicator
    figsize : tuple
        Figure size

    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    fig, ax = plt.subplots(figsize=figsize)

    # Stratify by median
    median_score = df[score_col].median()
    df = df.copy()
    df['group'] = (df[score_col] > median_score).map({True: 'High', False: 'Low'})

    # Plot KM curves
    kmf = KaplanMeierFitter()

    colors = {'Low': '#2ecc71', 'High': '#e74c3c'}

    for group in ['Low', 'High']:
        mask = df['group'] == group
        kmf.fit(
            df.loc[mask, time_col],
            df.loc[mask, event_col],
            label=f'{group} {score_col}'
        )
        kmf.plot_survival_function(ax=ax, color=colors[group])

    # Log-rank test
    high_mask = df['group'] == 'High'
    result = logrank_test(
        df.loc[high_mask, time_col],
        df.loc[~high_mask, time_col],
        df.loc[high_mask, event_col],
        df.loc[~high_mask, event_col]
    )

    ax.set_title(f'Survival by {score_col}\nLog-rank p = {result.p_value:.4f}')
    ax.set_xlabel('Time')
    ax.set_ylabel('Survival Probability')

    plt.tight_layout()

    return fig
