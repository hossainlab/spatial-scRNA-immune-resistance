"""
Spatial Feature Extraction

Extract features from spatial transcriptomics data for resistance prediction.
"""

import numpy as np
import pandas as pd


class SpatialFeatureExtractor:
    """
    Extract spatial features for resistance prediction.

    Extracts features from:
    - Cell type abundances
    - Spatial statistics
    - Niche compositions
    - Gene expression signatures

    Examples
    --------
    >>> extractor = SpatialFeatureExtractor()
    >>> features = extractor.extract_all_features(adata)
    """

    def __init__(self):
        pass

    def extract_abundance_features(self, adata):
        """
        Extract cell type abundance statistics.

        Parameters
        ----------
        adata : AnnData
            Spatial data with abundance columns in obs

        Returns
        -------
        dict
            Abundance feature dictionary
        """
        features = {}

        # Get abundance columns
        abundance_cols = [c for c in adata.obs.columns if c.startswith('abundance_')]

        for col in abundance_cols:
            ct = col.replace('abundance_', '')
            values = adata.obs[col].values

            # Basic statistics
            features[f'{ct}_mean'] = np.mean(values)
            features[f'{ct}_std'] = np.std(values)
            features[f'{ct}_max'] = np.max(values)
            features[f'{ct}_min'] = np.min(values)
            features[f'{ct}_median'] = np.median(values)
            features[f'{ct}_q25'] = np.percentile(values, 25)
            features[f'{ct}_q75'] = np.percentile(values, 75)

            # Spatial concentration (Gini coefficient)
            features[f'{ct}_gini'] = self._gini_coefficient(values)

            # Fraction of spots with this cell type
            features[f'{ct}_presence'] = np.mean(values > 0)

        return features

    def extract_spatial_stats(self, adata):
        """
        Extract spatial statistics features.

        Parameters
        ----------
        adata : AnnData
            Spatial data with computed spatial statistics

        Returns
        -------
        dict
            Spatial statistics features
        """
        features = {}

        # Neighborhood enrichment results
        if 'nhood_enrichment' in adata.uns:
            try:
                nhood = adata.uns['nhood_enrichment']['zscore']

                if isinstance(nhood, (pd.DataFrame, np.ndarray)):
                    if isinstance(nhood, pd.DataFrame):
                        cell_types = nhood.index.tolist()
                    else:
                        cell_types = [f'type_{i}' for i in range(nhood.shape[0])]

                    # Extract key interactions
                    for i, ct1 in enumerate(cell_types):
                        for j, ct2 in enumerate(cell_types):
                            if i <= j:  # Upper triangle
                                if isinstance(nhood, pd.DataFrame):
                                    val = nhood.iloc[i, j]
                                else:
                                    val = nhood[i, j]
                                features[f'nhood_{ct1}_{ct2}'] = val
            except Exception:
                pass

        # Spatial autocorrelation (Moran's I)
        if 'moranI' in adata.uns:
            try:
                moran = adata.uns['moranI']
                if isinstance(moran, pd.DataFrame):
                    for gene in moran.index:
                        if 'I' in moran.columns:
                            features[f'moran_{gene}'] = moran.loc[gene, 'I']
            except Exception:
                pass

        # Ripley's statistics if available
        if 'ripley' in adata.uns:
            try:
                ripley = adata.uns['ripley']
                for key in ripley:
                    features[f'ripley_{key}'] = ripley[key]
            except Exception:
                pass

        return features

    def extract_niche_features(self, adata):
        """
        Extract features from spatial niches.

        Parameters
        ----------
        adata : AnnData
            Spatial data with niche annotations

        Returns
        -------
        dict
            Niche features
        """
        features = {}

        if 'spatial_niche' in adata.obs.columns:
            niche_counts = adata.obs['spatial_niche'].value_counts(normalize=True)

            for niche, prop in niche_counts.items():
                features[f'niche_{niche}_proportion'] = prop

            # Niche diversity (Shannon entropy)
            proportions = niche_counts.values
            proportions = proportions[proportions > 0]
            entropy = -np.sum(proportions * np.log(proportions))
            features['niche_diversity'] = entropy

            # Number of niches
            features['n_niches'] = len(niche_counts)

        return features

    def extract_expression_features(self, adata):
        """
        Extract gene expression signature features.

        Parameters
        ----------
        adata : AnnData
            Spatial data with expression scores

        Returns
        -------
        dict
            Expression features
        """
        features = {}

        # Look for score columns
        score_cols = [c for c in adata.obs.columns if c.endswith('_score')]

        for col in score_cols:
            values = adata.obs[col].values

            features[f'{col}_mean'] = np.mean(values)
            features[f'{col}_std'] = np.std(values)
            features[f'{col}_max'] = np.max(values)
            features[f'{col}_q75'] = np.percentile(values, 75)

        return features

    def extract_all_features(self, adata):
        """
        Extract all spatial features.

        Parameters
        ----------
        adata : AnnData
            Spatial AnnData object

        Returns
        -------
        dict
            All extracted features
        """
        features = {}

        # Combine all feature types
        features.update(self.extract_abundance_features(adata))
        features.update(self.extract_spatial_stats(adata))
        features.update(self.extract_niche_features(adata))
        features.update(self.extract_expression_features(adata))

        # Add sample-level info
        features['n_spots'] = adata.n_obs
        features['n_genes'] = adata.n_vars

        return features

    @staticmethod
    def _gini_coefficient(values):
        """
        Calculate Gini coefficient for spatial concentration.

        Parameters
        ----------
        values : array-like
            Values to compute Gini coefficient for

        Returns
        -------
        float
            Gini coefficient (0 = perfect equality, 1 = perfect inequality)
        """
        values = np.array(values).flatten()
        values = values[~np.isnan(values)]

        if len(values) == 0 or np.sum(values) == 0:
            return 0

        values = np.sort(values)
        n = len(values)
        index = np.arange(1, n + 1)

        return (2 * np.sum(index * values) / (n * np.sum(values))) - (n + 1) / n
