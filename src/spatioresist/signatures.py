"""
Resistance Signatures

Pre-defined gene signatures for immunotherapy resistance.
"""

import numpy as np
import scanpy as sc


class ResistanceSignatures:
    """
    Collection of resistance-related gene signatures.

    Contains curated gene signatures for:
    - Innate resistance
    - Acquired resistance
    - Immune exclusion
    - T cell exhaustion
    - Macrophage polarization
    """

    # Innate resistance signature
    INNATE_RESISTANCE = [
        'AXL', 'TWIST1', 'WNT5A', 'MITF', 'SERPINF1',
        'VIM', 'ZEB1', 'CDH2', 'SNAI1', 'SNAI2'
    ]

    # Acquired resistance signature
    ACQUIRED_RESISTANCE = [
        'SERPINF1', 'VEGFA', 'MYC', 'HIF1A', 'LDHA',
        'SLC2A1', 'PKM', 'ENO1', 'GAPDH', 'PGK1'
    ]

    # Immune exclusion signature
    IMMUNE_EXCLUSION = [
        'TGFB1', 'FAP', 'POSTN', 'COL1A1', 'COL3A1',
        'FN1', 'VIM', 'ACTA2', 'THY1', 'PDGFRB'
    ]

    # T cell exhaustion markers
    T_CELL_EXHAUSTION = [
        'PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2',
        'TOX', 'TOX2', 'ENTPD1', 'ITGAE', 'CXCL13'
    ]

    # Cytotoxic T cell signature
    T_CELL_CYTOTOXIC = [
        'GZMA', 'GZMB', 'GZMK', 'PRF1', 'IFNG',
        'TNF', 'NKG7', 'GNLY', 'FASLG', 'CST7'
    ]

    # M1 macrophage (anti-tumor)
    MACROPHAGE_M1 = [
        'CD80', 'CD86', 'IL1B', 'TNF', 'NOS2',
        'IL6', 'IL12A', 'IL12B', 'CXCL9', 'CXCL10'
    ]

    # M2 macrophage (pro-tumor)
    MACROPHAGE_M2 = [
        'CD163', 'MRC1', 'IL10', 'ARG1', 'TGFB1',
        'CCL18', 'CCL22', 'MSR1', 'CD209', 'VEGFA'
    ]

    # Interferon gamma response
    IFNG_RESPONSE = [
        'IDO1', 'CXCL9', 'CXCL10', 'CXCL11', 'STAT1',
        'IRF1', 'GBP1', 'GBP2', 'GBP4', 'GBP5'
    ]

    # Tumor proliferation
    PROLIFERATION = [
        'MKI67', 'TOP2A', 'PCNA', 'MCM2', 'MCM3',
        'MCM4', 'MCM5', 'MCM6', 'MCM7', 'CDK1'
    ]

    @classmethod
    def get_all_signatures(cls):
        """
        Get all available signatures.

        Returns
        -------
        dict
            Dictionary of signature names to gene lists
        """
        return {
            'innate_resistance': cls.INNATE_RESISTANCE,
            'acquired_resistance': cls.ACQUIRED_RESISTANCE,
            'immune_exclusion': cls.IMMUNE_EXCLUSION,
            't_cell_exhaustion': cls.T_CELL_EXHAUSTION,
            't_cell_cytotoxic': cls.T_CELL_CYTOTOXIC,
            'macrophage_m1': cls.MACROPHAGE_M1,
            'macrophage_m2': cls.MACROPHAGE_M2,
            'ifng_response': cls.IFNG_RESPONSE,
            'proliferation': cls.PROLIFERATION
        }

    @classmethod
    def score_signatures(cls, adata, signatures=None, prefix='sig'):
        """
        Score all signatures in an AnnData object.

        Parameters
        ----------
        adata : AnnData
            Annotated data matrix
        signatures : dict, optional
            Custom signatures to use. If None, uses all default signatures.
        prefix : str
            Prefix for score column names

        Returns
        -------
        AnnData
            Data with signature scores added to obs
        """
        if signatures is None:
            signatures = cls.get_all_signatures()

        for name, genes in signatures.items():
            # Filter to genes present in data
            present_genes = [g for g in genes if g in adata.var_names]

            if len(present_genes) > 0:
                sc.tl.score_genes(
                    adata,
                    gene_list=present_genes,
                    score_name=f'{prefix}_{name}'
                )
            else:
                adata.obs[f'{prefix}_{name}'] = 0

        return adata

    @classmethod
    def score_bulk(cls, expression_df, signature_genes):
        """
        Score bulk RNA-seq samples using a gene signature.

        Uses mean z-score method.

        Parameters
        ----------
        expression_df : pd.DataFrame
            Expression matrix (genes x samples)
        signature_genes : list
            Gene list for scoring

        Returns
        -------
        pd.Series
            Score for each sample
        """
        import pandas as pd

        # Filter to present genes
        present = [g for g in signature_genes if g in expression_df.index]

        if len(present) == 0:
            return pd.Series(0, index=expression_df.columns)

        # Z-score normalization
        expr_subset = expression_df.loc[present]
        z_scores = (expr_subset - expr_subset.mean()) / expr_subset.std()

        # Mean z-score across genes
        return z_scores.mean(axis=0)
