"""
SpatioResist: Spatial Resistance Prediction Framework

A machine learning framework for predicting immunotherapy resistance
from spatial transcriptomics data.

Author: Spatial-scRNA-seq Immunotherapy Resistance Atlas Project
"""

from .model import SpatioResist
from .features import SpatialFeatureExtractor
from .signatures import ResistanceSignatures
from .visualization import plot_resistance_scores, plot_spatial_risk

__version__ = "0.1.0"

__all__ = [
    "SpatioResist",
    "SpatialFeatureExtractor",
    "ResistanceSignatures",
    "plot_resistance_scores",
    "plot_spatial_risk",
]
