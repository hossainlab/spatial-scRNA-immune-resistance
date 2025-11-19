"""
SpatioResist Model

Core machine learning model for spatial resistance prediction.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import joblib

from .features import SpatialFeatureExtractor


class SpatioResist:
    """
    SpatioResist: Spatial Resistance Prediction Model.

    Predicts immunotherapy resistance from spatial transcriptomics features
    including cell type abundances, spatial statistics, and niche compositions.

    Parameters
    ----------
    model_type : str
        Model type: 'rf' (Random Forest), 'gb' (Gradient Boosting), 'lr' (Logistic Regression)
    random_state : int
        Random seed for reproducibility
    **kwargs
        Additional parameters passed to the underlying model

    Examples
    --------
    >>> from spatioresist import SpatioResist
    >>> model = SpatioResist(model_type='rf')
    >>> X = model.extract_features(adata_list)
    >>> model.fit(X, y)
    >>> risk_scores = model.predict_proba(X_new)
    """

    def __init__(self, model_type='rf', random_state=42, **kwargs):
        self.model_type = model_type
        self.random_state = random_state
        self.feature_extractor = SpatialFeatureExtractor()
        self.scaler = StandardScaler()
        self.feature_names = None
        self.is_fitted = False

        # Initialize model
        if model_type == 'rf':
            self.model = RandomForestClassifier(
                n_estimators=kwargs.get('n_estimators', 100),
                max_depth=kwargs.get('max_depth', 10),
                random_state=random_state,
                n_jobs=-1
            )
        elif model_type == 'gb':
            self.model = GradientBoostingClassifier(
                n_estimators=kwargs.get('n_estimators', 100),
                max_depth=kwargs.get('max_depth', 5),
                random_state=random_state
            )
        elif model_type == 'lr':
            self.model = LogisticRegression(
                random_state=random_state,
                max_iter=1000
            )
        else:
            raise ValueError(f"Unknown model_type: {model_type}")

    def extract_features(self, adata_list):
        """
        Extract spatial features from multiple samples.

        Parameters
        ----------
        adata_list : list of AnnData
            List of deconvolved spatial AnnData objects

        Returns
        -------
        pd.DataFrame
            Feature matrix (samples x features)
        """
        feature_dicts = []

        for adata in adata_list:
            features = self.feature_extractor.extract_all_features(adata)
            feature_dicts.append(features)

        df = pd.DataFrame(feature_dicts)
        self.feature_names = df.columns.tolist()

        return df

    def fit(self, X, y):
        """
        Train the resistance prediction model.

        Parameters
        ----------
        X : pd.DataFrame or np.ndarray
            Feature matrix (samples x features)
        y : array-like
            Response labels (0=responder, 1=non-responder)

        Returns
        -------
        self
            Fitted model
        """
        # Store feature names if not already set
        if isinstance(X, pd.DataFrame) and self.feature_names is None:
            self.feature_names = X.columns.tolist()

        # Handle missing values
        X_clean = pd.DataFrame(X).fillna(0).values if isinstance(X, pd.DataFrame) else np.nan_to_num(X)

        # Scale features
        X_scaled = self.scaler.fit_transform(X_clean)

        # Train model
        self.model.fit(X_scaled, y)
        self.is_fitted = True

        return self

    def predict(self, X):
        """
        Predict resistance class.

        Parameters
        ----------
        X : pd.DataFrame or np.ndarray
            Feature matrix

        Returns
        -------
        np.ndarray
            Predicted classes (0=responder, 1=non-responder)
        """
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")

        X_clean = pd.DataFrame(X).fillna(0).values if isinstance(X, pd.DataFrame) else np.nan_to_num(X)
        X_scaled = self.scaler.transform(X_clean)

        return self.model.predict(X_scaled)

    def predict_proba(self, X):
        """
        Predict resistance probability.

        Parameters
        ----------
        X : pd.DataFrame or np.ndarray
            Feature matrix

        Returns
        -------
        np.ndarray
            Predicted probabilities of resistance
        """
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")

        X_clean = pd.DataFrame(X).fillna(0).values if isinstance(X, pd.DataFrame) else np.nan_to_num(X)
        X_scaled = self.scaler.transform(X_clean)

        return self.model.predict_proba(X_scaled)[:, 1]

    def score_sample(self, adata):
        """
        Generate resistance risk score for a single spatial sample.

        Parameters
        ----------
        adata : AnnData
            Deconvolved spatial data

        Returns
        -------
        float
            Resistance risk score (0-1)
        dict
            Feature contributions
        """
        if not self.is_fitted:
            raise ValueError("Model not fitted. Call fit() first.")

        # Extract features
        features = self.feature_extractor.extract_all_features(adata)
        X = pd.DataFrame([features])

        # Ensure correct feature order
        if self.feature_names:
            missing = set(self.feature_names) - set(X.columns)
            for col in missing:
                X[col] = 0
            X = X[self.feature_names]

        # Get prediction
        risk_score = self.predict_proba(X)[0]

        # Get feature contributions
        importance = self.get_feature_importance()

        contributions = {}
        if importance is not None:
            for feat in importance.head(10).index:
                if feat in features:
                    contributions[feat] = {
                        'value': features[feat],
                        'importance': importance[feat]
                    }

        return risk_score, contributions

    def get_feature_importance(self):
        """
        Get feature importance scores.

        Returns
        -------
        pd.Series
            Feature importance values sorted by importance
        """
        if not self.is_fitted:
            return None

        if hasattr(self.model, 'feature_importances_'):
            importance = self.model.feature_importances_
        elif hasattr(self.model, 'coef_'):
            importance = np.abs(self.model.coef_[0])
        else:
            return None

        if self.feature_names is None:
            return pd.Series(importance).sort_values(ascending=False)

        return pd.Series(importance, index=self.feature_names).sort_values(ascending=False)

    def save(self, path):
        """
        Save model to disk.

        Parameters
        ----------
        path : str or Path
            Output path
        """
        joblib.dump({
            'model': self.model,
            'scaler': self.scaler,
            'feature_names': self.feature_names,
            'model_type': self.model_type,
            'random_state': self.random_state,
            'is_fitted': self.is_fitted
        }, path)

    @classmethod
    def load(cls, path):
        """
        Load model from disk.

        Parameters
        ----------
        path : str or Path
            Model path

        Returns
        -------
        SpatioResist
            Loaded model
        """
        data = joblib.load(path)

        instance = cls(
            model_type=data['model_type'],
            random_state=data['random_state']
        )
        instance.model = data['model']
        instance.scaler = data['scaler']
        instance.feature_names = data['feature_names']
        instance.is_fitted = data['is_fitted']

        return instance

    def __repr__(self):
        status = "fitted" if self.is_fitted else "not fitted"
        return f"SpatioResist(model_type='{self.model_type}', {status})"
