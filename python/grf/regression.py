import numpy as np
from ._grf_python import regression_train, regression_predict, regression_predict_oob


class RegressionForest:
    def __init__(
        self,
        n_estimators=100,
        mtry=None,
        min_node_size=5,
        honesty=True,
        honesty_fraction=0.5,
        honesty_prune_leaves=True,
        alpha=0.05,
        imbalance_penalty=0,
        ci_group_size=2,
        sample_fraction=0.5,
        clusters=None,
        samples_per_cluster=0,
        compute_oob_predictions=False,
        num_threads=0,
        seed=42,
    ):
        self.n_estimators = n_estimators
        self.mtry = mtry
        self.min_node_size = min_node_size
        self.honesty = honesty
        self.honesty_fraction = honesty_fraction
        self.honesty_prune_leaves = honesty_prune_leaves
        self.alpha = alpha
        self.imbalance_penalty = imbalance_penalty
        self.ci_group_size = ci_group_size
        self.sample_fraction = sample_fraction
        self.clusters = clusters if clusters is not None else []
        self.samples_per_cluster = samples_per_cluster
        self.compute_oob_predictions = compute_oob_predictions
        self.num_threads = num_threads
        self.seed = seed
        self.forest = None

    def fit(self, X, y):
        if not isinstance(X, np.ndarray) or not isinstance(y, np.ndarray):
            raise ValueError("X and y must be numpy arrays")

        if X.ndim != 2:
            raise ValueError("X must be a 2-dimensional array")

        if y.ndim != 1:
            raise ValueError("y must be a 1-dimensional array")

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y must have the same number of samples")

        # Combine X and y into a single matrix
        train_matrix = np.column_stack((X, y))

        # Set mtry if not specified
        if self.mtry is None:
            self.mtry = max(np.ceil(np.sqrt(X.shape[1])), 1)

        self.forest = regression_train(
            train_matrix,
            X.shape[1],  # outcome_index (last column)
            0,  # sample_weight_index (not used, set to 0)
            False,  # use_sample_weights
            int(self.mtry),
            self.n_estimators,
            self.min_node_size,
            self.sample_fraction,
            self.honesty,
            self.honesty_fraction,
            self.honesty_prune_leaves,
            self.ci_group_size,
            self.alpha,
            self.imbalance_penalty,
            self.clusters,
            self.samples_per_cluster,
            self.compute_oob_predictions,
            self.num_threads,
            self.seed,
        )
        return self

    def predict(self, X):
        if self.forest is None:
            raise ValueError("Model not fitted. Call 'fit' first.")

        if not isinstance(X, np.ndarray):
            raise ValueError("X must be a numpy array")

        if X.ndim != 2:
            raise ValueError("X must be a 2-dimensional array")

        # The forest is nested inside the 'forest' key
        forest_data = self.forest["forest"]

        predictions = regression_predict(
            forest_data,
            np.column_stack((X, np.zeros(X.shape[0]))),  # Add dummy column for outcome
            X.shape[1],  # outcome_index
            X,
            self.num_threads,
            False,  # estimate_variance
        )
        return predictions["predictions"]

    def predict_oob(self):
        if self.forest is None:
            raise ValueError("Model not fitted. Call 'fit' first.")

        # The forest is nested inside the 'forest' key
        forest_data = self.forest["forest"]

        predictions = regression_predict_oob(
            forest_data, self.num_threads, False  # estimate_variance
        )

        return predictions["predictions"]
