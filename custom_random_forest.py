import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from multiprocessing import Pool


class RandomForestClassifierCustom(BaseEstimator):
    """
    Custom RF classifier that supports parallel fit and predict
    """
    def __init__(
        self,
        random_state,
        n_estimators=10,
        min_samples_leaf=1,
        max_depth=None,
        max_features=None,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.min_samples_leaf = min_samples_leaf
        self.trees = []
        self.feat_ids_by_tree = []

    def _fit(self, args):
        X, y, n_estimators = args
        n_samples, n_features = X.shape
        np.random.seed(self.random_state + n_estimators)
        feat_ids = np.random.choice(n_features, self.max_features, replace=False)

        sample_indices = np.random.randint(0, n_samples, size=n_samples)
        X_sampled = X[sample_indices][:, feat_ids].reshape(X.shape[0], -1)
        y_sampled = y[sample_indices]

        tree = DecisionTreeClassifier(
            max_depth=self.max_depth,
            max_features=self.max_features,
            random_state=self.random_state,
        )
        tree.fit(X_sampled, y_sampled)
        return tree, feat_ids

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))
        args = [(X, y, n) for n in range(self.n_estimators)]
        with Pool(n_jobs) as pool:
            self.trees, self.feat_ids_by_tree = map(
                list, zip(*pool.map(self._fit, args))
            )
        return self

    def _predict_proba(self, args):
        tree, feat_id, X = args
        proba = tree.predict_proba(X[:, feat_id])
        return proba

    def predict_proba(self, X, n_jobs=1):
        args = [
            (tree, feat_ids, X)
            for tree, feat_ids in zip(self.trees, self.feat_ids_by_tree)
        ]
        with Pool(n_jobs) as pool:
            probas = pool.map(self._predict_proba, args)
        return np.mean(probas, axis=0)

    def predict(self, X, n_jobs=1):
        proba = self.predict_proba(X, n_jobs)
        prediction = np.argmax(proba, axis=1)
        return prediction
