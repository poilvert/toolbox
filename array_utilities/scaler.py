#!/usr/bin/env python

# Authors: Nicolas Poilvert, <nicolas.poilvert@gmail.com>
# Licence: BSD 3-clause

import numpy as np

class OnlineScaler(object):

    def __init__(self, d):
        """initializes an Online Scaler object

        Parameters
        ----------

        `d`: int
            dimensionality of the feature vectors

        attributes
        ----------

        `r_ns`: int64
            running number of sample feature vectors

        `r_mn`: 1D array of float64
            running mean feature vector

        `r_va`: 1D array of float64
            running vector of feature-wise variance

        methods
        -------

        `partial_fit`:
            updates the values of `r_ns`, `r_mn` and `r_va`
            using the new data in the batch

        `transform`:
            returns a batch, where all the feature vectors have
            been translated uniformly by the current mean feature
            vector and normalized by the current standard deviation
            (i.e. the square root of the current variance) feature
            wise.

        conventions
        -----------
        'r' means 'running' like in 'running average'
        'ns' stands for 'number of samples'
        'mn' stands for 'mean'
        'va' stands for 'variance'
        'd' is the dimensionality of the feature vector
        """

        self.r_ns = np.int64(0)
        self.r_mn = np.zeros(d, dtype=np.float64)
        self.r_va = np.zeros(d, dtype=np.float64)
        self.d = np.int64(d)

    def partial_fit(self, X):

        # -- X has to be 2D
        if X.ndim != 2:
            raise ValueError('training batch must be 2D')

        # -- makes sure that the batch X has a second
        #    dimension equal to d
        if X.shape[1] != self.d:
            raise ValueError('wrong number of features')

        # -- makes sure that the batch has at least one
        #    sample
        if X.shape[0] < 1:
            raise ValueError('batch is too small')

        # -- extract number of samples in batch and
        #    computes the sum of X and X**2 of the batch
        b_ns = np.int64(X.shape[0])
        b_su = np.float64(X.sum(axis=0))
        b_sq = np.float64((X**2).sum(axis=0))

        # -- updating the values for the number of samples
        #    the mean vector X and the variance. 'c' stands
        #    for 'current'
        c_ns = self.r_ns
        c_mn = self.r_mn
        c_va = self.r_va

        self.r_ns = c_ns + b_ns
        self.r_mn = 1. / self.r_ns * (c_ns * c_mn + b_su)
        self.r_va = 1. / self.r_ns * (c_ns * (c_va + c_mn**2) + b_sq) \
                    - self.r_mn**2

        return self

    def transform(self, X):

        _mean = (self.r_mn[np.newaxis, :]).astype(X.dtype)
        _std = np.sqrt((self.r_va[np.newaxis, :]).astype(X.dtype))

        out = X - _mean

        if (_std == 0).any():
            import warnings
            warnings.warn("Zero standard deviation detected!")
            _std = 1
        out /= _std

        return out

    def fit_transform(self, X):

        self.partial_fit(X)
        out = self.transform(X)

        return out
