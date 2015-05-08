#!/usr/bin/env python

"""
Small test suite for online scaler
"""

import numpy as np
from numpy.testing import assert_allclose
from pytest import raises

from scaler import OnlineScaler

DTYPE = np.float32
RTOL = 1e-6
ATOL = 1e-5

def test_wrong_feature_dimension():

    n_samples = 1000
    n_features = 100
    data = np.random.randn(n_samples, n_features).astype(DTYPE)

    scaler = OnlineScaler(n_features + 2)
    raises(ValueError, scaler.partial_fit, data)

def test_one_batch_fit():

    n_samples = 1000
    n_features = 100
    data = np.random.randn(n_samples, n_features).astype(DTYPE)

    # -- reference outputs
    ref_mean = data.mean(axis=0)
    ref_std = data.std(axis=0)
    ref_out = (data - ref_mean[np.newaxis]) / ref_std[np.newaxis, :]

    # -- calculation with OnlineScaler
    scaler = OnlineScaler(n_features)
    scaler.partial_fit(data)
    out_mean = scaler.r_mn
    out_std = np.sqrt(scaler.r_va)
    out = scaler.transform(data)

    # -- tests
    assert_allclose(out_mean, ref_mean, rtol=RTOL, atol=ATOL)
    assert_allclose(out_std, ref_std, rtol=RTOL, atol=ATOL)
    assert_allclose(out, ref_out, rtol=RTOL, atol=ATOL)

def test_one_batch_fit_transform():

    n_samples = 1000
    n_features = 100
    data = np.random.randn(n_samples, n_features).astype(DTYPE)

    # -- reference outputs
    ref_mean = data.mean(axis=0)
    ref_std = data.std(axis=0)
    ref_out = (data - ref_mean[np.newaxis]) / ref_std[np.newaxis, :]

    # -- calculation with OnlineScaler
    scaler = OnlineScaler(n_features)
    out = scaler.fit_transform(data)

    # -- tests
    assert_allclose(out, ref_out, rtol=RTOL, atol=ATOL)

def test_multiple_batch_fit_for_mean_and_std():

    n_samples = 1000
    n_features = 100
    n_batches = 5
    data = np.random.randn(n_batches * n_samples, n_features).astype(DTYPE)

    # -- reference outputs
    ref_mean = data.mean(axis=0)
    ref_std = data.std(axis=0)

    # -- calculation with OnlineScaler
    scaler = OnlineScaler(n_features)
    for i in xrange(n_batches):
        scaler.partial_fit(data[(i * n_samples):((i + 1) * n_samples), :])
    out_mean = scaler.r_mn
    out_std = np.sqrt(scaler.r_va)

    # -- tests
    assert_allclose(out_mean, ref_mean, rtol=RTOL, atol=ATOL)
    assert_allclose(out_std, ref_std, rtol=RTOL, atol=ATOL)
