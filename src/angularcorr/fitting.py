"""Binomial-likelihood fitting of exponential response coefficients."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import xlogy

from angularcorr.models import CoefficientFit


def fit_response_coefficient(
    path_lengths: np.ndarray,
    trials: np.ndarray,
    successes: np.ndarray,
    *,
    response: str,
    minimum_path_length: float = 0.05,
) -> CoefficientFit:
    paths = np.asarray(path_lengths, dtype=float)
    n = np.asarray(trials, dtype=np.int64)
    k = np.asarray(successes, dtype=np.int64)
    if minimum_path_length < 0:
        raise ValueError("Minimum path length cannot be negative.")
    valid = (n > 0) & (paths > minimum_path_length) & (k > 0) & (k < n)
    paths, n, k = paths[valid], n[valid], k[valid]
    if paths.size == 0:
        raise ValueError("No populated angular bins are available for fitting.")
    if np.any(paths <= 0) or np.any(k < 0) or np.any(k > n):
        raise ValueError("Invalid path lengths or binomial event counts.")
    if int(np.sum(k)) == 0 or int(np.sum(k)) == int(np.sum(n)):
        raise ValueError("Coefficient fitting requires both successes and failures.")

    def negative_log_likelihood(mu: float) -> float:
        probability = -np.expm1(-mu * paths)
        log_failure = -mu * paths
        value = xlogy(k, probability) + (n - k) * log_failure
        return -float(np.sum(value))

    fitted = minimize_scalar(
        negative_log_likelihood,
        method="bounded",
        bounds=(np.finfo(float).eps, 10.0),
        options={"xatol": 1e-13},
    )
    if not fitted.success or not np.isfinite(fitted.fun):
        raise RuntimeError(f"Coefficient fit failed: {fitted.message}")

    mu = float(fitted.x)
    survival = np.exp(-mu * paths)
    probability = 1.0 - survival
    derivative = paths * survival
    fisher_information = float(
        np.sum(n * derivative**2 / (probability * (1.0 - probability)))
    )
    standard_error = fisher_information**-0.5
    z_95 = 1.959963984540054
    return CoefficientFit(
        response=response,
        mu=mu,
        standard_error=standard_error,
        confidence_low=max(0.0, mu - z_95 * standard_error),
        confidence_high=mu + z_95 * standard_error,
        negative_log_likelihood=float(fitted.fun),
        bins_used=int(paths.size),
        events_used=int(np.sum(n)),
        successes=int(np.sum(k)),
    )
