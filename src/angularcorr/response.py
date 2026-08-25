"""Directional detector-response functions."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def efficiency(mu: float, path_length: ArrayLike) -> NDArray[np.float64]:
    """Return ``1 - exp(-mu*x)`` for a positive linear coefficient."""

    if mu <= 0 or not np.isfinite(mu):
        raise ValueError("The linear response coefficient must be positive and finite.")
    paths = np.asarray(path_length, dtype=float)
    if np.any(~np.isfinite(paths)) or np.any(paths < 0):
        raise ValueError("Path lengths must be finite and non-negative.")
    return -np.expm1(-mu * paths)
