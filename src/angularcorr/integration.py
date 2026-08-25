"""Numerical angular averages for the SNP annular-detector formalism."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.integrate import quad

from angularcorr.geometry import AnnularGeometry
from angularcorr.models import Coefficients, IntegrationResult, MomentErrors, Moments
from angularcorr.response import efficiency


def _moment_values(
    geometry: AnnularGeometry,
    coefficients: Coefficients,
    theta: np.ndarray,
) -> tuple[np.ndarray, ...]:
    path = geometry.path_length(theta)
    peak = efficiency(coefficients.mu_peak, path)
    total = efficiency(coefficients.mu_total, path)
    return peak, total, peak * total, total**2, peak**2


def _integrate_upper_half(
    geometry: AnnularGeometry,
    function: Callable[[float], float],
    *,
    epsabs: float,
    epsrel: float,
) -> tuple[float, float]:
    """Integrate a symmetric angular average over the accepted upper half."""

    first, first_error = quad(
        function,
        geometry.theta_min,
        geometry.theta_cap,
        epsabs=epsabs,
        epsrel=epsrel,
    )
    second, second_error = quad(
        function,
        geometry.theta_cap,
        np.pi / 2.0,
        epsabs=epsabs,
        epsrel=epsrel,
    )
    return first + second, first_error + second_error


def calculate_moments_quad(
    geometry: AnnularGeometry,
    coefficients: Coefficients,
    *,
    epsabs: float = 1e-12,
    epsrel: float = 1e-10,
) -> IntegrationResult:
    """Calculate moments with piecewise adaptive SciPy quadrature."""

    if epsabs <= 0 or epsrel <= 0:
        raise ValueError("Integration tolerances must be positive.")

    def component(index: int) -> Callable[[float], float]:
        def integrand(theta: float) -> float:
            values = _moment_values(
                geometry, coefficients, np.asarray(theta, dtype=float)
            )
            return float(values[index] * np.sin(theta))

        return integrand

    results = [
        _integrate_upper_half(
            geometry, component(index), epsabs=epsabs, epsrel=epsrel
        )
        for index in range(5)
    ]
    values = [item[0] for item in results]
    errors = [item[1] for item in results]
    return IntegrationResult(
        method="scipy_quad_symmetry_reduced",
        moments=Moments(*values),
        errors=MomentErrors(*errors),
    )


def calculate_moments_midpoint(
    geometry: AnnularGeometry,
    coefficients: Coefficients,
    *,
    intervals: int = 1_000_000,
) -> IntegrationResult:
    """Independently integrate the complete accepted angular interval."""

    if intervals < 2:
        raise ValueError("Midpoint integration needs at least two intervals.")
    width = (geometry.theta_max - geometry.theta_min) / intervals
    theta = geometry.theta_min + (np.arange(intervals) + 0.5) * width
    weight = np.sin(theta) * 0.5 * width
    values = _moment_values(geometry, coefficients, theta)
    moments = Moments(*(float(np.sum(value * weight)) for value in values))
    return IntegrationResult(method=f"midpoint_full_interval_{intervals}", moments=moments)
