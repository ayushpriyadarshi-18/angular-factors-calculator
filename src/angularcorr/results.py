"""Algebraic calculation of angular-correlation factors."""

from __future__ import annotations

from angularcorr.models import CorrelationFactors, Moments


def calculate_factors(moments: Moments) -> CorrelationFactors:
    if moments.mean_peak <= 0 or moments.mean_total <= 0:
        raise ValueError("Mean peak and total efficiencies must be positive.")
    return CorrelationFactors(
        w_l1=moments.mean_peak_total
        / (moments.mean_peak * moments.mean_total),
        w_l2=moments.mean_total_squared / moments.mean_total**2,
        w_g=moments.mean_peak_squared / moments.mean_peak**2,
    )
