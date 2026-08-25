"""End-to-end scientific workflows."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from angularcorr.binning import AnalysisSettings, BinnedEvents, bin_events
from angularcorr.fitting import fit_response_coefficient
from angularcorr.geometry import AnnularGeometry
from angularcorr.integration import calculate_moments_midpoint, calculate_moments_quad
from angularcorr.models import Coefficients, MassCoefficients
from angularcorr.results import calculate_factors
from angularcorr.root_io import read_events


def _bin_summary(bins: BinnedEvents) -> dict[str, Any]:
    return {
        "bins": bins.bins,
        "emitted_events": bins.emitted_events,
        "accepted_events": bins.accepted_events,
        "total_interaction_events": bins.total_interaction_events,
        "full_energy_events": bins.full_energy_events,
    }


def analyze_root(
    path: str | Path,
    geometry: AnnularGeometry,
    settings: AnalysisSettings = AnalysisSettings(),
    *,
    density: float,
    midpoint_intervals: int = 1_000_000,
) -> tuple[dict[str, Any], BinnedEvents]:
    if density <= 0:
        raise ValueError("Density must be positive.")
    inspection, events = read_events(path)
    binned = bin_events(events, geometry, settings)
    paths = geometry.path_length(binned.centers_rad)
    total_fit = fit_response_coefficient(
        paths, binned.n_accepted, binned.n_total, response="total_interaction"
    )
    peak_fit = fit_response_coefficient(
        paths, binned.n_accepted, binned.n_peak, response="full_energy_peak"
    )
    coefficients = Coefficients(mu_total=total_fit.mu, mu_peak=peak_fit.mu)
    mass_coefficients = MassCoefficients.from_linear(coefficients, density)
    quad_result = calculate_moments_quad(geometry, coefficients)
    midpoint_result = calculate_moments_midpoint(
        geometry, coefficients, intervals=midpoint_intervals
    )
    factors = calculate_factors(quad_result.moments)
    midpoint_factors = calculate_factors(midpoint_result.moments)
    expected_accepted = geometry.fraction_of_full_sphere * binned.emitted_events
    result = {
        "schema_version": 1,
        "root": inspection.to_dict(),
        "geometry": geometry.to_dict(),
        "settings": asdict(settings),
        "event_summary": {
            **_bin_summary(binned),
            "expected_accepted_events": expected_accepted,
            "acceptance_residual_events": binned.accepted_events - expected_accepted,
            "theta_consistency_failures": events.theta_consistency_failures,
        },
        "fits": {
            "total": total_fit.to_dict(),
            "peak": peak_fit.to_dict(),
        },
        "material": {"density_g_cm3": density},
        "coefficients": {
            "linear_cm_inverse": coefficients.to_dict(),
            "mass_cm2_g": mass_coefficients.to_dict(),
        },
        "integration": {
            "quad": quad_result.to_dict(),
            "midpoint": midpoint_result.to_dict(),
            "maximum_moment_difference": max(
                abs(a - b)
                for a, b in zip(
                    quad_result.moments.to_dict().values(),
                    midpoint_result.moments.to_dict().values(),
                    strict=True,
                )
            ),
        },
        "factors": factors.to_dict(),
        "midpoint_factors": midpoint_factors.to_dict(),
    }
    return result, binned


def calculate_from_mass_coefficients(
    geometry: AnnularGeometry,
    mass_coefficients: MassCoefficients,
    *,
    density: float,
    midpoint_intervals: int = 1_000_000,
) -> dict[str, Any]:
    coefficients = mass_coefficients.to_linear(density)
    quad_result = calculate_moments_quad(geometry, coefficients)
    midpoint_result = calculate_moments_midpoint(
        geometry, coefficients, intervals=midpoint_intervals
    )
    return {
        "schema_version": 1,
        "geometry": geometry.to_dict(),
        "material": {"density_g_cm3": density},
        "coefficients": {
            "linear_cm_inverse": coefficients.to_dict(),
            "mass_cm2_g": mass_coefficients.to_dict(),
        },
        "integration": {
            "quad": quad_result.to_dict(),
            "midpoint": midpoint_result.to_dict(),
        },
        "factors": calculate_factors(quad_result.moments).to_dict(),
    }


def binned_events_as_columns(bins: BinnedEvents) -> dict[str, np.ndarray]:
    return {
        "theta_lower_deg": np.degrees(bins.lower_edges_rad),
        "theta_upper_deg": np.degrees(bins.upper_edges_rad),
        "theta_center_deg": np.degrees(bins.centers_rad),
        "n_accepted": bins.n_accepted,
        "n_total": bins.n_total,
        "n_peak": bins.n_peak,
    }
