"""Geometrical acceptance, event selections, and angular binning."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from angularcorr.geometry import AnnularGeometry
from angularcorr.root_io import EventArrays


@dataclass(frozen=True)
class AnalysisSettings:
    photon_energy_keV: float = 511.0
    energy_half_window_keV: float = 0.01
    angular_bin_width_deg: float = 0.25

    def __post_init__(self) -> None:
        if self.photon_energy_keV <= 0:
            raise ValueError("Photon energy must be positive.")
        if self.energy_half_window_keV <= 0:
            raise ValueError("Full-energy half-window must be positive.")
        if self.angular_bin_width_deg <= 0:
            raise ValueError("Angular-bin width must be positive.")


@dataclass(frozen=True)
class BinnedEvents:
    lower_edges_rad: np.ndarray
    upper_edges_rad: np.ndarray
    centers_rad: np.ndarray
    n_accepted: np.ndarray
    n_total: np.ndarray
    n_peak: np.ndarray
    emitted_events: int
    accepted_events: int
    total_interaction_events: int
    full_energy_events: int

    @property
    def bins(self) -> int:
        return int(self.centers_rad.size)


def _angular_edges(geometry: AnnularGeometry, width_deg: float) -> np.ndarray:
    width = np.radians(width_deg)
    span = np.pi / 2.0 - geometry.theta_min
    complete_bins = int(np.floor(span / width))
    edges = geometry.theta_min + np.arange(complete_bins + 1, dtype=float) * width
    tolerance = 64.0 * np.finfo(float).eps
    if np.pi / 2.0 - edges[-1] > tolerance:
        edges = np.append(edges, np.pi / 2.0)
    else:
        edges[-1] = np.pi / 2.0
    return edges


def bin_events(
    events: EventArrays,
    geometry: AnnularGeometry,
    settings: AnalysisSettings = AnalysisSettings(),
) -> BinnedEvents:
    folded_theta = np.minimum(events.theta_rad, np.pi - events.theta_rad)
    accepted = folded_theta >= geometry.theta_min
    full_energy = (
        np.abs(events.energy_crystal_keV - settings.photon_energy_keV)
        <= settings.energy_half_window_keV
    )
    edges = _angular_edges(geometry, settings.angular_bin_width_deg)
    accepted_theta = folded_theta[accepted]
    n_accepted, _ = np.histogram(accepted_theta, bins=edges)
    n_total, _ = np.histogram(
        folded_theta[accepted & events.primary_interacted], bins=edges
    )
    n_peak, _ = np.histogram(folded_theta[accepted & full_energy], bins=edges)
    if np.any(n_total > n_accepted) or np.any(n_peak > n_accepted):
        raise ArithmeticError("A response count exceeds its bin denominator.")
    return BinnedEvents(
        lower_edges_rad=edges[:-1],
        upper_edges_rad=edges[1:],
        centers_rad=(edges[:-1] + edges[1:]) * 0.5,
        n_accepted=n_accepted,
        n_total=n_total,
        n_peak=n_peak,
        emitted_events=events.entries,
        accepted_events=int(np.count_nonzero(accepted)),
        total_interaction_events=int(np.count_nonzero(accepted & events.primary_interacted)),
        full_energy_events=int(np.count_nonzero(accepted & full_energy)),
    )
