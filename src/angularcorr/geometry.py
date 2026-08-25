"""Geometry and path length for a centered bare annular detector."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class AnnularGeometry:
    """Dimensions of an annular cylinder, expressed in centimetres."""

    r1: float
    r2: float
    h: float

    def __post_init__(self) -> None:
        if self.r1 <= 0:
            raise ValueError("Inner radius r1 must be positive.")
        if self.r2 <= self.r1:
            raise ValueError("Outer radius r2 must be greater than r1.")
        if self.h <= 0:
            raise ValueError("Half-length h must be positive.")

    @property
    def theta_min(self) -> float:
        return float(np.arctan(self.r1 / self.h))

    @property
    def theta_cap(self) -> float:
        return float(np.arctan(self.r2 / self.h))

    @property
    def theta_max(self) -> float:
        return float(np.pi - self.theta_min)

    @property
    def solid_angle_sr(self) -> float:
        return float(4.0 * np.pi * np.cos(self.theta_min))

    @property
    def fraction_of_full_sphere(self) -> float:
        return float(np.cos(self.theta_min))

    def accepts_theta(self, theta: ArrayLike) -> NDArray[np.bool_]:
        values = np.asarray(theta, dtype=float)
        return (values >= self.theta_min) & (values <= self.theta_max)

    def path_length(self, theta: ArrayLike) -> NDArray[np.float64]:
        """NaI path length using the SNP paper's symmetric angle definition.

        Angles must lie inside the active angular interval. The returned array
        has the same shape as the supplied values, including for scalar input.
        """

        values = np.asarray(theta, dtype=float)
        if np.any(~np.isfinite(values)):
            raise ValueError("Angles must be finite.")
        if np.any(~self.accepts_theta(values)):
            raise ValueError("Angles must lie inside the detector acceptance.")

        theta_bar = np.minimum(values, np.pi - values)
        first_region = theta_bar < self.theta_cap
        with np.errstate(divide="raise", invalid="raise"):
            region_1 = self.h / np.cos(theta_bar) - self.r1 / np.sin(theta_bar)
            region_2 = (self.r2 - self.r1) / np.sin(theta_bar)
        result = np.where(first_region, region_1, region_2)

        roundoff_tolerance = 64.0 * np.finfo(float).eps * max(self.r2, self.h)
        if np.any(result < -roundoff_tolerance):
            raise ArithmeticError("Calculated a negative detector path length.")
        return np.maximum(result, 0.0)

    def to_dict(self) -> dict[str, float]:
        result = asdict(self)
        result.update(
            theta_min_rad=self.theta_min,
            theta_cap_rad=self.theta_cap,
            theta_max_rad=self.theta_max,
            theta_min_deg=float(np.degrees(self.theta_min)),
            theta_cap_deg=float(np.degrees(self.theta_cap)),
            theta_max_deg=float(np.degrees(self.theta_max)),
            solid_angle_sr=self.solid_angle_sr,
            fraction_of_full_sphere=self.fraction_of_full_sphere,
        )
        return result
