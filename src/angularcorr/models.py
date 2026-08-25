"""Validated data objects shared by the scientific workflows."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class Coefficients:
    """Linear response coefficients in inverse centimetres."""

    mu_total: float
    mu_peak: float

    def __post_init__(self) -> None:
        if self.mu_total <= 0 or self.mu_peak <= 0:
            raise ValueError("Response coefficients must be positive.")

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class MassCoefficients:
    """Mass response coefficients in square centimetres per gram."""

    mu_over_rho_total: float
    mu_over_rho_peak: float

    def __post_init__(self) -> None:
        if self.mu_over_rho_total <= 0 or self.mu_over_rho_peak <= 0:
            raise ValueError("Mass response coefficients must be positive.")

    def to_linear(self, density: float) -> Coefficients:
        if density <= 0:
            raise ValueError("Density must be positive.")
        return Coefficients(
            mu_total=self.mu_over_rho_total * density,
            mu_peak=self.mu_over_rho_peak * density,
        )

    @classmethod
    def from_linear(cls, coefficients: Coefficients, density: float) -> "MassCoefficients":
        if density <= 0:
            raise ValueError("Density must be positive.")
        return cls(
            mu_over_rho_total=coefficients.mu_total / density,
            mu_over_rho_peak=coefficients.mu_peak / density,
        )

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class Moments:
    """Five angular averages needed by the three correlation factors."""

    mean_peak: float
    mean_total: float
    mean_peak_total: float
    mean_total_squared: float
    mean_peak_squared: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class MomentErrors:
    """Absolute quadrature error estimates corresponding to ``Moments``."""

    mean_peak: float
    mean_total: float
    mean_peak_total: float
    mean_total_squared: float
    mean_peak_squared: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class CorrelationFactors:
    w_l1: float
    w_l2: float
    w_g: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class IntegrationResult:
    method: str
    moments: Moments
    errors: MomentErrors | None = None

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "method": self.method,
            "moments": self.moments.to_dict(),
        }
        if self.errors is not None:
            result["errors"] = self.errors.to_dict()
        return result


@dataclass(frozen=True)
class CoefficientFit:
    response: str
    mu: float
    standard_error: float
    confidence_low: float
    confidence_high: float
    negative_log_likelihood: float
    bins_used: int
    events_used: int
    successes: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
