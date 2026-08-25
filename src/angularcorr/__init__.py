"""Scientific tools for annular-detector angular-correlation factors."""

from angularcorr.geometry import AnnularGeometry
from angularcorr.integration import calculate_moments_midpoint, calculate_moments_quad
from angularcorr.models import Coefficients, CorrelationFactors, MassCoefficients, Moments
from angularcorr.results import calculate_factors

__all__ = [
    "AnnularGeometry",
    "Coefficients",
    "CorrelationFactors",
    "MassCoefficients",
    "Moments",
    "calculate_factors",
    "calculate_moments_midpoint",
    "calculate_moments_quad",
]

__version__ = "0.1.0.dev0"
