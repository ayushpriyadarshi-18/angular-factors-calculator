import unittest

from angularcorr.geometry import AnnularGeometry
from angularcorr.integration import calculate_moments_midpoint, calculate_moments_quad
from angularcorr.models import Coefficients
from angularcorr.results import calculate_factors


class IntegrationTests(unittest.TestCase):
    def test_quad_and_full_interval_midpoint_agree(self) -> None:
        geometry = AnnularGeometry(r1=5.08, r2=10.16, h=7.62)
        coefficients = Coefficients(mu_total=0.3437988, mu_peak=0.1952923)
        adaptive = calculate_moments_quad(geometry, coefficients)
        midpoint = calculate_moments_midpoint(
            geometry, coefficients, intervals=400_000
        )
        for quad_value, midpoint_value in zip(
            adaptive.moments.to_dict().values(),
            midpoint.moments.to_dict().values(),
            strict=True,
        ):
            self.assertAlmostEqual(quad_value, midpoint_value, places=10)

    def test_factors_are_finite_and_greater_than_one(self) -> None:
        geometry = AnnularGeometry(r1=5.08, r2=10.16, h=7.62)
        coefficients = Coefficients(mu_total=0.3437988, mu_peak=0.1952923)
        factors = calculate_factors(calculate_moments_quad(geometry, coefficients).moments)
        self.assertGreater(factors.w_l1, 1.0)
        self.assertGreater(factors.w_l2, 1.0)
        self.assertGreater(factors.w_g, 1.0)


if __name__ == "__main__":
    unittest.main()
