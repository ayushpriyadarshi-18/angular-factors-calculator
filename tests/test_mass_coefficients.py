import unittest

from angularcorr.models import Coefficients, MassCoefficients


class MassCoefficientTests(unittest.TestCase):
    def test_converts_mass_coefficients_to_linear_and_back(self) -> None:
        mass = MassCoefficients(
            mu_over_rho_total=0.0936771,
            mu_over_rho_peak=0.0532134,
        )
        linear = mass.to_linear(3.67)
        self.assertAlmostEqual(linear.mu_total, 0.3437950, places=6)
        self.assertAlmostEqual(linear.mu_peak, 0.1952932, places=6)
        reconstructed = MassCoefficients.from_linear(linear, 3.67)
        self.assertEqual(reconstructed, mass)

    def test_rejects_nonpositive_values(self) -> None:
        with self.assertRaises(ValueError):
            MassCoefficients(mu_over_rho_total=0, mu_over_rho_peak=0.05)
        with self.assertRaises(ValueError):
            MassCoefficients.from_linear(Coefficients(0.3, 0.2), 0)


if __name__ == "__main__":
    unittest.main()
