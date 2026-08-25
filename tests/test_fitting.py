import unittest

import numpy as np

from angularcorr.fitting import fit_response_coefficient


class CoefficientFittingTests(unittest.TestCase):
    def test_recovers_coefficient_from_exact_expected_counts(self) -> None:
        paths = np.linspace(0.1, 8.0, 120)
        trials = np.full(paths.shape, 1_000_000, dtype=np.int64)
        expected_mu = 0.27
        probabilities = 1.0 - np.exp(-expected_mu * paths)
        successes = np.rint(trials * probabilities).astype(np.int64)
        fitted = fit_response_coefficient(
            paths, trials, successes, response="synthetic"
        )
        self.assertAlmostEqual(fitted.mu, expected_mu, places=6)
        self.assertGreater(fitted.standard_error, 0)


if __name__ == "__main__":
    unittest.main()
