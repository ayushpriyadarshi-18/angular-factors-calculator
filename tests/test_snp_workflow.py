import os
import unittest

from angularcorr.geometry import AnnularGeometry
from angularcorr.workflows import analyze_root


class SnpWorkflowTests(unittest.TestCase):
    @unittest.skipUnless(
        os.environ.get("ANGULARCORR_REFERENCE_ROOT"),
        "ANGULARCORR_REFERENCE_ROOT is not set",
    )
    def test_complete_snp_regression(self) -> None:
        result, _ = analyze_root(
            os.environ["ANGULARCORR_REFERENCE_ROOT"],
            AnnularGeometry(r1=5.08, r2=10.16, h=7.62),
            density=3.67,
            midpoint_intervals=200_000,
        )
        summary = result["event_summary"]
        self.assertEqual(summary["emitted_events"], 1_000_000)
        self.assertEqual(summary["accepted_events"], 831_618)
        self.assertEqual(summary["total_interaction_events"], 649_484)
        self.assertEqual(summary["full_energy_events"], 494_474)

        coefficients = result["coefficients"]["linear_cm_inverse"]
        self.assertAlmostEqual(coefficients["mu_total"], 0.3437951, places=6)
        self.assertAlmostEqual(coefficients["mu_peak"], 0.1952932, places=6)
        mass_coefficients = result["coefficients"]["mass_cm2_g"]
        self.assertAlmostEqual(
            mass_coefficients["mu_over_rho_total"], 0.0936771, places=6
        )
        self.assertAlmostEqual(
            mass_coefficients["mu_over_rho_peak"], 0.0532134, places=6
        )
        self.assertEqual(result["material"]["density_g_cm3"], 3.67)

        moments = result["integration"]["quad"]["moments"]
        self.assertAlmostEqual(moments["mean_peak"], 0.4953246801, places=7)
        self.assertAlmostEqual(moments["mean_total"], 0.6497108638, places=7)
        self.assertAlmostEqual(moments["mean_peak_total"], 0.4057288992, places=7)
        self.assertAlmostEqual(moments["mean_total_squared"], 0.5292696853, places=7)
        self.assertAlmostEqual(moments["mean_peak_squared"], 0.3114398906, places=7)

        factors = result["factors"]
        self.assertAlmostEqual(factors["w_l1"], 1.260741, places=6)
        self.assertAlmostEqual(factors["w_l2"], 1.253825, places=6)
        self.assertAlmostEqual(factors["w_g"], 1.269388, places=6)
        self.assertLess(result["integration"]["maximum_moment_difference"], 1e-10)


if __name__ == "__main__":
    unittest.main()
