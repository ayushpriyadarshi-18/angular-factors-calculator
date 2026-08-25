import os
import unittest

from angularcorr.binning import AnalysisSettings, bin_events
from angularcorr.geometry import AnnularGeometry
from angularcorr.root_io import read_events


class ReferenceRootTests(unittest.TestCase):
    @unittest.skipUnless(
        os.environ.get("ANGULARCORR_REFERENCE_ROOT"),
        "ANGULARCORR_REFERENCE_ROOT is not set",
    )
    def test_reference_file_acceptance_and_schema(self) -> None:
        path = os.environ["ANGULARCORR_REFERENCE_ROOT"]
        inspection, events = read_events(path)
        self.assertTrue(inspection.compatible)
        self.assertEqual(inspection.entries, 1_000_000)
        self.assertEqual(events.theta_consistency_failures, 0)
        binned = bin_events(
            events,
            AnnularGeometry(r1=5.08, r2=10.16, h=7.62),
            AnalysisSettings(),
        )
        self.assertEqual(binned.accepted_events, 831_618)
        self.assertEqual(binned.total_interaction_events, 649_484)
        self.assertEqual(binned.full_energy_events, 494_474)


if __name__ == "__main__":
    unittest.main()
