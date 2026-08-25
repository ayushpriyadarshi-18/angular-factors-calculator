import math
import unittest

import numpy as np

from angularcorr.geometry import AnnularGeometry


class AnnularGeometryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.geometry = AnnularGeometry(r1=5.08, r2=10.16, h=7.62)

    def test_published_angular_limits_and_solid_angle(self) -> None:
        self.assertAlmostEqual(math.degrees(self.geometry.theta_min), 33.6900675, places=6)
        self.assertAlmostEqual(math.degrees(self.geometry.theta_cap), 53.1301024, places=6)
        self.assertAlmostEqual(math.degrees(self.geometry.theta_max), 146.3099325, places=6)
        self.assertAlmostEqual(self.geometry.solid_angle_sr, 10.4558524, places=6)
        self.assertAlmostEqual(self.geometry.fraction_of_full_sphere, 0.8320502943, places=9)

    def test_path_is_centrosymmetric(self) -> None:
        theta = np.linspace(self.geometry.theta_min, np.pi / 2, 101)
        np.testing.assert_allclose(
            self.geometry.path_length(theta),
            self.geometry.path_length(np.pi - theta),
            rtol=0,
            atol=1e-13,
        )

    def test_path_is_continuous_at_cap(self) -> None:
        value = float(self.geometry.path_length(self.geometry.theta_cap))
        expected = (self.geometry.r2 - self.geometry.r1) / np.sin(
            self.geometry.theta_cap
        )
        self.assertAlmostEqual(value, expected, places=13)

    def test_invalid_dimensions_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            AnnularGeometry(r1=0, r2=10, h=5)
        with self.assertRaises(ValueError):
            AnnularGeometry(r1=5, r2=5, h=5)
        with self.assertRaises(ValueError):
            AnnularGeometry(r1=5, r2=10, h=-1)


if __name__ == "__main__":
    unittest.main()
