import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rayleigh_gans.form_factors import (
    F_disk,
    F_disk_q,
    F_sphere,
    F_sphere_q,
    F_thick_spherical_shell_q,
    F_thin_spherical_shell_q,
    P_from_F,
    P_gaussian_coil,
    P_gaussian_coil_q,
    P_guinier_q,
    Rg_from_rod,
    Rg_from_sphere_radius,
    q_from_theta,
    q_from_theta_deg,
)


class FormFactorTests(unittest.TestCase):
    def test_angle_wrappers_match_q_space_functions(self):
        theta = np.linspace(0.0, 1.2, 9)
        nm = 1.33
        wavelength = 532.0
        radius = 175.0
        rg = 250.0
        q = q_from_theta(theta, nm, wavelength)

        np.testing.assert_allclose(F_sphere(theta, radius, nm, wavelength), F_sphere_q(q, radius))
        np.testing.assert_allclose(F_disk(theta, radius, nm, wavelength), F_disk_q(q, radius))
        np.testing.assert_allclose(
            P_gaussian_coil(theta, rg, nm, wavelength),
            P_gaussian_coil_q(q, rg),
        )

    def test_degree_wrapper_matches_radians(self):
        theta_deg = np.array([0.0, 15.0, 35.0, 70.0])
        nm = 1.4
        wavelength = 488.0

        np.testing.assert_allclose(
            q_from_theta_deg(theta_deg, nm, wavelength),
            q_from_theta(np.deg2rad(theta_deg), nm, wavelength),
        )

    def test_form_factors_are_normalized_at_zero_q(self):
        q0 = np.array([0.0])

        self.assertAlmostEqual(float(F_sphere_q(q0, 25.0)[0]), 1.0)
        self.assertAlmostEqual(float(F_thin_spherical_shell_q(q0, 25.0)[0]), 1.0)
        self.assertAlmostEqual(float(F_thick_spherical_shell_q(q0, 10.0, 25.0)[0]), 1.0)
        self.assertAlmostEqual(float(F_disk_q(q0, 25.0)[0]), 1.0)
        self.assertAlmostEqual(float(P_gaussian_coil_q(q0, 30.0)[0]), 1.0)
        self.assertAlmostEqual(float(P_guinier_q(q0, 30.0)[0]), 1.0)

    def test_invalid_shell_geometry_raises(self):
        with self.assertRaises(ValueError):
            F_thick_spherical_shell_q(np.array([0.0, 0.1]), 12.0, 12.0)

    def test_intensity_helper_handles_complex_amplitudes(self):
        amplitude = np.array([1 + 2j, -3j, 4.0])
        expected = np.array([5.0, 9.0, 16.0])
        np.testing.assert_allclose(P_from_F(amplitude), expected)

    def test_geometry_helpers_match_known_relations(self):
        self.assertAlmostEqual(Rg_from_sphere_radius(10.0), np.sqrt(3.0 / 5.0) * 10.0)
        self.assertAlmostEqual(Rg_from_rod(12.0), np.sqrt(12.0))


if __name__ == "__main__":
    unittest.main()
