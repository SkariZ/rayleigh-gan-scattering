import numpy as np
from scipy import special

# ============================================================
# Wavevector definitions
# ============================================================

def k_medium(nm: float, wavelength: float) -> float:
    """
    Wave number in the medium:
        k = 2π n / λ

    wavelength and all lengths must share the same units (e.g. nm).
    """
    return 2.0 * np.pi * nm / wavelength


def q_from_theta(theta_rad, nm: float, wavelength: float):
    """
    Scattering vector magnitude (elastic scattering):

        q = 2 k sin(theta / 2)
          = (4π n / λ) sin(theta / 2)

    theta_rad MUST be in radians.
    """
    return 2.0 * k_medium(nm, wavelength) * np.sin(theta_rad / 2.0)


def q_from_theta_deg(theta_deg, nm: float, wavelength: float):
    """
    Explicit convenience wrapper for degrees -> radians.
    """
    return q_from_theta(np.deg2rad(theta_deg), nm, wavelength)


# ============================================================
# Numerically stable helpers
# ============================================================

def _sinc(x):
    """
    sin(x)/x with a stable limit at x -> 0.
    """
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    m = np.abs(x) > 1e-12
    out[m] = np.sin(x[m]) / x[m]
    return out


def _two_j1_over_x(x):
    """
    2*J1(x)/x with a stable limit at x -> 0.
    """
    x = np.asarray(x, dtype=float)
    out = np.ones_like(x)
    m = np.abs(x) > 1e-12
    out[m] = 2.0 * special.j1(x[m]) / x[m]
    return out


# ============================================================
# Form factors — AMPLITUDES (normalized: F(0)=1)
# ============================================================

def F_sphere_q(q, radius):
    """
    Sphere amplitude form factor:

        F(q) = 3 [sin(qR) - qR cos(qR)] / (qR)^3
    """
    x = np.asarray(q) * radius
    out = np.ones_like(x, dtype=float)

    m = np.abs(x) > 1e-8
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / (xm**3)
    return out


def F_thin_spherical_shell_q(q, radius):
    """
    Thin spherical surface ("skin") amplitude:

        F(q) = sin(qR)/(qR)
    """
    x = np.asarray(q) * radius
    return _sinc(x)


def F_thick_spherical_shell_q(q, r_inner, r_outer):
    """
    Finite-thickness spherical shell (difference of spheres),
    normalized by shell volume.
    """
    if r_outer <= r_inner:
        raise ValueError("r_outer must be larger than r_inner")

    Vo = r_outer**3
    Vi = r_inner**3

    Fo = F_sphere_q(q, r_outer)
    Fi = F_sphere_q(q, r_inner)

    return (Vo * Fo - Vi * Fi) / (Vo - Vi)


def F_disk_q(q_perp, radius):
    """
    Uniform 2D disk amplitude:

        F(q_perp) = 2 J1(q_perp R) / (q_perp R)

    q_perp is the component of q perpendicular to the disk.
    """
    x = np.asarray(q_perp) * radius
    return _two_j1_over_x(x)


# ============================================================
# Form factors — INTENSITIES
# ============================================================

def P_from_F(F):
    """
    Intensity form factor:
        P(q) = |F(q)|^2
    """
    F = np.asarray(F)
    return np.abs(F) ** 2


# ============================================================
# Polymer & aggregate models (INTENSITY, normalized)
# ============================================================

def P_gaussian_coil_q(q, Rg):
    """
    Debye form factor for an ideal Gaussian polymer coil:

        P(q) = 2 [exp(-u) + u - 1] / u^2
        u = (q Rg)^2
    """
    u = (np.asarray(q) * Rg) ** 2
    out = np.ones_like(u, dtype=float)

    m = u > 1e-12
    um = u[m]
    out[m] = 2.0 * (np.exp(-um) + um - 1.0) / (um**2)
    return out


def P_guinier_q(q, Rg):
    """
    Guinier approximation (small-q limit):

        P(q) ≈ exp(-q^2 Rg^2 / 3)

    Valid for q*Rg ≲ 1.
    """
    return np.exp(-(np.asarray(q) * Rg) ** 2 / 3.0)


# ============================================================
# θ-based convenience wrappers (radians)
# ============================================================

def F_sphere(theta_rad, radius, nm, wavelength):
    q = q_from_theta(theta_rad, nm, wavelength)
    return F_sphere_q(q, radius)


def F_thin_spherical_shell(theta_rad, radius, nm, wavelength):
    q = q_from_theta(theta_rad, nm, wavelength)
    return F_thin_spherical_shell_q(q, radius)


def F_thick_spherical_shell(theta_rad, r_inner, r_outer, nm, wavelength):
    q = q_from_theta(theta_rad, nm, wavelength)
    return F_thick_spherical_shell_q(q, r_inner, r_outer)


def F_disk(theta_rad, radius, nm, wavelength):
    """
    Uses q_perp ≈ |q| as a simple default.
    """
    q = q_from_theta(theta_rad, nm, wavelength)
    return F_disk_q(q, radius)


def P_gaussian_coil(theta_rad, Rg, nm, wavelength):
    q = q_from_theta(theta_rad, nm, wavelength)
    return P_gaussian_coil_q(q, Rg)


def P_guinier(theta_rad, Rg, nm, wavelength):
    q = q_from_theta(theta_rad, nm, wavelength)
    return P_guinier_q(q, Rg)


# ============================================================
# Geometry helpers
# ============================================================

def Rg_from_sphere_radius(radius):
    """
    Radius of gyration of a solid sphere:
        Rg = sqrt(3/5) R
    """
    return np.sqrt(3.0 / 5.0) * radius


def Rg_from_rod(length):
    """
    Radius of gyration of a thin rod:
        Rg = L / sqrt(12)
    """
    return length / np.sqrt(12.0)