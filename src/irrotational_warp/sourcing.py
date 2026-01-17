"""
Simple sourcing models for plausibility studies.

These are TOY models for comparing geometric "required stress-energy" against
parameterized positive-energy sources. NOT for engineering claims.
"""

from __future__ import annotations

import numpy as np


class GaussianShellSource:
    """
    Gaussian shell source: energy density peaked on a spherical shell.
    
    ρ_source(r) = A exp(-(r-r0)²/(2σ²))
    
    This represents a simplified charged shell or plasma torus.
    """
    
    def __init__(self, *, amplitude: float, r0: float, sigma: float):
        """
        Parameters
        ----------
        amplitude : float
            Peak energy density
        r0 : float
            Shell radius
        sigma : float
            Shell thickness (Gaussian width)
        """
        self.amplitude = amplitude
        self.r0 = r0
        self.sigma = sigma
    
    def evaluate(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """Evaluate source energy density at given points."""
        r = np.sqrt(x * x + y * y + z * z)
        return self.amplitude * np.exp(-((r - self.r0) ** 2) / (2 * self.sigma ** 2))
    
    def total_energy(self, r_max: float, n_shells: int = 1000) -> float:
        """
        Compute total integrated energy via spherical shells.
        
        E = ∫ ρ(r) 4πr² dr
        """
        r_vals = np.linspace(0, r_max, n_shells)
        dr = r_vals[1] - r_vals[0]
        rho_vals = self.amplitude * np.exp(-((r_vals - self.r0) ** 2) / (2 * self.sigma ** 2))
        integrand = rho_vals * 4 * np.pi * r_vals ** 2
        return float(np.sum(integrand) * dr)


class UniformDiskSource:
    """
    Uniform disk source: constant energy density within a disk.
    
    ρ_source = A  for r_cyl < r_disk and |z| < h/2
    ρ_source = 0  otherwise
    
    This represents a simplified matter disk or coil geometry.
    """
    
    def __init__(self, *, amplitude: float, r_disk: float, thickness: float):
        """
        Parameters
        ----------
        amplitude : float
            Energy density within disk
        r_disk : float
            Disk radius (cylindrical)
        thickness : float
            Disk thickness (in z direction)
        """
        self.amplitude = amplitude
        self.r_disk = r_disk
        self.thickness = thickness
    
    def evaluate(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """Evaluate source energy density at given points."""
        r_cyl = np.sqrt(x * x + y * y)
        mask = (r_cyl <= self.r_disk) & (np.abs(z) <= self.thickness / 2)
        rho = np.zeros_like(x)
        rho[mask] = self.amplitude
        return rho
    
    def total_energy(self) -> float:
        """Compute total integrated energy (analytic)."""
        volume = np.pi * self.r_disk ** 2 * self.thickness
        return self.amplitude * volume


class SmoothToroidalSource:
    """
    Smooth toroidal (ring) source: Gaussian in both major and minor radii.
    
    ρ_source(R, z) = A exp(-(R-R0)²/(2σ_R²)) exp(-z²/(2σ_z²))
    
    where R = sqrt(x² + y²) is cylindrical radius.
    
    This represents a plasma torus or distributed coil system.
    """
    
    def __init__(
        self,
        *,
        amplitude: float,
        major_radius: float,
        major_width: float,
        minor_width: float,
    ):
        """
        Parameters
        ----------
        amplitude : float
            Peak energy density
        major_radius : float
            Major radius R0 (ring center)
        major_width : float
            Major radius Gaussian width σ_R
        minor_width : float
            Minor radius Gaussian width σ_z
        """
        self.amplitude = amplitude
        self.major_radius = major_radius
        self.major_width = major_width
        self.minor_width = minor_width
    
    def evaluate(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """Evaluate source energy density at given points."""
        r_cyl = np.sqrt(x * x + y * y)
        radial_factor = np.exp(-((r_cyl - self.major_radius) ** 2) / (2 * self.major_width ** 2))
        vertical_factor = np.exp(-(z ** 2) / (2 * self.minor_width ** 2))
        return self.amplitude * radial_factor * vertical_factor
    
    def total_energy(self, r_max: float, z_max: float, nr: int = 200, nz: int = 200) -> float:
        """
        Compute total integrated energy via cylindrical grid.
        
        E = ∫∫ ρ(R,z) 2πR dR dz
        """
        r_vals = np.linspace(0, r_max, nr)
        z_vals = np.linspace(-z_max, z_max, nz)
        dr = r_vals[1] - r_vals[0]
        dz = z_vals[1] - z_vals[0]
        
        R, Z = np.meshgrid(r_vals, z_vals, indexing="ij")
        rho = self.evaluate(R, np.zeros_like(R), Z)
        integrand = rho * 2 * np.pi * R
        
        return float(np.sum(integrand) * dr * dz)


def plausibility_ratio(
    required_energy: float,
    source_energy: float,
) -> float:
    """
    Compute plausibility ratio: source_energy / |required_energy|.
    
    Ratio >> 1: source may be plausible (excess capacity)
    Ratio ~ 1: marginal plausibility
    Ratio << 1: source insufficient
    
    This is a VERY crude metric and does not account for:
    - Spatial distribution mismatch
    - Stress-energy tensor components (not just ρ)
    - Dynamic effects / backreaction
    """
    return source_energy / abs(required_energy)
