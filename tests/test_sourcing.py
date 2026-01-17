"""Test sourcing models."""

import numpy as np
import pytest

from irrotational_warp.sourcing import (
    GaussianShellSource,
    UniformDiskSource,
    SmoothToroidalSource,
    plausibility_ratio,
)


def test_gaussian_shell_basic():
    """Test GaussianShellSource evaluation."""
    source = GaussianShellSource(amplitude=1.0, r0=5.0, sigma=1.0)
    
    # Peak at r = r0
    x = np.array([5.0, 0.0, 0.0])
    y = np.array([0.0, 0.0, 0.0])
    z = np.array([0.0, 0.0, 0.0])
    rho = source.evaluate(x, y, z)
    assert rho[0] == pytest.approx(1.0, abs=1e-10)
    
    # Zero far from shell
    x_far = np.array([100.0])
    rho_far = source.evaluate(x_far, np.zeros(1), np.zeros(1))
    assert rho_far[0] < 1e-10


def test_gaussian_shell_energy():
    """Test GaussianShellSource energy integration."""
    source = GaussianShellSource(amplitude=1.0, r0=5.0, sigma=1.0)
    energy = source.total_energy(r_max=20.0, n_shells=2000)
    
    # Should be positive and finite
    assert energy > 0
    assert np.isfinite(energy)
    
    # Analytic check: E ≈ A * 4π * r0² * sqrt(2π) * sigma for narrow shell
    expected = 1.0 * 4 * np.pi * 5.0**2 * np.sqrt(2 * np.pi) * 1.0
    assert energy == pytest.approx(expected, rel=0.1)


def test_uniform_disk_evaluation():
    """Test UniformDiskSource evaluation."""
    source = UniformDiskSource(amplitude=2.0, r_disk=3.0, thickness=1.0)
    
    # Inside disk
    x_in = np.array([1.0, 0.0])
    y_in = np.array([0.0, 0.0])
    z_in = np.array([0.0, 0.3])
    rho_in = source.evaluate(x_in, y_in, z_in)
    assert np.all(rho_in == 2.0)
    
    # Outside disk
    x_out = np.array([5.0, 0.0])
    y_out = np.array([0.0, 0.0])
    z_out = np.array([0.0, 2.0])
    rho_out = source.evaluate(x_out, y_out, z_out)
    assert np.all(rho_out == 0.0)


def test_uniform_disk_energy():
    """Test UniformDiskSource energy (analytic)."""
    source = UniformDiskSource(amplitude=2.0, r_disk=3.0, thickness=1.0)
    energy = source.total_energy()
    
    # V = π r² h
    expected = 2.0 * np.pi * 3.0**2 * 1.0
    assert energy == pytest.approx(expected, abs=1e-10)


def test_toroidal_source_basic():
    """Test SmoothToroidalSource evaluation."""
    source = SmoothToroidalSource(
        amplitude=1.0,
        major_radius=5.0,
        major_width=1.0,
        minor_width=0.5,
    )
    
    # Peak at R=R0, z=0
    x = np.array([5.0])
    y = np.array([0.0])
    z = np.array([0.0])
    rho = source.evaluate(x, y, z)
    assert rho[0] == pytest.approx(1.0, abs=1e-10)
    
    # Decays away from ring
    x_off = np.array([10.0])
    rho_off = source.evaluate(x_off, np.zeros(1), np.zeros(1))
    assert rho_off[0] < 0.1


def test_toroidal_energy():
    """Test SmoothToroidalSource energy integration."""
    source = SmoothToroidalSource(
        amplitude=1.0,
        major_radius=5.0,
        major_width=1.0,
        minor_width=0.5,
    )
    energy = source.total_energy(r_max=15.0, z_max=5.0, nr=300, nz=300)
    
    assert energy > 0
    assert np.isfinite(energy)


def test_plausibility_ratio():
    """Test plausibility ratio calculation."""
    # Excess capacity
    ratio = plausibility_ratio(required_energy=-10.0, source_energy=100.0)
    assert ratio == pytest.approx(10.0)
    
    # Marginal
    ratio = plausibility_ratio(required_energy=-50.0, source_energy=50.0)
    assert ratio == pytest.approx(1.0)
    
    # Insufficient
    ratio = plausibility_ratio(required_energy=-100.0, source_energy=10.0)
    assert ratio == pytest.approx(0.1)
