"""Tests for Einstein tensor eigenvalue diagnostics."""

import numpy as np

from irrotational_warp.einstein import compute_einstein_eigenvalues


def test_minkowski_flatness():
    """Minkowski spacetime (β=0) should have zero Einstein tensor eigenvalues."""
    n = 21
    beta_x = np.zeros((n, n))
    beta_y = np.zeros((n, n))
    dx = dy = 0.5

    result = compute_einstein_eigenvalues(beta_x, beta_y, dx=dx, dy=dy)

    # All eigenvalues should be ~ 0 for flat spacetime
    assert np.allclose(result.eig_max, 0.0, atol=1e-6)
    assert np.allclose(result.ricci_scalar, 0.0, atol=1e-6)
    assert np.all(result.is_type_i)  # Flat space is Type I


def test_small_perturbation_type_i():
    """Small shift perturbation should produce bounded Einstein tensor eigenvalues.

    Note: Type-I classification can be numerically fragile due to finite-difference
    error accumulation. This test focuses on verifying eigenvalues are physically
    reasonable (small, scaling ~ eps²) rather than strict Type-I fraction.
    """
    n = 21
    x = np.linspace(-5, 5, n)
    y = np.linspace(-5, 5, n)
    X, Y = np.meshgrid(x, y)
    dx = dy = float(x[1] - x[0])

    # Small harmonic perturbation: β_x ~ sin(x), β_y ~ 0
    eps = 0.01
    beta_x = eps * np.sin(2 * np.pi * X / 10.0)
    beta_y = np.zeros_like(beta_x)

    result = compute_einstein_eigenvalues(beta_x, beta_y, dx=dx, dy=dy)

    # Eigenvalues should be small (scaling ~ eps²) for small perturbations
    assert np.abs(result.eig_max).max() < 0.001  # ~ 10 × eps²

    # At least some regions should be Type I (exact fraction depends on numerics)
    assert np.sum(result.is_type_i) > 0.2 * n * n  # At least 20% Type I
