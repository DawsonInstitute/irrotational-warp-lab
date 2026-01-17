"""Tests for tail correction module."""

import numpy as np

from irrotational_warp.tail import (
    compute_radial_average_z0,
    fit_power_law_decay,
    extrapolate_tail_integral_2d,
    compute_tail_correction,
)


def test_radial_average_uniform_field():
    """Radial average of a uniform field should be constant."""
    n = 51
    x = np.linspace(-10, 10, n)
    y = np.linspace(-10, 10, n)
    field = np.ones((n, n)) * 5.0

    r_bins, field_avg = compute_radial_average_z0(field, x, y, n_bins=20)

    # All bins should be ~ 5.0 (allow edge effects)
    assert np.allclose(field_avg, 5.0, atol=0.1)


def test_power_law_fitting():
    """Fit should recover known exponent for synthetic 1/r^4 decay."""
    r = np.linspace(5, 20, 50)
    amplitude_true = 100.0
    exponent_true = 4.0
    rho = amplitude_true / r**exponent_true

    exponent_fit, amplitude_fit, rms = fit_power_law_decay(
        r, rho, r_min=5.0, r_max=20.0
    )

    # Should recover exponent and amplitude within a few percent
    assert np.abs(exponent_fit - exponent_true) < 0.1
    assert np.abs(amplitude_fit - amplitude_true) / amplitude_true < 0.05
    assert rms < 0.01  # Good fit in log-space


def test_tail_integral_convergence():
    """Tail integral should converge for n > 2, diverge for n <= 2."""
    r_start = 10.0
    amplitude = 50.0

    # Convergent case: n = 4 > 2
    tail_pos, tail_neg = extrapolate_tail_integral_2d(
        amplitude=amplitude, exponent=4.0, r_start=r_start
    )
    assert tail_pos > 0  # Should have finite positive tail
    assert tail_neg == 0.0

    # Divergent case: n = 1.5 < 2
    tail_pos_div, tail_neg_div = extrapolate_tail_integral_2d(
        amplitude=amplitude, exponent=1.5, r_start=r_start
    )
    assert tail_pos_div == 0.0  # No correction for divergent case
    assert tail_neg_div == 0.0


def test_full_tail_correction_pipeline():
    """Full pipeline on synthetic field with known decay."""
    n = 101
    x = np.linspace(-20, 20, n)
    y = np.linspace(-20, 20, n)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2 + 1e-12)

    # Create field with 1/r^4 decay in far-field
    field = 200.0 / R**4

    result = compute_tail_correction(field, x, y, fit_r_min=12.0, n_bins=30)

    # Check fitted exponent is near 4
    assert 3.5 < result.exponent < 4.5

    # Tail integral should be positive and finite
    assert result.tail_integral_pos > 0
    assert result.tail_integral_neg == 0.0  # Positive field

    # Uncertainty should be reasonable (< tail integral)
    assert result.tail_uncertainty < result.tail_integral_pos
