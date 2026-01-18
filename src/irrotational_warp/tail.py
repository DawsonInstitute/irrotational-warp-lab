"""Tail correction for finite-box error control.

Addresses the issue that energy integrals on finite grids miss the far-field tail.
For potentials decaying as ~1/r^n, we can fit the radial average ⟨ρ⟩(r) in the
far-field and analytically extrapolate to infinity.

Reference: Standard practice in numerical relativity for finite-box corrections.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TailCorrectionResult:
    """Results from tail fitting and extrapolation.

    Attributes:
        r_bins: Radial bin centers [n_bins]
        rho_avg: Angle-averaged density ⟨ρ⟩(r) [n_bins]
        fit_r_min: Inner radius of fit region (start of far-field)
        fit_r_max: Outer radius of fit region (edge of grid)
        exponent: Fitted power-law exponent n (ρ ~ 1/r^n)
        amplitude: Fitted amplitude A (ρ ~ A/r^n)
        tail_integral_pos: Extrapolated E⁺_tail from R to infinity
        tail_integral_neg: Extrapolated |E⁻_tail| from R to infinity
        tail_uncertainty: Estimated uncertainty in tail integrals
        fit_residual_rms: RMS residual of fit in log-log space
    """

    r_bins: np.ndarray
    rho_avg: np.ndarray
    fit_r_min: float
    fit_r_max: float
    exponent: float
    amplitude: float
    tail_integral_pos: float
    tail_integral_neg: float
    tail_uncertainty: float
    fit_residual_rms: float


def compute_radial_average_z0(
    field: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    n_bins: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute angle-averaged radial profile ⟨field⟩(r) for z=0 slice.

    Args:
        field: 2D field to average [ny, nx]
        x: 1D x coordinates [nx]
        y: 1D y coordinates [ny]
        n_bins: Number of radial bins

    Returns:
        r_bins: Radial bin centers
        field_avg: Angle-averaged field values at each r
    """
    ny, nx = field.shape
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)

    r_max = np.max(R)
    r_edges = np.linspace(0, r_max, n_bins + 1)
    r_bins = 0.5 * (r_edges[:-1] + r_edges[1:])

    field_avg = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (R >= r_edges[i]) & (R < r_edges[i + 1])
        if np.sum(mask) > 0:
            field_avg[i] = np.mean(field[mask])
        else:
            field_avg[i] = 0.0

    return r_bins, field_avg


def fit_power_law_decay(
    r: np.ndarray,
    rho: np.ndarray,
    r_min: float,
    r_max: float,
) -> tuple[float, float, float]:
    """Fit ρ(r) ~ A/r^n in the range [r_min, r_max].

    Uses log-log linear regression: log(|ρ|) = log(|A|) - n*log(r)

    Args:
        r: Radial bin centers
        rho: Radial average of field
        r_min: Start of fit region (far-field onset)
        r_max: End of fit region (grid boundary)

    Returns:
        exponent: Power-law exponent n
        amplitude: Amplitude A
        rms_residual: RMS residual in log-log space
    """
    # Select fit region
    mask = (r >= r_min) & (r <= r_max) & (np.abs(rho) > 1e-12)
    if np.sum(mask) < 3:
        # Not enough points for fit
        return 0.0, 0.0, np.inf

    r_fit = r[mask]
    rho_fit = rho[mask]

    # Log-log linear fit
    log_r = np.log(r_fit)
    log_rho_abs = np.log(np.abs(rho_fit))

    # Linear regression: y = a + b*x → log|ρ| = log|A| - n*log(r)
    coeffs = np.polyfit(log_r, log_rho_abs, deg=1)
    exponent = -coeffs[0]  # -slope = n
    amplitude = np.exp(coeffs[1])  # intercept = log|A|

    # Sign of amplitude from average sign of rho in fit region
    if np.mean(rho_fit) < 0:
        amplitude = -amplitude

    # Compute residuals
    log_rho_pred = coeffs[1] + coeffs[0] * log_r
    residuals = log_rho_abs - log_rho_pred
    rms_residual = float(np.sqrt(np.mean(residuals**2)))

    return float(exponent), float(amplitude), rms_residual


def extrapolate_tail_integral_2d(
    amplitude: float,
    exponent: float,
    r_start: float,
) -> tuple[float, float]:
    """Analytically integrate A/r^n from r_start to infinity (2D area integral).

    For 2D slice: dA = r*dr*dθ → ∫ρ dA = 2π ∫_R^∞ (A/r^n) * r dr = 2πA ∫_R^∞ r^(1-n) dr

    Convergence requires n > 2 for 2D (r^(1-n) integrable as r→∞).

    Returns:
        integral_pos: Tail contribution if A > 0
        integral_neg: Tail contribution if A < 0 (as positive magnitude)
    """
    # Check convergence
    if exponent <= 2.0:
        # Integral diverges or grows; no tail correction
        return 0.0, 0.0

    # ∫_R^∞ r^(1-n) dr = [r^(2-n)/(2-n)]_R^∞ = -R^(2-n)/(2-n) for n > 2
    tail_1d = -r_start ** (2.0 - exponent) / (2.0 - exponent)
    tail = 2.0 * np.pi * amplitude * tail_1d

    if tail > 0:
        return float(tail), 0.0
    else:
        return 0.0, float(-tail)


def estimate_tail_uncertainty(
    amplitude: float,
    exponent: float,
    r_start: float,
    rms_residual: float,
) -> float:
    """Estimate uncertainty in tail integral from fit residuals.

    Uses linear propagation of RMS residual in log-space to amplitude uncertainty,
    then propagates to integral uncertainty.

    Returns:
        uncertainty: Estimated 1-sigma uncertainty in tail integral
    """
    # Amplitude uncertainty from log-space RMS (crude approximation)
    delta_log_A = rms_residual
    delta_A = np.abs(amplitude) * delta_log_A

    # Propagate to integral (linear approximation)
    if exponent <= 2.0:
        return 0.0

    tail_1d = -r_start ** (2.0 - exponent) / (2.0 - exponent)
    delta_tail = 2.0 * np.pi * delta_A * np.abs(tail_1d)

    return float(delta_tail)


def compute_tail_correction(
    field: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    fit_r_min: float | None = None,
    n_bins: int = 50,
) -> TailCorrectionResult:
    """Compute full tail correction pipeline for a 2D z=0 slice field.

    Args:
        field: 2D field (e.g., rho_adm) [ny, nx]
        x: 1D x coordinates [nx]
        y: 1D y coordinates [ny]
        fit_r_min: Inner radius for fit region (default: 0.7 * r_max)
        n_bins: Number of radial bins for averaging

    Returns:
        TailCorrectionResult with fitted parameters and tail integrals
    """
    # 1. Radial averaging
    r_bins, rho_avg = compute_radial_average_z0(field, x, y, n_bins=n_bins)

    # 2. Auto-select fit region if not provided
    r_max = r_bins[-1]
    if fit_r_min is None:
        fit_r_min = 0.7 * r_max  # Start fit at 70% of domain radius

    fit_r_max = r_max

    # 3. Fit power-law decay
    exponent, amplitude, rms_residual = fit_power_law_decay(
        r_bins, rho_avg, r_min=fit_r_min, r_max=fit_r_max
    )

    # 4. Extrapolate tail
    tail_pos, tail_neg = extrapolate_tail_integral_2d(amplitude, exponent, r_start=fit_r_max)

    # 5. Estimate uncertainty
    uncertainty = estimate_tail_uncertainty(amplitude, exponent, fit_r_max, rms_residual)

    return TailCorrectionResult(
        r_bins=r_bins,
        rho_avg=rho_avg,
        fit_r_min=fit_r_min,
        fit_r_max=fit_r_max,
        exponent=exponent,
        amplitude=amplitude,
        tail_integral_pos=tail_pos,
        tail_integral_neg=tail_neg,
        tail_uncertainty=uncertainty,
        fit_residual_rms=rms_residual,
    )
