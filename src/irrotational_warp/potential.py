from __future__ import annotations

import numpy as np


def tanh_wall(r: np.ndarray, rho: float, sigma: float) -> np.ndarray:
    """Smooth top-hat-like wall function f(r/rho) in [0, 1].

    f(0) ~ 1 (inside), f(r>>rho) ~ 0 (outside)
    """
    xi = r / rho
    return 0.5 * (1.0 + np.tanh(sigma * (1.0 - xi)))


def phi_dipole_cartesian(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    rho: float,
    sigma: float,
    v: float,
    eps: float = 1e-12,
) -> np.ndarray:
    """Rodal-like axisymmetric dipole potential oriented along +x.

    Φ = v * rho * f(r/rho) * cosθ, with cosθ = x/r.

    Parameters are geometric; v is dimensionless (v/c).
    """
    r = np.sqrt(x * x + y * y + z * z)
    f = tanh_wall(r, rho=rho, sigma=sigma)
    costheta = x / (r + eps)
    return v * rho * f * costheta


def g_rodal_exact(r: np.ndarray, rho: float, sigma: float) -> np.ndarray:
    """Exact g(r) function from Rodal et al (2025) Eq. 824.
    
    g(r) = [ 2rσ sinh(ρσ) + cosh(ρσ) ( ln cosh(σ(r-ρ)) - ln cosh(σ(r+ρ)) ) ] 
           / [ 2rσ sinh(ρσ) ]
           
    (Simplified from Eq 824 by using 2sinh(x/2)cosh(x/2)=sinh(x))
    """
    rho_sigma = rho * sigma
    sinh_rs = np.sinh(rho_sigma)
    cosh_rs = np.cosh(rho_sigma)
    
    # Numerator term 1: 2 * r * sigma * sinh(rho*sigma)
    term1 = 2 * r * sigma * sinh_rs
    
    # Numerator term 2: cosh(rho*sigma) * (ln_cosh_minus - ln_cosh_plus)
    # Use logaddexp for stable log(cosh(x)) = logaddexp(x, -x) - log(2)
    val_minus = sigma * (r - rho)
    val_plus  = sigma * (r + rho)
    
    log_cosh_minus = np.logaddexp(val_minus, -val_minus)
    log_cosh_plus  = np.logaddexp(val_plus, -val_plus)
    # The -log(2) terms cancel out in the difference
    
    log_diff = log_cosh_minus - log_cosh_plus
    term2 = cosh_rs * log_diff
    
    numerator = term1 + term2
    denominator = 2 * r * sigma * sinh_rs
    
    # Avoid division by zero at r=0
    # As r->0, g(r) -> 0.
    valid_mask = r > 1e-12
    res = np.zeros_like(r)
    res[valid_mask] = numerator[valid_mask] / denominator[valid_mask]
    
    return res


def phi_exact_rodal(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    rho: float,
    sigma: float,
    v: float,
    eps: float = 1e-12,
) -> np.ndarray:
    """Exact Rodal (2025) scalar potential.
    
    Φ = v * r * g(r) * cosθ
    """
    r = np.sqrt(x * x + y * y + z * z)
    g = g_rodal_exact(r, rho, sigma)
    costheta = x / (r + eps)
    return v * r * g * costheta
