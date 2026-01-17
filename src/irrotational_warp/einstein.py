"""Einstein tensor eigenvalue diagnostics for invariant energy conditions.

This module implements Track A (direct metric → Einstein tensor) for 2D z=0 slices.
Computes:
- 4D metric g_μν from ADM variables (unit lapse, flat spatial metric)
- Christoffel symbols Γ^α_μν via finite differences
- Ricci tensor R_μν and scalar R
- Mixed Einstein tensor G^μ_ν
- Eigenvalues of G^μ_ν for Type-I classification and proper energy density ρ_p

Reference: Rodal (2025), Hawking & Ellis (1973) for Type classification.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fd import central_diff_2d


@dataclass(frozen=True)
class EinsteinResult:
    """Einstein tensor diagnostics on a 2D slice.

    Attributes:
        eig_max: Maximum real part of G^μ_ν eigenvalues (≈ proper energy density ρ_p)
        eig_all: All four eigenvalues of G^μ_ν at each grid point [4, ny, nx]
        ricci_scalar: Ricci scalar R [ny, nx]
        is_type_i: Boolean mask where all eigenvalues are real (Type I spacetime) [ny, nx]
    """

    eig_max: np.ndarray  # [ny, nx]
    eig_all: np.ndarray  # [4, ny, nx]
    ricci_scalar: np.ndarray  # [ny, nx]
    is_type_i: np.ndarray  # [ny, nx] boolean


def compute_metric_z0(
    beta_x: np.ndarray,
    beta_y: np.ndarray,
) -> np.ndarray:
    """Construct 4D covariant metric g_μν on z=0 slice.

    For unit lapse α=1, flat spatial metric γ_ij=δ_ij, shift β^i=(βx, βy, 0):
        g_00 = -α² + β^i β_i = -1 + (βx² + βy²)
        g_0i = β_i
        g_ij = γ_ij = δ_ij

    Returns:
        g_cov: 4D metric [4, 4, ny, nx] in (t,x,y,z) coords
    """
    ny, nx = beta_x.shape
    g = np.zeros((4, 4, ny, nx))

    # g_00 = -1 + β²
    g[0, 0] = -1.0 + beta_x**2 + beta_y**2

    # g_0i = β_i (mixed time-space components)
    g[0, 1] = beta_x
    g[1, 0] = beta_x
    g[0, 2] = beta_y
    g[2, 0] = beta_y

    # g_ij = δ_ij (flat spatial metric on z=0 slice)
    g[1, 1] = 1.0
    g[2, 2] = 1.0
    g[3, 3] = 1.0

    return g


def compute_metric_inverse(g_cov: np.ndarray) -> np.ndarray:
    """Invert 4D metric at each grid point.

    Args:
        g_cov: [4, 4, ny, nx]

    Returns:
        g_con: [4, 4, ny, nx] contravariant metric
    """
    ny, nx = g_cov.shape[2], g_cov.shape[3]
    g_con = np.zeros_like(g_cov)

    # Invert 4x4 matrix at each (i,j)
    for i in range(ny):
        for j in range(nx):
            g_con[:, :, i, j] = np.linalg.inv(g_cov[:, :, i, j])

    return g_con


def compute_christoffel(
    g_cov: np.ndarray,
    g_con: np.ndarray,
    dx: float,
    dy: float,
) -> np.ndarray:
    """Compute Christoffel symbols Γ^α_μν via finite differences.

    Γ^α_μν = (1/2) g^{ασ} (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)

    For 2D slice (z=0), only spatial derivatives in x,y (indices 1,2) are nonzero.
    Assumes stationarity (∂_t g_μν = 0) and z-independence (∂_z g_μν = 0).

    Args:
        g_cov: [4, 4, ny, nx]
        g_con: [4, 4, ny, nx]
        dx, dy: grid spacings

    Returns:
        Gamma: [4, 4, 4, ny, nx] Christoffel symbols Γ^α_μν
    """
    ny, nx = g_cov.shape[2], g_cov.shape[3]
    Gamma = np.zeros((4, 4, 4, ny, nx))

    # Compute ∂_μ g_αβ for spatial derivatives (μ=1,2 only)
    dg_dx = np.zeros_like(g_cov)  # ∂_x g_αβ
    dg_dy = np.zeros_like(g_cov)  # ∂_y g_αβ

    for alpha in range(4):
        for beta in range(4):
            dg_dx[alpha, beta], dg_dy[alpha, beta] = central_diff_2d(
                g_cov[alpha, beta], dx=dx, dy=dy
            )

    # Γ^α_μν = (1/2) g^{ασ} (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    # Loop over indices; use Einstein summation for contraction
    for mu in range(4):
        for nu in range(4):
            for alpha in range(4):
                # Spatial derivatives only (sigma=1,2 for x,y)
                for sigma in range(4):
                    # Derivative terms
                    if mu == 1:  # ∂_x
                        d_mu_g_sigma_nu = dg_dx[sigma, nu]
                    elif mu == 2:  # ∂_y
                        d_mu_g_sigma_nu = dg_dy[sigma, nu]
                    else:
                        d_mu_g_sigma_nu = 0.0

                    if nu == 1:
                        d_nu_g_sigma_mu = dg_dx[sigma, mu]
                    elif nu == 2:
                        d_nu_g_sigma_mu = dg_dy[sigma, mu]
                    else:
                        d_nu_g_sigma_mu = 0.0

                    if sigma == 1:
                        d_sigma_g_mu_nu = dg_dx[mu, nu]
                    elif sigma == 2:
                        d_sigma_g_mu_nu = dg_dy[mu, nu]
                    else:
                        d_sigma_g_mu_nu = 0.0

                    term = d_mu_g_sigma_nu + d_nu_g_sigma_mu - d_sigma_g_mu_nu
                    Gamma[alpha, mu, nu] += 0.5 * g_con[alpha, sigma] * term

    return Gamma


def compute_ricci_tensor(
    g_con: np.ndarray,
    Gamma: np.ndarray,
    dx: float,
    dy: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute Ricci tensor R_μν and scalar R.

    R_μν = ∂_σ Γ^σ_μν - ∂_ν Γ^σ_μσ + Γ^σ_ρσ Γ^ρ_μν - Γ^σ_ρν Γ^ρ_μσ

    This is the most expensive part; uses finite differences for ∂Γ terms.

    Returns:
        R_munu: Ricci tensor [4, 4, ny, nx]
        R: Ricci scalar [ny, nx]
    """
    ny, nx = Gamma.shape[3], Gamma.shape[4]
    R_munu = np.zeros((4, 4, ny, nx))

    # Compute ∂_σ Γ^σ_μν and ∂_ν Γ^σ_μσ
    for mu in range(4):
        for nu in range(4):
            # ∂_σ Γ^σ_μν summed over σ=1,2 (x,y derivatives)
            dGamma_term1 = 0.0
            for sigma in range(4):
                if sigma == 1:  # ∂_x Γ^σ_μν
                    dG, _ = central_diff_2d(Gamma[sigma, mu, nu], dx=dx, dy=dy)
                    dGamma_term1 += dG
                elif sigma == 2:  # ∂_y Γ^σ_μν
                    _, dG = central_diff_2d(Gamma[sigma, mu, nu], dx=dx, dy=dy)
                    dGamma_term1 += dG

            # ∂_ν Γ^σ_μσ summed over σ
            dGamma_term2 = 0.0
            for sigma in range(4):
                if nu == 1:
                    dG, _ = central_diff_2d(Gamma[sigma, mu, sigma], dx=dx, dy=dy)
                    dGamma_term2 += dG
                elif nu == 2:
                    _, dG = central_diff_2d(Gamma[sigma, mu, sigma], dx=dx, dy=dy)
                    dGamma_term2 += dG

            # Quadratic Γ terms: Γ^σ_ρσ Γ^ρ_μν - Γ^σ_ρν Γ^ρ_μσ
            quad_term = 0.0
            for sigma in range(4):
                for rho in range(4):
                    quad_term += (
                        Gamma[sigma, rho, sigma] * Gamma[rho, mu, nu]
                        - Gamma[sigma, rho, nu] * Gamma[rho, mu, sigma]
                    )

            R_munu[mu, nu] = dGamma_term1 - dGamma_term2 + quad_term

    # Ricci scalar R = g^{μν} R_μν
    R = np.zeros((ny, nx))
    for mu in range(4):
        for nu in range(4):
            R += g_con[mu, nu] * R_munu[mu, nu]

    return R_munu, R


def compute_einstein_tensor(
    g_con: np.ndarray,
    R_munu: np.ndarray,
    R: np.ndarray,
) -> np.ndarray:
    """Compute mixed Einstein tensor G^μ_ν.

    G^μ_ν = R^μ_ν - (1/2) δ^μ_ν R
    where R^μ_ν = g^{μσ} R_σν

    Returns:
        G_mixed: [4, 4, ny, nx] mixed Einstein tensor
    """
    ny, nx = R.shape
    # Raise first index of Ricci: R^μ_ν = g^{μσ} R_σν
    R_mixed = np.zeros((4, 4, ny, nx))
    for mu in range(4):
        for nu in range(4):
            for sigma in range(4):
                R_mixed[mu, nu] += g_con[mu, sigma] * R_munu[sigma, nu]

    # G^μ_ν = R^μ_ν - (1/2) δ^μ_ν R
    G_mixed = np.zeros_like(R_mixed)
    for mu in range(4):
        for nu in range(4):
            delta_mu_nu = 1.0 if mu == nu else 0.0
            G_mixed[mu, nu] = R_mixed[mu, nu] - 0.5 * delta_mu_nu * R

    return G_mixed


def compute_einstein_eigenvalues(
    beta_x: np.ndarray,
    beta_y: np.ndarray,
    dx: float,
    dy: float,
) -> EinsteinResult:
    """Compute Einstein tensor eigenvalues on z=0 slice.

    Full pipeline: metric → Christoffel → Ricci → Einstein → eigenvalues.

    This is expensive (O(n²) finite-diff ops on 4D tensors); use modest grid sizes.

    Returns:
        EinsteinResult with eigenvalues and Type-I classification.
    """
    ny, nx = beta_x.shape

    # 1. Construct metric
    g_cov = compute_metric_z0(beta_x, beta_y)
    g_con = compute_metric_inverse(g_cov)

    # 2. Christoffel symbols
    Gamma = compute_christoffel(g_cov, g_con, dx=dx, dy=dy)

    # 3. Ricci tensor and scalar
    R_munu, R = compute_ricci_tensor(g_con, Gamma, dx=dx, dy=dy)

    # 4. Mixed Einstein tensor
    G_mixed = compute_einstein_tensor(g_con, R_munu, R)

    # 5. Eigenvalues at each grid point
    eig_all = np.zeros((4, ny, nx), dtype=complex)
    eig_max = np.zeros((ny, nx))
    is_type_i = np.zeros((ny, nx), dtype=bool)

    for i in range(ny):
        for j in range(nx):
            # Extract 4x4 matrix G^μ_ν at (i,j)
            G_mat = G_mixed[:, :, i, j]
            eigs = np.linalg.eigvals(G_mat)
            eig_all[:, i, j] = eigs

            # Type I: all eigenvalues real (imaginary parts ~ 0 within tolerance)
            is_type_i[i, j] = np.all(np.abs(eigs.imag) < 1e-8)

            # Proper energy density ≈ max real part (Rodal convention)
            eig_max[i, j] = float(np.max(eigs.real))

    return EinsteinResult(
        eig_max=eig_max,
        eig_all=eig_all,
        ricci_scalar=R,
        is_type_i=is_type_i,
    )
