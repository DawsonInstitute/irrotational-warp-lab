
import numpy as np
import sys
from irrotational_warp.potential import phi_exact_rodal
import time
import time

def compute_energy_axisymmetric(rho, sigma, v, x_range, y_range, nx, ny):
    """
    Compute ADM energy density assuming axisymmetry around x-axis.
    y represents the cylindrical radius (rho_cyl).
    """
    x = np.linspace(x_range[0], x_range[1], nx)
    # y starts from slightly > 0 to avoid 1/y singularity, or handle carefully
    # paper uses spherical grid, we use cylindrical for integration convenience
    y = np.linspace(y_range[0], y_range[1], ny)
    
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    X, Y = np.meshgrid(x, y, indexing='xy')
    # Z is 0 on this slice
    Z = np.zeros_like(X)
    
    print("Computing Potential...")
    phi = phi_exact_rodal(X, Y, Z, rho=rho, sigma=sigma, v=v)
    
    print("Computing Derivatives...")
    # gradients
    # np.gradient returns [d/dy, d/dx] for 2D array with indexing='xy' (actually typically axis 0 is y, axis 1 is x)
    # if indexing='xy', shape is (ny, nx). axis 0 is y.
    
    grad_y, grad_x = np.gradient(phi, dy, dx, edge_order=2)
    
    # Second derivatives
    # Phi_xx
    _, phi_xx = np.gradient(grad_x, dy, dx, edge_order=2)
    
    # Phi_yy
    phi_yy, _ = np.gradient(grad_y, dy, dx, edge_order=2)
    
    # Phi_xy
    _, phi_xy = np.gradient(grad_y, dy, dx, edge_order=2)
    
    # K components
    # Shift beta = - grad Phi.
    # K_ij = - D_i D_j Phi  (for flat space)
    # Actually K_ij = 1/2 (D_i beta_j + D_j beta_i).
    # Since beta_i = - partial_i Phi, K_ij = - partial_i partial_j Phi.
    
    K_xx = -phi_xx
    K_yy = -phi_yy
    K_xy = -phi_xy
    
    # K_zz = - (1/y) * partial_y Phi
    # Handle y=0 limit: use K_yy
    with np.errstate(divide='ignore', invalid='ignore'):
         K_zz = - grad_y / Y
         
    # Fix singularity at y=0 (row 0 if y[0]==0)
    # If using linspace(0, ...), the first row is problematic.
    # We can replace the first row with K_yy[0, :]
    if y[0] == 0:
        K_zz[0, :] = K_yy[0, :]
    
    print("Computing Energy Density...")
    # Trace K = K_xx + K_yy + K_zz
    K_trace = K_xx + K_yy + K_zz
    
    # K_squared = sum K_ij K^ij
    # sum sum K_ij K_ij = K_xx^2 + K_yy^2 + K_zz^2 + 2*K_xy^2  (since K_yx = K_xy)
    # Note: K_xz = K_yz = 0 due to symmetry on the slice y>0, z=0
    K_sq = K_xx**2 + K_yy**2 + K_zz**2 + 2 * K_xy**2
    
    # rho = (K^2 - K_ab K^ab) / 16pi
    # Wait, paper Eq. 2290 says K^2 - K_ab K^ab = 2 lambda_H = -2 kappa rho_p
    # So rho_p = - (K^2 - K_ab K^ab) / (2 kappa).
    # Wait.
    # K^2 - K_ab K^ab = R (Scalar curvature of slice) + ...
    # Standard ADM constraint (alpha=1, flat space):
    # 16 pi rho = R + K^2 - Tr(K^2)
    # Since R=0 (flat space), 16 pi rho = (Tr K)^2 - Tr(K^2).
    # same as paper.
    
    rho_adm = (K_trace**2 - K_sq) / (16.0 * np.pi)
    
    # Integration
    # dV = 2 * pi * y * dx * dy
    dV = 2 * np.pi * Y * dx * dy
    
    E_density = rho_adm * dV
    
    # Summing up
    E_total = np.sum(E_density)
    
    # Separate Positive and Negative
    pos_mask = E_density > 0
    neg_mask = E_density < 0
    
    E_pos = np.sum(E_density[pos_mask])
    E_neg = np.sum(E_density[neg_mask])
    
    # Compute integrals for specific radii R < R_max
    # We use R = sqrt(x^2 + y^2)
    R_grid = np.sqrt(X**2 + Y**2)
    
    def get_E_at_R(limit_R):
        mask = R_grid <= limit_R
        e_dens_sub = E_density[mask]
        e_pos_sub = np.sum(e_dens_sub[e_dens_sub > 0])
        e_neg_sub = np.sum(e_dens_sub[e_dens_sub < 0])
        return e_pos_sub, e_neg_sub

    return E_pos, E_neg, E_total, X, Y, rho_adm, get_E_at_R

def main():
    rho = 5.0
    sigma = 4.0
    v = 1.0
    
    # Window: 12*rho = 60.
    limit = 12 * rho * 1.1 
    
    nx = 2400 
    ny = 1200
    
    print(f"Running Exact Rodal Simulation")
    print(f"Parameters: rho={rho}, sigma={sigma}, v={v}")
    print(f"Grid: {nx}x{ny}, limit={limit}")
    
    start_t = time.time()
    E_pos, E_neg, E_total, X, Y, rho_field, get_E_at_R = compute_energy_axisymmetric(
        rho, sigma, v, (-limit, limit), (0, limit), nx, ny
    )
    end_t = time.time()
    
    print(f"Computation Time: {end_t - start_t:.2f}s")
    
    # Tail Correction Analysis
    R1 = 8 * rho
    R2 = 12 * rho
    
    E_pos_1, E_neg_1 = get_E_at_R(R1)
    E_pos_2, E_neg_2 = get_E_at_R(R2)
    
    print("-" * 30)
    print(f"Energy at R1={R1:.1f} (8rho):")
    print(f"  E_pos: {E_pos_1:.6e}")
    print(f"  E_neg: {E_neg_1:.6e}")
    
    print(f"Energy at R2={R2:.1f} (12rho):")
    print(f"  E_pos: {E_pos_2:.6e}")
    print(f"  E_neg: {E_neg_2:.6e}")
    print(f"  Ratio: {E_pos_2/abs(E_neg_2):.4f} (Paper ~1.07)")
    print("-" * 30)
    
    # Two-point extrapolation
    # E_inf = E(R2) + (R1/(R2-R1)) * (E(R2) - E(R1))
    # Factor = R1 / (R2 - R1) = 40 / 20 = 2.0
    factor = R1 / (R2 - R1)
    
    E_pos_inf = E_pos_2 + factor * (E_pos_2 - E_pos_1)
    E_neg_inf = E_neg_2 + factor * (E_neg_2 - E_neg_1)
    
    print("Tail Extrapolated (Infinity):")
    print(f"  E_pos_inf: {E_pos_inf:.6e}")
    print(f"  E_neg_inf: {E_neg_inf:.6e}")
    
    E_net_inf = E_pos_inf + E_neg_inf
    E_abs_inf = E_pos_inf + abs(E_neg_inf)
    
    print("-" * 30)
    print(f"Net Energy (Inf): {E_net_inf:.6e}")
    print(f"Net / Abs (Inf): {abs(E_net_inf) / E_abs_inf:.6%} (Paper ~0.04%)")

if __name__ == "__main__":
    main()
