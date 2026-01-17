import numpy as np

from irrotational_warp.adm import compute_slice_z0
from irrotational_warp.sweep import sweep_sigma_z0


def test_v_zero_gives_zero_density():
    res = compute_slice_z0(rho=10.0, sigma=5.0, v=0.0, extent=10.0, n=81)
    assert np.allclose(res.rho_adm, 0.0, atol=1e-10)


def test_y_symmetry_on_z0_slice():
    res = compute_slice_z0(rho=10.0, sigma=5.0, v=1.2, extent=10.0, n=101)
    flipped = np.flip(res.rho_adm, axis=0)
    assert np.allclose(res.rho_adm, flipped, atol=1e-6)


def test_sweep_has_expected_length_and_zero_v():
    """Test sigma sweep on z=0 slice.

    Note: uses n=41 (not 81) to avoid test timeouts.
    Full production sweeps in scripts/ can use higher resolution.
    """
    sigmas = np.linspace(1.0, 3.0, 5)
    pts = sweep_sigma_z0(rho=10.0, v=0.0, sigma_values=sigmas, extent=10.0, n=41)
    assert len(pts) == 5
    assert all(abs(p.e_pos) < 1e-10 and abs(p.e_neg) < 1e-10 and abs(p.e_net) < 1e-10 for p in pts)
