"""Tests for 2D (sigma, v) parameter sweep."""
import numpy as np
from irrotational_warp.sweep import sweep_2d_z0


def test_sweep_2d_output_structure():
    """Verify that sweep_2d_z0 returns expected structure."""
    sigma_values = np.array([2.0, 4.0])
    v_values = np.array([1.0, 1.5])
    
    points = sweep_2d_z0(
        rho=10.0,
        sigma_values=sigma_values,
        v_values=v_values,
        extent=10.0,
        n=31,
    )
    
    # Should have 2 × 2 = 4 points
    assert len(points) == 4
    
    # Each point should have all required fields
    for p in points:
        assert hasattr(p, "sigma")
        assert hasattr(p, "v")
        assert hasattr(p, "e_pos")
        assert hasattr(p, "e_neg")
        assert hasattr(p, "e_net")
        assert hasattr(p, "neg_fraction")
        
        # Energy values should be finite
        assert np.isfinite(p.e_pos)
        assert np.isfinite(p.e_neg)
        assert np.isfinite(p.e_net)
        
        # Positive energy should be positive
        assert p.e_pos > 0.0
        
        # Negative energy is stored as MAGNITUDE (positive value)
        assert p.e_neg > 0.0
        
        # Net energy should equal e_pos - e_neg
        assert np.isclose(p.e_net, p.e_pos - p.e_neg, rtol=1e-10)
        
        # Neg fraction should be in [0, 1]
        assert 0.0 <= p.neg_fraction <= 1.0


def test_sweep_2d_parameter_coverage():
    """Verify that all (sigma, v) pairs are covered."""
    sigma_values = np.array([2.0, 3.0, 4.0])
    v_values = np.array([1.0, 1.5])
    
    points = sweep_2d_z0(
        rho=10.0,
        sigma_values=sigma_values,
        v_values=v_values,
        extent=10.0,
        n=31,
    )
    
    # Extract unique (sigma, v) pairs
    pairs = {(p.sigma, p.v) for p in points}
    
    # Should have exactly 3 × 2 = 6 unique pairs
    assert len(pairs) == 6
    
    # Verify all expected pairs are present
    for sigma in sigma_values:
        for v in v_values:
            assert (sigma, v) in pairs
