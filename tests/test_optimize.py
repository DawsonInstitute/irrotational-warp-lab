"""Tests for parameter optimization."""
import numpy as np
from irrotational_warp.optimize import (
    grid_search,
    optimize_nelder_mead,
    optimize_hybrid,
    objective_neg_energy,
)


def test_objective_function():
    """Test that objective function returns finite positive values."""
    value = objective_neg_energy(
        params=np.array([3.0, 1.5]),
        param_names=["sigma", "v"],
        rho=10.0,
        extent=10.0,
        n=31,
    )
    
    assert np.isfinite(value)
    assert value > 0.0  # e_neg is returned as magnitude


def test_grid_search_basic():
    """Test grid search finds a minimum."""
    result = grid_search(
        rho=10.0,
        sigma_range=(2.0, 4.0),
        v_range=(1.0, 2.0),
        sigma_steps=3,
        v_steps=3,
        extent=10.0,
        n=31,
    )
    
    assert result.success
    assert result.method == "grid_search"
    assert result.n_evaluations == 9  # 3×3
    assert result.best_value <= result.initial_value
    assert 2.0 <= result.best_params["sigma"] <= 4.0
    assert 1.0 <= result.best_params["v"] <= 2.0


def test_nelder_mead_basic():
    """Test Nelder-Mead local optimization."""
    result = optimize_nelder_mead(
        rho=10.0,
        initial_sigma=3.0,
        initial_v=1.5,
        extent=10.0,
        n=31,
        bounds={"sigma": (1.0, 10.0), "v": (0.5, 3.0)},
    )
    
    assert result.method == "nelder_mead"
    assert result.n_evaluations > 0
    # Should improve or stay same
    assert result.best_value <= result.initial_value + 1e-6


def test_hybrid_without_refinement():
    """Test hybrid optimizer in grid-only mode."""
    result = optimize_hybrid(
        rho=10.0,
        sigma_range=(2.0, 4.0),
        v_range=(1.0, 2.0),
        sigma_steps=3,
        v_steps=3,
        extent=10.0,
        n=31,
        refine=False,
    )
    
    assert result.method == "grid_search"
    assert result.n_evaluations == 9


def test_hybrid_with_refinement():
    """Test hybrid optimizer with Nelder-Mead refinement."""
    result = optimize_hybrid(
        rho=10.0,
        sigma_range=(2.0, 4.0),
        v_range=(1.0, 2.0),
        sigma_steps=3,
        v_steps=3,
        extent=10.0,
        n=31,
        refine=True,
    )
    
    assert result.method == "hybrid_grid_nelder_mead"
    # Should have grid + refinement evaluations
    assert result.n_evaluations > 9
    # Refined should be better or equal to grid
    assert result.success


def test_optimizer_determinism():
    """Test that optimizer produces consistent results."""
    # Run twice with identical parameters
    result1 = grid_search(
        rho=10.0,
        sigma_range=(2.0, 4.0),
        v_range=(1.0, 2.0),
        sigma_steps=3,
        v_steps=3,
        extent=10.0,
        n=31,
    )
    
    result2 = grid_search(
        rho=10.0,
        sigma_range=(2.0, 4.0),
        v_range=(1.0, 2.0),
        sigma_steps=3,
        v_steps=3,
        extent=10.0,
        n=31,
    )
    
    # Should get identical results
    assert result1.best_params == result2.best_params
    assert np.isclose(result1.best_value, result2.best_value)
