"""Optimization tools for finding minimal negative energy configurations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import minimize

from .adm import compute_slice_z0, integrate_signed


@dataclass(frozen=True)
class OptimizationResult:
    """Result from parameter optimization."""
    best_params: dict[str, float]
    best_value: float
    initial_params: dict[str, float]
    initial_value: float
    n_evaluations: int
    success: bool
    method: str
    message: str


def objective_neg_energy(
    params: np.ndarray,
    param_names: list[str],
    rho: float,
    extent: float,
    n: int,
) -> float:
    """Objective function: minimize |E⁻| (magnitude of negative energy integral).
    
    Args:
        params: Array of parameter values [sigma, v, ...]
        param_names: Names corresponding to params (e.g., ["sigma", "v"])
        rho: Fixed bubble radius
        extent: Grid half-width
        n: Grid resolution
    
    Returns:
        Magnitude of negative energy integral (to minimize)
    """
    param_dict = dict(zip(param_names, params))
    
    res = compute_slice_z0(
        rho=rho,
        sigma=param_dict.get("sigma", 3.0),
        v=param_dict.get("v", 1.5),
        extent=extent,
        n=n,
    )
    
    _, e_neg, _ = integrate_signed(res.rho_adm, dx=res.dx, dy=res.dy)
    
    # Return magnitude (e_neg is already positive from integrate_signed)
    return float(e_neg)


def grid_search(
    *,
    rho: float,
    sigma_range: tuple[float, float],
    v_range: tuple[float, float],
    sigma_steps: int,
    v_steps: int,
    extent: float,
    n: int,
) -> OptimizationResult:
    """Exhaustive grid search over parameter space.
    
    Args:
        rho: Fixed bubble radius
        sigma_range: (min, max) for sigma parameter
        v_range: (min, max) for v parameter
        sigma_steps: Number of sigma grid points
        v_steps: Number of v grid points
        extent: Grid half-width
        n: Spatial grid resolution
    
    Returns:
        OptimizationResult with best parameters from grid
    """
    sigma_vals = np.linspace(sigma_range[0], sigma_range[1], sigma_steps)
    v_vals = np.linspace(v_range[0], v_range[1], v_steps)
    
    best_value = float('inf')
    best_sigma = sigma_vals[0]
    best_v = v_vals[0]
    
    n_evals = 0
    
    for sigma in sigma_vals:
        for v in v_vals:
            value = objective_neg_energy(
                params=np.array([sigma, v]),
                param_names=["sigma", "v"],
                rho=rho,
                extent=extent,
                n=n,
            )
            n_evals += 1
            
            if value < best_value:
                best_value = value
                best_sigma = sigma
                best_v = v
    
    return OptimizationResult(
        best_params={"sigma": best_sigma, "v": best_v},
        best_value=best_value,
        initial_params={"sigma": sigma_vals[0], "v": v_vals[0]},
        initial_value=objective_neg_energy(
            params=np.array([sigma_vals[0], v_vals[0]]),
            param_names=["sigma", "v"],
            rho=rho,
            extent=extent,
            n=n,
        ),
        n_evaluations=n_evals,
        success=True,
        method="grid_search",
        message=f"Evaluated {n_evals} grid points",
    )


def optimize_nelder_mead(
    *,
    rho: float,
    initial_sigma: float,
    initial_v: float,
    extent: float,
    n: int,
    bounds: dict[str, tuple[float, float]] | None = None,
) -> OptimizationResult:
    """Local optimization using Nelder-Mead simplex method.
    
    Args:
        rho: Fixed bubble radius
        initial_sigma: Starting sigma value
        initial_v: Starting v value
        extent: Grid half-width
        n: Spatial grid resolution
        bounds: Optional bounds dict {"sigma": (min, max), "v": (min, max)}
    
    Returns:
        OptimizationResult with optimized parameters
    """
    x0 = np.array([initial_sigma, initial_v])
    param_names = ["sigma", "v"]
    
    # Nelder-Mead doesn't support bounds directly, so we use a penalty method
    if bounds is not None:
        sigma_bounds = bounds.get("sigma", (0.1, 100.0))
        v_bounds = bounds.get("v", (0.01, 10.0))
        
        def bounded_objective(params):
            sigma, v = params
            # Add large penalty if outside bounds
            penalty = 0.0
            if sigma < sigma_bounds[0] or sigma > sigma_bounds[1]:
                penalty += 1e10
            if v < v_bounds[0] or v > v_bounds[1]:
                penalty += 1e10
            
            if penalty > 0:
                return penalty
            
            return objective_neg_energy(params, param_names, rho, extent, n)
        
        obj_func = bounded_objective
    else:
        obj_func = lambda params: objective_neg_energy(params, param_names, rho, extent, n)
    
    # Track number of evaluations
    n_evals = [0]
    def counting_objective(params):
        n_evals[0] += 1
        return obj_func(params)
    
    initial_value = obj_func(x0)
    
    result = minimize(
        counting_objective,
        x0,
        method='Nelder-Mead',
        options={'maxiter': 100, 'xatol': 1e-4, 'fatol': 1e-6}
    )
    
    return OptimizationResult(
        best_params={"sigma": result.x[0], "v": result.x[1]},
        best_value=result.fun,
        initial_params={"sigma": initial_sigma, "v": initial_v},
        initial_value=initial_value,
        n_evaluations=n_evals[0],
        success=result.success,
        method="nelder_mead",
        message=result.message,
    )


def optimize_hybrid(
    *,
    rho: float,
    sigma_range: tuple[float, float],
    v_range: tuple[float, float],
    sigma_steps: int,
    v_steps: int,
    extent: float,
    n: int,
    refine: bool = True,
) -> OptimizationResult:
    """Hybrid optimization: grid search followed by local refinement.
    
    Args:
        rho: Fixed bubble radius
        sigma_range: (min, max) for sigma parameter
        v_range: (min, max) for v parameter
        sigma_steps: Number of sigma grid points for initial search
        v_steps: Number of v grid points for initial search
        extent: Grid half-width
        n: Spatial grid resolution
        refine: If True, apply Nelder-Mead refinement after grid search
    
    Returns:
        OptimizationResult with best parameters (refined if refine=True)
    """
    # Phase 1: Grid search
    grid_result = grid_search(
        rho=rho,
        sigma_range=sigma_range,
        v_range=v_range,
        sigma_steps=sigma_steps,
        v_steps=v_steps,
        extent=extent,
        n=n,
    )
    
    if not refine:
        return grid_result
    
    # Phase 2: Local refinement from grid optimum
    refined_result = optimize_nelder_mead(
        rho=rho,
        initial_sigma=grid_result.best_params["sigma"],
        initial_v=grid_result.best_params["v"],
        extent=extent,
        n=n,
        bounds={"sigma": sigma_range, "v": v_range},
    )
    
    # Combine results
    return OptimizationResult(
        best_params=refined_result.best_params,
        best_value=refined_result.best_value,
        initial_params=grid_result.initial_params,
        initial_value=grid_result.initial_value,
        n_evaluations=grid_result.n_evaluations + refined_result.n_evaluations,
        success=refined_result.success,
        method="hybrid_grid_nelder_mead",
        message=f"Grid: {grid_result.n_evaluations} evals, Refinement: {refined_result.n_evaluations} evals",
    )
