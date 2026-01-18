#!/usr/bin/env python3
"""Test Bayesian optimization against grid+NM method.

This script validates the Bayesian optimization implementation by:
1. Running both methods on the same problem
2. Verifying reproducibility (seeded runs)
3. Comparing convergence efficiency
"""
from __future__ import annotations

import json
from pathlib import Path

from irrotational_warp.optimize import optimize_bayesian, optimize_hybrid


def test_bayesian_reproducibility():
    """Test that seeded Bayesian runs are reproducible."""
    print("=" * 60)
    print("Test 1: Bayesian Reproducibility (Seeded Runs)")
    print("=" * 60)
    
    params = {
        "rho": 10.0,
        "sigma_range": (2.0, 8.0),
        "v_range": (0.8, 2.0),
        "extent": 20.0,
        "n": 41,  # Small for speed
        "n_calls": 20,
        "n_initial_points": 5,
        "random_state": 42,
    }
    
    print(f"Running Bayesian optimization with seed={params['random_state']}...")
    result1 = optimize_bayesian(**params)
    
    print("Re-running with same seed...")
    result2 = optimize_bayesian(**params)
    
    # Check reproducibility
    sigma_match = abs(result1.best_params["sigma"] - result2.best_params["sigma"]) < 1e-10
    v_match = abs(result1.best_params["v"] - result2.best_params["v"]) < 1e-10
    value_match = abs(result1.best_value - result2.best_value) < 1e-10
    
    print(f"\nRun 1: σ={result1.best_params['sigma']:.6f}, v={result1.best_params['v']:.6f}, |E⁻|={result1.best_value:.6e}")
    print(f"Run 2: σ={result2.best_params['sigma']:.6f}, v={result2.best_params['v']:.6f}, |E⁻|={result2.best_value:.6e}")
    print(f"\nReproducibility: {'✓ PASS' if (sigma_match and v_match and value_match) else '✗ FAIL'}")
    
    if sigma_match and v_match and value_match:
        print("  Results are bit-identical (expected for seeded runs)")
    else:
        print("  WARNING: Results differ despite identical seed!")
    
    return result1


def test_bayesian_vs_hybrid():
    """Compare Bayesian optimization efficiency against grid+NM."""
    print("\n" + "=" * 60)
    print("Test 2: Bayesian vs Hybrid (Convergence Efficiency)")
    print("=" * 60)
    
    common_params = {
        "rho": 10.0,
        "sigma_range": (2.0, 8.0),
        "v_range": (0.8, 2.0),
        "extent": 20.0,
        "n": 41,
    }
    
    # Hybrid (grid + Nelder-Mead)
    print("\nRunning hybrid optimization (5×5 grid + NM refinement)...")
    hybrid_result = optimize_hybrid(
        **common_params,
        sigma_steps=5,
        v_steps=5,
        refine=True,
    )
    
    # Bayesian (GP)
    print("Running Bayesian optimization (30 calls, 8 initial)...")
    bayes_result = optimize_bayesian(
        **common_params,
        n_calls=30,
        n_initial_points=8,
        random_state=123,
    )
    
    print(f"\n{'Method':<20} {'Best |E⁻|':<15} {'σ':<10} {'v':<10} {'Evals':<10}")
    print("-" * 65)
    print(f"{'Hybrid (Grid+NM)':<20} {hybrid_result.best_value:<15.6e} "
          f"{hybrid_result.best_params['sigma']:<10.4f} "
          f"{hybrid_result.best_params['v']:<10.4f} "
          f"{hybrid_result.n_evaluations:<10}")
    print(f"{'Bayesian (GP)':<20} {bayes_result.best_value:<15.6e} "
          f"{bayes_result.best_params['sigma']:<10.4f} "
          f"{bayes_result.best_params['v']:<10.4f} "
          f"{bayes_result.n_evaluations:<10}")
    
    # Analysis
    improvement_ratio = hybrid_result.best_value / bayes_result.best_value
    eval_ratio = hybrid_result.n_evaluations / bayes_result.n_evaluations
    
    print("\nAnalysis:")
    print(f"  Value ratio (hybrid/bayes): {improvement_ratio:.3f}")
    print(f"  Eval ratio (hybrid/bayes):  {eval_ratio:.2f}x")
    
    if bayes_result.best_value < hybrid_result.best_value:
        print(f"  → Bayesian found better optimum ({((1 - improvement_ratio) * 100):.1f}% lower |E⁻|)")
    else:
        print(f"  → Hybrid found better optimum ({((improvement_ratio - 1) * 100):.1f}% lower |E⁻|)")
    
    if bayes_result.n_evaluations < hybrid_result.n_evaluations:
        print(f"  → Bayesian used fewer evaluations ({eval_ratio:.1f}x reduction)")
    else:
        print("  → Hybrid used fewer evaluations")
    
    return hybrid_result, bayes_result


def test_bayesian_bounds():
    """Verify Bayesian optimization respects bounds."""
    print("\n" + "=" * 60)
    print("Test 3: Bounds Verification")
    print("=" * 60)
    
    sigma_min, sigma_max = 3.0, 6.0
    v_min, v_max = 1.0, 1.5
    
    print(f"Constraints: σ ∈ [{sigma_min}, {sigma_max}], v ∈ [{v_min}, {v_max}]")
    
    result = optimize_bayesian(
        rho=10.0,
        sigma_range=(sigma_min, sigma_max),
        v_range=(v_min, v_max),
        extent=20.0,
        n=41,
        n_calls=15,
        n_initial_points=5,
        random_state=456,
    )
    
    sigma_ok = sigma_min <= result.best_params["sigma"] <= sigma_max
    v_ok = v_min <= result.best_params["v"] <= v_max
    
    print("\nBest parameters:")
    print(f"  σ = {result.best_params['sigma']:.4f}  {'✓' if sigma_ok else '✗ OUT OF BOUNDS'}")
    print(f"  v = {result.best_params['v']:.4f}  {'✓' if v_ok else '✗ OUT OF BOUNDS'}")
    
    if sigma_ok and v_ok:
        print("\n✓ PASS: All parameters within specified bounds")
    else:
        print("\n✗ FAIL: Parameters violated bounds!")
    
    return result


def save_comparison_results(hybrid_result, bayes_result):
    """Save comparison results to JSON."""
    output_path = Path("results/optimization/bayesian_comparison.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    data = {
        "hybrid": {
            "best_params": hybrid_result.best_params,
            "best_value": hybrid_result.best_value,
            "n_evaluations": hybrid_result.n_evaluations,
            "method": hybrid_result.method,
        },
        "bayesian": {
            "best_params": bayes_result.best_params,
            "best_value": bayes_result.best_value,
            "n_evaluations": bayes_result.n_evaluations,
            "method": bayes_result.method,
        },
        "comparison": {
            "value_ratio_hybrid_over_bayes": hybrid_result.best_value / bayes_result.best_value,
            "eval_ratio_hybrid_over_bayes": hybrid_result.n_evaluations / bayes_result.n_evaluations,
        },
    }
    
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    
    print(f"\nSaved comparison to {output_path}")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("BAYESIAN OPTIMIZATION TEST SUITE")
    print("=" * 60)
    
    try:
        # Test 1: Reproducibility
        bayes_result = test_bayesian_reproducibility()
        
        # Test 2: Compare with hybrid
        hybrid_result, bayes_result = test_bayesian_vs_hybrid()
        
        # Test 3: Bounds
        test_bayesian_bounds()
        
        # Save comparison
        save_comparison_results(hybrid_result, bayes_result)
        
        print("\n" + "=" * 60)
        print("ALL TESTS COMPLETE")
        print("=" * 60)
        
    except ImportError as e:
        if "scikit-optimize" in str(e) or "skopt" in str(e):
            print("\n✗ ERROR: scikit-optimize not installed")
            print("Install with: pip install scikit-optimize")
        else:
            raise
