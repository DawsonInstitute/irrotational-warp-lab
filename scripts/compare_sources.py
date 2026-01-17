#!/usr/bin/env python3
"""
Compare geometric "required" negative energy against toy positive-energy sources.

This is a PLAUSIBILITY study, NOT an engineering design.

Usage:
    python scripts/compare_sources.py --rho 5 --sigma 4 --v 1 --n 101 \
        --out results/sourcing/comparison.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from irrotational_warp.sourcing import (
    GaussianShellSource,
    UniformDiskSource,
    SmoothToroidalSource,
    plausibility_ratio,
)
from reproduce_rodal_exact import compute_energy_axisymmetric, _git_info


def compare_sources(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    nx: int,
    ny: int,
    backend: str,
    dtype: str,
):
    """Compare required energy against toy sources."""
    
    print("=" * 60)
    print("SOURCING COMPARISON STUDY")
    print("=" * 60)
    print(f"Warp parameters: rho={rho}, sigma={sigma}, v={v}")
    print(f"Grid: {nx}×{ny} (axisym)")
    print("=" * 60)
    
    # Compute required energy from warp geometry
    print("\nComputing warp energy requirement...")
    run = compute_energy_axisymmetric(
        rho=rho,
        sigma=sigma,
        v=v,
        extent=extent,
        nx=nx,
        ny=ny,
        backend=backend,
        dtype=dtype,
    )
    
    required_neg = run["energies"]["e_neg"]
    required_pos = run["energies"]["e_pos"]
    
    print(f"Required negative energy: {required_neg:.6e}")
    print(f"Required positive energy: {required_pos:.6e}")
    
    # Define toy sources
    sources = {
        "gaussian_shell_r5": GaussianShellSource(amplitude=1.0, r0=rho, sigma=rho/5),
        "gaussian_shell_r10": GaussianShellSource(amplitude=1.0, r0=rho, sigma=rho/10),
        "uniform_disk_thin": UniformDiskSource(amplitude=1.0, r_disk=rho, thickness=0.5),
        "uniform_disk_thick": UniformDiskSource(amplitude=1.0, r_disk=rho, thickness=2.0),
        "toroidal_narrow": SmoothToroidalSource(
            amplitude=1.0,
            major_radius=rho,
            major_width=rho/5,
            minor_width=rho/10,
        ),
        "toroidal_wide": SmoothToroidalSource(
            amplitude=1.0,
            major_radius=rho,
            major_width=rho/2,
            minor_width=rho/4,
        ),
    }
    
    # Evaluate source energies
    print("\nEvaluating toy source energies...")
    source_results = {}
    
    for name, source in sources.items():
        if isinstance(source, GaussianShellSource):
            energy = source.total_energy(r_max=extent, n_shells=1000)
        elif isinstance(source, UniformDiskSource):
            energy = source.total_energy()
        elif isinstance(source, SmoothToroidalSource):
            energy = source.total_energy(r_max=extent, z_max=extent/2, nr=200, nz=200)
        else:
            energy = 0.0
        
        ratio = plausibility_ratio(required_neg, energy)
        
        source_results[name] = {
            "type": source.__class__.__name__,
            "energy": energy,
            "plausibility_ratio": ratio,
            "params": {
                k: v for k, v in source.__dict__.items()
            },
        }
        
        print(f"  {name:25s}: E={energy:.6e}, ratio={ratio:.3f}")
    
    # Summary
    print("\n" + "=" * 60)
    print("PLAUSIBILITY ASSESSMENT")
    print("=" * 60)
    
    ratios = [r["plausibility_ratio"] for r in source_results.values()]
    best_name = max(source_results.keys(), key=lambda k: source_results[k]["plausibility_ratio"])
    worst_name = min(source_results.keys(), key=lambda k: source_results[k]["plausibility_ratio"])
    
    print(f"Ratio range: [{min(ratios):.3f}, {max(ratios):.3f}]")
    print(f"Best source: {best_name} (ratio={source_results[best_name]['plausibility_ratio']:.3f})")
    print(f"Worst source: {worst_name} (ratio={source_results[worst_name]['plausibility_ratio']:.3f})")
    
    print("\nInterpretation:")
    if max(ratios) > 10:
        print("  ✓ Some sources have excess capacity (ratio >> 1)")
    elif max(ratios) > 1:
        print("  ~ Some sources are marginally plausible (ratio ~ 1)")
    else:
        print("  ✗ All toy sources insufficient (ratio << 1)")
    
    print("\nCAVEAT: This is a CRUDE energy budget comparison.")
    print("It does NOT account for:")
    print("  - Spatial distribution mismatch")
    print("  - Full stress-energy tensor components")
    print("  - Dynamic effects / backreaction")
    print("  - Quantum field theory corrections")
    print("=" * 60)
    
    return {
        "warp": {
            "params": {"rho": rho, "sigma": sigma, "v": v},
            "energies": run["energies"],
        },
        "sources": source_results,
        "summary": {
            "ratio_min": min(ratios),
            "ratio_max": max(ratios),
            "best_source": best_name,
            "worst_source": worst_name,
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare warp energy vs sources")
    parser.add_argument("--backend", default="numpy", choices=["numpy", "cupy"])
    parser.add_argument("--dtype", default="float64", choices=["float64", "float32"])
    parser.add_argument("--rho", type=float, default=5.0)
    parser.add_argument("--sigma", type=float, default=4.0)
    parser.add_argument("--v", type=float, default=1.0)
    parser.add_argument("--extent-mult", type=float, default=12.0)
    parser.add_argument("--extent-pad", type=float, default=1.1)
    parser.add_argument("--nx", type=int, default=1200)
    parser.add_argument("--ny", type=int, default=600)
    parser.add_argument("--out", type=str, required=True)
    
    args = parser.parse_args(argv)
    
    extent = args.extent_mult * args.rho * args.extent_pad
    
    results = compare_sources(
        rho=args.rho,
        sigma=args.sigma,
        v=args.v,
        extent=extent,
        nx=args.nx,
        ny=args.ny,
        backend=args.backend,
        dtype=args.dtype,
    )
    
    payload = {
        "meta": {
            "timestamp_utc": np.datetime64("now", "s").astype(str),
            "git": _git_info(),
            "python": sys.version,
        },
        "comparison": results,
    }
    
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    
    print(f"\nResults written to: {args.out}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
