#!/usr/bin/env python3
"""
Profile 3D Hessian computation to identify optimization opportunities.

Usage:
    python scripts/profile_3d.py --n 80 --profile-out results/profiling/profile_n80.prof
"""

import argparse
import cProfile
import pstats
import sys
import time
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from reproduce_rodal_exact import compute_energy_3d


def profile_3d_run(n: int, rho: float, sigma: float, v: float, backend: str, dtype: str):
    """Profile a single 3D run."""
    extent = 12.0 * rho * 1.1
    
    print("=" * 60)
    print(f"PROFILING 3D RUN (n={n})")
    print("=" * 60)
    
    # Warm-up run
    print("Warm-up run...")
    t0 = time.perf_counter()
    _ = compute_energy_3d(
        rho=rho, sigma=sigma, v=v, extent=extent, n=n,
        backend=backend, dtype=dtype
    )
    warmup_time = time.perf_counter() - t0
    print(f"Warm-up completed in {warmup_time:.2f}s")
    
    # Profiled run
    print("\nProfiled run...")
    profiler = cProfile.Profile()
    profiler.enable()
    
    t0 = time.perf_counter()
    result = compute_energy_3d(
        rho=rho, sigma=sigma, v=v, extent=extent, n=n,
        backend=backend, dtype=dtype
    )
    elapsed = time.perf_counter() - t0
    
    profiler.disable()
    
    print(f"Profiled run completed in {elapsed:.2f}s")
    print("=" * 60)
    
    return profiler, result, elapsed


def analyze_profile(profiler, n_funcs: int = 20):
    """Analyze and print profiling results."""
    print("\nTOP FUNCTIONS BY CUMULATIVE TIME:")
    print("-" * 60)
    
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.strip_dirs()
    ps.sort_stats('cumulative')
    ps.print_stats(n_funcs)
    
    output = s.getvalue()
    print(output)
    
    print("\nTOP FUNCTIONS BY TOTAL TIME:")
    print("-" * 60)
    
    s = StringIO()
    ps = pstats.Stats(profiler, stream=s)
    ps.strip_dirs()
    ps.sort_stats('tottime')
    ps.print_stats(n_funcs)
    
    output = s.getvalue()
    print(output)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Profile 3D Hessian computation")
    parser.add_argument("--n", type=int, default=80)
    parser.add_argument("--backend", default="numpy", choices=["numpy", "cupy"])
    parser.add_argument("--dtype", default="float64", choices=["float64", "float32"])
    parser.add_argument("--rho", type=float, default=5.0)
    parser.add_argument("--sigma", type=float, default=4.0)
    parser.add_argument("--v", type=float, default=1.0)
    parser.add_argument("--profile-out", type=str, default=None,
                       help="Save profile data to file")
    parser.add_argument("--n-funcs", type=int, default=20,
                       help="Number of functions to show in report")
    
    args = parser.parse_args(argv)
    
    profiler, result, elapsed = profile_3d_run(
        n=args.n,
        rho=args.rho,
        sigma=args.sigma,
        v=args.v,
        backend=args.backend,
        dtype=args.dtype,
    )
    
    analyze_profile(profiler, n_funcs=args.n_funcs)
    
    # Performance summary
    n_points = args.n ** 3
    points_per_sec = n_points / elapsed
    
    print("\n" + "=" * 60)
    print("PERFORMANCE SUMMARY")
    print("=" * 60)
    print(f"Grid points: {n_points:,}")
    print(f"Total time: {elapsed:.2f}s")
    print(f"Throughput: {points_per_sec:,.0f} points/second")
    print(f"Time per point: {elapsed/n_points*1e6:.2f} μs")
    
    energies = result['energies']
    print("\nEnergy results:")
    print(f"  E+: {energies['e_pos']:.6e}")
    print(f"  E-: {energies['e_neg']:.6e}")
    print(f"  neg_fraction: {energies['neg_fraction']:.6f}")
    print("=" * 60)
    
    # Save profile data
    if args.profile_out:
        Path(args.profile_out).parent.mkdir(parents=True, exist_ok=True)
        profiler.dump_stats(args.profile_out)
        print(f"\nProfile data saved to: {args.profile_out}")
        print(f"View with: python -m pstats {args.profile_out}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
