#!/usr/bin/env python3
"""Generate publication-quality figures for the paper.

This script regenerates all figures referenced in docs/paper/main.tex
from results data with consistent styling.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Set publication style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.figsize': (3.5, 2.5),  # Single column width
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})


def load_json(path: Path) -> dict:
    """Load JSON data file."""
    with open(path) as f:
        return json.load(f)


def make_convergence_3d(out_path: Path):
    """Generate 3D convergence study figure."""
    print("  Loading convergence data...")
    
    # Check if convergence study exists
    data_path = Path("results/experiments/convergence/study_3d.json")
    if not data_path.exists():
        print(f"  ⚠ Data not found: {data_path}")
        print("  Run: python scripts/convergence_study_3d.py")
        return
    
    data = load_json(data_path)
    
    # Extract data
    resolutions = [p['n'] for p in data['points']]
    tail_imbalances = [p['tail_imbalance_pct'] for p in data['points']]
    ratios_r2 = [p['ratio_at_r2'] for p in data['points']]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.5))
    
    # Panel 1: Tail imbalance
    ax1.plot(resolutions, tail_imbalances, 'o-', label='Numerical', linewidth=2)
    ax1.axhline(0.04, color='red', linestyle='--', label="Rodal's ~0.04%")
    ax1.set_xlabel('Grid resolution $n$')
    ax1.set_ylabel('Tail imbalance (%)')
    ax1.set_title('(a) Convergence to Literature Value')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Ratio at R2
    ax2.plot(resolutions, ratios_r2, 's-', color='C1', linewidth=2, label='Numerical')
    ax2.axhline(1.07, color='red', linestyle='--', label="Rodal's ~1.07")
    ax2.set_xlabel('Grid resolution $n$')
    ax2.set_ylabel('Energy ratio at $R_2$')
    ax2.set_title('(b) Ratio $E^+/|E^-|$ at $r=2\\rho$')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"  ✓ Saved: {out_path}")


def make_superluminal(out_path: Path):
    """Generate superluminal velocity sweep figure."""
    print("  Loading superluminal sweep data...")
    
    data_path = Path("results/experiments/superluminal/sweep_v1_to_3.json")
    if not data_path.exists():
        print(f"  ⚠ Data not found: {data_path}")
        print("  Run: python scripts/sweep_superluminal.py ...")
        return
    
    data = load_json(data_path)
    
    velocities = [p['v'] for p in data['points']]
    e_pos = [p['e_pos'] for p in data['points']]
    e_neg = [p['e_neg'] for p in data['points']]
    tail_pct = [p['tail_imbalance_pct'] for p in data['points']]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.5))
    
    # Panel 1: Energy scaling
    ax1.plot(velocities, e_pos, 'o-', label='$E^+$', linewidth=2)
    ax1.plot(velocities, e_neg, 's-', label='$|E^-|$', linewidth=2)
    
    # Show v^2 scaling reference
    v_ref = np.array(velocities)
    e_ref_scale = e_pos[0] * (v_ref / velocities[0])**2
    ax1.plot(velocities, e_ref_scale, 'k--', alpha=0.5, label='$\\propto v^2$')
    
    ax1.set_xlabel('Velocity parameter $v$')
    ax1.set_ylabel('Energy (dimensionless)')
    ax1.set_title('(a) Energy Scaling')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # Panel 2: Tail imbalance (should be constant)
    ax2.plot(velocities, tail_pct, 'o-', color='C2', linewidth=2)
    ax2.axhline(np.mean(tail_pct), color='gray', linestyle='--', 
                label=f'Mean: {np.mean(tail_pct):.3f}%')
    ax2.set_xlabel('Velocity parameter $v$')
    ax2.set_ylabel('Tail imbalance (%)')
    ax2.set_title('(b) Velocity-Independent Imbalance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"  ✓ Saved: {out_path}")


def make_optimization(out_path: Path):
    """Generate optimization comparison figure."""
    print("  Loading optimization comparison data...")
    
    data_path = Path("results/optimization/bayesian_comparison.json")
    if not data_path.exists():
        print(f"  ⚠ Data not found: {data_path}")
        print("  Run: python scripts/test_bayesian_optimization.py")
        return
    
    data = load_json(data_path)
    
    methods = ['Hybrid\n(Grid+NM)', 'Bayesian\n(GP)']
    evals = [data['hybrid']['n_evaluations'], data['bayesian']['n_evaluations']]
    values = [data['hybrid']['best_value'], data['bayesian']['best_value']]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.5))
    
    # Panel 1: Number of evaluations
    bars1 = ax1.bar(methods, evals, color=['C0', 'C1'], alpha=0.7)
    ax1.set_ylabel('Number of evaluations')
    ax1.set_title('(a) Computational Efficiency')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, val in zip(bars1, evals):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(val)}', ha='center', va='bottom')
    
    # Panel 2: Best objective value
    bars2 = ax2.bar(methods, values, color=['C0', 'C1'], alpha=0.7)
    ax2.set_ylabel('Best $|E^-|$ (dimensionless)')
    ax2.set_title('(b) Optimization Quality')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, val in zip(bars2, values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.3f}', ha='center', va='bottom')
    
    # Add efficiency annotation
    efficiency_gain = evals[0] / evals[1]
    fig.text(0.5, 0.02, 
             f'Bayesian: {efficiency_gain:.1f}× fewer evaluations for equivalent quality',
             ha='center', fontsize=9, style='italic')
    
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(out_path)
    print(f"  ✓ Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate paper figures")
    parser.add_argument('--figure', type=str, required=True,
                       choices=['convergence_3d', 'superluminal', 'optimization', 'all'],
                       help='Which figure to generate')
    parser.add_argument('--out', type=Path, default=None,
                       help='Output path (default: docs/paper/figures/<figure>.pdf)')
    
    args = parser.parse_args()
    
    if args.out is None:
        args.out = Path(f"docs/paper/figures/{args.figure}.pdf")
    
    args.out.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating figure: {args.figure}")
    
    if args.figure == 'convergence_3d' or args.figure == 'all':
        make_convergence_3d(args.out if args.figure != 'all' else 
                           Path("docs/paper/figures/convergence_3d.pdf"))
    
    if args.figure == 'superluminal' or args.figure == 'all':
        make_superluminal(args.out if args.figure != 'all' else 
                         Path("docs/paper/figures/superluminal_sweep.pdf"))
    
    if args.figure == 'optimization' or args.figure == 'all':
        make_optimization(args.out if args.figure != 'all' else 
                         Path("docs/paper/figures/optimization_comparison.pdf"))
    
    print("✓ Figure generation complete")


if __name__ == "__main__":
    main()
