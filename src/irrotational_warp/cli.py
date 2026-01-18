import argparse
from pathlib import Path

import numpy as np

from .viz import plot_slice, plot_heatmap_2d, plot_sweep
from .io import write_summary_json, get_git_info
from .sweep import sweep_sigma_z0, sweep_2d_z0
from .optimize import optimize_hybrid, optimize_bayesian


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="irrotational_warp")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_plot = sub.add_parser("plot-slice", help="Render a 2D z=0 slice of rho_adm")
    p_plot.add_argument("--rho", type=float, default=10.0, help="Bubble radius (geometric units)")
    p_plot.add_argument("--sigma", type=float, default=5.0, help="Wall sharpness (1/thickness)")
    p_plot.add_argument("--v", type=float, default=1.5, help="Dimensionless v/c")
    p_plot.add_argument("--extent", type=float, default=20.0, help="Half-width of plot domain")
    p_plot.add_argument("--n", type=int, default=301, help="Grid resolution per axis")
    p_plot.add_argument("--einstein", action="store_true", help="Compute Einstein tensor eigenvalues (slower)")
    p_plot.add_argument("--tail-correction", action="store_true", help="Compute tail extrapolation for far-field decay")
    p_plot.add_argument("--out", type=Path, default=Path("results/slice.png"))
    p_plot.add_argument("--json-out", type=Path, default=Path("results/summary.json"))

    p_sweep = sub.add_parser("sweep", help="Sweep sigma values for the z=0 slice diagnostic")
    p_sweep.add_argument("--rho", type=float, default=10.0, help="Bubble radius (geometric units)")
    p_sweep.add_argument("--sigma-min", type=float, default=1.0)
    p_sweep.add_argument("--sigma-max", type=float, default=10.0)
    p_sweep.add_argument("--sigma-steps", type=int, default=10)
    p_sweep.add_argument("--v", type=float, default=1.5, help="Dimensionless v/c")
    p_sweep.add_argument("--extent", type=float, default=20.0, help="Half-width of plot domain")
    p_sweep.add_argument("--n", type=int, default=201, help="Grid resolution per axis")
    p_sweep.add_argument("--out", type=Path, default=Path("results/sweep.json"))
    p_sweep.add_argument(
        "--out-plot",
        type=Path,
        default=None,
        help="Optional plot output (E+, |E-|, Enet, neg_fraction vs sigma)",
    )

    p_sweep2d = sub.add_parser("sweep-2d", help="2D parameter sweep over (sigma, v) with heatmap visualization")
    p_sweep2d.add_argument("--rho", type=float, default=10.0, help="Bubble radius (geometric units)")
    p_sweep2d.add_argument("--sigma-min", type=float, default=1.0)
    p_sweep2d.add_argument("--sigma-max", type=float, default=10.0)
    p_sweep2d.add_argument("--sigma-steps", type=int, default=5)
    p_sweep2d.add_argument("--v-min", type=float, default=0.5)
    p_sweep2d.add_argument("--v-max", type=float, default=2.0)
    p_sweep2d.add_argument("--v-steps", type=int, default=5)
    p_sweep2d.add_argument("--extent", type=float, default=20.0, help="Half-width of plot domain")
    p_sweep2d.add_argument("--n", type=int, default=101, help="Grid resolution per axis")
    p_sweep2d.add_argument("--out-json", type=Path, default=Path("results/sweep_2d.json"))
    p_sweep2d.add_argument("--out-plot", type=Path, default=Path("results/sweep_2d_heatmap.png"))

    p_opt = sub.add_parser("optimize", help="Find parameters minimizing negative energy magnitude |E⁻|")
    p_opt.add_argument("--rho", type=float, default=10.0, help="Bubble radius (geometric units)")
    p_opt.add_argument("--sigma-min", type=float, default=1.0)
    p_opt.add_argument("--sigma-max", type=float, default=10.0)
    p_opt.add_argument("--sigma-steps", type=int, default=5, help="Grid search resolution (sigma)")
    p_opt.add_argument("--v-min", type=float, default=0.5)
    p_opt.add_argument("--v-max", type=float, default=2.0)
    p_opt.add_argument("--v-steps", type=int, default=5, help="Grid search resolution (v)")
    p_opt.add_argument("--extent", type=float, default=20.0, help="Half-width of plot domain")
    p_opt.add_argument("--n", type=int, default=71, help="Grid resolution per axis")
    p_opt.add_argument(
        "--method",
        type=str,
        default="hybrid",
        choices=["grid", "hybrid", "bayes"],
        help="Optimization method: grid (exhaustive), hybrid (grid+Nelder-Mead), bayes (Bayesian/GP)",
    )
    p_opt.add_argument("--refine", action="store_true", help="[hybrid only] Apply Nelder-Mead refinement after grid search")
    p_opt.add_argument("--n-calls", type=int, default=50, help="[bayes only] Total number of evaluations")
    p_opt.add_argument("--n-initial", type=int, default=10, help="[bayes only] Random evaluations before GP fitting")
    p_opt.add_argument("--random-state", type=int, default=None, help="[bayes only] Random seed for reproducibility")
    p_opt.add_argument("--out", type=Path, default=Path("results/optimization.json"))

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "plot-slice":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        plot_slice(
            rho=args.rho,
            sigma=args.sigma,
            v=args.v,
            extent=args.extent,
            n=args.n,
            use_einstein=args.einstein,
            use_tail_correction=args.tail_correction,
            out_path=args.out,
            json_out_path=args.json_out,
        )
        return 0

    if args.cmd == "sweep":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        if args.out_plot is not None:
            args.out_plot.parent.mkdir(parents=True, exist_ok=True)
        sigma_values = np.linspace(args.sigma_min, args.sigma_max, args.sigma_steps)
        points = sweep_sigma_z0(
            rho=args.rho,
            v=args.v,
            sigma_values=sigma_values,
            extent=args.extent,
            n=args.n,
        )
        write_summary_json(
            args.out,
            {
                "git": get_git_info(),
                "params": {
                    "rho": args.rho,
                    "v": args.v,
                    "extent": args.extent,
                    "n": args.n,
                    "sigma_min": args.sigma_min,
                    "sigma_max": args.sigma_max,
                    "sigma_steps": args.sigma_steps,
                },
                "points": [p.__dict__ for p in points],
            },
        )

        if args.out_plot is not None:
            plot_sweep(points, output_path=str(args.out_plot), show=False)
            print(f"  Saved sweep plot to {args.out_plot}")
            print(f"  Saved sweep JSON to {args.out}")
        return 0

    if args.cmd == "sweep-2d":
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_plot.parent.mkdir(parents=True, exist_ok=True)

        sigma_values = np.linspace(args.sigma_min, args.sigma_max, args.sigma_steps)
        v_values = np.linspace(args.v_min, args.v_max, args.v_steps)

        print(f"Running 2D sweep over {args.sigma_steps} × {args.v_steps} = {args.sigma_steps * args.v_steps} points...")
        points = sweep_2d_z0(
            rho=args.rho,
            sigma_values=sigma_values,
            v_values=v_values,
            extent=args.extent,
            n=args.n,
        )

        # Save JSON output
        write_summary_json(
            args.out_json,
            {
                "git": get_git_info(),
                "params": {
                    "rho": args.rho,
                    "extent": args.extent,
                    "n": args.n,
                    "sigma_min": args.sigma_min,
                    "sigma_max": args.sigma_max,
                    "sigma_steps": args.sigma_steps,
                    "v_min": args.v_min,
                    "v_max": args.v_max,
                    "v_steps": args.v_steps,
                },
                "points": [p.__dict__ for p in points],
            },
        )

        # Generate heatmap visualization
        plot_heatmap_2d(points, output_path=str(args.out_plot))
        print(f"  Saved heatmap to {args.out_plot}")
        print(f"  Saved JSON to {args.out_json}")
        return 0

    if args.cmd == "optimize":
        args.out.parent.mkdir(parents=True, exist_ok=True)

        print(f"Running optimization over sigma ∈ [{args.sigma_min}, {args.sigma_max}], v ∈ [{args.v_min}, {args.v_max}]")
        
        if args.method == "bayes":
            print(f"  Bayesian optimization (GP): {args.n_calls} total calls ({args.n_initial} initial random)")
            if args.random_state is not None:
                print(f"  Random seed: {args.random_state} (reproducible)")
            
            result = optimize_bayesian(
                rho=args.rho,
                sigma_range=(args.sigma_min, args.sigma_max),
                v_range=(args.v_min, args.v_max),
                extent=args.extent,
                n=args.n,
                n_calls=args.n_calls,
                n_initial_points=args.n_initial,
                random_state=args.random_state,
            )
        elif args.method == "grid":
            print(f"  Grid search only: {args.sigma_steps} × {args.v_steps} = {args.sigma_steps * args.v_steps} evaluations")
            result = optimize_hybrid(
                rho=args.rho,
                sigma_range=(args.sigma_min, args.sigma_max),
                v_range=(args.v_min, args.v_max),
                sigma_steps=args.sigma_steps,
                v_steps=args.v_steps,
                extent=args.extent,
                n=args.n,
                refine=False,
            )
        else:  # hybrid
            if args.refine:
                print(f"  Grid search: {args.sigma_steps} × {args.v_steps} = {args.sigma_steps * args.v_steps} evaluations")
                print("  Then: Nelder-Mead local refinement")
            else:
                print(f"  Grid search only: {args.sigma_steps} × {args.v_steps} = {args.sigma_steps * args.v_steps} evaluations")
            
            result = optimize_hybrid(
                rho=args.rho,
                sigma_range=(args.sigma_min, args.sigma_max),
                v_range=(args.v_min, args.v_max),
                sigma_steps=args.sigma_steps,
                v_steps=args.v_steps,
                extent=args.extent,
                n=args.n,
                refine=args.refine,
            )

        # Get git info for provenance
        git_info = get_git_info()

        # Save result
        output_data = {
            "git": git_info,
            "params": {
                "rho": args.rho,
                "extent": args.extent,
                "n": args.n,
                "sigma_range": [args.sigma_min, args.sigma_max],
                "v_range": [args.v_min, args.v_max],
                "method": args.method,
            },
            "optimization": {
                "best_params": result.best_params,
                "best_value": result.best_value,
                "initial_params": result.initial_params,
                "initial_value": result.initial_value,
                "n_evaluations": result.n_evaluations,
                "success": result.success,
                "method": result.method,
                "message": result.message,
            },
        }
        
        # Add method-specific params
        if args.method == "bayes":
            output_data["params"]["n_calls"] = args.n_calls
            output_data["params"]["n_initial"] = args.n_initial
            if args.random_state is not None:
                output_data["params"]["random_state"] = args.random_state
        else:
            output_data["params"]["sigma_steps"] = args.sigma_steps
            output_data["params"]["v_steps"] = args.v_steps
            if args.method == "hybrid":
                output_data["params"]["refine"] = args.refine
        
        write_summary_json(args.out, output_data)

        print("\n✓ Optimization complete:")
        print(f"  Best |E⁻|: {result.best_value:.6e}")
        print(f"  Best σ: {result.best_params['sigma']:.4f}")
        print(f"  Best v: {result.best_params['v']:.4f}")
        print(f"  Evaluations: {result.n_evaluations}")
        print(f"  Method: {result.method}")
        print(f"  Saved to {args.out}")
        return 0

    raise ValueError(f"Unknown command: {args.cmd}")
