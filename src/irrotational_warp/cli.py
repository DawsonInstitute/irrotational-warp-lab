import argparse
from pathlib import Path

import numpy as np

from .viz import plot_slice
from .io import write_summary_json
from .sweep import sweep_sigma_z0


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
        return 0

    raise ValueError(f"Unknown command: {args.cmd}")
