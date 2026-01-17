from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .adm import compute_slice_z0, integrate_signed


@dataclass(frozen=True)
class SweepPoint:
    sigma: float
    e_pos: float
    e_neg: float
    e_net: float
    neg_fraction: float


@dataclass(frozen=True)
class SweepPoint2D:
    """Result for a single (sigma, v) point in 2D parameter sweep."""
    sigma: float
    v: float
    e_pos: float
    e_neg: float
    e_net: float
    neg_fraction: float


def sweep_sigma_z0(
    *,
    rho: float,
    v: float,
    sigma_values: np.ndarray,
    extent: float,
    n: int,
) -> list[SweepPoint]:
    """Sweep sigma for the z=0 2D slice diagnostic.

    This uses the current 2D slice approximation in `compute_slice_z0` and returns
    area-integrated quantities.
    """
    points: list[SweepPoint] = []
    for sigma in sigma_values:
        res = compute_slice_z0(rho=rho, sigma=float(sigma), v=v, extent=extent, n=n)
        e_pos, e_neg, e_net = integrate_signed(res.rho_adm, dx=res.dx, dy=res.dy)
        denom = e_pos + e_neg
        neg_fraction = (e_neg / denom) if denom > 0.0 else 0.0
        points.append(
            SweepPoint(
                sigma=float(sigma),
                e_pos=float(e_pos),
                e_neg=float(e_neg),
                e_net=float(e_net),
                neg_fraction=float(neg_fraction),
            )
        )
    return points


def sweep_2d_z0(
    *,
    rho: float,
    sigma_values: np.ndarray,
    v_values: np.ndarray,
    extent: float,
    n: int,
) -> list[SweepPoint2D]:
    """2D parameter sweep over (sigma, v) for z=0 slice diagnostic.

    Computes energy integrals at each (sigma, v) grid point for optimization
    and heatmap visualization.

    Returns:
        List of SweepPoint2D results (length = len(sigma_values) × len(v_values))
    """
    points: list[SweepPoint2D] = []
    total = len(sigma_values) * len(v_values)
    count = 0

    for sigma in sigma_values:
        for v in v_values:
            res = compute_slice_z0(rho=rho, sigma=float(sigma), v=float(v), extent=extent, n=n)
            e_pos, e_neg, e_net = integrate_signed(res.rho_adm, dx=res.dx, dy=res.dy)
            denom = e_pos + e_neg
            neg_fraction = (e_neg / denom) if denom > 0.0 else 0.0

            points.append(
                SweepPoint2D(
                    sigma=float(sigma),
                    v=float(v),
                    e_pos=float(e_pos),
                    e_neg=float(e_neg),
                    e_net=float(e_net),
                    neg_fraction=float(neg_fraction),
                )
            )

            count += 1
            if count % 10 == 0:
                print(f"  Progress: {count}/{total} points computed")

    return points
