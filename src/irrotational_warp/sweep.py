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
