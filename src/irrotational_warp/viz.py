from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .adm import compute_slice_z0, integrate_signed
from .einstein import compute_einstein_eigenvalues
from .io import write_summary_json
from .sweep import SweepPoint, SweepPoint2D


def plot_slice(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    n: int,
    use_einstein: bool = False,
    use_tail_correction: bool = False,
    out_path: Path,
    json_out_path: Path,
) -> None:
    res = compute_slice_z0(rho=rho, sigma=sigma, v=v, extent=extent, n=n, tail_correction=use_tail_correction)

    epos, eneg, enet = integrate_signed(res.rho_adm, dx=res.dx, dy=res.dy)

    # Optionally compute Einstein tensor eigenvalues
    einstein_result = None
    if use_einstein:
        einstein_result = compute_einstein_eigenvalues(res.beta_x, res.beta_y, dx=res.dx, dy=res.dy)

    # Choose diagnostic to plot
    if use_einstein and einstein_result is not None:
        field = einstein_result.eig_max
        title = f"Einstein tensor max eigenvalue (z=0) | v={v}, rho={rho}, sigma={sigma}"
        label = r"$\rho_p$ (max eig of $G^\mu_\nu$)"
    else:
        field = res.rho_adm
        title = f"ADM energy density (z=0 slice) | v={v}, rho={rho}, sigma={sigma}"
        label = r"$\rho_{\rm ADM}$ (geometric units)"

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    vlim = float(np.nanmax(np.abs(field)))
    im = ax.imshow(
        field,
        origin="lower",
        extent=[res.x[0], res.x[-1], res.y[0], res.y[-1]],
        cmap="coolwarm",
        vmin=-vlim,
        vmax=vlim,
        interpolation="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label=label)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    # Build JSON output
    json_data = {
        "params": {
            "rho": rho,
            "sigma": sigma,
            "v": v,
            "extent": extent,
            "n": n,
            "use_einstein": use_einstein,
            "use_tail_correction": use_tail_correction,
        },
        "integrals_2d": {"E_pos": epos, "E_neg": eneg, "E_net": enet},
        "rho_adm": res.rho_adm,
    }
    if res.tail_correction is not None:
        json_data["tail_correction"] = {
            "exponent": res.tail_correction.exponent,
            "amplitude": res.tail_correction.amplitude,
            "fit_r_min": res.tail_correction.fit_r_min,
            "fit_r_max": res.tail_correction.fit_r_max,
            "tail_integral_pos": res.tail_correction.tail_integral_pos,
            "tail_integral_neg": res.tail_correction.tail_integral_neg,
            "tail_uncertainty": res.tail_correction.tail_uncertainty,
            "fit_residual_rms": res.tail_correction.fit_residual_rms,
            "E_pos_corrected": epos + res.tail_correction.tail_integral_pos,
            "E_neg_corrected": eneg + res.tail_correction.tail_integral_neg,
            "E_net_corrected": (epos + res.tail_correction.tail_integral_pos)
            - (eneg + res.tail_correction.tail_integral_neg),
        }
    if einstein_result is not None:
        json_data["einstein"] = {
            "eig_max": einstein_result.eig_max,
            "ricci_scalar": einstein_result.ricci_scalar,
            "type_i_fraction": float(einstein_result.is_type_i.sum()) / einstein_result.is_type_i.size,
        }

    write_summary_json(json_out_path, json_data)


def plot_sweep(
    points: list[SweepPoint],
    output_path: str | None = None,
    show: bool = False,
) -> None:
    """Create a multi-panel plot for a sigma sweep.

    Shows e_pos, e_neg, e_net, neg_fraction vs sigma.
    """
    sigmas = np.array([p.sigma for p in points])
    e_pos = np.array([p.e_pos for p in points])
    e_neg = np.array([abs(p.e_neg) for p in points])  # Show magnitude
    e_net = np.array([p.e_net for p in points])
    neg_frac = np.array([p.neg_fraction for p in points])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].plot(sigmas, e_pos, "o-", color="blue")
    axes[0, 0].set_xlabel("σ")
    axes[0, 0].set_ylabel("E⁺")
    axes[0, 0].set_title("Positive Energy Integral")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(sigmas, e_neg, "o-", color="red")
    axes[0, 1].set_xlabel("σ")
    axes[0, 1].set_ylabel("|E⁻|")
    axes[0, 1].set_title("Negative Energy Integral (magnitude)")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(sigmas, e_net, "o-", color="green")
    axes[1, 0].axhline(0, color="black", linestyle="--", linewidth=0.8)
    axes[1, 0].set_xlabel("σ")
    axes[1, 0].set_ylabel("E_net")
    axes[1, 0].set_title("Net Energy (E⁺ - |E⁻|)")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(sigmas, neg_frac, "o-", color="purple")
    axes[1, 1].set_xlabel("σ")
    axes[1, 1].set_ylabel("Neg Fraction")
    axes[1, 1].set_ylim([0, 1])
    axes[1, 1].set_title("|E⁻| / (E⁺ + |E⁻|)")
    axes[1, 1].grid(True, alpha=0.3)

    fig.tight_layout()

    if show:
        plt.show()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_heatmap_2d(
    points: list[SweepPoint2D],
    output_path: str | None = None,
    show: bool = False,
) -> None:
    """Create heatmap visualizations for 2D (sigma, v) parameter sweep.

    Shows:
    - E^- (negative energy density integral)
    - E^+ (positive energy density integral)
    - neg_fraction = |E^-| / (E^+ + |E^-|)
    """
    # Extract unique sigma and v values
    sigma_vals = sorted({p.sigma for p in points})
    v_vals = sorted({p.v for p in points})

    n_sigma = len(sigma_vals)
    n_v = len(v_vals)

    # Build 2D grids
    e_neg_grid = np.full((n_v, n_sigma), np.nan)
    e_pos_grid = np.full((n_v, n_sigma), np.nan)
    neg_frac_grid = np.full((n_v, n_sigma), np.nan)

    for pt in points:
        i_sigma = sigma_vals.index(pt.sigma)
        i_v = v_vals.index(pt.v)
        e_neg_grid[i_v, i_sigma] = abs(pt.e_neg)  # Show magnitude
        e_pos_grid[i_v, i_sigma] = pt.e_pos
        neg_frac_grid[i_v, i_sigma] = pt.neg_fraction

    # Create 3-panel heatmap
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    extent = [sigma_vals[0], sigma_vals[-1], v_vals[0], v_vals[-1]]

    # E^- magnitude
    im0 = axes[0].imshow(
        e_neg_grid,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap="viridis",
        interpolation="nearest",
    )
    axes[0].set_xlabel("σ")
    axes[0].set_ylabel("v/c")
    axes[0].set_title("|E⁻| (area integral)")
    fig.colorbar(im0, ax=axes[0])

    # E^+
    im1 = axes[1].imshow(
        e_pos_grid,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap="plasma",
        interpolation="nearest",
    )
    axes[1].set_xlabel("σ")
    axes[1].set_ylabel("v/c")
    axes[1].set_title("E⁺ (area integral)")
    fig.colorbar(im1, ax=axes[1])

    # Neg fraction
    im2 = axes[2].imshow(
        neg_frac_grid,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap="RdYlBu_r",
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
    )
    axes[2].set_xlabel("σ")
    axes[2].set_ylabel("v/c")
    axes[2].set_title("Neg Fraction |E⁻|/(E⁺+|E⁻|)")
    fig.colorbar(im2, ax=axes[2])

    fig.tight_layout()

    if show:
        plt.show()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
