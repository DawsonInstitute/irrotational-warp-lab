from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .adm import compute_slice_z0, integrate_signed
from .einstein import compute_einstein_eigenvalues
from .io import write_summary_json


def plot_slice(
    *,
    rho: float,
    sigma: float,
    v: float,
    extent: float,
    n: int,
    use_einstein: bool = False,
    out_path: Path,
    json_out_path: Path,
) -> None:
    res = compute_slice_z0(rho=rho, sigma=sigma, v=v, extent=extent, n=n)

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
        "params": {"rho": rho, "sigma": sigma, "v": v, "extent": extent, "n": n, "use_einstein": use_einstein},
        "integrals_2d": {"E_pos": epos, "E_neg": eneg, "E_net": enet},
        "rho_adm": res.rho_adm,
    }
    if einstein_result is not None:
        json_data["einstein"] = {
            "eig_max": einstein_result.eig_max,
            "ricci_scalar": einstein_result.ricci_scalar,
            "type_i_fraction": float(einstein_result.is_type_i.sum()) / einstein_result.is_type_i.size,
        }

    write_summary_json(json_out_path, json_data)
