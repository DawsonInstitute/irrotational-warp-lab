# Notes / Caveats

## Diagnostics we compute

### Fast diagnostics (2D z=0 slice)
- `rho_adm` is computed from the 3+1 Hamiltonian constraint assuming **flat spatial slices** and **unit lapse**.
- This is a useful fast diagnostic but is **observer-dependent**; it is not the Rodal/McMonigal invariant "proper energy density".
- Grid resolution: can use n=101 or higher; scales well.

### Invariant diagnostics (M2 - Einstein tensor eigenvalues)
- `eig_max` (max eigenvalue of G^μ_ν) approximates **proper energy density ρ_p** (Rodal convention).
- `type_i_fraction` reports fraction of grid with real eigenvalues (Type I spacetime per Hawking-Ellis).
- **Numerically expensive**: use n≤51 for 2D slices; finite-difference errors accumulate through Christoffel → Ricci → Einstein chain.
- Enable with `--einstein` flag in `plot-slice` command.

## Units

All quantities are currently in geometric units with $G=c=1$.

- `v` is dimensionless and represents $v/c$.
- Global integrals computed from `rho_adm` are in “geometric” energy units.

## Numerics

- Finite differences are used; grid resolution and domain size materially affect peak values.
- Use convergence checks before interpreting magnitude.

### Tail correction (M4 - finite-box error control)
- `tail_correction` flag enables power-law fit to radial average ⟨ρ⟩(r) in far-field.
- Fits ρ ~ A/r^n and analytically integrates tail from R to infinity (convergent for n > 2).
- Output includes `E_pos_corrected`, `E_neg_corrected` with tail contribution added.
- Useful for validating that finite-grid energies are converged; typically tail ~ 0.1-1% of grid integral for well-resolved fields.

## Next priorities

- **M5 (multi-parameter optimization)**: 2D heatmaps (sigma vs v), Bayesian search for minimal E⁻.
- **M1 extension**: Full 3D volume integration (currently 2D z=0 slice approximation).
- **Validation (M6)**: Reproduce Rodal/McMonigal parameter regimes and compare diagnostics.
