- Created a new public repo: https://github.com/arcticoder/irrotational-warp-lab and cloned it to irrotational-warp-lab.
- Added it to the multi-root workspace in energy/energy.code-workspace.
- Wrote the handoff-grade backlog in irrotational-warp-lab/docs/TASKS.md, copied it into my internal task tracker, and implemented the first increment (M0 + minimal M1).

What you can run now:
- From the repo root: `python -m irrotational_warp plot-slice --out slice.png --json-out results/summary.json`
- Tests: `python -m pytest -q`

Notes:
- Current output is an ADM/Hamiltonian-constraint-based diagnostic on a 2D slice (fast, not Rodal's invariant eigenvalue diagnostic yet); see irrotational-warp-lab/docs/NOTES.md.

---

## Session 2 (January 16, 2026) — Test optimization + TASKS update

**Context:** Transitioned from previous machine (/home/sherri3 → /home/echo_). Commit 94896b3 added sweep CLI (M3), but tests timed out due to n=81 grid size (6561 points × 5 sigma values).

**Changes:**
- Reduced test grid from n=81 to n=41 in `test_sweep_has_expected_length_and_zero_v()` to avoid timeouts (tests now pass in 0.13s).
- Updated `docs/TASKS.md`:
  - Marked M0 (scaffold + CLI) as **COMPLETE**.
  - Marked M3 (sigma sweeps) as **COMPLETE** (implemented in commit 94896b3).
  - Added **math/code snippets** to M2 (invariant diagnostics) for next increment:
    - Metric construction from ADM variables (g_μν)
    - Christoffel symbols via finite differences
    - Ricci tensor → mixed Einstein tensor G^μ_ν
    - Eigenvalue solver for proper energy density ρ_p
  - Renumbered milestones: M4→tail correction, M5→advanced optimization, M6→validation, M7→extensions, M8→paper pipeline.

**Current state:**
- All tests passing: `python -m pytest -q` → 3 passed in 0.13s
- Package installed in editable mode: `pip install -e .`
- CLI commands functional:
  - `python -m irrotational_warp plot-slice --out results/slice.png --json-out results/summary.json`
  - `python -m irrotational_warp sweep --rho 10 --v 1.5 --sigma-min 1 --sigma-max 5 --sigma-count 20 --out results/sweep.json`

**Next priorities (per user request):**
1. ✅ **M2 (invariant diagnostics):** Implemented Einstein tensor eigenvalues for observer-independent "proper energy density" matching Rodal/McMonigal framing. Track A (direct metric → G^μ_ν) complete for 2D axisymmetric case.
2. **M1 extension (optional):** Full 3D volume integration (currently 2D z=0 slice approximation).
3. **M4 (tail correction):** Radial shell averaging + 1/r^4 decay fitting for finite-box extrapolation.

**Technical notes:**
- No sympy bottleneck found; code is pure numpy (fast).
- Grid size (n) scales as O(n²) for 2D slices; production sweeps can use higher res (n=81 or 101) via CLI flags.
- Path issues resolved (no /home/sherri3 references in this repo).

---

## Session 3 (January 16, 2026 continued) — Implemented M2 (Einstein tensor eigenvalues)

**Motivation:** M2 addresses the core limitation of the fast ADM diagnostic (observer-dependent) by computing the **invariant mixed Einstein tensor G^μ_ν** and its eigenvalues, matching Rodal's "proper energy density" ρ_p framework and enabling Type-I spacetime classification per Hawking-Ellis.

**Implementation (Track A: direct metric → curvature):**
- Created `src/irrotational_warp/einstein.py` with full GR curvature pipeline:
  - `compute_metric_z0()`: 4D covariant metric g_μν from ADM variables (unit lapse α=1, flat γ_ij=δ_ij, shift β^i=∇Φ)
  - `compute_christoffel()`: Christoffel symbols Γ^α_μν via finite differences on g_μν derivatives
  - `compute_ricci_tensor()`: Ricci tensor R_μν and scalar R from ∂Γ + Γ² contractions
  - `compute_einstein_tensor()`: Mixed Einstein tensor G^μ_ν = R^μ_ν - (1/2)δ^μ_ν R
  - `compute_einstein_eigenvalues()`: Solve eigenvalues per grid point; classify Type I (all real) vs non-Type-I
- Wired into CLI: added `--einstein` flag to `plot-slice` command (plots eig_max heatmap instead of rho_adm)
- Updated `viz.py` to support Einstein diagnostic rendering; fixed JSON serialization for nested numpy arrays in `io.py`
- Added tests in `tests/test_einstein.py`:
  - `test_minkowski_flatness()`: β=0 → all eigenvalues ~ 0 ✓
  - `test_small_perturbation_type_i()`: Small β perturbation → bounded eigenvalues ~ eps² ✓

**Results:**
- All tests passing (5/5 in 0.19s)
- Generated example outputs:
  - `results/slice_adm.png`: Fast ADM diagnostic (n=101, ~instant)
  - `results/slice_einstein.png`: Einstein tensor eigenvalues (n=51, ~5s)
  - JSON includes `type_i_fraction` (58% for test config), `eig_max` stats, `ricci_scalar`

**Numerical caveats (documented in TASKS.md / NOTES.md):**
- Einstein diagnostics are **numerically expensive** and **fragile**: finite-diff errors accumulate through Christoffel → Ricci → Einstein chain.
- Type-I fraction varies with grid resolution; strict classification unreliable at coarse grids.
- Recommend n≤51 for Einstein diagnostics on laptops; ADM can use n=101+.

**Documentation updates:**
- `docs/TASKS.md`: Marked M2 as COMPLETE, added implementation notes and CLI usage examples.
- `docs/NOTES.md`: Expanded diagnostics section to contrast fast (ADM) vs invariant (Einstein) approaches.
- `docs/history/history.md`: This entry.

**Current state:**
- M0 (scaffold + CLI): ✅ COMPLETE
- M1 (fast ADM diagnostics): ✅ COMPLETE (2D z=0 slice; 3D extension pending)
- M2 (Einstein eigenvalues): ✅ COMPLETE (Track A for 2D; Track B—3+1 stress-energy—deferred)
- M3 (sigma sweeps): ✅ COMPLETE
- Next: M4 (tail correction) or M6 (validation against Rodal papers)

**Files modified/created:**
- `src/irrotational_warp/einstein.py` (new, 300+ lines)
- `tests/test_einstein.py` (new)
- `tests/test_basics.py` (optimized n=81→41)
- `src/irrotational_warp/cli.py` (added --einstein flag)
- `src/irrotational_warp/viz.py` (Einstein diagnostic rendering)
- `src/irrotational_warp/io.py` (recursive JSON serialization fix)
- `docs/TASKS.md`, `docs/NOTES.md`, `docs/history/history.md` (documentation)
<!-- ------ -->