# Contributing

Thanks for considering a contribution! This repo is set up as a research codebase with a reproducible paper pipeline.

## Development setup

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

## Common workflows

- Run tests:

```bash
make test
```

- Run tests with coverage:

```bash
make test-cov
```

- Lint and format:

```bash
make lint
make format
```

- Build the paper (requires a LaTeX installation with `pdflatex` + `bibtex`):

```bash
make all
```

## What to include in a PR

- A short description of the scientific/numerical motivation.
- Reproduction steps (command line + expected artifact paths in `results/` and/or `papers/figures/`).
- Tests for bug fixes or numerics changes when feasible.

## Style

- Keep changes focused and avoid unrelated refactors.
- Prefer small, deterministic tests (fixed seeds, small grids).
- Use `ruff` for linting/formatting.

## Reporting results

When adding new experiments or scripts:
- Ensure there is an explicit `--out` argument (JSON/PNG/PDF) and do not default to writing in the repo root.
- Store larger experiment outputs under `results/experiments/`.
