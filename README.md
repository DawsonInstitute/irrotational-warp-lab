# Irrotational Warp Lab

Research code for exploring irrotational (curl-free) shift-vector warp metrics and energy-condition diagnostics.

## Quickstart

```bash
cd /home/echo_/Code/asciimath/irrotational-warp-lab

python -m venv .venv
. .venv/bin/activate
python -m pip install -e ".[dev]"

python -m irrotational_warp plot-slice --out results/slice.png --json-out results/summary.json
pytest -q
```

See [docs/TASKS.md](docs/TASKS.md) for the research plan and milestones.
