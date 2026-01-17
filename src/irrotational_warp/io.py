from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np


def _to_jsonable(obj: Any) -> Any:
    """Recursively convert numpy types to JSON-serializable primitives."""
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return {
            "shape": list(obj.shape),
            "min": float(np.min(obj)),
            "max": float(np.max(obj)),
            "mean": float(np.mean(obj)),
        }
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(item) for item in obj]
    return obj


def write_summary_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cooked = _to_jsonable(payload)
    path.write_text(json.dumps(cooked, indent=2, sort_keys=True) + "\n")
