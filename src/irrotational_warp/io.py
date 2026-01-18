from __future__ import annotations

import json
import subprocess
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


def get_git_sha() -> str | None:
    """Get current git commit SHA.
    
    Returns:
        Git SHA string if available, None otherwise
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return None


def get_git_info() -> dict[str, str | None]:
    """Get git repository information for provenance tracking.
    
    Returns:
        Dictionary with git metadata (sha, branch, status)
    """
    info = {
        "sha": None,
        "branch": None,
        "dirty": None,
    }
    
    try:
        # Get SHA
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        info["sha"] = result.stdout.strip()
        
        # Get branch
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        info["branch"] = result.stdout.strip()
        
        # Check if dirty
        result = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True,
            text=True,
            check=True,
            timeout=5,
        )
        info["dirty"] = "yes" if result.stdout.strip() else "no"
        
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    return info
