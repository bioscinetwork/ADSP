#!/usr/bin/env python3
"""
Docking Automation Suite — Vina Output Parser
===============================================
Parses Vina multi-model PDBQT output and log files.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from models import VinaResult


def count_models(pdbqt_path: Path) -> int:
    """Count the number of MODEL records in a multi-model PDBQT file."""
    if not pdbqt_path.exists():
        return 0
    content = pdbqt_path.read_text(encoding="utf-8", errors="replace")
    count = content.count("MODEL")
    return count if count > 0 else (1 if "ATOM" in content else 0)


def extract_affinities(pdbqt_path: Path) -> List[VinaResult]:
    """Extract binding affinities from a Vina output PDBQT.

    Parses REMARK VINA RESULT lines:
        REMARK VINA RESULT:    -7.3      0.000      0.000

    Returns:
        List of VinaResult objects with pose number, affinity, RMSD bounds.
    """
    results: List[VinaResult] = []

    if not pdbqt_path.exists():
        return results

    content = pdbqt_path.read_text(encoding="utf-8", errors="replace")
    pattern = re.compile(
        r"REMARK\s+VINA\s+RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    )

    pose = 0
    for line in content.splitlines():
        match = pattern.match(line.strip())
        if match:
            pose += 1
            results.append(VinaResult(
                pose=pose,
                binding_affinity=float(match.group(1)),
                rmsd_lower_bound=float(match.group(2)),
                rmsd_upper_bound=float(match.group(3)),
            ))

    return results


def parse_vina_log(log_path: Path) -> Dict[str, str]:
    """Parse a Vina stdout/log file for parameters used.

    Returns dict with keys like 'exhaustiveness', 'num_modes', etc.
    """
    params: Dict[str, str] = {}

    if not log_path.exists():
        return params

    content = log_path.read_text(encoding="utf-8", errors="replace")

    # Common parameter patterns from Vina output
    patterns = {
        "exhaustiveness": r"exhaustiveness\s*[=:]\s*(\d+)",
        "num_modes": r"num_modes\s*[=:]\s*(\d+)",
        "energy_range": r"energy_range\s*[=:]\s*([\d.]+)",
        "cpu": r"cpu\s*[=:]\s*(\d+)",
        "seed": r"seed\s*[=:]\s*(\d+)",
        "center_x": r"center_x\s*[=:]\s*([-\d.]+)",
        "center_y": r"center_y\s*[=:]\s*([-\d.]+)",
        "center_z": r"center_z\s*[=:]\s*([-\d.]+)",
        "size_x": r"size_x\s*[=:]\s*([-\d.]+)",
        "size_y": r"size_y\s*[=:]\s*([-\d.]+)",
        "size_z": r"size_z\s*[=:]\s*([-\d.]+)",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            params[key] = match.group(1)

    return params


def format_vina_results_table(results: List[VinaResult]) -> str:
    """Format Vina results as a human-readable table.

    Returns:
        Formatted string for display/logging.
    """
    if not results:
        return "  No results found."

    lines = [
        "  Pose  Affinity(kcal/mol)  RMSD_lb  RMSD_ub",
        "  " + "-" * 48,
    ]
    for r in results:
        lines.append(
            f"  {r.pose:4d}  {r.binding_affinity:17.1f}  "
            f"{r.rmsd_lower_bound:7.3f}  {r.rmsd_upper_bound:7.3f}"
        )
    return "\n".join(lines)
