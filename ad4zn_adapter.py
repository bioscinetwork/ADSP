"""Isolated AutoDock4Zn reference-workflow adapter.

This module deliberately does not import or copy the legacy Python 2 scripts.
It performs read-only availability and input checks and reports when the
reference workflow cannot be executed in the current ADSP process.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import shutil
import zipfile
from typing import List, Optional

from ad4_compatibility import AD4ZN_PROFILE, extract_atom_types


@dataclass(frozen=True)
class AD4ZnPreflight:
    status: str
    errors: List[str]
    warnings: List[str]
    archive_path: Optional[str]
    parameter_member: str
    receptor_atom_types: List[str]
    has_zinc: bool
    dependencies: dict

    def to_dict(self):
        return asdict(self)


def preflight_ad4zn(receptor_pdbqt: Path | str,
                    archive_path: Path | str | None = None,
                    mgltools_root: Path | str | None = None) -> AD4ZnPreflight:
    """Check specialized workflow prerequisites without executing legacy code."""
    root = Path(__file__).resolve().parent
    archive = Path(archive_path) if archive_path else root / "AutoDock4Zn-Pipeline-main.zip"
    member = "AutoDock4Zn-Pipeline-main/AD4Zn.dat"
    errors: List[str] = []
    warnings: List[str] = []
    if not archive.is_file():
        errors.append(f"AutoDock4Zn reference archive not found: {archive}")
    elif not zipfile.is_zipfile(archive):
        errors.append(f"AutoDock4Zn reference archive is not a valid ZIP: {archive}")
    else:
        with zipfile.ZipFile(archive) as package:
            names = set(package.namelist())
            if member not in names:
                errors.append(f"Specialized parameter member not found: {member}")
            required = [
                "zinc_pseudo.py", "prepare_gpf4zn.py", "prepare_dpf42.py",
                "prepare_receptor4.py", "prepare_ligand4.py",
            ]
            for filename in required:
                if f"AutoDock4Zn-Pipeline-main/{filename}" not in names:
                    errors.append(f"Reference workflow member not found: {filename}")

    receptor = Path(receptor_pdbqt)
    atom_types = extract_atom_types(receptor)
    has_zinc = any(value.upper() in {"ZN", "ZINC"} for value in atom_types)
    if not has_zinc:
        errors.append("SPECIALIZED_WORKFLOW_NOT_APPLICABLE: no Zn atom type was detected.")

    mgl_root = Path(mgltools_root) if mgltools_root else None
    python2 = None
    if mgl_root:
        python2 = mgl_root / "python.exe"
    else:
        python2 = Path(shutil.which("python2") or "")
    dependencies = {
        "legacy_python": str(python2) if python2 and python2.is_file() else None,
        "autodock4": str(root / "bin" / "autodock4.exe") if (root / "bin" / "autodock4.exe").is_file() else None,
        "autogrid4": str(root / "bin" / "autogrid4.exe") if (root / "bin" / "autogrid4.exe").is_file() else None,
    }
    if not dependencies["legacy_python"]:
        errors.append("SPECIALIZED_WORKFLOW_UNAVAILABLE: Python 2/MGLTools runtime was not found.")
    if dependencies["autodock4"] is None or dependencies["autogrid4"] is None:
        errors.append("SPECIALIZED_WORKFLOW_UNAVAILABLE: AutoDock4/AutoGrid4 executable is missing.")
    if not errors:
        warnings.append("Reference workflow remains isolated; no ADSP-native AutoDock4Zn execution is implemented.")
    return AD4ZnPreflight(
        status="AVAILABLE" if not errors else "SPECIALIZED_WORKFLOW_UNAVAILABLE",
        errors=errors,
        warnings=warnings,
        archive_path=str(archive) if archive.exists() else None,
        parameter_member=member,
        receptor_atom_types=atom_types,
        has_zinc=has_zinc,
        dependencies=dependencies,
    )


def run_ad4zn(*args, **kwargs):
    """Fail explicitly until an isolated legacy execution environment is supplied."""
    raise RuntimeError(
        "SPECIALIZED_WORKFLOW_UNAVAILABLE: AutoDock4Zn reference scripts require "
        "an isolated Python 2/MGLTools runtime; standard AutoDock4 fallback is disabled."
    )

