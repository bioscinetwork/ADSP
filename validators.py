#!/usr/bin/env python3
"""
Docking Automation Suite — Validation Module
==============================================
Input, runtime, and post-processing validators for docking jobs.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger("docking_automation.validators")


def inspect_pdbqt_charges(path: Path, expected_total: float | None = None) -> dict:
    """Report PDBQT charge-field completeness and numerical total.

    Numeric completeness is a file-format check; it does not validate the
    scientific correctness of a charge model.
    """
    import math
    path = Path(path)
    records = [line for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
               if line.startswith(("ATOM", "HETATM"))]
    values, missing, malformed = [], [], []
    for index, line in enumerate(records, 1):
        raw = line[70:76].strip() if len(line) >= 76 else ""
        try:
            value = float(raw)
            if not math.isfinite(value):
                raise ValueError("non-finite")
            values.append(value)
        except (ValueError, TypeError):
            if not raw:
                missing.append(index)
            else:
                malformed.append(index)
    total = sum(values) if len(values) == len(records) else None
    return {
        "path": str(path), "total_atoms": len(records), "atoms_with_charge": len(values),
        "atoms_without_charge": missing, "atoms_with_invalid_charge": malformed,
        "completeness_fraction": (len(values) / len(records)) if records else 0.0,
        "calculated_total_charge": total,
        "expected_total_charge": expected_total,
        "expected_charge_status": "NOT_DECLARED" if expected_total is None else
            ("MATCH" if total is not None and abs(total - expected_total) <= 1e-3 else "MISMATCH"),
        "scientific_correctness_established": False,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Input Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate_pdbqt(path: Path) -> Tuple[bool, List[str]]:
    """Validate a PDBQT file for basic structural integrity.

    Checks:
    - File exists and is readable
    - File is non-empty
    - Contains ATOM or HETATM records

    Returns:
        (is_valid, list of error/warning messages)
    """
    errors: List[str] = []

    if not path.exists():
        return False, [f"File does not exist: {path}"]
    if not path.is_file():
        return False, [f"Not a file: {path}"]

    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return False, [f"Cannot read file {path}: {e}"]

    if not content.strip():
        return False, [f"File is empty: {path}"]

    # Check fixed-column ATOM/HETATM records.  PDBQT is not a whitespace-only
    # format: charge occupies columns 71-76 and the AutoDock atom type occupies
    # columns 78-79 (1-based).  Reject malformed producers explicitly.
    has_atoms = False
    for line_no, line in enumerate(content.splitlines(), 1):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        has_atoms = True
        prefix = f"{path}:{line_no}"
        if len(line) < 78:
            errors.append(f"{prefix}: ATOM/HETATM record is shorter than PDBQT columns 1-78")
            continue
        try:
            int(line[6:11].strip())
        except ValueError:
            errors.append(f"{prefix}: invalid atom serial in columns 7-11")
        if not line[12:16].strip():
            errors.append(f"{prefix}: blank atom name in columns 13-16")
        if not line[17:20].strip():
            errors.append(f"{prefix}: blank residue name in columns 18-20")
        for label, start, end in (("x", 30, 38), ("y", 38, 46), ("z", 46, 54)):
            try:
                float(line[start:end])
            except ValueError:
                errors.append(f"{prefix}: invalid {label} coordinate in columns {start + 1}-{end}")
        try:
            charge = float(line[70:76].strip())
            if not __import__("math").isfinite(charge):
                raise ValueError
        except ValueError:
            errors.append(f"{prefix}: invalid partial charge in columns 71-76")
        if not line[77:79].strip():
            errors.append(f"{prefix}: blank AutoDock atom type in columns 78-79")

    if not has_atoms:
        errors.append(f"No ATOM/HETATM records found in: {path}")
        return False, errors

    return not errors, errors


def validate_pdbqt_multimodel(path: Path) -> Tuple[bool, int, List[str]]:
    """Validate a multi-model PDBQT file (Vina output).

    Returns:
        (is_valid, model_count, list of messages)
    """
    errors: List[str] = []
    valid, base_errors = validate_pdbqt(path)
    if not valid:
        return False, 0, base_errors

    content = path.read_text(encoding="utf-8", errors="replace")
    model_count = content.count("MODEL")
    endmodel_count = content.count("ENDMDL")

    if model_count == 0:
        # Single model PDBQT (no MODEL/ENDMDL markers) — still valid
        return True, 1, []

    if model_count != endmodel_count:
        errors.append(
            f"MODEL/ENDMDL mismatch: {model_count} MODEL vs "
            f"{endmodel_count} ENDMDL in {path}"
        )
        return False, model_count, errors

    return True, model_count, errors


def validate_vina_config(path: Path) -> Tuple[bool, List[str]]:
    """Validate a Vina configuration file (config.txt).

    Required parameters: center_x, center_y, center_z, size_x, size_y, size_z

    Returns:
        (is_valid, list of error messages)
    """
    errors: List[str] = []

    if not path.exists():
        return False, [f"Config file does not exist: {path}"]

    try:
        content = path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return False, [f"Cannot read config file {path}: {e}"]

    required_params = ["center_x", "center_y", "center_z",
                       "size_x", "size_y", "size_z"]
    found_params = set()

    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for param in required_params:
            if line.startswith(param):
                # Validate it has a value
                match = re.match(rf"{param}\s*=\s*([\d.\-+eE]+)", line)
                if match:
                    found_params.add(param)

    missing = set(required_params) - found_params
    if missing:
        errors.append(
            f"Missing required parameters in {path.name}: {', '.join(sorted(missing))}"
        )
        return False, errors

    return True, []


def parse_vina_config(path: Path) -> dict:
    """Parse a Vina config file and return a dictionary of parameters.

    Returns dict with keys like 'center_x', 'size_x', 'exhaustiveness', etc.
    """
    params = {}
    if not path.exists():
        return params

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"(\w+)\s*=\s*(.+)", line)
        if match:
            key = match.group(1).strip()
            value = match.group(2).strip()
            # Try to convert to numeric
            try:
                if "." in value:
                    params[key] = float(value)
                else:
                    params[key] = int(value)
            except ValueError:
                params[key] = value

    return params


# ═══════════════════════════════════════════════════════════════════════════════
# Runtime Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate_vina_output(output_path: Path) -> Tuple[bool, int, List[str]]:
    """Validate Vina docking output.

    Returns:
        (is_valid, model_count, list of messages)
    """
    if not output_path.exists():
        return False, 0, [f"Output file does not exist: {output_path}"]
    if output_path.stat().st_size == 0:
        return False, 0, [f"Output file is empty: {output_path}"]

    return validate_pdbqt_multimodel(output_path)


def validate_dlg(dlg_path: Path) -> Tuple[bool, List[str]]:
    """Validate an AutoDock4 DLG output file.

    Checks:
    - File exists and is non-empty
    - Contains docking result markers

    Returns:
        (is_valid, list of messages)
    """
    if not dlg_path.exists():
        return False, [f"DLG file does not exist: {dlg_path}"]
    if dlg_path.stat().st_size == 0:
        return False, [f"DLG file is empty: {dlg_path}"]

    try:
        content = dlg_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return False, [f"Cannot read DLG file: {e}"]

    # Check for key DLG markers
    has_docked = "DOCKED:" in content
    has_energy = "Estimated Free Energy of Binding" in content

    messages: List[str] = []
    if not has_docked:
        messages.append("No DOCKED: blocks found in DLG")
    if not has_energy:
        messages.append("No binding energy found in DLG")

    is_valid = has_docked and has_energy
    return is_valid, messages


def validate_grid_maps(receptor_dir: Path, receptor_stem: str) -> Tuple[bool, List[str]]:
    """Validate that AutoGrid4 map files exist for a receptor.

    Checks for .maps.fld file and at least one .map file.

    Returns:
        (is_valid, list of messages)
    """
    messages: List[str] = []

    fld_files = list(receptor_dir.glob(f"{receptor_stem}*.maps.fld"))
    if not fld_files:
        fld_files = list(receptor_dir.glob("*.maps.fld"))
    if not fld_files:
        messages.append(f"No .maps.fld file found in {receptor_dir}")
        return False, messages

    map_files = list(receptor_dir.glob("*.map"))
    if not map_files:
        messages.append(f"No .map files found in {receptor_dir}")
        return False, messages

    return True, messages


# ═══════════════════════════════════════════════════════════════════════════════
# Post-processing Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate_split_output(poses_dir: Path, expected_count: int) -> Tuple[bool, List[str]]:
    """Validate vina_split output.

    Checks that the number of split PDBQT files matches the expected model count.

    Returns:
        (is_valid, list of messages)
    """
    if not poses_dir.exists():
        return False, [f"Poses directory does not exist: {poses_dir}"]

    pose_files = sorted(poses_dir.glob("*.pdbqt"))
    actual_count = len(pose_files)

    if actual_count == 0:
        return False, [f"No split pose files found in {poses_dir}"]

    messages: List[str] = []
    if actual_count != expected_count:
        messages.append(
            f"Split file count mismatch: expected {expected_count}, "
            f"found {actual_count} in {poses_dir}"
        )
        return False, messages

    return True, messages


# ═══════════════════════════════════════════════════════════════════════════════
# Experimental Re-docking RMSD Validation Engine
# ═══════════════════════════════════════════════════════════════════════════════

def compute_pose_validation_rmsd(
    reference_path: Path,
    pose_file_path: Path,
    pose_index: int = 1,
) -> Tuple[float, float, int]:
    """Calculate (positional_rmsd, shape_kabsch_rmsd, matched_atoms) between
    an experimental reference ligand (SDF, MOL, MOL2, PDB, PDBQT, CIF) and a docked pose
    (from Vina PDBQT or AD4 DLG).

    Returns:
        (positional_rmsd, shape_kabsch_rmsd, matched_atom_count)
    """
    ref_p = Path(reference_path)
    pose_p = Path(pose_file_path)
    if not ref_p.is_file():
        raise FileNotFoundError(f"Reference ligand file not found: {ref_p}")
    if not pose_p.is_file():
        raise FileNotFoundError(f"Pose file not found: {pose_p}")

    from complex_builder import _extract_pose_lines
    from dlg_extract import _autodock_to_element, RMSDEngine

    try:
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
        import numpy as np
        has_rdkit = True
    except ImportError:
        has_rdkit = False

    ref_ext = ref_p.suffix.lower()
    pose_lines = _extract_pose_lines(pose_p, pose_index)
    if not pose_lines:
        raise ValueError(f"Pose {pose_index} coordinates could not be extracted from {pose_p}")

    pose_atoms = [l for l in pose_lines if l.startswith(("ATOM", "HETATM"))]
    if not pose_atoms:
        raise ValueError(f"No ATOM/HETATM lines found for pose {pose_index} in {pose_p}")

    # Build PDB block for pose with valid element symbols
    pose_pdb_lines = []
    for l in pose_atoms:
        nm = l[12:16].strip()
        raw = l[77:79].strip() if len(l) >= 79 else ""
        elem = _autodock_to_element(raw, nm)
        pose_pdb_lines.append(f"{l[:66].ljust(76)}{elem:>2s}\n")

    # ── Tier 1: OpenBabel obrms Engine (Hungarian Symmetry & Automorphism Matching) ──
    try:
        import shutil
        import subprocess
        import tempfile
        from executables import resolve_executable

        obrms_exe = resolve_executable(None, "obrms.exe")
        if not obrms_exe.is_file():
            w = shutil.which("obrms") or shutil.which("obrms.exe")
            if w:
                obrms_exe = Path(w)

        if obrms_exe.is_file():
            # Write single pose to a temporary file
            with tempfile.NamedTemporaryFile(suffix=".pdbqt", mode="w", delete=False) as tmp:
                tmp.writelines(pose_lines)
                tmp_pose_path = Path(tmp.name)

            try:
                # 1. In-situ / Positional RMSD (no rotation/translation)
                r_pos = subprocess.run([str(obrms_exe), str(ref_p), str(tmp_pose_path)], capture_output=True, text=True, timeout=15)
                # 2. Minimized / Shape RMSD (-m Hungarian / Kabsch alignment)
                r_min = subprocess.run([str(obrms_exe), "-m", str(ref_p), str(tmp_pose_path)], capture_output=True, text=True, timeout=15)

                pos_val = None
                for line in r_pos.stdout.splitlines():
                    if "RMSD" in line:
                        parts = line.strip().split()
                        if parts and parts[-1].lower() != "inf":
                            try:
                                pos_val = float(parts[-1])
                            except ValueError:
                                pass

                min_val = None
                for line in r_min.stdout.splitlines():
                    if "RMSD" in line:
                        parts = line.strip().split()
                        if parts and parts[-1].lower() != "inf":
                            try:
                                min_val = float(parts[-1])
                            except ValueError:
                                pass

                if pos_val is not None:
                    n_heavy = sum(1 for l in pose_atoms if not l[12:16].strip().startswith("H"))
                    if min_val is None:
                        min_val = pos_val
                    return round(pos_val, 3), round(min_val, 3), n_heavy
            finally:
                try:
                    tmp_pose_path.unlink()
                except Exception:
                    pass
    except Exception as e:
        logger.debug(f"obrms validation engine note: {e}")

    # ── Tier 2: RDKit MCS Fallback (Handles non-identical scaffolds / partial matches) ──

    if has_rdkit:
        ref_mol = None
        if ref_ext in (".sdf", ".mol"):
            suppl = Chem.SDMolSupplier(str(ref_p), removeHs=False)
            if suppl and len(suppl) > 0:
                ref_mol = suppl[0]
        elif ref_ext == ".mol2":
            ref_mol = Chem.MolFromMol2File(str(ref_p), removeHs=False)
        elif ref_ext in (".pdb", ".pdbqt"):
            ref_lines = [l for l in ref_p.read_text(encoding="utf-8", errors="replace").splitlines() if l.startswith(("ATOM", "HETATM"))]
            pdb_conv = []
            for l in ref_lines:
                nm = l[12:16].strip()
                raw = l[77:79].strip() if len(l) >= 79 else ""
                elem = _autodock_to_element(raw, nm)
                pdb_conv.append(f"{l[:66].ljust(76)}{elem:>2s}\n")
            ref_mol = Chem.MolFromPDBBlock("".join(pdb_conv), removeHs=False)

        if ref_mol is not None:
            ref_hvy = Chem.RemoveHs(ref_mol)
            pose_mol = Chem.MolFromPDBBlock("".join(pose_pdb_lines), removeHs=False)
            if pose_mol is not None:
                pose_hvy = Chem.RemoveHs(pose_mol)
                mcs = rdFMCS.FindMCS(
                    [ref_hvy, pose_hvy],
                    timeout=5,
                    matchValences=False,
                    bondCompare=rdFMCS.BondCompare.CompareAny,
                )
                if not mcs.canceled and mcs.numAtoms > 0:
                    common = Chem.MolFromSmarts(mcs.smartsString)
                    rm = ref_hvy.GetSubstructMatches(common)
                    pm = pose_hvy.GetSubstructMatches(common)
                    if rm and pm:
                        mapping = [
                            (ri, pi) for ri, pi in zip(rm[0], pm[0])
                            if ref_hvy.GetAtomWithIdx(ri).GetSymbol() == pose_hvy.GetAtomWithIdx(pi).GetSymbol()
                        ]
                        if mapping:
                            rc = ref_hvy.GetConformer()
                            pc = pose_hvy.GetConformer()
                            P = np.array([list(rc.GetAtomPosition(i)) for i, _ in mapping])
                            Q = np.array([list(pc.GetAtomPosition(j)) for _, j in mapping])
                            pos_rmsd = float(np.sqrt(((P - Q) ** 2).sum() / len(P)))
                            kab_rmsd = RMSDEngine.kabsch_rmsd(P, Q)
                            return round(pos_rmsd, 3), round(kab_rmsd, 3), len(mapping)

    # Never substitute atom-order coordinate comparison for a failed chemical
    # mapping. A partial/incorrect mapping can produce a plausible but invalid
    # validation number.
    raise ValueError(
        "RMSD mapping failed: no chemically defensible atom correspondence was found; "
        "coordinate-order fallback is disabled"
    )
