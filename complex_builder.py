#!/usr/bin/env python3
"""
AutoDock Suite Pro -- PDBQT Complex Builder
==========================================
Generates merged receptor + docked-ligand PDBQT complex files for each pose.

These files can be opened directly in:
  - UCSF Chimera / ChimeraX
  - VMD
  - Discovery Studio Visualizer
  - OpenBabel GUI
  - PyMol (File > Open, manually)
  - Any molecular viewer that reads PDBQT / PDB format

No external dependencies -- pure Python 3.

Replaces the PyMol PSE/PML generation which required the PyMol Python API to be
installed and importable. Complex PDBQT files work without any extra software.

Usage:
    from complex_builder import build_complex_pdbqt, build_job_complexes

Author: AutoDock Suite Pro Development Team
Version: 0.3.0
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger("docking_automation.complex_builder")


# ==============================================================================
# Core Complex Builder
# ==============================================================================

def _extract_pose_lines(output_pdbqt: Path, pose_index: int) -> List[str]:
    """Extract ATOM/HETATM lines for pose N from a multi-model PDBQT.

    Handles Vina (MODEL/ENDMDL) and AutoDock4 (DOCKED: prefix) formats.
    Returns raw lines for the requested pose, or empty list if not found.
    """
    try:
        text = output_pdbqt.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        logger.warning(f"Cannot read pose file {output_pdbqt}: {e}")
        return []

    lines = text.splitlines(keepends=True)

    # AutoDock4 DLG format: DOCKED: prefix (check FIRST, because DLG files also have un-prefixed MODEL lines in clustering summary)
    if any(l.startswith("DOCKED:") for l in lines):
        current_run = 0
        in_target = False
        pose_lines_ad4: List[str] = []
        for line in lines:
            if not line.startswith("DOCKED:"):
                continue
            stripped = line[8:] if line.startswith("DOCKED: ") else line[7:]
            if stripped.startswith("MODEL"):
                current_run += 1
                in_target = (current_run == pose_index)
                continue
            if stripped.startswith("ENDMDL"):
                if in_target:
                    return pose_lines_ad4
                in_target = False
                continue
            if in_target and stripped.startswith(("ATOM", "HETATM")):
                pose_lines_ad4.append(stripped)
        return pose_lines_ad4

    # Vina format: MODEL / ENDMDL blocks
    if any(l.startswith("MODEL") for l in lines):
        current_model = 0
        in_target = False
        pose_lines: List[str] = []
        for line in lines:
            if line.startswith("MODEL"):
                current_model += 1
                in_target = (current_model == pose_index)
                continue
            if line.startswith("ENDMDL"):
                if in_target:
                    return pose_lines
                in_target = False
                continue
            if in_target and line.startswith(("ATOM", "HETATM")):
                pose_lines.append(line)
        return pose_lines

    # Fallback: single-model PDBQT (already-split pose file)
    if pose_index == 1:
        return [l for l in lines if l.startswith(("ATOM", "HETATM"))]

    return []


def build_complex_pdbqt(
    receptor_path: Path,
    ligand_output_pdbqt: Path,
    output_path: Path,
    pose_index: int = 1,
    receptor_label: str = "",
    ligand_label: str = "",
) -> Optional[Path]:
    """Build a merged receptor+ligand PDBQT complex file for a single docked pose.

    The output contains:
      1. REMARK header with provenance info
      2. Receptor ATOM/HETATM records (unchanged)
      3. TER separator
      4. Ligand ATOM/HETATM records for the requested pose

    Args:
        receptor_path:        Rigid receptor PDBQT path.
        ligand_output_pdbqt:  Docked output PDBQT (multi-model from Vina or AD4).
        output_path:          Where to write the merged complex file.
        pose_index:           1-based pose index (default 1 = best pose).
        receptor_label:       Human-readable label for REMARK header.
        ligand_label:         Human-readable label for REMARK header.

    Returns:
        Path to the written file, or None on failure.
    """
    if not receptor_path or not receptor_path.exists():
        logger.warning(f"Receptor file not found: {receptor_path}")
        return None
    if not ligand_output_pdbqt or not ligand_output_pdbqt.exists():
        logger.warning(f"Ligand output file not found: {ligand_output_pdbqt}")
        return None

    # Receptor lines
    try:
        rec_text = receptor_path.read_text(encoding="utf-8", errors="replace")
        rec_lines = [
            l for l in rec_text.splitlines(keepends=True)
            if l.startswith(("ATOM", "HETATM", "TER"))
        ]
    except Exception as e:
        logger.warning(f"Cannot read receptor {receptor_path}: {e}")
        return None

    # Ligand pose lines
    lig_lines = _extract_pose_lines(ligand_output_pdbqt, pose_index)
    if not lig_lines:
        logger.warning(f"Pose {pose_index} not found in {ligand_output_pdbqt}. Complex not created.")
        return None

    # Write merged complex
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8", newline="\n") as f:
            f.write("REMARK AutoDock Suite Pro -- Docked Complex\n")
            f.write(f"REMARK Receptor : {receptor_label or receptor_path.name}\n")
            f.write(f"REMARK Ligand   : {ligand_label or ligand_output_pdbqt.stem}\n")
            f.write(f"REMARK Pose     : {pose_index}\n")
            f.write(f"REMARK Source   : {ligand_output_pdbqt}\n")
            f.write("REMARK\n")
            f.write("REMARK Open with: UCSF Chimera, ChimeraX, VMD,\n")
            f.write("REMARK           Discovery Studio, PyMol (File > Open),\n")
            f.write("REMARK           or any PDBQT/PDB-compatible viewer.\n")
            f.write("REMARK\n")
            f.writelines(rec_lines)
            if rec_lines and not rec_lines[-1].rstrip().endswith("TER"):
                f.write("TER\n")
            f.write(f"REMARK Ligand Pose {pose_index} -- BEGIN\n")
            f.writelines(lig_lines)
            f.write("END\n")

        from validators import validate_pdbqt
        valid, messages = validate_pdbqt(output_path)
        if not valid:
            logger.error("Generated complex PDBQT is invalid: %s", "; ".join(messages))
            return None

        logger.info(f"Complex written: {output_path}")
        return output_path

    except Exception as e:
        logger.warning(f"Failed to write complex {output_path}: {e}")
        return None


def convert_complex_to_pdb(pdbqt_path: Path, pdb_path: Path) -> bool:
    """Convert a merged complex PDBQT file to standard PDB format using OpenBabel.
    
    Tries OpenBabel Python bindings first, then the obabel CLI executable, and finally
    falls back to pure-Python element sanitization.
    """
    # 1. Try Python OpenBabel / Pybel
    try:
        from openbabel import pybel
        mols = list(pybel.readfile("pdbqt", str(pdbqt_path)))
        if mols:
            mols[0].write("pdb", str(pdb_path), overwrite=True)
            if pdb_path.is_file() and pdb_path.stat().st_size > 0:
                logger.debug(f"Converted complex to PDB via pybel: {pdb_path}")
                return True
    except Exception as e:
        logger.debug(f"Pybel conversion note: {e}")

    # 2. Try obabel command-line binary
    try:
        import shutil
        import subprocess
        from executables import resolve_executable
        ob_exe = resolve_executable(None, "obabel.exe")
        if not ob_exe.is_file():
            w = shutil.which("obabel") or shutil.which("obabel.exe")
            if w:
                ob_exe = Path(w)
        if ob_exe.is_file():
            cmd = [str(ob_exe), "-ipdbqt", str(pdbqt_path), "-opdb", "-O", str(pdb_path)]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if res.returncode == 0 and pdb_path.is_file() and pdb_path.stat().st_size > 0:
                logger.debug(f"Converted complex to PDB via obabel binary: {pdb_path}")
                return True
    except Exception as e:
        logger.debug(f"obabel CLI conversion note: {e}")

    # 3. Clean Python PDB converter fallback (strips PDBQT charge and atom-type columns)
    try:
        from dlg_extract import _autodock_to_element
        text = pdbqt_path.read_text(encoding="utf-8", errors="replace")
        out_lines = []
        for line in text.splitlines(keepends=True):
            if line.startswith(("ATOM", "HETATM")):
                nm = line[12:16].strip()
                raw_elem = line[77:79].strip() if len(line) >= 79 else ""
                elem = _autodock_to_element(raw_elem, nm)
                out_lines.append(f"{line[:66].ljust(76)}{elem:>2s}\n")
            elif line.startswith(("REMARK", "TER", "END")):
                out_lines.append(line)
        pdb_path.write_text("".join(out_lines), encoding="utf-8")
        return True
    except Exception as e:
        logger.warning(f"Fallback PDB conversion failed: {e}")
        return False


def build_job_complexes(
    receptor_path: Path,
    ligand_output_pdbqt: Path,
    output_dir: Path,
    n_poses: int = 3,
    receptor_label: str = "",
    ligand_label: str = "",
) -> List[Path]:
    """Build complex PDB and PDBQT files for the top N poses of a docking job.

    Creates standard .pdb complexes (via OpenBabel) as well as .pdbqt files for universal
    compatibility across Chimera, PyMOL, VMD, Discovery Studio, and web viewers.

    Args:
        receptor_path:        Rigid receptor PDBQT path.
        ligand_output_pdbqt:  Docked output PDBQT (multi-model).
        output_dir:           Directory to write complex_pose_1.pdb, etc.
        n_poses:              Number of top poses to extract (default 3).
        receptor_label:       Human-readable label for REMARK.
        ligand_label:         Human-readable label for REMARK.

    Returns:
        List of Paths for successfully written complex files (.pdb preferred).
    """
    written: List[Path] = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for i in range(1, n_poses + 1):
        out_pdbqt = output_dir / f"complex_pose_{i}.pdbqt"
        out_pdb = output_dir / f"complex_pose_{i}.pdb"

        res_pdbqt = build_complex_pdbqt(
            receptor_path=receptor_path,
            ligand_output_pdbqt=ligand_output_pdbqt,
            output_path=out_pdbqt,
            pose_index=i,
            receptor_label=receptor_label,
            ligand_label=ligand_label,
        )
        if res_pdbqt:
            convert_complex_to_pdb(out_pdbqt, out_pdb)
            if out_pdb.is_file() and out_pdb.stat().st_size > 0:
                written.append(out_pdb)
            else:
                written.append(res_pdbqt)

    if written:
        logger.info(f"Built {len(written)} complex file(s) in {output_dir}")
    return written

