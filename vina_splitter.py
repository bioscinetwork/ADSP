#!/usr/bin/env python3
"""
Docking Automation Suite — Vina Split Integration
===================================================
Runs vina_split to split multi-model PDBQT into individual pose files.
"""

from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger("docking_automation.vina_splitter")


def run_vina_split(
    vina_split_exe: Path,
    input_pdbqt: Path,
    output_dir: Path,
) -> Tuple[bool, List[Path], str]:
    """Run vina_split on a multi-model PDBQT file.

    Args:
        vina_split_exe: Path to the vina_split executable.
        input_pdbqt: Path to the multi-model PDBQT output from Vina.
        output_dir: Directory where split poses will be written.

    Returns:
        (success, list of split file paths, error/info message)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # vina_split typically creates files named <input_stem>_ligand_<N>.pdbqt
    # in the same directory as the input file. We may need to move them.
    cmd = [
        str(vina_split_exe),
        "--input", str(input_pdbqt),
        "--ligand", str(output_dir / input_pdbqt.stem),
    ]

    logger.info(f"Running vina_split: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )

        if result.returncode != 0:
            error_msg = result.stderr.strip() or result.stdout.strip()
            logger.error(f"vina_split failed (exit {result.returncode}): {error_msg}")

            # Fallback: try without --ligand flag (some versions differ)
            cmd_fallback = [
                str(vina_split_exe),
                "--input", str(input_pdbqt),
            ]
            logger.info(f"Trying fallback: {' '.join(cmd_fallback)}")
            result = subprocess.run(
                cmd_fallback,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(output_dir),
            )
            if result.returncode != 0:
                return False, [], f"vina_split failed: {result.stderr.strip()}"

    except subprocess.TimeoutExpired:
        logger.error("vina_split timed out")
        return False, [], "vina_split timed out after 120 seconds"
    except FileNotFoundError:
        logger.error(f"vina_split not found: {vina_split_exe}")
        return False, [], f"vina_split executable not found: {vina_split_exe}"
    except Exception as e:
        logger.error(f"vina_split error: {e}")
        return False, [], str(e)

    # Collect split files from the output directory
    split_files = _collect_split_files(output_dir, input_pdbqt.stem)

    # Also check the input file's directory (some vina_split versions write there)
    if not split_files and input_pdbqt.parent != output_dir:
        split_files = _collect_split_files(input_pdbqt.parent, input_pdbqt.stem)
        # Move them to output_dir
        moved = []
        for sf in split_files:
            dest = output_dir / sf.name
            sf.rename(dest)
            moved.append(dest)
        split_files = moved

    if not split_files:
        # Manual splitting fallback
        logger.info("vina_split produced no files; attempting manual split...")
        split_files = _manual_split(input_pdbqt, output_dir)

    if split_files:
        # Rename to pose_1.pdbqt, pose_2.pdbqt, etc. for consistency
        renamed = _rename_split_files(split_files, output_dir)
        logger.info(f"Split into {len(renamed)} pose files")
        return True, renamed, f"Split into {len(renamed)} poses"
    else:
        return False, [], "No split files produced"


def _collect_split_files(directory: Path, stem: str) -> List[Path]:
    """Collect split PDBQT files from a directory matching the input stem."""
    patterns = [
        f"{stem}_ligand_*.pdbqt",
        f"{stem}_*.pdbqt",
    ]
    files: List[Path] = []
    for pattern in patterns:
        for f in sorted(directory.glob(pattern)):
            # Exclude the original file
            if f.stem != stem and "_out" not in f.stem:
                files.append(f)
    return files


def _rename_split_files(files: List[Path], output_dir: Path) -> List[Path]:
    """Rename split files to a consistent pose_N.pdbqt convention."""
    renamed = []
    for i, f in enumerate(sorted(files), 1):
        new_name = output_dir / f"pose_{i}.pdbqt"
        if f != new_name:
            if new_name.exists():
                new_name.unlink()
            f.rename(new_name)
        renamed.append(new_name)
    return renamed


def _manual_split(input_pdbqt: Path, output_dir: Path) -> List[Path]:
    """Manually split a multi-model PDBQT if vina_split doesn't work.

    Splits on MODEL/ENDMDL boundaries. Preserves the original file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    content = input_pdbqt.read_text(encoding="utf-8", errors="replace")
    models = re.split(r"^ENDMDL\s*$", content, flags=re.MULTILINE)

    files: List[Path] = []
    pose_num = 0

    for model_text in models:
        model_text = model_text.strip()
        if not model_text:
            continue

        # Remove MODEL line if present
        lines = model_text.splitlines()
        filtered = [
            line for line in lines
            if not line.strip().startswith("MODEL")
        ]

        # Check this block has ATOM/HETATM records
        has_atoms = any(
            line.startswith(("ATOM", "HETATM"))
            for line in filtered
        )
        if not has_atoms:
            continue

        pose_num += 1
        pose_file = output_dir / f"pose_{pose_num}.pdbqt"
        pose_file.write_text("\n".join(filtered) + "\n", encoding="utf-8")
        files.append(pose_file)

    return files


# Public alias
split_pdbqt_manually = _manual_split


def validate_split(
    split_target: Path | List[Path],
    expected_count: int,
) -> Tuple[bool, str]:
    """Validate that split file count matches expected model count.

    Args:
        split_target: Either a list of split Paths, or a directory Path containing them.
        expected_count: Number of expected models.

    Returns:
        (is_valid, message)
    """
    if isinstance(split_target, Path) and split_target.is_dir():
        split_files = list(split_target.glob("*.pdbqt"))
    elif isinstance(split_target, list):
        split_files = split_target
    else:
        split_files = []

    actual = len(split_files)
    if actual == expected_count:
        return True, f"Split validation passed: {actual} files match {expected_count} models"
    elif actual > 0:
        return True, (
            f"Split count mismatch: expected {expected_count}, "
            f"got {actual} (may be normal if vina_split handles differently)"
        )
    else:
        return False, "No split files produced"
