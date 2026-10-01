#!/usr/bin/env python3
"""
AutoDock Suite Pro — AutoDock 4 Workflow
=========================================
Full orchestration of AutoDock 4 docking: GPF/DPF generation,
AutoGrid4 execution (with GLG grid log), AutoDock4 execution (DLG),
and DLG parsing via dlg_extract.py (BSNDVP™).

Folder Convention:
    receptors/RECEPTOR/rigid/receptor.pdbqt   ← pre-prepared PDBQT
    receptors/RECEPTOR/config.txt             ← grid center/size (Vina fmt)
    ligands/LIGAND.pdbqt                      ← pre-prepared PDBQT

Usage:
    python autodock4_workflow.py --config project_config.toml
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys

# Configure UTF-8 output on Windows consoles
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import shutil
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from config import ProjectConfig, load_config, validate_config
from ad4_compatibility import (
    STANDARD_AD4_PROFILE,
    AtomTypeCompatibilityError,
    extract_atom_types as _compat_extract_atom_types,
    preflight_job as preflight_ad4_job,
)
from executables import (
    get_autodock4_version, get_autogrid4_version, resolve_executable, validate_executable,
)
from job_manager import (
    ProgressTracker, build_job_queue, filter_jobs,
    load_job_status, run_jobs_parallel, save_job_status, validate_job_inputs,
)
from logging_utils import get_run_timestamp, setup_global_logger, setup_job_logger
from models import (
    SUITE_BANNER, AnalysisStatus, AutoDock4Metrics, CanonicalPose,
    ClusterInfo, DockingJob, DockingMode, DockingResult, Engine,
    ExecutionStatus, JobStatus, ProvenanceRecord, ResumeMode,
    ThermodynamicAnalysis, calculate_inhibition_constant, __version__,
)
from process_manager import compute_file_hash, process_manager, terminate_process

logger = logging.getLogger("docking_automation")


# ═══════════════════════════════════════════════════════════════════════════════
# GPF / DPF Template Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_gpf(
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    config: Optional[ProjectConfig] = None,
    grid_center: Optional[Dict[str, float]] = None,
    grid_size: Optional[Dict[str, float]] = None,
    flex_receptor_path: Optional[Path] = None,
    spacing: float = 0.375,
    smooth: float = 0.5,
    dielectric: float = -0.1465,
) -> Path:
    """Generate an AutoGrid4 Grid Parameter File (GPF).

    Args:
        receptor_path: Path to the rigid receptor PDBQT.
        ligand_path: Path to the ligand PDBQT.
        output_path: Where to write the GPF.
        config: Project configuration.
        grid_center: Optional {x, y, z} center override (Angstroms).
        grid_size: Optional {x, y, z} size override (Angstroms).
        flex_receptor_path: Optional flexible sidechain PDBQT for AD4 flexible docking.
        spacing: Grid spacing in Angstroms (default 0.375).
        smooth: Smoothing factor (default 0.5).
        dielectric: Dielectric constant (default -0.1465).

    Returns:
        Path to the generated GPF.
    """
    try:
        ensure_receptor_has_charges(receptor_path)
    except Exception:
        pass

    # Support single ligand path or list of ligand paths to compute union of atom types
    if isinstance(ligand_path, (list, tuple)):
        all_lig_types = []
        for lp in ligand_path:
            all_lig_types.extend(_extract_atom_types(Path(lp)))
        ligand_types = sorted(set(all_lig_types))
    else:
        ligand_types = _extract_atom_types(Path(ligand_path))
    receptor_types = _extract_atom_types(receptor_path)

    # In flexible docking, AutoGrid4 must compute maps for all moving atoms (ligand + flexres)
    if flex_receptor_path and flex_receptor_path.exists():
        flex_types = _extract_atom_types(flex_receptor_path)
        moving_types = sorted(set(ligand_types + flex_types))
    else:
        moving_types = sorted(set(ligand_types))

    moving_types_str = " ".join(moving_types)

    # Grid center and size (from config.txt or GUI in Angstroms)
    cx = grid_center.get("x", 0.0) if grid_center else 0.0
    cy = grid_center.get("y", 0.0) if grid_center else 0.0
    cz = grid_center.get("z", 0.0) if grid_center else 0.0
    sx = float(grid_size.get("x", 25.0)) if grid_size else 25.0
    sy = float(grid_size.get("y", 25.0)) if grid_size else 25.0
    sz = float(grid_size.get("z", 25.0)) if grid_size else 25.0

    grid_spacing = float(spacing) if spacing and float(spacing) > 0 else 0.375

    # AutoGrid4 npts = size_in_angstrom / spacing.
    # Must be an even integer clamped to [20, 126] for standard AutoGrid4.
    def _calc_npts(val_angstrom: float, sp: float, axis: str) -> int:
        raw_pts = (val_angstrom / sp) if sp > 0 else val_angstrom
        pts = int(round(raw_pts))
        if pts % 2 != 0:
            pts += 1
        clamped_pts = max(20, min(126, pts))
        if abs(clamped_pts - raw_pts) > 1e-4:
            actual_size = clamped_pts * sp
            logger.info(
                f"[AutoGrid4 Grid Transparency] Axis {axis}: Requested {val_angstrom:.2f} Å ({raw_pts:.2f} pts) "
                f"adjusted to {clamped_pts} pts ({actual_size:.2f} Å). Reason: AutoGrid4 requires an even integer in [20, 126]."
            )
        return clamped_pts

    npts_x = _calc_npts(sx, grid_spacing, "X")
    npts_y = _calc_npts(sy, grid_spacing, "Y")
    npts_z = _calc_npts(sz, grid_spacing, "Z")

    gpf_lines = [
        f"# AutoGrid4 Parameter File — AutoDock Suite Pro",
        f"# Requested grid: {sx:.2f}x{sy:.2f}x{sz:.2f} Å (spacing {grid_spacing:.3f} Å)",
        f"# Actual grid   : {npts_x}x{npts_y}x{npts_z} pts ({npts_x*grid_spacing:.2f}x{npts_y*grid_spacing:.2f}x{npts_z*grid_spacing:.2f} Å)",
        f"npts {npts_x} {npts_y} {npts_z}",
        f"gridfld {receptor_path.stem}.maps.fld",
        f"spacing {grid_spacing:.3f}",
        f"receptor_types {' '.join(sorted(set(receptor_types)))}",
        f"ligand_types {moving_types_str}",
        f"receptor {receptor_path.name}",
        f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}",
        f"smooth {float(smooth):.1f}",
    ]

    # Explicitly bind generated standard maps to the selected parameter asset.
    # AutoDock4Zn syntax remains isolated from this ordinary generator.
    try:
        from ad4_compatibility import PARAMETER_PROFILES
        profile_id = getattr(config, "ad4_parameter_profile", "ad4_standard_4.2") if config else "ad4_standard_4.2"
        profile = PARAMETER_PROFILES.get(profile_id)
        if profile and profile.parameter_file and profile.parameter_file.lower().endswith(".dat"):
            parameter_path = Path(profile.parameter_file)
            if not parameter_path.is_absolute() and profile.source_path:
                parameter_path = Path(profile.source_path) / parameter_path
            gpf_lines.append(f"parameter_file {parameter_path.resolve()}")
    except (OSError, ValueError):
        logger.warning("Could not resolve selected AutoDock4 parameter file for GPF.")

    # Individual map lines for each moving atom type
    for at in moving_types:
        gpf_lines.append(f"map {receptor_path.stem}.{at}.map")

    gpf_lines += [
        f"elecmap {receptor_path.stem}.e.map",
        f"dsolvmap {receptor_path.stem}.d.map",
        f"dielectric {float(dielectric):.4f}",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Ensure receptor PDBQT exists in the same directory as the GPF for AutoGrid4
    rec_target = output_path.parent / receptor_path.name
    if not rec_target.is_file() and receptor_path.is_file():
        shutil.copy2(receptor_path, rec_target)
    if rec_target.is_file():
        try:
            ensure_receptor_has_charges(rec_target)
        except Exception:
            pass
    output_path.write_text("\n".join(gpf_lines) + "\n", encoding="utf-8")
    return output_path


def generate_gpf_file(
    receptor_pdbqt: Union[Path, str],
    output_gpf: Union[Path, str],
    center: Union[Tuple[float, float, float], List[float], Dict[str, float]],
    size: Union[Tuple[float, float, float], List[float], Dict[str, float]],
    spacing: float = 0.375,
    ligand_pdbqt: Optional[Union[Path, str, List[Union[Path, str]]]] = None,
    ligand_types: Optional[List[str]] = None,
    flex_receptor_path: Optional[Union[Path, str]] = None,
    smooth: float = 0.5,
    dielectric: float = -0.1465,
) -> Path:
    """Generate an AutoGrid4 GPF file with flexible parameter input formats.

    Supports center/size as tuples, lists, or dicts, and handles automatic
    ligand atom typing or standard virtual screening fallback types.
    """
    rec_path = Path(receptor_pdbqt)
    out_path = Path(output_gpf)

    # Normalize center
    if isinstance(center, (list, tuple)) and len(center) >= 3:
        grid_center = {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])}
    elif isinstance(center, dict):
        grid_center = {k: float(v) for k, v in center.items()}
    else:
        grid_center = {"x": 0.0, "y": 0.0, "z": 0.0}

    # Normalize size
    if isinstance(size, (list, tuple)) and len(size) >= 3:
        grid_size = {"x": float(size[0]), "y": float(size[1]), "z": float(size[2])}
    elif isinstance(size, dict):
        grid_size = {k: float(v) for k, v in size.items()}
    else:
        grid_size = {"x": 25.0, "y": 25.0, "z": 25.0}

    # Determine ligand types
    types: List[str] = []
    if ligand_types:
        types = list(ligand_types)
    elif ligand_pdbqt:
        if isinstance(ligand_pdbqt, (list, tuple)):
            for lp in ligand_pdbqt:
                types.extend(_extract_atom_types(Path(lp)))
        else:
            types.extend(_extract_atom_types(Path(ligand_pdbqt)))

    # If no ligand types found or specified, try scanning ligands/ or fallback to standard screening types
    if not types:
        search_dirs = [
            out_path.parent.parent / "ligands",
            out_path.parent / "ligands",
            Path("ligands"),
        ]
        found_ligands = []
        for sdir in search_dirs:
            if sdir.is_dir():
                found_ligands = list(sdir.glob("*.pdbqt"))
                if found_ligands:
                    break
        if found_ligands:
            for lp in found_ligands[:50]:
                types.extend(_extract_atom_types(lp))

    if not types:
        # Standard AutoDock 4 virtual screening atom types
        types = ["A", "Br", "C", "Cl", "F", "HD", "I", "N", "NA", "OA", "P", "S", "SA"]

    types = sorted(set(types))

    flex_path = Path(flex_receptor_path) if flex_receptor_path else None

    # Receptor types
    try:
        ensure_receptor_has_charges(rec_path)
    except Exception:
        pass

    receptor_types = _extract_atom_types(rec_path)
    if not receptor_types:
        receptor_types = ["C", "A", "N", "O", "S", "HD", "H", "NA", "OA", "SA"]

    if flex_path and flex_path.exists():
        flex_types = _extract_atom_types(flex_path)
        moving_types = sorted(set(types + flex_types))
    else:
        moving_types = sorted(set(types))

    moving_types_str = " ".join(moving_types)

    grid_spacing = float(spacing) if spacing and float(spacing) > 0 else 0.375

    def _calc_npts(val_angstrom: float, sp: float) -> int:
        pts = int(round(val_angstrom / sp)) if sp > 0 else int(val_angstrom)
        if pts % 2 != 0:
            pts += 1
        return max(20, min(126, pts))

    npts_x = _calc_npts(grid_size.get("x", 25.0), grid_spacing)
    npts_y = _calc_npts(grid_size.get("y", 25.0), grid_spacing)
    npts_z = _calc_npts(grid_size.get("z", 25.0), grid_spacing)

    cx = grid_center.get("x", 0.0)
    cy = grid_center.get("y", 0.0)
    cz = grid_center.get("z", 0.0)

    gpf_lines = [
        f"npts {npts_x} {npts_y} {npts_z}",
        f"gridfld {rec_path.stem}.maps.fld",
        f"spacing {grid_spacing:.3f}",
        f"receptor_types {' '.join(sorted(set(receptor_types)))}",
        f"ligand_types {moving_types_str}",
        f"receptor {rec_path.name}",
        f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}",
        f"smooth {float(smooth):.1f}",
    ]

    for at in moving_types:
        gpf_lines.append(f"map {rec_path.stem}.{at}.map")

    gpf_lines += [
        f"elecmap {rec_path.stem}.e.map",
        f"dsolvmap {rec_path.stem}.d.map",
        f"dielectric {float(dielectric):.4f}",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    rec_target = out_path.parent / rec_path.name
    if not rec_target.is_file() and rec_path.is_file():
        try:
            shutil.copy2(rec_path, rec_target)
        except Exception:
            pass

    out_path.write_text("\n".join(gpf_lines) + "\n", encoding="utf-8")
    return out_path


def _calculate_ligand_center(ligand_path: Path) -> Tuple[float, float, float]:
    """Calculate the geometric center of the ligand's coordinates (the 'about' center in AutoDock 4).

    AutoDock 4 requires the 'about' keyword to specify the small molecule's own center of rotation
    (in its local coordinates), NOT the receptor grid center. Rotating around the receptor grid
    center shifts the molecule outside the grid box during docking and produces massive positive
    energy penalties (+5e+08 kcal/mol).
    """
    xs, ys, zs = [], [], []
    root_xs, root_ys, root_zs = [], [], []
    in_root = False
    try:
        for line in ligand_path.read_text(encoding="utf-8", errors="replace").splitlines():
            line_s = line.strip()
            if line_s == "ROOT":
                in_root = True
                continue
            elif line_s == "ENDROOT":
                in_root = False
                continue
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                    if in_root:
                        root_xs.append(x)
                        root_ys.append(y)
                        root_zs.append(z)
                except (ValueError, IndexError):
                    pass
    except Exception:
        pass
    if root_xs:
        return (sum(root_xs) / len(root_xs), sum(root_ys) / len(root_ys), sum(root_zs) / len(root_zs))
    if xs:
        return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    return (0.0, 0.0, 0.0)


def generate_dpf(
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    config: ProjectConfig,
    grid_center: Optional[Dict[str, float]] = None,
    flex_receptor_path: Optional[Path] = None,
    algorithm: str = "LGA",
) -> Path:
    """Generate an AutoDock4 Docking Parameter File (DPF).

    Supports:
        - LGA (Lamarckian Genetic Algorithm): GA + Pseudo-Solis-Wets local search
        - GA (Genetic Algorithm): standard GA without local search
        - LS (Local Search): pure local search optimization
        - Flexible receptor sidechains via 'flexres' directive
    """
    ligand_types = _extract_atom_types(ligand_path)
    if flex_receptor_path and flex_receptor_path.exists():
        flex_types = _extract_atom_types(flex_receptor_path)
        moving_types = sorted(set(ligand_types + flex_types))
    else:
        moving_types = sorted(set(ligand_types))

    moving_types_str = " ".join(moving_types)
    n_torsions = _count_torsions(ligand_path)
    about_x, about_y, about_z = _calculate_ligand_center(ligand_path)

    # Seed
    seed1, seed2 = "pid", "time"
    if config.ad4_seed.upper() != "AUTO":
        parts = config.ad4_seed.split()
        if len(parts) >= 2:
            seed1, seed2 = parts[0], parts[1]
        else:
            seed1, seed2 = parts[0], parts[0]

    dpf_lines = [
        f"autodock_parameter_version 4.2",
        f"outlev {getattr(config, 'outlev', 1)}",
        f"intelec",
        f"seed {seed1} {seed2}",
        f"ligand_types {moving_types_str}",
        f"fld {receptor_path.stem}.maps.fld",
    ]
    try:
        from ad4_compatibility import PARAMETER_PROFILES
        profile_id = getattr(config, "ad4_parameter_profile", "ad4_standard_4.2")
        profile = PARAMETER_PROFILES.get(profile_id)
        if profile and profile.parameter_file and profile.parameter_file.lower().endswith(".dat"):
            parameter_path = Path(profile.parameter_file)
            if not parameter_path.is_absolute() and profile.source_path:
                parameter_path = Path(profile.source_path) / parameter_path
            dpf_lines.append(f"parameter_file {parameter_path.resolve()}")
    except (OSError, ValueError):
        logger.warning("Could not resolve selected AutoDock4 parameter file for DPF.")

    # Map references for all moving atom types
    for at in moving_types:
        dpf_lines.append(f"map {receptor_path.stem}.{at}.map")

    dpf_lines += [
        f"elecmap {receptor_path.stem}.e.map",
        f"desolvmap {receptor_path.stem}.d.map",
        f"move {ligand_path.name}",
        f"about {about_x:.3f} {about_y:.3f} {about_z:.3f}",
        f"tran0 random",
        f"quaternion0 random",
        f"dihe0 random",
        f"torsdof {n_torsions}",
    ]

    # Flexible residue support (AutoDock 4 flexres)
    if flex_receptor_path and flex_receptor_path.exists():
        dpf_lines.append(f"flexres {flex_receptor_path.name}")

    dpf_lines.append("")

    # Algorithm selection
    alg = algorithm.upper()
    if alg in ("LGA", "GA"):
        dpf_lines += [
            f"ga_pop_size {config.ga_pop_size}",
            f"ga_num_evals {config.ga_num_evals}",
            f"ga_num_generations {config.ga_num_generations}",
            f"ga_elitism {getattr(config, 'ga_elitism', 1)}",
            f"ga_mutation_rate {getattr(config, 'ga_mutation_rate', 0.02)}",
            f"ga_crossover_rate {getattr(config, 'ga_crossover_rate', 0.8)}",
            f"ga_window_size {getattr(config, 'ga_window_size', 10)}",
            f"ga_cauchy_alpha {getattr(config, 'ga_cauchy_alpha', 0.0)}",
            f"ga_cauchy_beta {getattr(config, 'ga_cauchy_beta', 1.0)}",
            f"set_ga",
            f"",
        ]
        if alg == "LGA":
            # Lamarckian Genetic Algorithm with Pseudo-Solis-Wets Local Search
            dpf_lines += [
                f"sw_max_its {getattr(config, 'sw_max_its', 300)}",
                f"sw_max_succ {getattr(config, 'sw_max_succ', 4)}",
                f"sw_max_fail {getattr(config, 'sw_max_fail', 4)}",
                f"sw_rho {getattr(config, 'sw_rho', 1.0)}",
                f"sw_lb_rho {getattr(config, 'sw_lb_rho', 0.01)}",
                f"ls_search_freq {getattr(config, 'ls_search_freq', 0.06)}",
                f"set_psw1",
                f"",
            ]
        else:
            # Pure GA without local search
            dpf_lines += [
                f"ls_search_freq 0.0",
                f"",
            ]
    elif alg == "LS":
        # Pure Local Search
        dpf_lines += [
            f"sw_max_its {getattr(config, 'sw_max_its', 300)}",
            f"sw_max_succ {getattr(config, 'sw_max_succ', 4)}",
            f"sw_max_fail {getattr(config, 'sw_max_fail', 4)}",
            f"sw_rho {getattr(config, 'sw_rho', 1.0)}",
            f"sw_lb_rho {getattr(config, 'sw_lb_rho', 0.01)}",
            f"do_local_search 1",
            f"",
        ]

    dpf_lines += [
        f"unbound_model {getattr(config, 'unbound_model', 'bound')}",
        f"ga_run {config.ga_run}",
        f"",
        f"rmstol {config.rmstol}",
        f"analysis",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(dpf_lines) + "\n", encoding="utf-8")
    return output_path


def _extract_atom_types(pdbqt_path: Path) -> List[str]:
    """Extract AutoDock atom types through the central compatibility subsystem."""
    return _compat_extract_atom_types(pdbqt_path)


def _count_torsions(pdbqt_path: Path) -> int:
    """Count active torsions from PDBQT TORSDOF/REMARK lines."""
    if not pdbqt_path.exists():
        return 0
    for line in pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("TORSDOF"):
            try:
                return int(line.split()[1])
            except (IndexError, ValueError):
                pass
    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# Grid Center Extraction from Vina Config
# ═══════════════════════════════════════════════════════════════════════════════

def extract_grid_from_config(config_path: Optional[Path]) -> Dict[str, Dict[str, float]]:
    """Extract grid center and size from a Vina-style config.txt.

    Returns dict with 'center' and 'size' sub-dicts.
    """
    result: Dict[str, Dict[str, float]] = {
        "center": {"x": 0.0, "y": 0.0, "z": 0.0},
        "size": {"x": 60.0, "y": 60.0, "z": 60.0},
    }
    if not config_path or not config_path.exists():
        return result

    from validators import parse_vina_config
    params = parse_vina_config(config_path)

    for axis in ["x", "y", "z"]:
        if f"center_{axis}" in params:
            result["center"][axis] = float(params[f"center_{axis}"])
        if f"size_{axis}" in params:
            # Convert Angstroms to grid points (0.375 spacing)
            result["size"][axis] = float(params[f"size_{axis}"])

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Single Job Execution
# ═══════════════════════════════════════════════════════════════════════════════

def run_ad4_job(job: DockingJob, config: ProjectConfig) -> DockingJob:
    """Execute a single AutoDock4 docking job.

    Steps:
        1. Validate inputs
        2. Setup output directory
        3. Extract grid parameters from config.txt
        4. Generate GPF (Grid Parameter File)
        5. Generate DPF (Docking Parameter File)
        6. Copy receptor & ligand to output dir
        7. Run AutoGrid4  → produces .map files + .glg grid log
        8. Run AutoDock4  → produces .dlg docking log
        9. Validate DLG output
        10. Record results (gpf_path, glg_path, dpf_path, dlg_path on job)
    """
    job_log = setup_job_logger(config.log_directory, "AD4", job.job_id)
    if process_manager.is_cancelled():
        job.status = JobStatus.CANCELLED
        job.execution_status = ExecutionStatus.CANCELLED
        job_log.warning("Job cancelled before execution started.")
        return job

    job.status = JobStatus.RUNNING
    job.execution_status = ExecutionStatus.RUNNING

    # Step 1: Validate inputs
    input_errors = validate_job_inputs(job)
    if input_errors:
        for err in input_errors:
            job_log.error(err)
        job.errors.extend(input_errors)
        job.status = JobStatus.INVALID_INPUT
        job.execution_status = ExecutionStatus.FAILED
        return job

    # Validate atom types before emitting parameter files or starting AutoGrid4.
    compatibility_files = [("receptor", job.receptor_path), ("ligand", job.ligand_path)]
    if job.flex_receptor_path and job.flex_receptor_path.exists():
        compatibility_files.append(("flex_receptor", job.flex_receptor_path))
    try:
        compatibility = preflight_ad4_job(
            compatibility_files,
            profile_id=getattr(config, "ad4_parameter_profile", STANDARD_AD4_PROFILE.profile_id),
            engine="AUTODOCK4",
            required_map_types=sorted({
                atom_type for _, path in compatibility_files
                for atom_type in _extract_atom_types(path)
            }),
        )
        job.ad4_results["compatibility"] = compatibility
        job_log.info(
            "AD4 compatibility preflight passed: "
            f"profile={compatibility['parameter_profile']['profile_id']}; "
            f"metal atoms={len(compatibility['metals']['detected'])}; "
            "specialized metal workflow=not selected"
        )
    except AtomTypeCompatibilityError as exc:
        job.ad4_results["compatibility"] = exc.provenance
        job.errors.extend(exc.errors)
        for error in exc.errors:
            job_log.error(error)
        job.status = JobStatus.INVALID_INPUT
        job.execution_status = ExecutionStatus.FAILED
        return job

    # Step 2: Setup output directories
    output_dir = job.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 3: Extract grid parameters
    grid_params = extract_grid_from_config(job.config_path)

    # Determine receptor naming prefix (ensures flexible docking outputs never collide with rigid ones)
    is_flex = (job.docking_mode == DockingMode.FLEXIBLE)
    rec_prefix = f"{job.receptor_name}_flex" if is_flex and not job.receptor_name.lower().endswith("_flex") else job.receptor_name

    # Step 4: Generate GPF
    gpf_path = output_dir / f"{rec_prefix}.gpf"
    generate_gpf(
        job.receptor_path, job.ligand_path, gpf_path,
        config, grid_params["center"], grid_params["size"],
        flex_receptor_path=job.flex_receptor_path,
    )
    job.gpf_path = gpf_path
    job_log.info(f"Generated GPF: {gpf_path}")

    # Step 5: Generate DPF
    algorithm = getattr(config, "ad4_algorithm", "LGA")
    dpf_path = output_dir / f"{rec_prefix}_{job.ligand_name}.dpf"
    generate_dpf(
        job.receptor_path, job.ligand_path, dpf_path,
        config, grid_params["center"],
        flex_receptor_path=job.flex_receptor_path,
        algorithm=algorithm,
    )
    job.dpf_path = dpf_path
    job_log.info(f"Generated DPF: {dpf_path} (algorithm: {algorithm})")

    # Step 6: Copy receptor and ligand to output dir
    shutil.copy2(job.receptor_path, output_dir / job.receptor_path.name)
    shutil.copy2(job.ligand_path, output_dir / job.ligand_path.name)
    if job.flex_receptor_path and job.flex_receptor_path.exists():
        shutil.copy2(job.flex_receptor_path, output_dir / job.flex_receptor_path.name)

    # Step 7: Run AutoGrid4 → .map files + .glg grid log
    glg_path = output_dir / f"{rec_prefix}.glg"
    job.glg_path = glg_path
    if not config.reuse_existing_maps:
        success = _run_autogrid4(
            config.autogrid4_executable, gpf_path, output_dir, job, job_log
        )
        if not success:
            return job
        job_log.info(f"Grid log (GLG): {glg_path}")
    else:
        job_log.info("Reusing existing grid maps (skipping AutoGrid4)")

    # Step 8: Run AutoDock4 → .dlg docking log
    dlg_path = output_dir / f"{rec_prefix}_{job.ligand_name}.dlg"
    job.dlg_path = dlg_path

    success = _run_autodock4(
        config.autodock4_executable, dpf_path, dlg_path, output_dir, job, job_log
    )
    if not success:
        return job

    job.execution_status = ExecutionStatus.SUCCESS

    # Step 9: Validate DLG
    from validators import validate_dlg
    valid, messages = validate_dlg(dlg_path)
    if not valid:
        job_log.error(f"DLG validation failed: {messages}")
        job.errors.extend(messages)
        job.status = JobStatus.PARSING_FAILED
        job.analysis_status = AnalysisStatus.FAILED
        return job

    job_log.info(f"DLG validated successfully: {dlg_path}")

    # Clean up transient AutoGrid map files to eliminate disk bloat (saves 50-100MB per run)
    if not getattr(config, "reuse_existing_maps", False):
        for pattern in ("*.map", "*.maps.fld", "*.maps.xyz"):
            for map_file in output_dir.glob(pattern):
                try:
                    map_file.unlink(missing_ok=True)
                except Exception as e:
                    job_log.debug(f"Could not remove transient map file {map_file}: {e}")

    # Step 10: Record results
    job.ad4_results["dlg_path"] = str(dlg_path)
    job.ad4_results["gpf_path"] = str(gpf_path)
    job.ad4_results["glg_path"] = str(glg_path)
    job.ad4_results["dpf_path"] = str(dpf_path)
    try:
        from dlg_extract import DLGParser
        parsed_dlg = DLGParser().parse(dlg_path)
        if parsed_dlg and parsed_dlg.poses:
            max_runs = getattr(config, "ga_run", None)
            if max_runs and len(parsed_dlg.poses) > max_runs:
                dlg_poses = parsed_dlg.poses[:max_runs]
            else:
                dlg_poses = parsed_dlg.poses

            # 1. Backward-compatible ad4_results dict
            if not isinstance(job.ad4_results, dict):
                job.ad4_results = {}
            job.ad4_results["poses"] = [
                {
                    "rank": p.rank,
                    "mode": p.rank,
                    "pose": p.rank,
                    "binding_affinity": p.binding_energy,
                    "binding_energy": p.binding_energy,
                    "ki_nM": p.ki_nM,
                    "ki_raw": p.ki_raw,
                    "ki_unit": p.ki_unit,
                    "intermol_energy": p.intermol_energy,
                    "internal_energy": p.internal_energy,
                    "torsional_energy": p.torsional_energy,
                    "unbound_energy": p.unbound_energy,
                    "cluster_id": p.cluster_id,
                    "cluster_size": p.cluster_size,
                    "cluster_rmsd": p.cluster_rmsd,
                    "rmsd_from_ref": p.rmsd_from_ref,
                }
                for p in dlg_poses
            ]
            job.ad4_results["best_energy"] = dlg_poses[0].binding_energy
            job.obtained_modes = len(dlg_poses)

            # 2. Canonical DockingResult model
            canonical_poses: List[CanonicalPose] = []
            for p in dlg_poses:
                ad4_met = AutoDock4Metrics(
                    intermolecular_energy=p.intermol_energy,
                    internal_energy=p.internal_energy,
                    torsional_energy=p.torsional_energy,
                    unbound_energy=p.unbound_energy,
                    vdW_hbond_desolvation_energy=p.vdW_hbond_desolvation_energy,
                    electrostatic_energy=p.electrostatic_energy,
                    dlg_source=p.dlg_source,
                    cluster_id=p.cluster_id,
                    cluster_size=p.cluster_size,
                    cluster_population=(p.cluster_size / parsed_dlg.num_runs * 100.0) if (parsed_dlg.num_runs and p.cluster_size) else None,
                    cluster_rmsd=p.cluster_rmsd,
                    rmsd_from_reference=p.rmsd_from_ref if p.rmsd_from_ref != 0.0 else None,
                )
                val_rmsd = p.validation_rmsd_conf or p.validation_rmsd_pos
                if val_rmsd is None and p.rmsd_from_ref != 0.0:
                    val_rmsd = p.rmsd_from_ref

                ki_fmt = f"{p.ki_raw:.2f} {p.ki_unit}" if p.ki_raw is not None else None
                canonical_poses.append(CanonicalPose(
                    rank=p.rank,
                    run_number=p.run_number if p.run_number else p.rank,
                    binding_energy=p.binding_energy,
                    estimated_ki_nM=p.ki_nM,
                    estimated_ki_formatted=ki_fmt,
                    validation_rmsd=val_rmsd,
                    ad4_metrics=ad4_met,
                    source_model=p.rank,
                    source_run=p.run_number if p.run_number else p.rank,
                    source_file=str(dlg_path),
                ))

            cluster_infos: List[ClusterInfo] = []
            for c in getattr(parsed_dlg, "clusters", []):
                pop_pct = (c.size / parsed_dlg.num_runs * 100.0) if (parsed_dlg.num_runs and c.size) else None
                cluster_infos.append(ClusterInfo(
                    cluster_id=c.cluster_id,
                    size=c.size,
                    population_percent=pop_pct,
                    lowest_energy=c.lowest_energy,
                    mean_energy=c.mean_energy,
                    cluster_rmsd=c.avg_rmsd,
                    representative_run=c.best_run,
                    representative_pose=c.best_run,
                    runs=list(getattr(c, "runs", [])),
                ))

            thermo = ThermodynamicAnalysis(
                info_entropy=getattr(parsed_dlg, "info_entropy", None),
                partition_function=getattr(parsed_dlg, "partition_function", None),
                stat_temperature=getattr(parsed_dlg, "stat_temperature", None),
                stat_free_energy=getattr(parsed_dlg, "stat_free_energy", None),
                stat_internal_energy=getattr(parsed_dlg, "stat_internal_energy", None),
                stat_entropy=getattr(parsed_dlg, "stat_entropy", None),
            )

            provenance = ProvenanceRecord(
                app_version=__version__,
                engine="AutoDock4",
                engine_version="4.2",
                docking_mode=job.docking_mode.value,
                receptor_file=str(job.receptor_path),
                receptor_hash=compute_file_hash(job.receptor_path) if job.receptor_path.exists() else "",
                ligand_file=str(job.ligand_path),
                ligand_hash=compute_file_hash(job.ligand_path) if job.ligand_path.exists() else "",
                config_file=str(job.config_path) if job.config_path else "",
                config_hash=compute_file_hash(job.config_path) if job.config_path and job.config_path.exists() else "",
                executable_path=str(config.autodock4_executable),
                output_files=[str(f) for f in [dlg_path, gpf_path, dpf_path, glg_path] if f and f.exists()],
            )

            job.canonical_result = DockingResult(
                job_id=job.job_id,
                engine=Engine.AUTODOCK4,
                docking_mode=job.docking_mode,
                best_binding_energy=parsed_dlg.poses[0].binding_energy,
                poses=canonical_poses,
                clusters=cluster_infos,
                thermodynamics=thermo,
                provenance=provenance,
                summary={
                    "total_runs": parsed_dlg.num_runs,
                    "total_poses": len(canonical_poses),
                    "total_clusters": len(cluster_infos),
                }
            )
    except Exception as e:
        job_log.warning(f"Could not extract DLG poses during job completion: {e}")

    job.status = JobStatus.SUCCESS
    job.analysis_status = AnalysisStatus.SUCCESS

    # ── Build receptor+ligand complex PDBQT files ─────────────────────────────
    # Generates viewable complex files (receptor + each docked pose) in complexes/
    # Works in UCSF Chimera, VMD, Discovery Studio, PyMol manually, etc.
    if job.dlg_path and job.dlg_path.exists() and job.receptor_path and job.receptor_path.exists():
        try:
            from complex_builder import build_job_complexes
            if job.canonical_result and job.canonical_result.poses:
                n_ad4_poses = len(job.canonical_result.poses)
            elif isinstance(job.ad4_results, dict):
                n_ad4_poses = len(job.ad4_results.get("poses", []))
            elif isinstance(job.ad4_results, list):
                n_ad4_poses = len(job.ad4_results)
            else:
                n_ad4_poses = 3
            n_poses_to_build = n_ad4_poses
            complex_dir = job.output_dir / "complexes"
            built = build_job_complexes(
                receptor_path=job.receptor_path,
                ligand_output_pdbqt=job.dlg_path,
                output_dir=complex_dir,
                n_poses=n_poses_to_build,
                receptor_label=job.receptor_name,
                ligand_label=job.ligand_name,
            )
            if built:
                job_log.info(f"Complex PDBQT files written to: {complex_dir}")
        except Exception as _cx_err:
            job_log.debug(f"Complex build skipped: {_cx_err}")

    job_log.info(f"Final status: {job.status.value}")
    return job


def _run_autogrid4(
    exe: Path, gpf_path: Path, work_dir: Path,
    job: DockingJob, job_log: logging.Logger,
) -> bool:
    """Run AutoGrid4 on the GPF file. Produces .map files + .glg grid log."""
    glg_name = f"{gpf_path.stem}.glg"
    cmd = [str(exe), "-p", gpf_path.name, "-l", glg_name]
    job_log.info(f"AutoGrid4 command: {' '.join(cmd)}")

    if process_manager.is_cancelled():
        job.status = JobStatus.CANCELLED
        job.execution_status = ExecutionStatus.CANCELLED
        job_log.warning("AutoGrid4 cancelled before start")
        return False

    proc = None
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(work_dir),
        )
        process_manager.register(proc.pid, job.job_id)

        while proc.poll() is None:
            if process_manager.is_cancelled():
                job_log.warning("Cancellation requested during AutoGrid4 execution")
                process_manager.terminate_process(proc.pid)
                job.status = JobStatus.CANCELLED
                job.execution_status = ExecutionStatus.CANCELLED
                return False
            time.sleep(0.2)

        stdout, stderr = proc.communicate(timeout=5)
        if proc.returncode != 0:
            error = (stderr or stdout or "").strip()
            job_log.error(f"AutoGrid4 failed (exit {proc.returncode}): {error}")
            job.errors.append(f"AutoGrid4 failed: {error}")
            job.status = JobStatus.EXECUTION_FAILED
            job.execution_status = ExecutionStatus.FAILED
            return False
        job_log.info("AutoGrid4 completed successfully")
        return True
    except subprocess.TimeoutExpired:
        job_log.error("AutoGrid4 timed out")
        job.errors.append("AutoGrid4 timed out")
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.TIMEOUT
        return False
    except Exception as e:
        job_log.error(f"AutoGrid4 error: {e}")
        job.errors.append(str(e))
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.FAILED
        return False
    finally:
        if proc and proc.pid:
            process_manager.unregister(proc.pid)


def _run_autodock4(
    exe: Path, dpf_path: Path, dlg_path: Path,
    work_dir: Path, job: DockingJob, job_log: logging.Logger,
    central_dlg_dir: Optional[Path] = None,
    central_dpf_dir: Optional[Path] = None,
    workspace_dlg_dir: Optional[Path] = None,
) -> bool:
    """Run AutoDock4 on the DPF file with real-time live run monitoring.

    Reads stdout continuously on a background thread to prevent OS pipe buffer
    exhaustion and deadlock, while monitoring DLG file progress in real time.
    """
    cmd = [str(exe), "-p", dpf_path.name, "-l", dlg_path.name]
    job_log.info(f"AutoDock4 command: {' '.join(cmd)}")

    if process_manager.is_cancelled():
        job.status = JobStatus.CANCELLED
        job.execution_status = ExecutionStatus.CANCELLED
        job_log.warning("AutoDock4 cancelled before start")
        return False

    print(f"\n  [AutoDock4] Launching docking for {job.receptor_name} + {job.ligand_name} ...")
    print(f"  [AutoDock4] DPF location : {dpf_path}")
    print(f"  [AutoDock4] DLG output   : {dlg_path}")
    sys.stdout.flush()

    start_time = time.time()
    proc = None
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(work_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        process_manager.register(proc.pid, job.job_id)

        stdout_lines: List[str] = []
        import threading

        def _drain_stdout() -> None:
            """Continuously drain stdout pipe to prevent OS buffer exhaustion deadlock."""
            try:
                if proc and proc.stdout:
                    for line in iter(proc.stdout.readline, ""):
                        stdout_lines.append(line)
            except Exception:
                pass
            finally:
                try:
                    if proc and proc.stdout:
                        proc.stdout.close()
                except Exception:
                    pass

        drainer = threading.Thread(target=_drain_stdout, daemon=True)
        drainer.start()

        last_reported_run = 0
        import re

        # Poll process and monitor DLG file as it gets populated
        while proc.poll() is None:
            if process_manager.is_cancelled():
                job_log.warning("Cancellation requested during AutoDock4 execution")
                process_manager.terminate_process(proc.pid)
                job.status = JobStatus.CANCELLED
                job.execution_status = ExecutionStatus.CANCELLED
                print(f"  [AutoDock4] ⏹ Cancelled by user.")
                return False

            time.sleep(0.5)
            if dlg_path.is_file():
                try:
                    lines = dlg_path.read_text(encoding="utf-8", errors="replace").splitlines()
                    for idx, line in enumerate(lines):
                        m_run = re.search(r"\[\s*Run\s+(\d+)\s+of\s+(\d+)", line)
                        if m_run:
                            run_idx = int(m_run.group(1))
                            total_runs = int(m_run.group(2))
                            if run_idx > last_reported_run:
                                energy_str = ""
                                for next_line in lines[idx:min(len(lines), idx + 25)]:
                                    if "Estimated Free Energy of Binding" in next_line:
                                        m_eng = re.search(r"=\s*([+-]?\d+\.?\d*)\s*kcal/mol", next_line)
                                        if m_eng:
                                            energy_str = f"  (ΔG = {m_eng.group(1)} kcal/mol)"
                                            break
                                print(f"  [AutoDock4]  ▶ Completed Run {run_idx}/{total_runs}{energy_str}")
                                sys.stdout.flush()
                                last_reported_run = run_idx
                except Exception:
                    pass

        drainer.join(timeout=5.0)
        stdout_rest = "".join(stdout_lines)
        job.elapsed_seconds = time.time() - start_time
        job.exit_code = proc.returncode
        job.stdout = stdout_rest
        job.stderr = ""

        if proc.returncode != 0:
            error = stdout_rest.strip() if stdout_rest else f"AutoDock4 failed with exit code {proc.returncode}"
            job_log.error(f"AutoDock4 failed (exit {proc.returncode}): {error}")
            job.errors.append(f"AutoDock4 failed: {error}")
            job.status = JobStatus.EXECUTION_FAILED
            job.execution_status = ExecutionStatus.FAILED
            print(f"  [AutoDock4] ✖ Execution failed (exit code {proc.returncode})")
            return False

        job_log.info(f"AutoDock4 completed in {job.elapsed_seconds:.1f}s")
        print(f"  [AutoDock4] ✔ Docking complete in {job.elapsed_seconds:.1f}s\n")
        sys.stdout.flush()
        return True
    except subprocess.TimeoutExpired:
        job.elapsed_seconds = time.time() - start_time
        job_log.error("AutoDock4 timed out")
        job.errors.append("AutoDock4 timed out")
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.TIMEOUT
        return False
    except Exception as e:
        job.elapsed_seconds = time.time() - start_time
        job_log.error(f"AutoDock4 error: {e}")
        job.errors.append(str(e))
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.FAILED
        return False
    finally:
        if proc and proc.pid:
            process_manager.unregister(proc.pid)


# ═══════════════════════════════════════════════════════════════════════════════
# Full AD4 Workflow
# ═══════════════════════════════════════════════════════════════════════════════

def run_ad4_workflow(
    config: ProjectConfig,
    resume_mode: ResumeMode = ResumeMode.RESUME,
    progress_callback: Optional[Any] = None,
) -> List[DockingJob]:
    """Run the complete AutoDock4 docking workflow.

    After all jobs complete, delegates DLG parsing to dlg_extract.py.
    """
    # Build job queue
    all_jobs = build_job_queue(config)
    if not all_jobs:
        logger.error("No docking jobs to execute.")
        return []

    # Load previous status with physical output verification
    previous_status = load_job_status(config.result_directory)
    jobs_to_run = filter_jobs(all_jobs, resume_mode, previous_status, validate_outputs=True)
    skipped = len(all_jobs) - len(jobs_to_run)

    # Display summary
    receptors = set(j.receptor_name for j in all_jobs)
    ligands = set(j.ligand_name for j in all_jobs)

    print()
    print("=" * 60)
    print(f"  {SUITE_BANNER}")
    print(f"  Engine: AutoDock 4.2")
    print("=" * 60)
    print(f"  Receptors   : {len(receptors)}")
    print(f"  Ligands     : {len(ligands)}")
    print(f"  Total jobs  : {len(all_jobs)}")
    if skipped:
        print(f"  Skipped     : {skipped} (previously completed)")
    print(f"  To execute  : {len(jobs_to_run)}")
    print(f"  GA runs     : {config.ga_run}")
    print(f"  Output files: GPF, GLG (grid log), DPF, DLG (docking log)")
    print("=" * 60)

    # Validate executables before execution
    if not validate_executable(config.autodock4_executable, "AutoDock 4"):
        config.autodock4_executable = resolve_executable(config.autodock4_executable, "autodock4.exe")
        if not validate_executable(config.autodock4_executable, "AutoDock 4"):
            err_msg = (
                f"AutoDock 4 executable not found or invalid: '{config.autodock4_executable}'.\n"
                f"Please verify that autodock4.exe is placed in the bin/ directory or configure its path in Settings."
            )
            logger.error(err_msg)
            print(f"\n[ERROR] {err_msg}\n")
            raise RuntimeError(err_msg)

    if not validate_executable(config.autogrid4_executable, "AutoGrid 4"):
        config.autogrid4_executable = resolve_executable(config.autogrid4_executable, "autogrid4.exe")
        if not validate_executable(config.autogrid4_executable, "AutoGrid 4"):
            err_msg = (
                f"AutoGrid 4 executable not found or invalid: '{config.autogrid4_executable}'.\n"
                f"Please verify that autogrid4.exe is placed in the bin/ directory or configure its path in Settings."
            )
            logger.error(err_msg)
            print(f"\n[ERROR] {err_msg}\n")
            raise RuntimeError(err_msg)

    # Dry-run mode
    if config.dry_run:
        print("\n  DRY-RUN MODE — No docking will be executed.")
        return all_jobs

    # Execute jobs (parallel or sequential)
    if config.max_workers > 1:
        run_jobs_parallel(
            jobs_to_run=jobs_to_run,
            all_jobs=all_jobs,
            config=config,
            runner_func=run_ad4_job,
            engine_name="AutoDock 4.2",
            skipped=skipped,
            progress_callback=progress_callback,
        )
    else:
        progress = ProgressTracker(len(all_jobs), "AutoDock 4.2")
        progress.completed = skipped

        for i, job in enumerate(jobs_to_run, 1):
            if process_manager.is_cancelled():
                logger.warning("AutoDock4 workflow cancelled by user.")
                break
            progress.display(skipped + i, job)
            job = run_ad4_job(job, config)
            progress.update(job)
            save_job_status(config.result_directory, all_jobs)
            if process_manager.is_cancelled():
                logger.warning("AutoDock4 workflow cancelled after job completion.")
                break

            mode_label = f"[{job.docking_mode.value}]"
            if job.status == JobStatus.SUCCESS:
                print(
                    f"  [OK] {job.job_id} {mode_label}: DLG generated"
                    + (f" | GLG: {job.glg_path.name}" if job.glg_path else "")
                )
            else:
                print(f"  [FAIL] {job.job_id} {mode_label}: {job.status.value}")
                for e in job.errors:
                    print(f"    -> {e}")

        progress.display_summary()

    if not process_manager.is_cancelled():
        # Post-processing: Invoke dlg_extract.py for batch analysis
        _run_dlg_analysis(config, all_jobs)

    _save_ad4_run_metadata(config, all_jobs)

    if not process_manager.is_cancelled():

        # Generate combined reports
        from reporting import generate_reports
        generate_reports(all_jobs, config)

    return all_jobs


def _save_ad4_run_metadata(config: ProjectConfig, jobs: List[DockingJob]) -> Path:
    """Persist additive AD4 typing/parameter provenance without changing reports."""
    import platform
    from ad4_compatibility import PARAMETER_PROFILES, STANDARD_AD4_PROFILE, build_job_provenance

    profile = PARAMETER_PROFILES.get(
        getattr(config, "ad4_parameter_profile", STANDARD_AD4_PROFILE.profile_id),
        STANDARD_AD4_PROFILE,
    )
    records = []
    for job in jobs:
        compatibility = job.ad4_results.get("compatibility") if isinstance(job.ad4_results, dict) else None
        if not compatibility:
            input_files = [("receptor", job.receptor_path), ("ligand", job.ligand_path)]
            if job.flex_receptor_path and job.flex_receptor_path.exists():
                input_files.append(("flex_receptor", job.flex_receptor_path))
            compatibility = build_job_provenance(input_files, profile, engine="AUTODOCK4")
        records.append({
            "job_id": job.job_id,
            "receptor": job.receptor_name,
            "ligand": job.ligand_name,
            "status": job.status.value,
            "compatibility": compatibility,
        })

    metadata = {
        "suite": "AutoDock Suite Pro",
        "suite_version": __version__,
        "timestamp": datetime.now().isoformat(),
        "operating_system": platform.platform(),
        "python_version": platform.python_version(),
        "engine": "AUTODOCK4",
        "autodock4_executable": str(config.autodock4_executable),
        "autogrid4_executable": str(config.autogrid4_executable),
        "parameter_profile": profile.to_dict(),
        "parameter_file_selected_by_adsp": profile.parameter_file,
        "specialized_metal_workflow_available": False,
        "specialized_metal_workflow_reference_available": True,
        "specialized_metal_workflow_status": "SPECIALIZED_WORKFLOW_UNAVAILABLE_UNTIL_LEGACY_DEPENDENCIES_ARE_PRESENT",
        "specialized_metal_workflow_selected": False,
        "scientific_status_levels": {
            "structural_preservation": True,
            "standard_ad4_typing": True,
            "standard_ad4_parameters": bool(profile.parameter_file or profile.runtime_parameter_resolution),
            "specialized_metal_parameters": False,
            "specialized_coordination_workflow": False,
        },
        "metal_support_semantics": "Standard atom-type recognition/parameters only; no specialized metal coordination treatment.",
        "jobs": records,
    }
    metadata_path = config.result_directory / "ad4_run_metadata.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    logger.info(f"AutoDock4 compatibility metadata saved: {metadata_path}")
    return metadata_path


run_autodock4_workflow = run_ad4_workflow


def _run_dlg_analysis(config: ProjectConfig, jobs: List[DockingJob]) -> None:
    """Run dlg_extract.py (BSNDVP) on all generated DLG files.

    Discovers DLG files directly from results/AD4 without duplicating files into
    multiple mirror directories, and outputs directly to BSNDVP_RESULTS/.
    """
    completed_jobs = [
        j for j in jobs
        if j.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)
        and (
            j.dlg_path
            or (isinstance(j.ad4_results, dict) and j.ad4_results.get("dlg_path"))
        )
    ]

    ad4_dir = config.result_directory / "AD4"
    has_dlgs = ad4_dir.is_dir() and any(ad4_dir.rglob("*.dlg"))
    if not completed_jobs and not has_dlgs:
        logger.info("No completed AD4 jobs or DLG files to analyze")
        return

    print()
    print("=" * 60)
    print("  BSNDVP™ DLG ANALYSIS")
    print("=" * 60)

    # Output directly to workspace BSNDVP_RESULTS (separate from results/ which only keeps VINA and AD4)
    root_dir = config.project_root if config.project_root and config.project_root.exists() else config.result_directory.parent
    output_dir = root_dir / "BSNDVP_RESULTS"
    output_dir.mkdir(parents=True, exist_ok=True)
    reports_dir = config.report_directory

    # Import and run DLGPipeline from dlg_extract.py
    try:
        if sys.platform == "win32":
            for _s in (sys.stdout, sys.stderr):
                if hasattr(_s, "reconfigure"):
                    try:
                        _s.reconfigure(encoding="utf-8", errors="replace")
                    except Exception:
                        pass

        import importlib.util

        # Multi-strategy path resolution: works in normal Python, PyInstaller frozen exe, and side-by-side deployments
        def _locate_dlg_extract() -> Path:
            candidates = [
                # Strategy 1: alongside __file__ (normal Python / development)
                Path(__file__).parent / "dlg_extract.py",
                # Strategy 2: PyInstaller _MEIPASS bundle root
                Path(getattr(sys, "_MEIPASS", Path(__file__).parent)) / "dlg_extract.py",
                # Strategy 3: next to the frozen executable
                Path(sys.executable).parent / "dlg_extract.py",
                # Strategy 4: _internal subfolder next to the exe (PyInstaller onedir layout)
                Path(sys.executable).parent / "_internal" / "dlg_extract.py",
            ]
            for candidate in candidates:
                if candidate.exists():
                    return candidate
            raise FileNotFoundError(
                f"dlg_extract.py not found. Searched:\n"
                + "\n".join(f"  {c}" for c in candidates)
                + f"\n\nRun manually: python dlg_extract.py --dlg \"{ad4_dir}\" --output \"{output_dir}\""
            )

        dlg_extract_path = _locate_dlg_extract()
        spec = importlib.util.spec_from_file_location("dlg_extract", str(dlg_extract_path))
        dlg_extract = importlib.util.module_from_spec(spec)
        sys.modules["dlg_extract"] = dlg_extract
        spec.loader.exec_module(dlg_extract)

        # Pass ad4_dir (results/AD4) directly — DLGPipeline searches recursively for *.dlg
        pipeline = dlg_extract.DLGPipeline(ad4_dir, output_dir, root_dir)
        pipeline.run()
        print(f"\n  ✔  BSNDVP™ analysis complete. Results saved to: {output_dir}\n")

        # Backfill parsed poses and canonical result to jobs if missing
        try:
            modified = False
            for j in completed_jobs:
                if not j.canonical_result:
                    j.ensure_canonical_result()
                    if j.canonical_result:
                        modified = True
            if modified:
                from job_manager import save_job_status
                save_job_status(config.result_directory, jobs)
        except Exception as ex:
            logger.debug(f"Canonical result backfilling skipped: {ex}")

        print(f"\n  [AutoDock4] ✔ BSNDVP™ Analysis complete!")
        print(f"  [AutoDock4] Results folder : {output_dir}")
        if config.project_root:
            print(f"  [AutoDock4] Workspace      : {config.project_root}")
        print(f"  [AutoDock4] Reports folder : {reports_dir}\n")
        logger.info(f"BSNDVP analysis complete: {output_dir}")
    except Exception as e:
        logger.error(f"DLG analysis failed: {e}", exc_info=True)
        # Preserve the raw DLG/engine result, but make downstream failure
        # explicit.  A successful executable run is not a successful whole
        # workflow when DLG analysis or report generation fails.
        for job in completed_jobs:
            job.analysis_status = AnalysisStatus.FAILED
            job.status = JobStatus.POST_PROCESSING_FAILED
            job.errors.append(f"DLG analysis/report generation failed: {e}")
        try:
            from job_manager import save_job_status
            save_job_status(config.result_directory, jobs)
        except Exception as save_err:
            logger.warning(f"Could not persist post-processing failure status: {save_err}")
        print(f"  [FAIL] DLG analysis failed: {e}")
        print(f"    You can run dlg_extract.py manually:")
        print(f"    python dlg_extract.py --dlg \"{ad4_dir}\" --output \"{output_dir}\"")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI Entry Point
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="autodock4_workflow.py",
        description=f"{SUITE_BANNER} — AutoDock 4 Workflow",
    )
    parser.add_argument(
        "--config", required=True, type=Path,
        help="Path to TOML configuration file"
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--version", action="version", version=SUITE_BANNER)

    args = parser.parse_args()

    config = load_config(args.config)
    config.engine = Engine.AUTODOCK4

    if args.dry_run:
        config.dry_run = True

    run_ts = get_run_timestamp()
    setup_global_logger(config.log_directory, run_ts)

    errors = validate_config(config)
    if errors:
        print(f"\n  CONFIGURATION ERRORS:")
        for err in errors:
            print(f"    [FAIL] {err}")
        sys.exit(1)

    if args.force:
        resume_mode = ResumeMode.FORCE
    elif args.rerun_failed:
        resume_mode = ResumeMode.RERUN_FAILED
    elif args.resume or config.resume:
        resume_mode = ResumeMode.RESUME
    else:
        resume_mode = ResumeMode.FORCE

    run_ad4_workflow(config, resume_mode)


if __name__ == "__main__":
    main()
