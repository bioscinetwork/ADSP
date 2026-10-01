#!/usr/bin/env python3
"""
AutoDock Suite Pro — Job Manager
==================================
Builds job queues, tracks status, handles resume/rerun, and displays progress.

Receptor Folder Convention (Option A):
    receptors/
        RECEPTOR_NAME/
            rigid/                      ← Rigid receptor PDBQT(s)
                receptor.pdbqt
            flex/                       ← Flexible receptor PDBQTs (optional)
                receptor_rigid.pdbqt   ← Rigid backbone part  (--receptor)
                receptor_flex.pdbqt    ← Flexible sidechain part (--flex)
            config.txt                  ← Vina grid config

    Mode is auto-detected per receptor:
        - rigid/ only          → RIGID docking
        - flex/ present        → FLEXIBLE docking
        - Global [docking] mode in config.toml overrides per-receptor detection
          only if set explicitly; per-receptor detection takes precedence.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from config import ProjectConfig
from models import (
    DockingJob, DockingMode, Engine, ExecutionStatus, JobStatus, ResumeMode,
    VinaResult, generate_job_id,
)
from validators import validate_pdbqt, validate_vina_config
from vina_parser import extract_affinities
from process_manager import compute_job_input_hashes, process_manager, validate_job_outputs

logger = logging.getLogger("docking_automation.job_manager")


# ═══════════════════════════════════════════════════════════════════════════════
# Receptor Discovery
# ═══════════════════════════════════════════════════════════════════════════════

def discover_receptors(config: ProjectConfig) -> List[Tuple[str, Path, DockingMode]]:
    """Discover receptor subfolders in the receptor directory.

    Expected structure:
        receptors/
        ├── receptor_01/
        │   ├── rigid/receptor_01.pdbqt        (rigid docking)
        │   ├── flex/receptor_01_rigid.pdbqt   (flexible docking: rigid backbone)
        │   ├── flex/receptor_01_flex.pdbqt    (flexible docking: flexible sidechains)
        │   └── config.txt
        └── receptor_02/
            ├── rigid/receptor_02.pdbqt
            └── config.txt

    Mode behavior:
        - If BOTH rigid/ and flex/ exist (and mode is AUTO or BOTH):
          The receptor is queued TWICE — once for RIGID docking and once for FLEXIBLE docking.
        - If only rigid/ exists: queued for RIGID docking.
        - If only flex/ exists: queued for FLEXIBLE docking.
        - If mode is set explicitly to RIGID or FLEXIBLE in config or CLI:
          Only that specific mode is queued.

    Returns:
        List of (receptor_name, receptor_dir_path, docking_mode).
    """
    receptors: List[Tuple[str, Path, DockingMode]] = []
    rec_dir = config.receptor_directory

    if not rec_dir.exists():
        logger.error(f"Receptor directory does not exist: {rec_dir}")
        return []

    for item in sorted(rec_dir.iterdir()):
        if not item.is_dir():
            continue

        name = item.name

        # Filter by selected receptors if specified
        if config.selected_receptors and name not in config.selected_receptors:
            continue

        # Check for rigid setup
        rigid_pdbqt = _find_rigid_pdbqt(item, name)
        has_rigid = rigid_pdbqt is not None

        # Check for flexible setup
        flex_dir = item / "flex"
        has_flex = False
        if flex_dir.is_dir():
            rigid_in_flex = _find_flex_rigid_pdbqt(flex_dir, name)
            flex_in_flex = _find_flex_sidechain_pdbqt(flex_dir, name)
            if rigid_in_flex is not None or flex_in_flex is not None:
                has_flex = True

        # Evaluate based on config.docking_mode (supported for both Vina and AutoDock4)
        target_mode = config.docking_mode

        queue_rigid = False
        queue_flex = False

        if target_mode == DockingMode.RIGID:
            if has_rigid:
                queue_rigid = True
            else:
                logger.warning(
                    f"Skipping receptor '{name}': mode=RIGID requested but no rigid PDBQT found "
                    f"in {item / 'rigid'} or {item}."
                )
        elif target_mode == DockingMode.FLEXIBLE:
            if has_flex:
                queue_flex = True
            else:
                logger.warning(
                    f"Skipping receptor '{name}': mode=FLEXIBLE requested but no valid flex/ files found "
                    f"in {item / 'flex'}."
                )
        else:
            # AUTO or BOTH: auto-detect from folders
            if has_rigid:
                queue_rigid = True
            if has_flex:
                queue_flex = True

            if not has_rigid and not has_flex:
                logger.warning(
                    f"Skipping receptor '{name}': no PDBQT found in "
                    f"{item / 'rigid'} or {item}, and no flex/ folder found."
                )
                continue

        if queue_rigid and queue_flex:
            logger.info(
                f"Receptor '{name}': discovered both rigid/ and flex/ setups — "
                f"will queue BOTH rigid and flexible docking."
            )
            receptors.append((name, item, DockingMode.RIGID))
            receptors.append((name, item, DockingMode.FLEXIBLE))
        elif queue_rigid:
            receptors.append((name, item, DockingMode.RIGID))
        elif queue_flex:
            receptors.append((name, item, DockingMode.FLEXIBLE))

    if not receptors:
        logger.warning(
            f"No receptors found in {rec_dir}. "
            f"Each receptor should be a subfolder containing either:\n"
            f"  - rigid/ with a .pdbqt file (rigid docking), OR\n"
            f"  - flex/ with *rigid*.pdbqt and *flex*.pdbqt files (flexible docking)."
        )

    return receptors



# ═══════════════════════════════════════════════════════════════════════════════
# File Finders
# ═══════════════════════════════════════════════════════════════════════════════

def _find_rigid_pdbqt(receptor_dir: Path, receptor_name: str) -> Optional[Path]:
    """Find the rigid receptor PDBQT.

    Search order:
        1. receptors/RECEPTOR/rigid/*.pdbqt   (preferred)
        2. receptors/RECEPTOR/*.pdbqt          (legacy fallback — not prefixed flex_/rigid_)
    """
    rigid_dir = receptor_dir / "rigid"
    if rigid_dir.is_dir():
        # Prefer exact name match first
        for candidate in [
            rigid_dir / f"{receptor_name}.pdbqt",
            rigid_dir / "receptor.pdbqt",
        ]:
            if candidate.exists():
                return candidate
        # Any PDBQT in the rigid/ folder
        pdbqts = list(rigid_dir.glob("*.pdbqt"))
        if pdbqts:
            return sorted(pdbqts)[0]

    # Legacy fallback: PDBQT files directly in the receptor folder
    for p in sorted(receptor_dir.glob("*.pdbqt")):
        stem = p.stem.lower()
        if not any(stem.startswith(pfx) for pfx in ("flex_", "rigid_", "flexible_")):
            return p

    return None


def _find_flex_rigid_pdbqt(flex_dir: Path, receptor_name: str) -> Optional[Path]:
    """Find the rigid backbone PDBQT inside the flex/ subfolder.

    This is the file passed to Vina's --receptor in flexible mode.
    Naming conventions accepted (in priority order):
        *rigid*.pdbqt
        *backbone*.pdbqt
        first *.pdbqt that does NOT contain 'flex' in its stem
    """
    if not flex_dir.is_dir():
        return None

    candidates = list(flex_dir.glob("*.pdbqt"))

    # Prefer files with 'rigid' or 'backbone' in the name
    for p in sorted(candidates):
        stem = p.stem.lower()
        if "rigid" in stem or "backbone" in stem:
            return p

    # Fallback: any file that does not contain 'flex' in the stem
    for p in sorted(candidates):
        if "flex" not in p.stem.lower():
            return p

    return None


def _find_flex_sidechain_pdbqt(flex_dir: Path, receptor_name: str) -> Optional[Path]:
    """Find the flexible sidechain PDBQT inside the flex/ subfolder.

    This is the file passed to Vina's --flex argument.
    Naming conventions accepted (in priority order):
        *flex*.pdbqt
        *sidechain*.pdbqt
        *flexible*.pdbqt
        Any file that does not match the rigid backbone finder
    """
    if not flex_dir.is_dir():
        return None

    candidates = list(flex_dir.glob("*.pdbqt"))

    # Prefer files with 'flex', 'sidechain', or 'flexible' in the name
    for p in sorted(candidates):
        stem = p.stem.lower()
        if "flex" in stem or "sidechain" in stem or "flexible" in stem:
            return p

    return None


def find_receptor_config(receptor_dir: Path, receptor_name: str = "") -> Optional[Path]:
    """Find the Vina config.txt in a receptor directory or parent."""
    candidates = [
        "config.txt", "vina_config.txt", "conf.txt", "grid.txt",
    ]
    if receptor_name:
        candidates.extend([
            f"{receptor_name}.txt",
            f"{receptor_name}_config.txt",
            f"{receptor_name}config.txt",
            f"{receptor_name}grid.txt",
        ])
    for name in candidates:
        p = receptor_dir / name
        if p.exists():
            return p
    # Also check parent folder
    if (receptor_dir.parent / "config.txt").exists():
        return receptor_dir.parent / "config.txt"
    return None


# ═══════════════════════════════════════════════════════════════════════════════
# Ligand Discovery
# ═══════════════════════════════════════════════════════════════════════════════

def discover_ligands(config: ProjectConfig) -> List[Tuple[str, Path]]:
    """Discover ligand PDBQT files in the ligand directory.

    Returns:
        List of (ligand_name, ligand_file_path).
    """
    ligands = []
    lig_dir = config.ligand_directory

    if not lig_dir.exists():
        logger.error(f"Ligand directory does not exist: {lig_dir}")
        return []

    for pdbqt in sorted(lig_dir.glob("*.pdbqt")):
        name = pdbqt.stem
        if config.selected_ligands and name not in config.selected_ligands:
            continue
        ligands.append((name, pdbqt))

    return ligands


# ═══════════════════════════════════════════════════════════════════════════════
# Job Queue Construction
# ═══════════════════════════════════════════════════════════════════════════════

def build_job_queue(config: ProjectConfig) -> List[DockingJob]:
    """Build the N×M Cartesian product of receptors × ligands.

    Each receptor is auto-detected as RIGID or FLEXIBLE based on its folder
    structure (presence of flex/ subfolder with PDBQT files).

    Returns:
        List of DockingJob objects ready for execution.
    """
    receptors = discover_receptors(config)
    ligands = discover_ligands(config)

    if not receptors:
        logger.error("No receptors found!")
        return []
    if not ligands:
        logger.error("No ligands found!")
        return []

    jobs: List[DockingJob] = []

    for rec_name, rec_dir, mode in receptors:
        flex_dir = rec_dir / "flex"

        # Resolve receptor paths based on detected mode
        if mode == DockingMode.FLEXIBLE:
            # --receptor = rigid backbone from flex/
            rigid_pdbqt = _find_flex_rigid_pdbqt(flex_dir, rec_name)
            # --flex = flexible sidechains from flex/
            flex_pdbqt = _find_flex_sidechain_pdbqt(flex_dir, rec_name)

            if rigid_pdbqt is None:
                logger.warning(
                    f"Receptor '{rec_name}': flex/ folder found but no rigid "
                    f"backbone PDBQT (*rigid*.pdbqt or *backbone*.pdbqt). "
                    f"Falling back to RIGID mode."
                )
                rigid_pdbqt = _find_rigid_pdbqt(rec_dir, rec_name)
                flex_pdbqt = None
                mode = DockingMode.RIGID

            if flex_pdbqt is None and mode == DockingMode.FLEXIBLE:
                logger.warning(
                    f"Receptor '{rec_name}': no flexible sidechain PDBQT found "
                    f"in flex/ (*flex*.pdbqt or *sidechain*.pdbqt). "
                    f"Falling back to RIGID mode."
                )
                rigid_pdbqt = _find_rigid_pdbqt(rec_dir, rec_name)
                flex_pdbqt = None
                mode = DockingMode.RIGID
        else:
            # RIGID mode: use rigid/ subfolder
            rigid_pdbqt = _find_rigid_pdbqt(rec_dir, rec_name)
            flex_pdbqt = None

        # Find Vina config
        config_txt = find_receptor_config(rec_dir, rec_name)

        mode_label = "FLEX" if mode == DockingMode.FLEXIBLE else "RIGID"
        logger.info(
            f"Receptor '{rec_name}': mode={mode_label}, "
            f"rigid={rigid_pdbqt}, flex={flex_pdbqt}, config={config_txt}"
        )

        for lig_name, lig_path in ligands:
            job_id = generate_job_id(rec_name, lig_name, docking_mode=mode,
                                     engine=config.engine)

            # Determine output directory
            engine_label = "VINA" if config.engine == Engine.VINA else "AD4"
            if mode == DockingMode.FLEXIBLE:
                rec_folder = rec_name if rec_name.lower().endswith("_flex") else f"{rec_name}_flex"
            else:
                rec_folder = rec_name
            output_dir = config.result_directory / engine_label / rec_folder / lig_name

            job = DockingJob(
                job_id=job_id,
                receptor_name=rec_name,
                ligand_name=lig_name,
                receptor_path=rigid_pdbqt,
                ligand_path=lig_path,
                config_path=config_txt,
                flex_receptor_path=flex_pdbqt,
                output_dir=output_dir,
                engine=config.engine,
                docking_mode=mode,
                requested_modes=(
                    config.num_modes if config.engine == Engine.VINA
                    else config.ad4_num_modes
                ),
            )
            jobs.append(job)

    return jobs


# ═══════════════════════════════════════════════════════════════════════════════
# Validate Job Inputs
# ═══════════════════════════════════════════════════════════════════════════════

def validate_job_inputs(job: DockingJob) -> List[str]:
    """Validate all input files for a job before execution.

    Returns list of error messages (empty if valid).
    """
    errors: List[str] = []

    # Receptor (rigid backbone)
    if not job.receptor_path:
        errors.append(f"[{job.job_id}] No receptor PDBQT found for {job.receptor_name}")
    elif not job.receptor_path.exists():
        errors.append(f"[{job.job_id}] Receptor PDBQT not found: {job.receptor_path}")
    else:
        valid, msgs = validate_pdbqt(job.receptor_path)
        if not valid:
            errors.extend(f"[{job.job_id}] Receptor: {m}" for m in msgs)

    # Ligand
    if not job.ligand_path:
        errors.append(f"[{job.job_id}] No ligand PDBQT found for {job.ligand_name}")
    elif not job.ligand_path.exists():
        errors.append(f"[{job.job_id}] Ligand PDBQT not found: {job.ligand_path}")
    else:
        valid, msgs = validate_pdbqt(job.ligand_path)
        if not valid:
            errors.extend(f"[{job.job_id}] Ligand: {m}" for m in msgs)

    # Vina config
    if job.engine == Engine.VINA:
        if not job.config_path:
            errors.append(
                f"[{job.job_id}] No config.txt found for receptor '{job.receptor_name}'. "
                f"Place a config.txt in receptors/{job.receptor_name}/ or run "
                f"'python main.py --prepare' to auto-generate one."
            )
        else:
            valid, msgs = validate_vina_config(job.config_path)
            if not valid:
                errors.extend(f"[{job.job_id}] Config: {m}" for m in msgs)

    # Flexible sidechain PDBQT
    if job.docking_mode == DockingMode.FLEXIBLE:
        if not job.flex_receptor_path:
            errors.append(
                f"[{job.job_id}] Flexible mode but no sidechain PDBQT found. "
                f"Place *flex*.pdbqt in receptors/{job.receptor_name}/flex/"
            )
        elif not job.flex_receptor_path.exists():
            errors.append(
                f"[{job.job_id}] Flexible sidechain PDBQT not found: {job.flex_receptor_path}"
            )
        else:
            valid, msgs = validate_pdbqt(job.flex_receptor_path)
            if not valid:
                errors.extend(f"[{job.job_id}] Flex receptor: {m}" for m in msgs)

    return errors


# ═══════════════════════════════════════════════════════════════════════════════
# Job Status Persistence (Resume Support)
# ═══════════════════════════════════════════════════════════════════════════════

STATUS_FILENAME = "job_status.json"


def load_job_status(status_dir: Path) -> Dict[str, Dict]:
    """Load job status from the JSON status file.

    Returns dict mapping job_id -> job data dict.
    """
    status_file = status_dir / STATUS_FILENAME
    if not status_file.exists():
        return {}

    try:
        data = json.loads(status_file.read_text(encoding="utf-8"))
        return {job["job_id"]: job for job in data.get("jobs", [])}
    except Exception as e:
        logger.warning(f"Could not load job status: {e}")
        return {}


def load_docking_jobs(status_dir: Path) -> List[DockingJob]:
    """Load and deserialize all DockingJob objects from the JSON status file."""
    status_file = status_dir / STATUS_FILENAME
    if not status_file.exists():
        return []

    try:
        data = json.loads(status_file.read_text(encoding="utf-8"))
        return [DockingJob.from_dict(j) for j in data.get("jobs", [])]
    except Exception as e:
        logger.warning(f"Could not load docking jobs: {e}")
        return []


def save_job_status(status_dir: Path, jobs: List[DockingJob]) -> None:
    """Save current job status to the JSON status file, merging with existing jobs by job_id."""
    status_dir.mkdir(parents=True, exist_ok=True)
    status_file = status_dir / STATUS_FILENAME

    merged_jobs: Dict[str, Dict] = {}
    if status_file.exists():
        try:
            prev_data = json.loads(status_file.read_text(encoding="utf-8"))
            for j_dict in prev_data.get("jobs", []):
                jid = j_dict.get("job_id")
                if jid:
                    merged_jobs[jid] = j_dict
        except Exception:
            pass

    for job in jobs:
        j_dict = job.to_dict()
        j_dict["input_hashes"] = compute_job_input_hashes(job)
        merged_jobs[job.job_id] = j_dict

    all_list = list(merged_jobs.values())
    data = {
        "schema_version": "2.0.0",
        "jobs": all_list,
        "total": len(all_list),
        "completed": sum(1 for j in all_list if j.get("status") in (
            JobStatus.SUCCESS.value, JobStatus.SUCCESS_WITH_WARNING.value,
            "SUCCESS", "SUCCESS_WITH_WARNING"
        )),
        "failed": sum(1 for j in all_list if j.get("status") in (
            JobStatus.EXECUTION_FAILED.value, JobStatus.PARSING_FAILED.value,
            JobStatus.SPLIT_FAILED.value, JobStatus.INVALID_INPUT.value,
            "EXECUTION_FAILED", "PARSING_FAILED", "SPLIT_FAILED", "INVALID_INPUT"
        )),
        "cancelled": sum(1 for j in all_list if j.get("status") in (
            JobStatus.CANCELLED.value, "CANCELLED"
        )),
    }

    status_file.write_text(
        json.dumps(data, indent=2, default=str),
        encoding="utf-8"
    )


def filter_jobs(
    jobs: List[DockingJob],
    mode: ResumeMode,
    previous_status: Dict[str, Dict],
    validate_outputs: bool = False,
) -> List[DockingJob]:
    """Filter jobs based on resume mode, output file validity, and input hashes.

    When validate_outputs is True:
      A job is ONLY skipped when:
        1. Its previous status is SUCCESS or SUCCESS_WITH_WARNING.
        2. Its physical output file(s) exist and are validated.
        3. Its input file hashes (receptor, ligand, config) match the previous run.

    Returns:
        Filtered list of jobs to execute.
    """
    if mode == ResumeMode.FORCE or not previous_status:
        return jobs

    filtered: List[DockingJob] = []
    for job in jobs:
        prev = previous_status.get(job.job_id)
        if prev is None:
            filtered.append(job)
            continue

        prev_status = prev.get("status", "PENDING")

        if mode in (ResumeMode.RESUME, ResumeMode.RERUN_FAILED):
            if prev_status in ("SUCCESS", "SUCCESS_WITH_WARNING"):
                test_job = DockingJob.from_dict(prev)
                if validate_outputs:
                    out_valid, reason = validate_job_outputs(test_job)
                    if not out_valid:
                        logger.warning(
                            f"Rerunning job '{job.job_id}': Previous checkpoint was SUCCESS, "
                            f"but underlying output is invalid/missing ({reason})."
                        )
                        filtered.append(job)
                        continue

                # Verify input hashes if previous record stored them
                prev_hashes = prev.get("input_hashes")
                if prev_hashes and isinstance(prev_hashes, dict):
                    curr_hashes = compute_job_input_hashes(job)
                    hash_mismatch = False
                    for key in ("receptor_hash", "ligand_hash", "config_hash", "flex_receptor_hash"):
                        if key in curr_hashes and key in prev_hashes:
                            if curr_hashes[key] != prev_hashes[key]:
                                logger.info(
                                    f"Rerunning job '{job.job_id}': Input {key} changed since previous run "
                                    f"(hash mismatch: {prev_hashes[key][:8]}... -> {curr_hashes[key][:8]}...)."
                                )
                                hash_mismatch = True
                                break
                    if hash_mismatch:
                        filtered.append(job)
                        continue

                # Output valid and hashes match -> safely skip
                job.status = JobStatus(prev_status)
                job.execution_status = ExecutionStatus.SUCCESS
                job.obtained_modes = prev.get("obtained_modes", 0)
                job.warnings = prev.get("warnings", [])
                job.elapsed_seconds = prev.get("elapsed_seconds", 0.0)
                job.output_pdbqt = Path(prev["output_pdbqt"]) if prev.get("output_pdbqt") else None
                job.output_log = Path(prev["output_log"]) if prev.get("output_log") else None
                job.gpf_path = Path(prev["gpf_path"]) if prev.get("gpf_path") else None
                job.glg_path = Path(prev["glg_path"]) if prev.get("glg_path") else None
                job.dpf_path = Path(prev["dpf_path"]) if prev.get("dpf_path") else None
                job.dlg_path = Path(prev["dlg_path"]) if prev.get("dlg_path") else None
                job.ad4_results = prev.get("ad4_results", {})
                job.canonical_result = test_job.canonical_result
                if prev.get("vina_results"):
                    job.vina_results = [
                        VinaResult(**r) if isinstance(r, dict) else r
                        for r in prev["vina_results"]
                    ]
                elif job.output_pdbqt and job.output_pdbqt.exists():
                    job.vina_results = extract_affinities(job.output_pdbqt)

                logger.info(
                    f"SKIPPED {job.job_id}: Reason: Existing completed result matches current "
                    f"receptor, ligand, engine and configuration hashes."
                )
            else:
                filtered.append(job)
        else:
            filtered.append(job)

    return filtered


# ═══════════════════════════════════════════════════════════════════════════════
# Progress Display
# ═══════════════════════════════════════════════════════════════════════════════

class ProgressTracker:
    """Tracks and displays docking progress."""

    def __init__(self, total_jobs: int, engine: str):
        self.total = total_jobs
        self.engine = engine
        self.completed = 0
        self.warnings = 0
        self.failed = 0
        self.start_time = time.time()

    def update(self, job: DockingJob) -> None:
        """Update counters after a job completes."""
        if job.status == JobStatus.SUCCESS:
            self.completed += 1
        elif job.status == JobStatus.SUCCESS_WITH_WARNING:
            self.completed += 1
            self.warnings += 1
        elif job.status in (JobStatus.EXECUTION_FAILED, JobStatus.PARSING_FAILED,
                            JobStatus.SPLIT_FAILED, JobStatus.INVALID_INPUT):
            self.failed += 1
        elif job.status == JobStatus.SKIPPED:
            pass

    def display(self, current_index: int, job: DockingJob) -> None:
        """Display progress information for the current job."""
        elapsed = time.time() - self.start_time
        remaining = self.total - current_index
        done = current_index - 1

        if done > 0:
            avg_time = elapsed / done
            eta_seconds = avg_time * remaining
            eta_str = _format_duration(eta_seconds)
        else:
            eta_str = "calculating..."

        mode_str = f" [{job.docking_mode.value}]" if hasattr(job, "docking_mode") else ""

        print()
        print("=" * 60)
        print(f"  AUTODOCK SUITE PRO — {self.engine}")
        print("=" * 60)
        print(f"  Current job : {current_index} / {self.total}")
        print(f"  Receptor    : {job.receptor_name}{mode_str}")
        print(f"  Ligand      : {job.ligand_name}")
        print(f"  Job ID      : {job.job_id}")
        print(f"  Completed   : {self.completed}")
        print(f"  Warnings    : {self.warnings}")
        print(f"  Failed      : {self.failed}")
        print(f"  Remaining   : {remaining}")
        print(f"  Elapsed     : {_format_duration(elapsed)}")
        print(f"  ETA (est.)  : {eta_str}")
        print("=" * 60)

    def display_summary(self) -> None:
        """Display final run summary."""
        elapsed = time.time() - self.start_time
        total_done = self.completed + self.failed

        print()
        print("=" * 60)
        print("  RUN COMPLETE")
        print("=" * 60)
        print(f"  Total jobs      : {self.total}")
        print(f"  Successful      : {self.completed}")
        print(f"  Warnings        : {self.warnings}")
        print(f"  Failed          : {self.failed}")
        print(f"  Skipped         : {self.total - total_done}")
        print(f"  Total time      : {_format_duration(elapsed)}")
        print("=" * 60)


def _format_duration(seconds: float) -> str:
    """Format seconds into HH:MM:SS."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def run_jobs_parallel(
    jobs_to_run: List[DockingJob],
    all_jobs: List[DockingJob],
    config: ProjectConfig,
    runner_func: Any,
    engine_name: str = "AutoDock",
    skipped: int = 0,
    progress_callback: Optional[Any] = None,
) -> List[DockingJob]:
    """Execute docking jobs concurrently using a thread pool.

    Subprocesses (vina.exe, autodock4.exe) execute on multiple CPU cores in parallel.
    Console reporting and job_status.json saving are synchronized via locks.
    """
    import concurrent.futures
    import threading

    lock = threading.Lock()
    progress = ProgressTracker(len(all_jobs), engine_name)
    progress.completed = skipped

    max_workers = max(1, config.max_workers)
    print()
    print("=" * 60)
    print(f"  PARALLEL DOCKING EXECUTION ({max_workers} WORKERS)")
    print("=" * 60)

    completed_counter = 0

    def _worker(job_item: DockingJob) -> DockingJob:
        nonlocal completed_counter
        if process_manager.is_cancelled():
            job_item.status = JobStatus.CANCELLED
            job_item.execution_status = ExecutionStatus.CANCELLED
            return job_item

        finished_job = runner_func(job_item, config)

        with lock:
            completed_counter += 1
            progress.update(finished_job)
            save_job_status(config.result_directory, all_jobs)

            current_idx = skipped + completed_counter
            if progress_callback:
                try:
                    progress_callback(current_idx, len(all_jobs), f"{finished_job.receptor_name} × {finished_job.ligand_name}")
                except Exception:
                    pass
            mode_label = f"[{finished_job.docking_mode.value}]"
            if finished_job.status == JobStatus.SUCCESS:
                affinity_str = (
                    f"best: {finished_job.vina_results[0].binding_affinity:.1f} kcal/mol"
                    if finished_job.vina_results else ""
                )
                print(f"  [{current_idx}/{len(all_jobs)}] [OK] {finished_job.job_id} {mode_label}: {finished_job.obtained_modes} poses  {affinity_str}")
            elif finished_job.status == JobStatus.SUCCESS_WITH_WARNING:
                print(f"  [{current_idx}/{len(all_jobs)}] [WARN] {finished_job.job_id} {mode_label}: {finished_job.obtained_modes}/{finished_job.requested_modes} poses")
            else:
                print(f"  [{current_idx}/{len(all_jobs)}] [FAIL] {finished_job.job_id} {mode_label}: {finished_job.status.value}")

        return finished_job

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_worker, j) for j in jobs_to_run]
        for f in concurrent.futures.as_completed(futures):
            try:
                f.result()
            except Exception as e:
                logger.error(f"Worker exception: {e}")

    progress.display_summary()
    return all_jobs

