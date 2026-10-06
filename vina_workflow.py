#!/usr/bin/env python3
"""
AutoDock Suite Pro — AutoDock Vina Workflow
============================================
Full orchestration of Vina docking: validation, execution, parsing,
splitting, retry, resume, and reporting.

Folder Convention:
    receptors/RECEPTOR/rigid/receptor.pdbqt         ← rigid docking
    receptors/RECEPTOR/flex/receptor_rigid.pdbqt    ← backbone (--receptor)
    receptors/RECEPTOR/flex/receptor_flex.pdbqt     ← sidechains (--flex)
    receptors/RECEPTOR/config.txt                   ← grid config

Usage:
    python vina_workflow.py --config project_config.toml
    python vina_workflow.py --config project_config.toml --dry-run
    python vina_workflow.py --config project_config.toml --resume
    python vina_workflow.py --config project_config.toml --rerun-failed
    python vina_workflow.py --config project_config.toml --force
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import subprocess
import sys

# Configure UTF-8 output on Windows consoles
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from config import ProjectConfig, load_config, validate_config
from executables import (
    get_cpu_count, get_vina_version, get_vina_split_version,
    resolve_executable, validate_executable,
)
from job_manager import (
    ProgressTracker, build_job_queue, filter_jobs,
    load_job_status, run_jobs_parallel, save_job_status, validate_job_inputs,
)
from logging_utils import get_run_timestamp, setup_global_logger, setup_job_logger
from models import (
    SUITE_BANNER, AnalysisStatus, AutoDock4Metrics, CanonicalPose, DockingJob,
    DockingMode, DockingResult, Engine, ExecutionStatus, JobStatus,
    ProvenanceRecord, ResumeMode, VinaMetrics, VinaResult,
    __version__,
)
from process_manager import compute_file_hash, process_manager
from validators import parse_vina_config, validate_vina_output
from vina_parser import count_models, extract_affinities, format_vina_results_table
from vina_splitter import run_vina_split, validate_split

logger = logging.getLogger("docking_automation")


# ═══════════════════════════════════════════════════════════════════════════════
# Single Job Execution
# ═══════════════════════════════════════════════════════════════════════════════

def run_vina_job(job: DockingJob, config: ProjectConfig) -> DockingJob:
    """Execute a single Vina docking job.

    Steps:
        1.  Validate inputs (receptor, ligand, config)
        2.  Determine rigid vs. flexible mode from job.docking_mode
        3.  Construct output paths and directories
        4.  Copy input files for reproducibility
        5.  Build Vina command
            - RIGID:    --receptor <rigid.pdbqt>
            - FLEXIBLE: --receptor <flex/rigid.pdbqt> --flex <flex/flex.pdbqt>
        6.  Launch vina.exe
        7.  Capture stdout/stderr and write permanent log
        8.  Check exit code
        9.  Validate output PDBQT
        10. Parse output (count models, extract affinities)
        11. Compare requested vs. obtained poses
        12. Run vina_split
        13. Validate split output
        14. Record final status

    Returns:
        Updated DockingJob with results and status.
    """
    job_log = setup_job_logger(config.log_directory, "VINA", job.job_id)
    job.status = JobStatus.RUNNING
    job.execution_status = ExecutionStatus.RUNNING

    # ── Step 1: Validate inputs ───────────────────────────────────────────────
    input_errors = validate_job_inputs(job)
    if input_errors:
        for err in input_errors:
            job_log.error(err)
        job.errors.extend(input_errors)
        job.status = JobStatus.INVALID_INPUT
        job.execution_status = ExecutionStatus.FAILED
        return job

    # ── Step 2: Determine mode ────────────────────────────────────────────────
    flex_mode = (
        job.docking_mode == DockingMode.FLEXIBLE
        and job.flex_receptor_path is not None
        and job.flex_receptor_path.exists()
    )
    mode_str = "FLEXIBLE" if flex_mode else "RIGID"
    job_log.info(f"Docking mode: {mode_str}")

    if job.docking_mode == DockingMode.FLEXIBLE and not flex_mode:
        job_log.warning(
            "Job was marked FLEXIBLE but flex sidechain PDBQT is missing or invalid. "
            "Falling back to RIGID mode."
        )

    # ── Step 3: Construct output paths ────────────────────────────────────────
    output_dir = job.output_dir or (config.result_directory / "VINA" / f"{job.receptor_name}_x_{job.ligand_name}")
    job.output_dir = output_dir
    input_dir = output_dir / "input"
    out_dir = output_dir / "output"
    poses_dir = output_dir / "poses"

    for d in [input_dir, out_dir, poses_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ── Step 4: Copy input files for reproducibility ──────────────────────────
    if job.receptor_path and job.receptor_path.exists():
        shutil.copy2(job.receptor_path, input_dir / job.receptor_path.name)
    if job.ligand_path and job.ligand_path.exists():
        shutil.copy2(job.ligand_path, input_dir / job.ligand_path.name)
    if job.config_path and job.config_path.exists():
        shutil.copy2(job.config_path, input_dir / job.config_path.name)
    if flex_mode and job.flex_receptor_path.exists():
        shutil.copy2(job.flex_receptor_path, input_dir / job.flex_receptor_path.name)

    output_pdbqt = out_dir / f"{job.ligand_name}_out.pdbqt"
    output_log = out_dir / f"{job.ligand_name}_log.txt"

    job.output_pdbqt = output_pdbqt
    job.output_log = output_log

    # ── Step 5: Build Vina command ────────────────────────────────────────────
    cmd = [str(config.vina_executable)]

    if flex_mode:
        # Flexible docking:
        #   --receptor = rigid backbone PDBQT
        #   --flex     = flexible sidechain PDBQT
        if job.flex_receptor_path:
            try:
                from prepare import sanitize_flex_pdbqt
                sanitize_flex_pdbqt(job.flex_receptor_path)
            except Exception:
                pass
        cmd.extend(["--receptor", str(job.receptor_path)])
        cmd.extend(["--flex", str(job.flex_receptor_path)])
        job_log.info(f"Flex receptor (rigid backbone): {job.receptor_path}")
        job_log.info(f"Flex receptor (sidechains):     {job.flex_receptor_path}")
    else:
        # Rigid docking:
        cmd.extend(["--receptor", str(job.receptor_path)])

    cmd.extend(["--ligand", str(job.ligand_path)])

    if job.config_path:
        cmd.extend(["--config", str(job.config_path)])

    cmd.extend(["--out", str(output_pdbqt)])

    # Override config values with global settings if not already in config file
    config_params: Dict = {}
    if job.config_path and job.config_path.exists():
        config_params = parse_vina_config(job.config_path)

    if "exhaustiveness" not in config_params:
        cmd.extend(["--exhaustiveness", str(config.exhaustiveness)])
    if "num_modes" not in config_params:
        cmd.extend(["--num_modes", str(config.num_modes)])
    if "energy_range" not in config_params:
        cmd.extend(["--energy_range", str(config.energy_range)])

    cpu_count = get_cpu_count(config.cpu)
    if "cpu" not in config_params:
        cmd.extend(["--cpu", str(cpu_count)])

    if config.seed.upper() != "AUTO" and "seed" not in config_params:
        cmd.extend(["--seed", config.seed])

    if getattr(config, "vina_scoring", "vina") in ("vinardo", "ad4") and "scoring" not in config_params:
        cmd.extend(["--scoring", config.vina_scoring])

    if getattr(config, "vina_min_rmsd", 1.0) != 1.0 and "min_rmsd" not in config_params:
        cmd.extend(["--min_rmsd", str(config.vina_min_rmsd)])

    job.command_executed = " ".join(cmd)
    job_log.info(f"Command: {job.command_executed}")

    print(f"\n  [Vina] Launching docking for {job.receptor_name} + {job.ligand_name} ...")
    print(f"  [Vina] Output PDBQT: {output_pdbqt}")
    sys.stdout.flush()

    # ── Steps 6-7: Launch Vina ────────────────────────────────────────────────
    start_time = time.time()
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        process_manager.register(job.job_id, proc)
        stdout_chunks = []
        try:
            for line in proc.stdout:
                if process_manager.is_cancelled():
                    break
                stdout_chunks.append(line)
                clean_line = line.strip()
                if clean_line:
                    if any(k in clean_line for k in ("Performing search", "%", "mode |", "-----", "affinity", "Refining", "Writing output")):
                        print(f"  [Vina] {clean_line}")
                        sys.stdout.flush()
            proc.wait(timeout=7200)
        finally:
            process_manager.unregister(job.job_id)

        if process_manager.is_cancelled():
            job.status = JobStatus.CANCELLED
            job.execution_status = ExecutionStatus.CANCELLED
            job.warnings.append("Vina execution cancelled by user")
            job_log.warning("Vina execution cancelled by user")
            return job

        job.exit_code = proc.returncode
        full_stdout = "".join(stdout_chunks)
        job.stdout = full_stdout
        job.stderr = ""
        job.elapsed_seconds = time.time() - start_time

        # Save stdout as permanent log
        output_log.write_text(full_stdout, encoding="utf-8")

        job_log.info(f"Exit code: {proc.returncode}")
        job_log.info(f"Elapsed: {job.elapsed_seconds:.1f}s")

        if proc.returncode == 0:
            print(f"  [Vina] ✔ Docking complete in {job.elapsed_seconds:.1f}s")
            sys.stdout.flush()

    except subprocess.TimeoutExpired:
        job.elapsed_seconds = time.time() - start_time
        job_log.error("Vina timed out after 2 hours")
        job.errors.append("Vina execution timed out")
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.TIMEOUT
        return job
    except Exception as e:
        job.elapsed_seconds = time.time() - start_time
        job_log.error(f"Vina execution error: {e}")
        job.errors.append(str(e))
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.FAILED
        return job

    # ── Step 8: Check exit status ─────────────────────────────────────────────
    if job.exit_code != 0:
        job_log.error(f"Vina exited with code {job.exit_code}")
        # Extract the meaningful error message from stdout/stderr, skipping the citation header
        out_lines = [l.strip() for l in (job.stdout or "").splitlines() if l.strip()]
        error_lines = [l for l in out_lines if any(k in l.lower() for k in ("error", "failed", "missing", "fatal", "cannot", "syntax"))]
        if error_lines:
            stderr_snippet = " ".join(error_lines[-3:])
        elif out_lines:
            stderr_snippet = out_lines[-1]
        else:
            stderr_snippet = ""

        if stderr_snippet:
            job_log.error(f"Error output:\n{stderr_snippet}")
            job.errors.append(f"Non-zero exit code: {job.exit_code} — {stderr_snippet}")
        else:
            job.errors.append(f"Non-zero exit code: {job.exit_code}")
        job.status = JobStatus.EXECUTION_FAILED
        job.execution_status = ExecutionStatus.FAILED
        return job

    job.execution_status = ExecutionStatus.SUCCESS

    # Copy output pose to central results/VINA directory
    central_vina_dir = config.result_directory / "VINA"
    if output_pdbqt.exists():
        try:
            central_vina_dir.mkdir(parents=True, exist_ok=True)
            dest_vina = central_vina_dir / f"{job.job_id}_out.pdbqt"
            shutil.copy2(output_pdbqt, dest_vina)
            print(f"  [Vina] ✔ Central Pose: {dest_vina}")
        except Exception:
            pass

    # ── Steps 9-10: Validate and parse output ─────────────────────────────────
    valid, model_count, messages = validate_vina_output(output_pdbqt)
    if not valid:
        job_log.error(f"Output validation failed: {messages}")
        job.errors.extend(messages)
        job.status = JobStatus.EXECUTION_FAILED
        job.analysis_status = AnalysisStatus.FAILED
        return job

    job.obtained_modes = model_count
    results = extract_affinities(output_pdbqt)
    job.vina_results = results

    job_log.info(f"Obtained {model_count} poses (requested {job.requested_modes})")
    job_log.info(format_vina_results_table(results))

    # ── Step 11: Compare requested vs. obtained poses ─────────────────────────
    if model_count == 0:
        job.status = JobStatus.EXECUTION_FAILED
        job.analysis_status = AnalysisStatus.FAILED
        job.errors.append("Zero poses obtained")
        return job
    elif model_count < job.requested_modes:
        warning_msg = (
            f"FEWER_MODES: Requested {job.requested_modes} modes, "
            f"obtained {model_count}. This is expected behavior — Vina may "
            f"return fewer modes based on search/output-selection criteria."
        )
        job.warnings.append(warning_msg)
        job_log.warning(warning_msg)
        job.analysis_status = AnalysisStatus.WARNING

    # ── Step 12: Run vina_split ───────────────────────────────────────────────
    split_success, split_files, split_msg = run_vina_split(
        config.vina_split_executable,
        output_pdbqt,
        poses_dir,
    )
    job_log.info(f"vina_split: {split_msg}")

    if split_success:
        job.split_poses = split_files
        # Step 13: Validate split output
        split_valid, split_val_msg = validate_split(split_files, model_count)
        job_log.info(f"Split validation: {split_val_msg}")
        if not split_valid:
            job.warnings.append(f"SPLIT_WARNING: {split_val_msg}")
    else:
        job.warnings.append(f"SPLIT_FAILED: {split_msg}")
        job_log.warning(f"Split failed: {split_msg}")
        # Split failure is a warning — the original multi-model PDBQT is preserved

    # ── Step 14: Record final status & Canonical Result ───────────────────────
    if job.analysis_status == AnalysisStatus.WARNING or job.warnings:
        job.status = JobStatus.SUCCESS_WITH_WARNING
    else:
        job.status = JobStatus.SUCCESS
        job.analysis_status = AnalysisStatus.SUCCESS

    # Build Authoritative Canonical Result Model
    canon = DockingResult(
        job_id=job.job_id,
        engine=Engine.VINA,
        docking_mode=job.docking_mode,
        success=(job.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)),
        summary={
            "engine": "AutoDock Vina",
            "docking_mode": job.docking_mode.value,
            "modes": len(job.vina_results),
            "best_affinity": job.vina_results[0].binding_affinity if job.vina_results else None,
        },
    )
    for r in job.vina_results:
        canon.poses.append(
            CanonicalPose(
                rank=r.pose,
                binding_score=r.binding_affinity,
                estimated_ki=None,
                estimated_ki_formatted="Not applicable — Vina does not report Ki",
                vina_metrics=VinaMetrics(
                    rmsd_lower_bound=r.rmsd_lower_bound,
                    rmsd_upper_bound=r.rmsd_upper_bound,
                ),
                source_model=r.pose,
                source_file=str(job.output_pdbqt) if job.output_pdbqt else None,
            )
        )
    canon.provenance = ProvenanceRecord(
        application_version=__version__,
        engine="AutoDock Vina",
        engine_version=get_vina_version(config.vina_executable),
        docking_mode=job.docking_mode.value,
        receptor_name=job.receptor_name,
        receptor_path=str(job.receptor_path) if job.receptor_path else None,
        receptor_hash=compute_file_hash(job.receptor_path),
        ligand_name=job.ligand_name,
        ligand_path=str(job.ligand_path) if job.ligand_path else None,
        ligand_hash=compute_file_hash(job.ligand_path),
        config_path=str(job.config_path) if job.config_path else None,
        config_hash=compute_file_hash(job.config_path),
        executable_path=str(config.vina_executable),
        executable_hash=compute_file_hash(config.vina_executable),
        elapsed_seconds=job.elapsed_seconds,
        output_files=[str(output_pdbqt)] if output_pdbqt else [],
        output_status=job.status.value,
        analysis_status=job.analysis_status.value,
    )
    job.canonical_result = canon

    # ── Step 15: Build receptor+ligand complex PDBQT files ────────────────────
    # These replace PyMol PSE/PML generation and work in any molecular viewer
    # (UCSF Chimera, VMD, Discovery Studio, PyMol manually, etc.)
    if job.output_pdbqt and job.output_pdbqt.exists() and job.receptor_path and job.receptor_path.exists():
        try:
            from complex_builder import build_job_complexes
            n_poses_to_build = len(job.vina_results) if job.vina_results else 0
            complex_dir = job.output_dir / "complexes"
            built = build_job_complexes(
                receptor_path=job.receptor_path,
                ligand_output_pdbqt=job.output_pdbqt,
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


# ═══════════════════════════════════════════════════════════════════════════════
# Retry Logic
# ═══════════════════════════════════════════════════════════════════════════════

def run_vina_with_retry(job: DockingJob, config: ProjectConfig) -> DockingJob:
    """Run a Vina job with optional retry on fewer-than-requested modes.

    If retry is enabled and the initial run produces fewer modes than requested,
    retries with increasing exhaustiveness values.

    Important:
    - Never claims increasing exhaustiveness guarantees more poses
    - Never fabricates poses
    - Preserves the original result
    - Records every retry attempt
    """
    job = run_vina_job(job, config)

    if not config.retry_enabled:
        return job

    # Only retry if we got poses but fewer than requested
    if (job.execution_status != ExecutionStatus.SUCCESS or
            job.obtained_modes >= job.requested_modes or
            job.obtained_modes == 0):
        return job

    original_modes = job.obtained_modes
    best_job = job

    for retry_exh in config.retry_exhaustiveness:
        job.retry_count += 1
        logger.info(
            f"Retry {job.retry_count} for {job.job_id}: "
            f"exhaustiveness {retry_exh} (previous: {original_modes} modes)"
        )

        job.retry_history.append({
            "attempt": job.retry_count,
            "exhaustiveness": retry_exh,
            "obtained_modes": job.obtained_modes,
        })

        original_exh = config.exhaustiveness
        config.exhaustiveness = retry_exh

        retry_job = DockingJob(
            job_id=job.job_id,
            receptor_name=job.receptor_name,
            ligand_name=job.ligand_name,
            receptor_path=job.receptor_path,
            ligand_path=job.ligand_path,
            config_path=job.config_path,
            flex_receptor_path=job.flex_receptor_path,
            output_dir=job.output_dir,
            engine=job.engine,
            docking_mode=job.docking_mode,
            requested_modes=job.requested_modes,
            retry_count=job.retry_count,
        )

        retry_job = run_vina_job(retry_job, config)
        config.exhaustiveness = original_exh

        if retry_job.obtained_modes > best_job.obtained_modes:
            best_job = retry_job
            best_job.retry_count = job.retry_count
            best_job.retry_history = job.retry_history

        if retry_job.obtained_modes >= job.requested_modes:
            break

    if best_job.obtained_modes < best_job.requested_modes:
        best_job.warnings.append(
            f"FEWER_MODES_FOUND: After {best_job.retry_count} retries, "
            f"best result is {best_job.obtained_modes} modes "
            f"(requested {best_job.requested_modes})."
        )

    return best_job


# ═══════════════════════════════════════════════════════════════════════════════
# Full Workflow Orchestration
# ═══════════════════════════════════════════════════════════════════════════════

def run_vina_workflow(
    config: ProjectConfig,
    resume_mode: ResumeMode = ResumeMode.RESUME,
    progress_callback: Optional[Any] = None,
) -> List[DockingJob]:
    """Run the complete Vina docking workflow.

    1. Build job queue (N receptors × M ligands)
    2. Validate all inputs
    3. Handle resume
    4. Execute jobs sequentially
    5. Generate reports

    Returns:
        List of all DockingJob objects with final status.
    """
    # Build job queue
    all_jobs = build_job_queue(config)
    if not all_jobs:
        logger.error("No docking jobs to execute.")
        return []

    logger.info(f"Built {len(all_jobs)} docking jobs")

    # Load previous status for resume
    previous_status = load_job_status(config.result_directory)

    # Filter jobs based on resume mode with physical output verification
    jobs_to_run = filter_jobs(all_jobs, resume_mode, previous_status, validate_outputs=True)
    skipped = len(all_jobs) - len(jobs_to_run)

    # Display job summary
    receptors = set(j.receptor_name for j in all_jobs)
    ligands = set(j.ligand_name for j in all_jobs)
    rigid_count = sum(1 for j in all_jobs if j.docking_mode == DockingMode.RIGID)
    flex_count = sum(1 for j in all_jobs if j.docking_mode == DockingMode.FLEXIBLE)

    print()
    print("=" * 60)
    print(f"  {SUITE_BANNER}")
    print(f"  Engine: AutoDock Vina")
    print("=" * 60)
    print(f"  Receptors   : {len(receptors)}")
    print(f"  Ligands     : {len(ligands)}")
    print(f"  Total jobs  : {len(all_jobs)}")
    print(f"    ↳ Rigid   : {rigid_count}")
    print(f"    ↳ Flexible: {flex_count}")
    if skipped:
        print(f"  Skipped     : {skipped} (previously completed)")
    print(f"  To execute  : {len(jobs_to_run)}")
    print("=" * 60)

    # Validate executables before execution
    if not validate_executable(config.vina_executable, "AutoDock Vina"):
        config.vina_executable = resolve_executable(config.vina_executable, "vina.exe")
        if not validate_executable(config.vina_executable, "AutoDock Vina"):
            err_msg = (
                f"AutoDock Vina executable not found or invalid: '{config.vina_executable}'.\n"
                f"Please verify that vina.exe is placed in the bin/ directory or configure its path in Settings."
            )
            logger.error(err_msg)
            print(f"\n[ERROR] {err_msg}\n")
            raise RuntimeError(err_msg)

    # Validate all inputs before starting
    all_errors = []
    for job in jobs_to_run:
        errors = validate_job_inputs(job)
        all_errors.extend(errors)

    if all_errors:
        print(f"\n  INPUT VALIDATION ERRORS ({len(all_errors)}):")
        for err in all_errors:
            print(f"    [FAIL] {err}")

    # Dry-run mode
    if config.dry_run:
        print("\n  DRY-RUN MODE — No docking will be executed.")
        print(f"\n  Commands that would be executed:")
        for i, job in enumerate(jobs_to_run, 1):
            mode_label = f"[{job.docking_mode.value}]"
            print(f"\n  [{i}/{len(jobs_to_run)}] {job.job_id} {mode_label}")
            if job.receptor_path:
                print(f"    Receptor (rigid): {job.receptor_path}")
            if job.docking_mode == DockingMode.FLEXIBLE and job.flex_receptor_path:
                print(f"    Receptor (flex) : {job.flex_receptor_path}")
            if job.ligand_path:
                print(f"    Ligand          : {job.ligand_path}")
            if job.config_path:
                print(f"    Config          : {job.config_path}")
        return all_jobs

    if all_errors:
        fatal = [e for e in all_errors if "not found" in e.lower() or "does not exist" in e.lower()]
        if fatal:
            logger.error("Fatal input errors detected. Fix before continuing.")
            return all_jobs

    # Execute jobs (parallel or sequential)
    runner = run_vina_with_retry if config.retry_enabled else run_vina_job

    if config.max_workers > 1:
        run_jobs_parallel(
            jobs_to_run=jobs_to_run,
            all_jobs=all_jobs,
            config=config,
            runner_func=runner,
            engine_name="AutoDock Vina",
            skipped=skipped,
            progress_callback=progress_callback,
        )
    else:
        progress = ProgressTracker(len(all_jobs), "AutoDock Vina")
        progress.completed = skipped

        for i, job in enumerate(jobs_to_run, 1):
            if process_manager.is_cancelled():
                logger.warning("Vina workflow cancelled by user. Halting remaining jobs.")
                print("\n  [Vina] ⏹ Workflow cancelled by user. Halting remaining jobs.")
                break
            progress.display(skipped + i, job)
            job = runner(job, config)
            progress.update(job)

            if progress_callback:
                try:
                    progress_callback(skipped + i, len(all_jobs), f"{job.receptor_name} × {job.ligand_name}")
                except Exception:
                    pass

            # Save status after each job (for resume capability)
            save_job_status(config.result_directory, all_jobs)

            # Report result
            mode_label = f"[{job.docking_mode.value}]"
            if job.status == JobStatus.SUCCESS:
                affinity_str = (
                    f"best affinity: {job.vina_results[0].binding_affinity:.1f} kcal/mol"
                    if job.vina_results else ""
                )
                print(f"  [OK] {job.job_id} {mode_label}: {job.obtained_modes} poses  {affinity_str}")
            elif job.status == JobStatus.SUCCESS_WITH_WARNING:
                print(f"  [WARN] {job.job_id} {mode_label}: {job.obtained_modes}/{job.requested_modes} poses")
                for w in job.warnings:
                    print(f"    -> {w}")
            else:
                print(f"  [FAIL] {job.job_id} {mode_label}: {job.status.value}")
                for e in job.errors:
                    print(f"    -> {e}")

        # Final summary
        progress.display_summary()

    # Generate reports
    from reporting import generate_reports
    generate_reports(all_jobs, config)

    # Save run metadata
    _save_run_metadata(config, all_jobs)

    return all_jobs


def _save_run_metadata(config: ProjectConfig, jobs: List[DockingJob]) -> None:
    """Save run metadata JSON for reproducibility."""
    import os
    import platform
    from ad4_compatibility import build_job_provenance

    vina_version = get_vina_version(config.vina_executable)
    vina_split_version = get_vina_split_version(config.vina_split_executable)

    rigid_jobs = [j for j in jobs if j.docking_mode == DockingMode.RIGID]
    flex_jobs = [j for j in jobs if j.docking_mode == DockingMode.FLEXIBLE]

    metadata = {
        "suite": "AutoDock Suite Pro",
        "suite_version": __version__,
        "timestamp": datetime.now().isoformat(),
        "operating_system": platform.platform(),
        "python_version": platform.python_version(),
        "vina_version": vina_version,
        "vina_split_version": vina_split_version,
        "vina_executable": str(config.vina_executable),
        "vina_split_executable": str(config.vina_split_executable),
        "project_root": str(config.project_root),
        "receptor_directory": str(config.receptor_directory),
        "ligand_directory": str(config.ligand_directory),
        "receptors": sorted(set(j.receptor_name for j in jobs)),
        "ligands": sorted(set(j.ligand_name for j in jobs)),
        "total_jobs": len(jobs),
        "rigid_jobs": len(rigid_jobs),
        "flexible_jobs": len(flex_jobs),
        "engine": config.engine.value,
        "exhaustiveness": config.exhaustiveness,
        "num_modes": config.num_modes,
        "energy_range": config.energy_range,
        "cpu": config.cpu,
        "seed": config.seed,
        "retry_enabled": config.retry_enabled,
        "retry_exhaustiveness": config.retry_exhaustiveness,
        "atom_type_provenance": [
            {
                "job_id": j.job_id,
                "compatibility": build_job_provenance(
                    [("receptor", j.receptor_path), ("ligand", j.ligand_path)]
                    + ([ ("flex_receptor", j.flex_receptor_path) ]
                       if j.flex_receptor_path and j.flex_receptor_path.exists() else []),
                    profile=None,
                    engine="VINA",
                ),
            }
            for j in jobs
        ],
        "specialized_metal_workflow_selected": False,
    }

    metadata_path = config.result_directory / "run_metadata.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(
        json.dumps(metadata, indent=2, default=str),
        encoding="utf-8"
    )
    logger.info(f"Run metadata saved: {metadata_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI Entry Point
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="vina_workflow.py",
        description=f"{SUITE_BANNER} — AutoDock Vina Workflow",
    )
    parser.add_argument(
        "--config", required=True, type=Path,
        help="Path to TOML configuration file"
    )
    parser.add_argument("--dry-run", action="store_true", help="Validate without executing")
    parser.add_argument("--resume", action="store_true", help="Skip previously completed jobs")
    parser.add_argument("--rerun-failed", action="store_true", help="Rerun failed jobs")
    parser.add_argument("--force", action="store_true", help="Force rerun all jobs")
    parser.add_argument("--version", action="version", version=SUITE_BANNER)

    args = parser.parse_args()

    config = load_config(args.config)
    config.engine = Engine.VINA

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

    run_vina_workflow(config, resume_mode)


if __name__ == "__main__":
    main()
