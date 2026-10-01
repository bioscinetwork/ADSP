#!/usr/bin/env python3
"""
AutoDock Suite Pro — Unified CLI Entry Point
=============================================
Command-line interface supporting Vina and AutoDock 4 workflows,
configuration generation, project scaffolding, file preparation,
validation, dry runs, reporting, and resuming.

Usage:
    python main.py --scaffold                           Create project folder structure
    python main.py --init                               Generate starter project_config.toml
    python main.py --prepare                            Auto-generate missing config.txt files (Vina)
    python main.py --prepare --engine AUTODOCK4         Auto-generate GPF stubs (AutoDock4)
    python main.py --validate --config project_config.toml
    python main.py --dry-run --config project_config.toml
    python main.py --config project_config.toml
    python main.py --config project_config.toml --engine AUTODOCK4
    python main.py --config project_config.toml --resume
    python main.py --config project_config.toml --rerun-failed
    python main.py --config project_config.toml --force
    python main.py --config project_config.toml --report-only
    python main.py --version
"""

from __future__ import annotations

import argparse
import logging
import platform
import sys

# Configure UTF-8 output on Windows consoles
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
from pathlib import Path
from typing import List, Optional

from executables import (
    configure_openbabel_env, get_autodock4_version, get_autogrid4_version,
    get_cpu_count, get_vina_split_version, get_vina_version,
    resolve_executable, validate_executable,
)

# Auto-configure OpenBabel environment (BABEL_LIBDIR, BABEL_DATADIR) immediately
configure_openbabel_env()

from config import (
    ProjectConfig, generate_template_config, load_config, validate_config,
)
from job_manager import (
    ProgressTracker, build_job_queue, discover_ligands,
    discover_receptors, filter_jobs, load_docking_jobs, load_job_status,
)
from logging_utils import get_run_timestamp, setup_global_logger
from models import (
    SUITE_BANNER, DockingMode, Engine, JobStatus, ResumeMode, __version__,
)
from reporting import generate_reports

logger = logging.getLogger("docking_automation")


# ═══════════════════════════════════════════════════════════════════════════════
# Helper Routines
# ═══════════════════════════════════════════════════════════════════════════════

def print_banner() -> None:
    """Print the suite title banner."""
    print()
    print("=" * 70)
    print(f"  {SUITE_BANNER}")
    print("  Cross-Engine Molecular Docking Automation Platform")
    print("=" * 70)


def print_version_info() -> None:
    """Display suite version, platform, and detected engine versions."""
    print_banner()
    print(f"  Python           : {platform.python_version()} ({platform.architecture()[0]})")
    print(f"  Operating System : {platform.platform()}")
    print(f"  CPU Logical Cores: {get_cpu_count()}")
    print("-" * 70)

    # Auto-detect executables across standard locations
    default_dir = Path("C:/Program Files (x86)/MGLTools-1.5.7")
    vina_path = resolve_executable(default_dir / "vina.exe", "vina.exe")
    vina_split_path = resolve_executable(default_dir / "vina_split.exe", "vina_split.exe")
    autogrid_path = resolve_executable(default_dir / "autogrid4.exe", "autogrid4.exe")
    autodock_path = resolve_executable(default_dir / "autodock4.exe", "autodock4.exe")

    vina_ver = get_vina_version(vina_path) if vina_path.exists() else "Not found"
    split_ver = get_vina_split_version(vina_split_path) if vina_split_path.exists() else "Not found"
    ag4_ver = get_autogrid4_version(autogrid_path) if autogrid_path.exists() else "Not found"
    ad4_ver = get_autodock4_version(autodock_path) if autodock_path.exists() else "Not found"

    print(f"  AutoDock Vina    : {vina_ver} (at {vina_path})")
    print(f"  Vina Split       : {split_ver} (at {vina_split_path})")
    print(f"  AutoGrid 4       : {ag4_ver} (at {autogrid_path})")
    print(f"  AutoDock 4       : {ad4_ver} (at {autodock_path})")
    print("=" * 70)


def perform_validation(config: ProjectConfig) -> bool:
    """Validate executables, directories, and configuration settings."""
    print()
    print("=" * 70)
    print("  CONFIGURATION & ENVIRONMENT VALIDATION")
    print("=" * 70)
    print(f"  Project Name   : {config.project_name or '(unnamed)'}")
    print(f"  Project Root   : {config.project_root}")
    print(f"  Engine         : {config.engine.value}")
    print(f"  Docking Mode   : {config.docking_mode.value}")
    print(f"  Receptors Dir  : {config.receptor_directory}")
    print(f"  Ligands Dir    : {config.ligand_directory}")
    print(f"  Results Dir    : {config.result_directory}")
    print(f"  Logs Dir       : {config.log_directory}")
    print(f"  Reports Dir    : {config.report_directory}")
    print("-" * 70)

    errors = validate_config(config)

    if errors:
        print("  [FAIL] VALIDATION FAILED:")
        for err in errors:
            print(f"    - {err}")
        print("=" * 70)
        return False

    print("  [OK] All configuration checks passed.")

    # Check input files presence
    rec_count = len(discover_receptors(config)) if config.receptor_directory.exists() else 0
    lig_count = len(discover_ligands(config)) if config.ligand_directory.exists() else 0
    print(f"  [OK] Receptors found : {rec_count}")
    print(f"  [OK] Ligands found   : {lig_count} (*.pdbqt)")
    if rec_count == 0:
        print(
            "    [WARN] No receptors found. Each receptor needs a subfolder with "
            "rigid/ containing a .pdbqt file.\n"
            "    Run: python main.py --scaffold  to create the folder structure."
        )
    if lig_count == 0:
        print("    [WARN] No ligand PDBQT files found in ligand directory.")
    print("=" * 70)
    return True


def perform_dry_run(config: ProjectConfig) -> None:
    """Simulate workflow execution without launching external docking tools."""
    print()
    print("=" * 70)
    print("  DRY-RUN SIMULATION")
    print("=" * 70)

    is_valid = perform_validation(config)
    if not is_valid:
        print("  Simulation aborted due to validation errors.")
        return

    jobs = build_job_queue(config)
    print(f"\n  Job Queue Construction:")
    print(f"  - Total Planned Jobs : {len(jobs)}")

    if not jobs:
        print("  [WARN] No jobs queued. Check receptor/ligand directories.")
        return

    print("\n  Sample Planned Jobs (first 5):")
    for i, job in enumerate(jobs[:5], 1):
        print(f"    [{i}] Job ID: {job.job_id}")
        rec_str = job.receptor_path.name if job.receptor_path else "(none)"
        lig_str = job.ligand_path.name if job.ligand_path else "(none)"
        mode_str = job.docking_mode.value
        print(f"        Receptor   : {rec_str}  [{mode_str}]")
        print(f"        Ligand     : {lig_str}")
        if config.engine == Engine.VINA:
            cfg_str = str(job.config_path) if job.config_path else "(none)"
            out_str = str(job.output_pdbqt) if job.output_pdbqt else "(will be set at runtime)"
            print(f"        Config     : {cfg_str}")
            if job.flex_receptor_path:
                flex_str = str(job.flex_receptor_path)
                print(f"        Flex PDBQT : {flex_str}")

    if len(jobs) > 5:
        print(f"    ... and {len(jobs) - 5} more jobs.")

    print("\n  Execution Strategy:")
    print(f"  - Mode           : {'Sequential (safe)' if config.sequential else 'Parallel'}")
    print(f"  - Resume Enabled : {config.resume}")
    print(f"  - Reporting      : CSV={config.export_csv}, Excel={config.export_excel}")
    print("=" * 70)


def perform_report_only(config: ProjectConfig) -> None:
    """Generate reports from existing job_status.json in the result directory."""
    print()
    print("=" * 70)
    print("  REPORT-ONLY MODE")
    print("=" * 70)

    jobs = load_docking_jobs(config.result_directory)
    if not jobs:
        print(f"  [FAIL] No existing job_status.json found in {config.result_directory}")
        print("    Cannot generate reports without previous docking runs.")
        return

    print(f"  Loaded {len(jobs)} jobs from {config.result_directory / 'job_status.json'}")
    generate_reports(jobs, config)


# ═══════════════════════════════════════════════════════════════════════════════
# Main Entry Point
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="main.py",
        description=f"{SUITE_BANNER} — Cross-Engine Docking Automation Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--config", type=Path, default=Path("project_config.toml"),
        help="Path to TOML configuration file (default: project_config.toml)"
    )
    parser.add_argument(
        "--engine", choices=["VINA", "AUTODOCK4", "BOTH"],
        help="Override engine selection in configuration (BOTH runs Vina then AutoDock4)"
    )
    parser.add_argument(
        "--mode", type=str,
        choices=["AUTO", "RIGID", "FLEXIBLE", "BOTH", "auto", "rigid", "flexible", "both"],
        help="Docking mode: AUTO (default: detects folders; runs both if both exist), RIGID, FLEXIBLE, or BOTH"
    )
    parser.add_argument(
        "--scaffold", action="store_true",
        help="Create the project folder structure (receptors/, ligands/, etc.)"
    )
    parser.add_argument(
        "--init", nargs="?", const="project_config.toml", type=str,
        help="Generate a starter project_config.toml template"
    )
    parser.add_argument(
        "--prepare", action="store_true",
        help=(
            "Auto-generate missing config.txt files (Vina) or GPF stubs (AutoDock4) "
            "for all receptors in the receptor directory"
        )
    )
    parser.add_argument(
        "--validate", action="store_true",
        help="Validate executables, paths, and config without running docking"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate and simulate job queue setup without running docking"
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume execution, skipping previously completed jobs"
    )
    parser.add_argument(
        "--rerun-failed", action="store_true",
        help="Rerun only jobs that previously failed or produced errors"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Force rerun of all jobs regardless of prior completion"
    )
    parser.add_argument(
        "--report-only", action="store_true",
        help="Generate CSV and Excel reports from an existing result directory"
    )
    parser.add_argument(
        "--receptor", type=str,
        help="Filter to a single receptor name"
    )
    parser.add_argument(
        "--ligand", type=str,
        help="Filter to a single ligand name"
    )
    parser.add_argument(
        "--parallel", "-p",
        nargs="?",
        const="AUTO",
        default=None,
        metavar="WORKERS",
        help="Run docking jobs in parallel using N workers (e.g. --parallel 4, or --parallel for AUTO)",
    )
    parser.add_argument(
        "--gui", action="store_true",
        help="Launch the AutoDock Suite Pro modern Desktop GUI",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose debug logging"
    )
    parser.add_argument(
        "--version", action="store_true",
        help="Show suite version and environment information"
    )

    # If run with no CLI flags on desktop, auto-launch GUI
    if len(sys.argv) == 1:
        try:
            from gui import launch_gui as _auto_launch_gui
            cfg_path = Path("project_config.toml")
            cfg = load_config(cfg_path) if cfg_path.exists() else None
            _auto_launch_gui(cfg)
            sys.exit(0)
        except Exception:
            pass  # Fall back to CLI parsing if headless or error

    args = parser.parse_args()

    # 1. Version action
    if args.version:
        print_version_info()
        sys.exit(0)

    # 2. Scaffold action (does not need a config file)
    if args.scaffold:
        print_banner()
        config_path = args.config
        if config_path.exists():
            config = load_config(config_path)
            project_root = config.project_root
        else:
            project_root = config_path.parent
        from prepare import create_scaffold
        create_scaffold(project_root)
        sys.exit(0)

    # 3. Init action
    if args.init:
        target_path = Path(args.init)
        if target_path.is_dir():
            target_path = target_path / "project_config.toml"
        print(f"Generating template configuration: {target_path}")
        generate_template_config(target_path)
        print("[OK] Template created. Edit the paths and parameters before running.")
        sys.exit(0)

    # 4. Load configuration
    config_path = args.config
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        print("Hint: Run with '--init' to create a starter configuration template.")
        print("Hint: Run with '--scaffold' to create the full project folder structure.")
        sys.exit(1)

    config = load_config(config_path)

    # Command-line overrides
    if args.engine:
        config.engine = Engine[args.engine]
    if args.mode:
        config.docking_mode = DockingMode[args.mode.upper()]
    if args.dry_run:
        config.dry_run = True
    if args.receptor:
        config.selected_receptors = [args.receptor]
    if args.ligand:
        config.selected_ligands = [args.ligand]

    if args.parallel is not None:
        if str(args.parallel).upper() == "AUTO":
            import os
            config.max_workers = max(1, (os.cpu_count() or 2) // 2)
        else:
            try:
                config.max_workers = max(1, int(args.parallel))
            except ValueError:
                config.max_workers = 2
        config.sequential = (config.max_workers == 1)

    if args.gui:
        from gui import launch_gui
        launch_gui(config)
        sys.exit(0)

    print_banner()

    # 5. Prepare action
    if args.prepare:
        from prepare import run_prepare
        run_prepare(config, args.engine)
        sys.exit(0)

    # 6. Validate action
    if args.validate:
        is_valid = perform_validation(config)
        sys.exit(0 if is_valid else 1)

    # Setup global logging
    run_ts = get_run_timestamp()
    log_level = logging.DEBUG if args.verbose else logging.INFO
    try:
        setup_global_logger(config.log_directory, run_ts, console_level=log_level)
    except OSError as e:
        print(f"\n  [ERROR] Error initializing log directory '{config.log_directory}': {e}")
        print("    Please check 'log_directory' and 'root' in your project configuration.")
        sys.exit(1)

    # 7. Dry-run action
    if args.dry_run:
        perform_dry_run(config)
        sys.exit(0)

    # 8. Report-only action
    if args.report_only:
        perform_report_only(config)
        sys.exit(0)

    # 9. Validate before execution
    errors = validate_config(config)
    if errors:
        print("\n  CONFIGURATION ERRORS:")
        for err in errors:
            print(f"    [FAIL] {err}")
        sys.exit(1)

    # Determine resume mode
    if args.force:
        resume_mode = ResumeMode.FORCE
    elif args.rerun_failed:
        resume_mode = ResumeMode.RERUN_FAILED
    elif args.resume or config.resume:
        resume_mode = ResumeMode.RESUME
    else:
        resume_mode = ResumeMode.FORCE

    # 10. Execute selected engine workflow via unified workflow service
    from workflow_service import execute_workflow
    execute_workflow(config, resume_mode=resume_mode)


if __name__ == "__main__":
    main()
