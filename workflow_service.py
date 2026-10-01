#!/usr/bin/env python3
"""
AutoDock Suite Pro — Unified Workflow Service
=============================================
Central workflow orchestrator shared identically by CLI (main.py) and PySide6 GUI (runner_tab.py).
Supports VINA, AUTODOCK4, and combined BOTH workflows with real-time cancellation,
verifiable checkpoints, and authoritative single-pass analysis.
"""

from __future__ import annotations

import copy
import logging
from typing import Any, Callable, List, Optional

from config import ProjectConfig
from models import DockingJob, Engine, ResumeMode
from process_manager import process_manager

logger = logging.getLogger("docking_automation.workflow_service")


def execute_workflow(
    config: ProjectConfig,
    resume_mode: ResumeMode = ResumeMode.RESUME,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
) -> List[DockingJob]:
    """Execute docking workflow for VINA, AUTODOCK4, or BOTH.

    Both CLI and GUI call this single authoritative service.
    """
    process_manager.reset()
    all_jobs: List[DockingJob] = []

    if config.engine == Engine.VINA:
        from vina_workflow import run_vina_workflow
        all_jobs = run_vina_workflow(
            config,
            resume_mode=resume_mode,
            progress_callback=progress_callback,
        ) or []

    elif config.engine == Engine.AUTODOCK4:
        from autodock4_workflow import run_autodock4_workflow
        all_jobs = run_autodock4_workflow(
            config,
            resume_mode=resume_mode,
            progress_callback=progress_callback,
        ) or []

    elif config.engine == Engine.BOTH:
        print()
        print("=" * 70)
        print("  DUAL-ENGINE WORKFLOW — AutoDock Vina + AutoDock 4.2.6 (sequential)")
        print("=" * 70)

        # Phase 1: Vina
        vina_config = copy.copy(config)
        vina_config.engine = Engine.VINA

        # Phase 2 config (needed for pre-count)
        ad4_config = copy.copy(config)
        ad4_config.engine = Engine.AUTODOCK4

        # Pre-count both job queues so progress spans both phases as a single total
        from job_manager import build_job_queue
        try:
            _vina_count = len(build_job_queue(vina_config))
        except Exception:
            _vina_count = 0
        try:
            _ad4_count = len(build_job_queue(ad4_config))
        except Exception:
            _ad4_count = 0
        _grand_total = max(1, _vina_count + _ad4_count)

        def _vina_progress(completed: int, total: int, label: str = "") -> None:
            if progress_callback:
                progress_callback(completed, _grand_total, f"[Vina] {label}")

        def _ad4_progress(completed: int, total: int, label: str = "") -> None:
            if progress_callback:
                progress_callback(_vina_count + completed, _grand_total, f"[AutoDock4] {label}")

        from vina_workflow import run_vina_workflow
        vina_jobs = run_vina_workflow(
            vina_config,
            resume_mode=resume_mode,
            progress_callback=_vina_progress if progress_callback else None,
        ) or []

        if process_manager.is_cancelled():
            logger.warning("Workflow cancelled after Vina phase. Halting AutoDock4 phase.")
            return vina_jobs

        # Phase 2: AutoDock4
        print("\n" + "=" * 70)
        print("  DUAL-ENGINE WORKFLOW — Phase 2 of 2: AutoDock 4.2.6")
        print("=" * 70)
        from autodock4_workflow import run_autodock4_workflow
        ad4_jobs = run_autodock4_workflow(
            ad4_config,
            resume_mode=resume_mode,
            progress_callback=_ad4_progress if progress_callback else None,
        ) or []

        all_jobs = (vina_jobs or []) + (ad4_jobs or [])

        # Phase 3: Combined Reporting
        if all_jobs:
            from reporting import generate_reports
            generate_reports(all_jobs, config)

        print()
        print("=" * 70)
        print("  ✔  DUAL-ENGINE RUN COMPLETE")
        print("=" * 70)

    else:
        logger.error(f"Unsupported engine: {config.engine}")

    return all_jobs
