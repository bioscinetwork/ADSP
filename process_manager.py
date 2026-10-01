#!/usr/bin/env python3
"""
AutoDock Suite Pro — Process Manager & Checkpoint Hash Verifier
===============================================================
Manages subprocess lifecycles, real-time process termination, cancellation propagation,
cryptographic input hashing (SHA-256), and verifiable resume checkpoints.
"""

from __future__ import annotations

import hashlib
import logging
import os
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from models import DockingJob, Engine

logger = logging.getLogger("docking_automation.process_manager")


# ═══════════════════════════════════════════════════════════════════════════════
# Cryptographic Input Hashing
# ═══════════════════════════════════════════════════════════════════════════════

def compute_file_hash(path: Optional[Path | str]) -> str:
    """Compute deterministic SHA-256 hash of a file's content.

    Returns empty string if file does not exist or cannot be read.
    """
    if not path:
        return ""
    p = Path(path)
    if not p.is_file():
        return ""
    hasher = hashlib.sha256()
    try:
        with open(p, "rb") as f:
            while chunk := f.read(65536):
                hasher.update(chunk)
        return hasher.hexdigest()
    except OSError as e:
        logger.warning(f"Could not hash file {p}: {e}")
        return ""


def compute_job_input_hashes(job: "DockingJob") -> Dict[str, str]:
    """Compute SHA-256 hashes of all input files for a docking job."""
    hashes: Dict[str, str] = {
        "receptor_hash": compute_file_hash(job.receptor_path),
        "ligand_hash": compute_file_hash(job.ligand_path),
        "config_hash": compute_file_hash(job.config_path),
    }
    if getattr(job, "flex_receptor_path", None):
        hashes["flex_receptor_hash"] = compute_file_hash(job.flex_receptor_path)
    return hashes


# ═══════════════════════════════════════════════════════════════════════════════
# Output Validation for Resuming
# ═══════════════════════════════════════════════════════════════════════════════

def validate_job_outputs(job: "DockingJob") -> Tuple[bool, str]:
    """Validate that the underlying outputs for a completed job actually exist and are non-empty.

    Returns:
        (is_valid, reason_if_invalid)
    """
    from models import Engine

    if job.engine == Engine.VINA:
        out_p = job.output_pdbqt
        if not out_p:
            return False, "Missing output_pdbqt path"
        p = Path(out_p)
        if not p.is_file():
            return False, f"Output file does not exist: {p}"
        if p.stat().st_size == 0:
            return False, f"Output file is empty: {p}"
        # Quick check for model record
        try:
            content = p.read_text(encoding="utf-8", errors="replace")[:2048]
            if "MODEL" not in content and "ATOM" not in content and "HETATM" not in content:
                return False, f"Output file does not contain valid coordinate models: {p}"
        except OSError as e:
            return False, f"Cannot read output file {p}: {e}"
        return True, ""

    elif job.engine == Engine.AUTODOCK4:
        dlg_p = job.dlg_path
        if not dlg_p and hasattr(job, "ad4_results") and isinstance(job.ad4_results, dict):
            dlg_str = job.ad4_results.get("dlg_path")
            if dlg_str:
                dlg_p = Path(dlg_str)
        if not dlg_p:
            return False, "Missing dlg_path for AutoDock4 job"
        p = Path(dlg_p)
        if not p.is_file():
            return False, f"AutoDock4 DLG file does not exist: {p}"
        if p.stat().st_size < 100:
            return False, f"AutoDock4 DLG file is incomplete or truncated: {p}"
        return True, ""

    return True, ""


# ═══════════════════════════════════════════════════════════════════════════════
# Process Termination & Lifecycle Management
# ═══════════════════════════════════════════════════════════════════════════════

def terminate_process(proc: Optional[subprocess.Popen], timeout: float = 3.0) -> None:
    """Safely terminate a subprocess and all child processes on Windows and Unix."""
    if proc is None or proc.poll() is not None:
        return

    pid = proc.pid
    logger.info(f"Terminating subprocess PID {pid}...")

    if sys.platform == "win32":
        try:
            # taskkill /F /T kills the process and all child processes it spawned
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(pid)],
                capture_output=True,
                timeout=timeout,
            )
        except Exception as e:
            logger.debug(f"taskkill failed: {e}; falling back to proc.kill()")
            try:
                proc.kill()
            except Exception:
                pass
    else:
        try:
            import signal
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except Exception:
            try:
                proc.kill()
            except Exception:
                pass

    try:
        proc.wait(timeout=2.0)
    except Exception:
        pass


class SubprocessManager:
    """Thread-safe registry for managing active docking subprocesses and cancellation."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._active_processes: Dict[str, subprocess.Popen] = {}
        self._cancel_event = threading.Event()

    def register(self, job_id: str, proc: subprocess.Popen) -> None:
        """Register an active subprocess for a job."""
        with self._lock:
            self._active_processes[job_id] = proc

    def unregister(self, job_id: str) -> None:
        """Unregister a finished subprocess for a job."""
        with self._lock:
            self._active_processes.pop(job_id, None)

    def request_cancellation(self) -> None:
        """Signal cancellation and terminate all currently registered docking subprocesses."""
        self._cancel_event.set()
        logger.warning("Cancellation requested! Terminating all active docking processes...")
        with self._lock:
            for job_id, proc in list(self._active_processes.items()):
                logger.info(f"Cancelling active subprocess for job {job_id} (PID {proc.pid})")
                terminate_process(proc)
            self._active_processes.clear()

    def is_cancelled(self) -> bool:
        """Check if cancellation has been requested."""
        return self._cancel_event.is_set()

    def reset(self) -> None:
        """Reset the cancellation state before a new run."""
        self._cancel_event.clear()
        with self._lock:
            self._active_processes.clear()


# Global singleton instance
process_manager = SubprocessManager()
