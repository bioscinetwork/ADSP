#!/usr/bin/env python3
"""
Docking Automation Suite — Logging Utilities
=============================================
Sets up global and per-job loggers with file + console output.
"""

from __future__ import annotations

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_global_logger(
    log_dir: Path,
    run_timestamp: Optional[str] = None,
    console_level: int = logging.INFO,
) -> logging.Logger:
    """Create and configure the global run logger.

    Writes to both a file and the console.

    Args:
        log_dir: Directory to store the log file.
        run_timestamp: Timestamp string for the log filename.
        console_level: Logging level for console output.

    Returns:
        Configured logger.
    """
    log_dir.mkdir(parents=True, exist_ok=True)

    if run_timestamp is None:
        run_timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    log_file = log_dir / f"run_{run_timestamp}.log"

    logger = logging.getLogger("docking_automation")
    logger.setLevel(logging.DEBUG)

    # Clear any existing handlers
    logger.handlers.clear()

    # File handler — captures everything (DEBUG+)
    file_handler = logging.FileHandler(str(log_file), encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(
        "[%(asctime)s] [%(levelname)-8s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logger.addHandler(file_handler)

    # Console handler — info and above
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_handler.setFormatter(logging.Formatter(
        "[%(levelname)-8s] %(message)s"
    ))
    logger.addHandler(console_handler)

    logger.info(f"Log file: {log_file}")
    return logger


def setup_job_logger(
    log_dir: Path,
    engine: str,
    job_id: str,
) -> logging.Logger:
    """Create a per-job logger that writes to a dedicated log file.

    Args:
        log_dir: Base log directory.
        engine: "VINA" or "AD4".
        job_id: The job identifier.

    Returns:
        Configured logger for this job.
    """
    job_log_dir = log_dir / engine
    job_log_dir.mkdir(parents=True, exist_ok=True)

    log_file = job_log_dir / f"{job_id}.log"

    logger_name = f"docking_automation.{engine}.{job_id}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers to avoid duplicates on resume
    logger.handlers.clear()

    file_handler = logging.FileHandler(str(log_file), encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(
        "[%(asctime)s] [%(levelname)-8s] [%(name)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logger.addHandler(file_handler)

    # Propagate to global logger
    logger.propagate = True

    return logger


def get_run_timestamp() -> str:
    """Generate a formatted timestamp for the current run."""
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
