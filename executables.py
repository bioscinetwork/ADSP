#!/usr/bin/env python3
"""
Docking Automation Suite — Executable Management
==================================================
Validates docking executables and extracts version information.
Supports self-contained standalone execution, dynamic directory detection,
and fallback across portable PyInstaller bundles and system installations.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger("docking_automation.executables")

# Candidate directories across diverse Windows PCs and installations
STANDARD_SYSTEM_DIRS = [
    Path("C:/Program Files (x86)/MGLTools-1.5.7"),
    Path("C:/Program Files/MGLTools-1.5.7"),
    Path("C:/ADT3"),
    Path("C:/Screening"),
    Path("C:/Program Files (x86)/The Scripps Research Institute/Vina"),
    Path("C:/Program Files/The Scripps Research Institute/Vina"),
    Path("C:/Program Files/Dockamon-PyRx/1.0/PyRx"),
]

DEFAULT_MGLTOOLS_DIR = Path("C:/Program Files (x86)/MGLTools-1.5.7")


def get_candidate_search_dirs() -> List[Path]:
    """Dynamically collects candidate binary directories in priority order.

    Priority:
      1. Portable app runtime dirs (next to sys.executable)
      2. PyInstaller onefile / onedir internal dirs (_internal/bin, _MEIPASS)
      3. Script workspace dirs (next to __file__, cwd/bin)
      4. Standard external Windows installation directories
    """
    candidates: List[Path] = []

    # 1. Executable-relative (compiled standalone app)
    try:
        exe_dir = Path(sys.executable).parent.resolve()
        candidates.append(exe_dir / "bin")
        candidates.append(exe_dir / "_internal" / "bin")
        candidates.append(exe_dir)
    except Exception:
        pass

    # 2. PyInstaller temporary folder if running as onefile
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        try:
            m_path = Path(meipass).resolve()
            candidates.append(m_path / "bin")
            candidates.append(m_path / "_internal" / "bin")
            candidates.append(m_path)
        except Exception:
            pass

    # 3. Script-relative and current working directory
    try:
        script_dir = Path(__file__).resolve().parent
        candidates.append(script_dir / "bin")
        candidates.append(script_dir / "_internal" / "bin")
        candidates.append(script_dir / "dist" / "AutoDockSuitePro" / "bin")
        candidates.append(script_dir / "dist" / "AutoDockSuitePro" / "_internal" / "bin")
    except Exception:
        pass

    try:
        cwd = Path.cwd().resolve()
        candidates.append(cwd / "bin")
        candidates.append(cwd / "_internal" / "bin")
        candidates.append(cwd / "dist" / "AutoDockSuitePro" / "bin")
        candidates.append(cwd / "dist" / "AutoDockSuitePro" / "_internal" / "bin")
    except Exception:
        pass

    # 4. Standard external system installations
    candidates.extend(STANDARD_SYSTEM_DIRS)

    # De-duplicate while preserving order
    seen = set()
    unique_candidates: List[Path] = []
    for d in candidates:
        norm = str(d).lower()
        if norm not in seen:
            seen.add(norm)
            unique_candidates.append(d)

    return unique_candidates


def resolve_executable(configured_path: Optional[Path | str], exe_name: str) -> Path:
    """Resolve an executable path with multi-PC automatic discovery fallback.

    1. Checks configured_path if valid, non-empty, non-directory, and existing.
    2. Searches dynamic candidate directories (bundled bin/, _internal/bin, etc.).
    3. Checks the system PATH environment variable.
    4. Falls back to the highest-priority candidate location (never empty or ".").

    Args:
        configured_path: User-provided path from TOML config or GUI.
        exe_name: Binary name (e.g. 'vina.exe', 'vina_split.exe', 'autogrid4.exe', 'autodock4.exe').

    Returns:
        Path to existing executable if found, or best candidate path.
    """
    # 1. Test explicitly configured path
    if configured_path:
        str_path = str(configured_path).strip()
        if str_path and str_path != ".":
            p = Path(str_path)
            try:
                if p.is_file():
                    return p.resolve()
            except OSError:
                pass

    # 2. Search candidate directories
    search_dirs = get_candidate_search_dirs()
    for directory in search_dirs:
        candidate = directory / exe_name
        try:
            if candidate.is_file():
                logger.debug(f"Auto-detected {exe_name} at {candidate}")
                return candidate.resolve()
        except OSError:
            continue

    # 3. Check system PATH
    which_path = shutil.which(exe_name)
    if which_path:
        p = Path(which_path)
        if p.is_file():
            return p.resolve()

    # 4. Fallback to bundled or default MGLTools location (guaranteed not "." or empty)
    if configured_path and str(configured_path).strip() and str(configured_path).strip() != ".":
        return Path(configured_path).resolve()

    # Best guess: bin/ next to executable
    try:
        best_candidate = Path(sys.executable).parent / "bin" / exe_name
        return best_candidate.resolve()
    except Exception:
        return DEFAULT_MGLTOOLS_DIR / exe_name


def validate_executable(exe_path: Path | str, name: str) -> bool:
    """Check that an executable exists and is a valid executable file.

    Args:
        exe_path: Path to the executable.
        name: Human-readable name for error messages.

    Returns:
        True if the executable exists and is a file.
    """
    if not exe_path or str(exe_path).strip() in (".", ""):
        logger.error(f"{name}: path is empty or not configured.")
        return False

    p = Path(exe_path)
    try:
        if not p.exists():
            logger.error(f"{name}: file not found at {p}")
            return False
        if not p.is_file():
            logger.error(f"{name}: path is a directory, not an executable file: {p}")
            return False
    except OSError as err:
        logger.error(f"{name}: OS error accessing {p}: {err}")
        return False

    return True


def get_vina_version(vina_path: Path | str) -> str:
    """Run `vina --version` and return the version string.

    Returns 'NOT FOUND' or 'UNKNOWN' if the version cannot be determined.
    """
    if not validate_executable(vina_path, "AutoDock Vina"):
        return "NOT FOUND"
    return _run_version_command(Path(vina_path), ["--version"], r"(?:AutoDock\s+)?Vina\s+(v?[\d.]+\b.*)")


def get_vina_split_version(vina_split_path: Path | str) -> str:
    """Attempt to get vina_split version. May not support --version."""
    if not validate_executable(vina_split_path, "Vina Split"):
        return "NOT FOUND"
    ver = _run_version_command(Path(vina_split_path), ["--version"], r"(v?[\d.]+)")
    if ver == "UNKNOWN":
        try:
            result = subprocess.run(
                [str(vina_split_path), "--help"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0 or "vina_split" in (result.stdout + result.stderr).lower():
                return "available"
        except Exception:
            pass
    return ver


def get_autodock4_version(autodock4_path: Path | str) -> str:
    """Run `autodock4 --version` and return the version string."""
    if not validate_executable(autodock4_path, "AutoDock 4"):
        return "NOT FOUND"
    return _run_version_command(
        Path(autodock4_path), ["--version"],
        r"AutoDock\s+([\d.]+)"
    )


def get_autogrid4_version(autogrid4_path: Path | str) -> str:
    """Run `autogrid4 --version` and return the version string."""
    if not validate_executable(autogrid4_path, "AutoGrid 4"):
        return "NOT FOUND"
    return _run_version_command(
        Path(autogrid4_path), ["--version"],
        r"AutoGrid\s+([\d.]+)"
    )


def get_obabel_version(obabel_path: Path | str) -> str:
    """Run `obabel -V` and return the OpenBabel version string."""
    if not validate_executable(obabel_path, "OpenBabel"):
        return "NOT FOUND"
    return _run_version_command(
        Path(obabel_path), ["-V"],
        r"Open\s+Babel\s+([\d.]+)"
    )


def configure_openbabel_env(base_dir: Optional[Path] = None) -> None:
    """Ensure OpenBabel data and plugin directories are set in the environment.

    In frozen PyInstaller bundles and on various Windows systems, OpenBabel requires
    BABEL_LIBDIR (pointing to directory with *.obf plugins) and BABEL_DATADIR (pointing
    to directory with atomtyp.txt, types.txt, etc.) to convert formats.
    """
    import os
    search_dirs = get_candidate_search_dirs()
    if base_dir:
        search_dirs.insert(0, Path(base_dir))

    # Also include Python package installation paths
    try:
        import openbabel as _ob_mod
        ob_pkg = Path(_ob_mod.__file__).parent.resolve()
        search_dirs.extend([
            ob_pkg / "bin",
            ob_pkg / "lib" / "openbabel" / "3.2.1",
            ob_pkg,
        ])
    except Exception:
        pass

    found_lib = False
    found_data = False

    for d in search_dirs:
        if not d.is_dir():
            continue

        # Look for plugins (*.obf files)
        if not found_lib and any(d.glob("*.obf")):
            os.environ["BABEL_LIBDIR"] = str(d.resolve())
            found_lib = True
            logger.debug(f"Configured BABEL_LIBDIR={d}")

        # Look for data directory
        data_candidates = [
            d / "data",
            d.parent / "share" / "openbabel" / "3.2.1",
            d / "share" / "openbabel" / "3.2.1",
        ]
        for data_candidate in data_candidates:
            if not found_data and data_candidate.is_dir() and (
                (data_candidate / "space-groups.txt").exists()
                or (data_candidate / "types.txt").exists()
                or (data_candidate / "elements.txt").exists()
                or (data_candidate / "atomtyp.txt").exists()
            ):
                os.environ["BABEL_DATADIR"] = str(data_candidate.resolve())
                found_data = True
                logger.debug(f"Configured BABEL_DATADIR={data_candidate}")
                break

        # Add to DLL search path on Windows
        if sys.platform == "win32":
            cur_path = os.environ.get("PATH", "")
            resolved_str = str(d.resolve())
            if resolved_str not in cur_path:
                os.environ["PATH"] = resolved_str + os.pathsep + cur_path
            if hasattr(os, "add_dll_directory"):
                try:
                    os.add_dll_directory(resolved_str)
                except Exception:
                    pass

    # Suppress non-fatal OpenBabel warnings (e.g. kekulize warnings on complex alkaloids)
    try:
        from openbabel import openbabel as _ob_core
        _ob_core.obErrorLog.SetOutputLevel(0)
    except Exception:
        pass


# Auto-configure OpenBabel environment on import
configure_openbabel_env()


def _run_version_command(exe_path: Path, args: list, pattern: str) -> str:
    """Run an executable with arguments and extract version using regex.

    Args:
        exe_path: Path to the executable.
        args: Command-line arguments.
        pattern: Regex pattern to extract version from combined stdout+stderr.

    Returns:
        Extracted version string or 'UNKNOWN'.
    """
    if not exe_path or not exe_path.is_file():
        return "NOT FOUND"

    try:
        result = subprocess.run(
            [str(exe_path)] + args,
            capture_output=True, text=True, timeout=5,
            stdin=subprocess.DEVNULL,
        )
        combined = result.stdout + "\n" + result.stderr
        match = re.search(pattern, combined, re.IGNORECASE)
        if match:
            return match.group(1).strip()
    except subprocess.TimeoutExpired:
        logger.warning(f"Timeout running {exe_path} {' '.join(args)}")
    except FileNotFoundError:
        logger.error(f"Executable not found: {exe_path}")
    except Exception as e:
        logger.warning(f"Error running {exe_path}: {e}")
    return "UNKNOWN"


def get_cpu_count(cpu_setting: str = "AUTO") -> int:
    """Resolve CPU setting to a concrete number.

    Args:
        cpu_setting: "AUTO" or a numeric string.

    Returns:
        Number of CPUs to use.
    """
    import os
    if cpu_setting.upper() == "AUTO":
        return os.cpu_count() or 1
    try:
        n = int(cpu_setting)
        return max(1, n)
    except ValueError:
        logger.warning(f"Invalid CPU setting '{cpu_setting}', using AUTO")
        return os.cpu_count() or 1
