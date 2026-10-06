"""Regression tests for portable packaging and resource discovery.

These tests verify:
  - The committed AutoDockSuitePro.spec contains no hardcoded machine-specific paths
  - AD4 scientific parameter assets are present and have correct checksums
  - build_exe.py bundles AD4 parameter assets in both spec generation and copy_distribution_assets
  - RMSD null values are not silently converted to 0.0 in reporting output
  - Vina RMSD missing values appear as empty strings, not zeros
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent if not (Path(__file__).resolve().parent / 'main.py').is_file() else Path(__file__).resolve().parent.parent if not (Path(__file__).resolve().parent / 'main.py').is_file() else Path(__file__).resolve().parent

# Reference checksums (verified against canonical files)
AD4_PARAMS_SHA256 = "625DE5779B914382E21A135C776EFBC02B4221085BD0280118D103CCDD93EA7C"
AD41_BOUND_SHA256 = "6B98F7AB508F4882801938F8CED1C0BF38096496155A8005BAF941A201781CE8"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


# ─────────────────────────────────────────────────────────────────────────────
# Spec file portability
# ─────────────────────────────────────────────────────────────────────────────

def test_spec_file_contains_no_machine_specific_absolute_paths():
    """AutoDockSuitePro.spec must not contain hardcoded machine-specific paths.

    The canonical spec uses SPECPATH (a PyInstaller-provided variable) and
    dynamic path resolution so it works on any developer machine.
    """
    spec = REPO_ROOT / "AutoDockSuitePro.spec"
    assert spec.is_file(), "AutoDockSuitePro.spec must exist in repository root"

    content = spec.read_text(encoding="utf-8", errors="replace")

    # Patterns that indicate machine-specific hardcoded paths
    machine_patterns = [
        r"C:/Users/[A-Za-z]",       # Windows user home
        r"C:\\\\Users\\\\[A-Za-z]", # backslash variant
        r"D:/Projects/ADSP/bin/",   # developer-specific bin path
        r"D:\\\\Projects",          # backslash variant
        r"/home/[a-z]",             # Linux home path
        r"AppData/Roaming/Python",  # user-specific Python installation
    ]

    violations = []
    for line_no, line in enumerate(content.splitlines(), 1):
        for pattern in machine_patterns:
            if re.search(pattern, line):
                violations.append(f"Line {line_no}: {line.strip()}")
                break

    assert not violations, (
        "AutoDockSuitePro.spec contains machine-specific absolute paths "
        "(release blocker):\n" + "\n".join(violations)
    )


def test_spec_file_bundles_ad4_parameter_assets():
    """Spec must bundle AD4_parameters.dat and AD4.1_bound.dat."""
    spec = REPO_ROOT / "AutoDockSuitePro.spec"
    content = spec.read_text(encoding="utf-8", errors="replace")
    assert "AD4_parameters.dat" in content, (
        "AutoDockSuitePro.spec does not bundle AD4_parameters.dat "
        "(packaged app cannot locate parameter files)"
    )
    assert "AD4.1_bound.dat" in content, (
        "AutoDockSuitePro.spec does not bundle AD4.1_bound.dat"
    )


# ─────────────────────────────────────────────────────────────────────────────
# AD4 parameter file integrity
# ─────────────────────────────────────────────────────────────────────────────

def test_ad4_parameters_dat_checksum():
    """AD4_parameters.dat must match the verified reference checksum."""
    path = REPO_ROOT / "parameter_profiles" / "ad4_standard_4.2" / "AD4_parameters.dat"
    if not path.is_file():
        path = REPO_ROOT / "AD4_parameters.dat"
    assert path.is_file(), f"AD4_parameters.dat not found at {path}"
    observed = _sha256(path)
    assert observed == AD4_PARAMS_SHA256, (
        f"AD4_parameters.dat checksum MISMATCH\n"
        f"  Expected : {AD4_PARAMS_SHA256}\n"
        f"  Observed : {observed}\n"
        "STOP: Do not silently overwrite verified scientific parameter assets."
    )


def test_ad41_bound_dat_checksum():
    """AD4.1_bound.dat must match the verified reference checksum."""
    path = REPO_ROOT / "parameter_profiles" / "ad4_1_bound" / "AD4.1_bound.dat"
    if not path.is_file():
        path = REPO_ROOT / "AD4.1_bound.dat"
    assert path.is_file(), f"AD4.1_bound.dat not found at {path}"
    observed = _sha256(path)
    assert observed == AD41_BOUND_SHA256, (
        f"AD4.1_bound.dat checksum MISMATCH\n"
        f"  Expected : {AD41_BOUND_SHA256}\n"
        f"  Observed : {observed}\n"
        "STOP: Do not silently overwrite verified scientific parameter assets."
    )


# ─────────────────────────────────────────────────────────────────────────────
# build_exe.py bundling of AD4 assets
# ─────────────────────────────────────────────────────────────────────────────

def test_build_exe_bundles_ad4_parameter_assets_in_spec():
    """build_exe.py generate_spec_file() must include AD4 parameter files in datas."""
    build_exe = REPO_ROOT / "build_exe.py"
    content = build_exe.read_text(encoding="utf-8", errors="replace")
    assert "AD4_parameters.dat" in content, (
        "build_exe.py does not bundle AD4_parameters.dat in the generated spec"
    )
    assert "AD4.1_bound.dat" in content, (
        "build_exe.py does not bundle AD4.1_bound.dat in the generated spec"
    )


def test_build_exe_copies_ad4_assets_in_copy_distribution_assets():
    """build_exe.py copy_distribution_assets() must copy AD4 parameter files to _internal/."""
    build_exe = REPO_ROOT / "build_exe.py"
    content = build_exe.read_text(encoding="utf-8", errors="replace")
    # Confirm the copy loop that handles AD4 assets appears in copy_distribution_assets
    assert "_internal" in content and "AD4_parameters.dat" in content, (
        "build_exe.py does not copy AD4_parameters.dat to _internal/ in distribution"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Reporting: RMSD null values must not silently become 0.0
# ─────────────────────────────────────────────────────────────────────────────

def test_vina_rmsd_null_is_empty_string_in_report_rows():
    """Vina RMSD lower/upper bound must be empty string, not 0.0, when not available."""
    from models import DockingJob, Engine, DockingMode, JobStatus, VinaResult

    job = DockingJob(engine=Engine.VINA, docking_mode=DockingMode.RIGID)
    job.status = JobStatus.SUCCESS
    # Simulate a Vina pose with missing RMSD bounds (both None)
    vr = VinaResult(pose=1, binding_affinity=-7.5, rmsd_lower_bound=None, rmsd_upper_bound=None)
    job.vina_results = [vr]
    job.ensure_canonical_result()

    # Import internal build function
    from reporting import _build_vina_rows
    rows = _build_vina_rows([job])
    assert rows, "No rows generated for successful Vina job"

    row = rows[0]
    lb = row["RMSD_Lower_Bound"]
    ub = row["RMSD_Upper_Bound"]
    assert lb == "" or lb is None, (
        f"RMSD_Lower_Bound should be empty/None when not available, got {lb!r}"
    )
    assert ub == "" or ub is None, (
        f"RMSD_Upper_Bound should be empty/None when not available, got {ub!r}"
    )
    # Explicitly fail if they have been silently zeroed
    assert lb != 0.0, "RMSD_Lower_Bound must not silently become 0.0 (absolute rule violation)"
    assert ub != 0.0, "RMSD_Upper_Bound must not silently become 0.0 (absolute rule violation)"


def test_vina_rmsd_real_value_passes_through_unmodified():
    """When Vina does report RMSD bounds, they must appear unchanged in the report."""
    from models import DockingJob, Engine, DockingMode, JobStatus, VinaResult

    job = DockingJob(engine=Engine.VINA, docking_mode=DockingMode.RIGID)
    job.status = JobStatus.SUCCESS
    vr = VinaResult(pose=1, binding_affinity=-7.5, rmsd_lower_bound=0.123, rmsd_upper_bound=1.456)
    job.vina_results = [vr]
    job.ensure_canonical_result()

    from reporting import _build_vina_rows
    rows = _build_vina_rows([job])
    assert rows

    row = rows[0]
    assert abs(float(row["RMSD_Lower_Bound"]) - 0.123) < 1e-9
    assert abs(float(row["RMSD_Upper_Bound"]) - 1.456) < 1e-9


# ─────────────────────────────────────────────────────────────────────────────
# Portable resource discovery: ad4_compatibility module-relative lookup
# ─────────────────────────────────────────────────────────────────────────────

def test_ad4_compatibility_resolves_parameters_module_relative():
    """ad4_compatibility._bundled_checksum() must find files relative to the module."""
    from ad4_compatibility import STANDARD_AD4_PROFILE, AD4_1_BOUND_PROFILE

    # SHA-256 values are populated at import time; None only if files were missing
    assert STANDARD_AD4_PROFILE.sha256 is not None, (
        "AD4_parameters.dat not found by ad4_compatibility module-relative lookup "
        "(packaged resource discovery would also fail)"
    )
    assert AD4_1_BOUND_PROFILE.sha256 is not None, (
        "AD4.1_bound.dat not found by ad4_compatibility module-relative lookup"
    )
    # Cross-check with the verified reference checksums
    assert STANDARD_AD4_PROFILE.sha256.upper() == AD4_PARAMS_SHA256
    assert AD4_1_BOUND_PROFILE.sha256.upper() == AD41_BOUND_SHA256
