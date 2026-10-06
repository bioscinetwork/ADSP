# Changelog

All notable changes to AutoDock Suite Platform (ADSP) are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.3.0] - 2026-10-06

### Added
- **Repository Reorganization**: Restructured repository into clean, professional GitHub architecture:
  - Curated benchmarks migrated to `benchmarks/` (`1CA2`, `1MBN`, `2NV6`) with dedicated metadata.
  - Structural inputs organized into `data/` (`receptors/`, `ligands/`, `reference/`).
  - Scientific parameter profiles isolated in `parameter_profiles/` with cryptographic manifests.
  - Specialized AutoDock4Zn pipeline quarantined in `external/autodock4zn/`.
  - Comprehensive automated test suite relocated into `tests/` (124 tests across unit, integration, and regression).
  - Documentation and developer guides organized under `docs/`.
  - Build and launcher scripts organized under `scripts/`.
- **Packaging Test Suite**: Dedicated packaging regression tests (`tests/test_packaging.py`) verifying spec file portability, absence of machine-private paths, AD4 parameter bundling, and distribution asset copying.
- **Standard Project Files**: Added `LICENSE` (GPL-2.0-or-later) and `pyproject.toml` configuration with build system metadata and entry points.

### Changed
- **Thermodynamic Semantic Corrections**:
  - Eliminated synthetic $K_i$ calculations from AutoDock Vina binding scores in GUI and reporting. Vina results strictly report empirical affinity (kcal/mol).
  - AutoDock4 $K_i$ and free energy values are strictly DLG-driven from AutoDock4's internal analysis.
  - Replaced legacy calculation in GUI with explicit `NotImplementedError` to prevent regressions.
- **DLG Reporting & Null Safety**:
  - Hardened `dlg_extract.py` and `reporting.py` against `NoneType` formatting crashes when optional DLG metrics (reference RMSD, internal energy, $K_i$) are absent.
  - Ensured null RMSD values render as empty strings / N/A rather than silently defaulting to `0.0` Å.
- **PDBQT Fixed-Column Enforcement**: Enforced strict fixed-column formatting checks (atom names, coords, partial charges, atom types) with Open Babel round-trip validation.
- **Portable Resource Discovery**: Made `ad4_compatibility.py`, `AutoDockSuitePro.spec`, and `build_exe.py` use portable, relative path resolution across both development and frozen PyInstaller environments.
- **Relocation Safety**: Verified standalone packaged executable can be moved across directories and drives with zero dependence on developer source directories.

### Fixed
- Fixed keyword argument mismatch in `autodock4_workflow.py` (`receptor_path`, `ligand_path`, `config_path` in `ProvenanceRecord`).
- Fixed `NoneType` format crash in `dlg_extract.py:_txt_report` when positional RMSD was present but conformational RMSD was absent.
- Removed machine-private paths from root `project_config.toml`.
- Added defensive output directory initialization in standalone workflow runners.

### Known Limitations
- AutoDock4Zn requires Python 2 / MGLTools legacy environment and is kept as an isolated reference asset.
- ADSP Kollman-style partial charge mode is an automated fallback approximation, not a replacement for quantum-derived ESP charges.
- 2NV6 INH-NAD adduct benchmark isolates crystallographic extraction from independent SDF without claiming false redocking validation.
