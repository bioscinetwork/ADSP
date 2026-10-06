# ADSP Automated Test Suite

This directory contains the automated test suites for AutoDock Suite Pro v0.3.0.

## Test Inventory

- `test_2nv6_benchmark.py`: Verification of 2NV6 / INH-NAD benchmark inputs and representation isolation.
- `test_ad4_compatibility.py`: AutoDock4 parameter profile parsing, metal typing, and CIF preparation tests.
- `test_canonical_upgrade.py`: Data model serialization and canonical result schema upgrade tests.
- `test_charge_diagnostics.py`: Formal and partial charge sanity checks and non-finite charge rejection.
- `test_cleanup_provenance.py`: Receptor and ligand cleanup provenance tracking and file safety.
- `test_dlg_reporting_nulls.py`: DLG parsing and null value formatting safety in reporting.
- `test_interaction_synchronization.py`: Pose selection synchronization across diagram, complex, and metrics.
- `test_level3_upgrade.py`: Workflow execution, engine selection, and output structure verification.
- `test_packaging.py`: Portable PyInstaller spec verification, parameter bundling, and relocation safety.
- `test_pdbqt_validation.py`: Fixed-column PDBQT syntax validation and Open Babel format compatibility.
- `test_results_semantics.py`: Engine-isolated thermodynamic semantics (no Vina thermodynamic inference).
- `test_rmsd_validation.py`: Hungarian automorphism and MCS graph-constrained RMSD calculations.
- `test_suite.py`: Comprehensive end-to-end and unit test suite across all ADSP subsystems.

## Running Tests

Execute with pytest from the project root:

```bash
pytest
```
