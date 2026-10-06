# AutoDock4Zn Specialized Pipeline (External Asset)

This directory isolates the specialized zinc coordination docking pipeline from ordinary AutoDock4 parameter profiles.

- **Archive**: `AutoDock4Zn-Pipeline-main.zip`
- **License**: GPLv3 (see archive member LICENSE)
- **Status**: External reference pipeline and coordination parameter asset (`AD4Zn.dat` with `TZ` pseudo-atom parameters).
- **Isolation Policy**:
  - AutoDock4Zn requires specialized pseudo-atom tetrahedral zinc coordination preparation (`prepare_gpf4zn.py`).
  - It is strictly quarantined from standard AutoGrid4 / AutoDock4 jobs in ADSP to prevent unvalidated coordination assumptions from affecting general metalloprotein docking.
