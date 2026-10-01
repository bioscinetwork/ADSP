RECEPTOR FOLDER — AutoDock Suite Pro
======================================
Place one subfolder per receptor here. Each subfolder must follow this layout:

  RECEPTOR_NAME/
    rigid/
      receptor.pdbqt     ← Your rigid receptor PDBQT (prepared externally)
    flex/                ← OPTIONAL: only needed for flexible docking
      receptor_rigid.pdbqt  ← Rigid backbone PDBQT  (passed to --receptor)
      receptor_flex.pdbqt   ← Flexible sidechain PDBQT (passed to --flex)
    config.txt           ← Vina grid config (run --prepare to auto-generate)

PREPARATION TOOLS (external — not part of this suite):
  • MGLTools / prepare_receptor4.py  — for receptor PDBQT preparation
  • MGLTools / prepare_flexreceptor4.py — for rigid/flex split
  • OpenBabel / obabel                — alternative preparation tool

AUTO-GENERATE CONFIG.TXT:
  Run: python main.py --prepare
  This will prompt for grid center (x, y, z) and box size for each receptor
  that is missing a config.txt.

FLEXIBLE DOCKING:
  The platform auto-detects FLEXIBLE mode when a 'flex/' subfolder exists
  and contains at least one PDBQT file.
