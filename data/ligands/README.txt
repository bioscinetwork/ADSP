LIGANDS FOLDER — AutoDock Suite Pro
=====================================
Place all prepared ligand PDBQT files directly in this folder.
Naming: compound_001.pdbqt, compound_002.pdbqt, etc.

PREPARATION TOOLS (external — not part of this suite):
  • MGLTools / prepare_ligand4.py  — for ligand PDBQT preparation
  • OpenBabel / obabel              — alternative preparation tool

The platform will automatically dock ALL *.pdbqt files in this folder
against ALL receptors in the receptors/ folder.

To select specific ligands, edit 'selected_ligands' in project_config.toml:
  [inputs]
  selected_ligands = ["compound_001", "compound_005"]
