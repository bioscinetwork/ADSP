"""
AutoDock Suite Pro — 2D Ligand Structure Depiction (depiction.py)
=================================================================
Renders a 2D structural diagram of a ligand from its PDBQT file, SDF, or SMILES.

Strategy (tried in order of reliability):
  1. In-memory render cache (instant return)
  2. Direct RDKit load if path is SDF/MOL/MOL2 or sibling exists in ligands/
  3. RDKit: parse SMILES from PDBQT REMARK line (fastest, most accurate)
  4. RDKit: convert via OpenBabel PDBQT → SDF → RDKit mol
  5. RDKit: parse directly from PDBQT atom coordinates (3D → 2D cleanup)
  6. Return None (caller shows placeholder text)

Handles aromatic alkaloids (e.g. Berberine, Palmatine) gracefully with
non-kekulized rendering fallback.

Author: AutoDock Suite Pro Development Team
Version: 0.3.0
"""

from __future__ import annotations

import logging
import re
import subprocess
from io import BytesIO
from pathlib import Path
from typing import Dict, Optional, Tuple

logger = logging.getLogger("docking_automation.depiction")

# In-memory PNG cache to prevent duplicate depiction renders
_RENDER_CACHE: Dict[str, bytes] = {}


def render_ligand_2d(
    pdbqt_path: Path | str,
    size: Tuple[int, int] = (300, 220),
    obabel_exe: Optional[Path | str] = None,
) -> Optional[bytes]:
    """Render a 2D structure diagram for a ligand file (PDBQT, SDF, MOL2).

    Args:
        pdbqt_path: Path to the ligand file.
        size:       (width, height) of output PNG in pixels.
        obabel_exe: Path to obabel executable (for conversion fallback).

    Returns:
        PNG image as raw bytes, or None if all methods fail.
    """
    path = Path(pdbqt_path)
    if not path.is_file():
        return None

    cache_key = f"{path.resolve()}_{size[0]}x{size[1]}_{path.stat().st_mtime}"
    if cache_key in _RENDER_CACHE:
        return _RENDER_CACHE[cache_key]

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.ERROR)
    except ImportError:
        logger.warning("RDKit not available — 2D depiction disabled.")
        return None

    mol = None

    # ── Strategy 1: Direct load from SDF/MOL2/SMILES or sibling reference ───
    mol = _mol_from_file_or_sibling(path)

    # ── Strategy 2: SMILES from PDBQT REMARK line ────────────────────────────
    if mol is None and path.suffix.lower() == ".pdbqt":
        mol = _mol_from_remark_smiles(path)

    # ── Strategy 3: OpenBabel PDBQT → SDF → RDKit ────────────────────────────
    if mol is None and obabel_exe:
        mol = _mol_via_obabel(path, obabel_exe)

    # ── Strategy 4: RDKit PDBQT direct parse (3D coordinates) ────────────────
    if mol is None:
        mol = _mol_from_pdbqt_direct(path)

    if mol is None:
        logger.debug(f"All depiction strategies failed for: {path.name}")
        return None

    # Generate 2D coordinates & Render PNG
    try:
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            pass

        try:
            # Try kekulized first
            img = Draw.MolToImage(mol, size=size, kekulize=True)
        except Exception:
            # Fallback for charged alkaloids and tricky aromatic systems (e.g. Berberine)
            img = Draw.MolToImage(mol, size=size, kekulize=False)

        buf = BytesIO()
        img.save(buf, format="PNG")
        raw_bytes = buf.getvalue()
        _RENDER_CACHE[cache_key] = raw_bytes
        return raw_bytes
    except Exception as e:
        logger.debug(f"RDKit Draw failed for {path.name}: {e}")
        return None


# ---------------------------------------------------------------------------
# Strategy Helpers
# ---------------------------------------------------------------------------

def _mol_from_file_or_sibling(path: Path):
    """Load molecule from original format or search for sibling SDF/MOL2."""
    try:
        from rdkit import Chem
        ext = path.suffix.lower()
        if ext in (".sdf", ".mol"):
            suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True)
            for m in suppl:
                if m:
                    return m
        elif ext == ".mol2":
            return Chem.MolFromMol2File(str(path), sanitize=True, removeHs=True)
        elif ext in (".smi", ".smiles"):
            txt = path.read_text(encoding="utf-8", errors="replace").strip()
            smi = txt.split()[0] if txt else ""
            if smi:
                return Chem.MolFromSmiles(smi)

        # If PDBQT, search for sibling SDF/MOL2 with the same stem
        stem = path.stem
        # Clean off common docking suffixes
        for suffix in ("_out", "_pose_01", "_run1", "_run01"):
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break

        dirs_to_check = [path.parent]
        if path.parent.name in ("DLG", "results", "BSNDVP_RESULTS", "EXTRACTED_POSES", "rigid", "flex"):
            dirs_to_check.append(path.parent.parent / "ligands")
            dirs_to_check.append(path.parent.parent)

        for d in dirs_to_check:
            if not d.is_dir():
                continue
            for cand_ext in (".sdf", ".mol2", ".mol", ".smi"):
                cand = d / f"{stem}{cand_ext}"
                if cand.is_file():
                    res = _mol_from_file_or_sibling(cand)
                    if res:
                        return res
    except Exception:
        pass
    return None


def _mol_from_remark_smiles(pdbqt_path: Path):
    """Extract SMILES from a REMARK line in the PDBQT (if present)."""
    try:
        from rdkit import Chem
        with open(pdbqt_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith("REMARK") and "SMILES" in stripped.upper():
                    m = re.search(r"SMILES\s*[:=]?\s*(\S+)", stripped, re.IGNORECASE)
                    if m:
                        smiles = m.group(1).strip()
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            return mol
    except Exception:
        pass
    return None


def _mol_via_obabel(pdbqt_path: Path, obabel_exe: Path | str):
    """Convert PDBQT → SDF via obabel subprocess, then parse with RDKit."""
    try:
        from rdkit import Chem
        result = subprocess.run(
            [str(obabel_exe), str(pdbqt_path), "-osdf", "--gen2D"],
            capture_output=True, text=True, timeout=4,
        )
        if result.returncode == 0 and result.stdout.strip():
            mol = Chem.MolFromMolBlock(result.stdout, sanitize=True, removeHs=True)
            if mol:
                return mol
    except Exception:
        pass
    return None


def _mol_from_pdbqt_direct(pdbqt_path: Path):
    """Parse atoms from PDBQT and build an RDKit mol with 3D→2D cleanup."""
    try:
        from rdkit import Chem

        # Build a minimal PDB block from ATOM/HETATM records
        pdb_lines = ["COMPND    LIGAND\n"]
        serial = 1
        with open(pdbqt_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("DOCKED:"):
                    line = line[7:].lstrip()
                if line.startswith(("ATOM", "HETATM")):
                    # Convert PDBQT ATOM type back to element symbol
                    atom_type = line[76:].strip().split()[0] if len(line) > 76 else ""
                    element = _ad_type_to_element(atom_type, line[12:16].strip())
                    # Write a clean PDB HETATM line
                    pdb_line = (
                        f"HETATM{serial:5d}  {line[12:16].strip():<3s} LIG A   1    "
                        f"{line[30:38]}{line[38:46]}{line[46:54]}"
                        f"  1.00  0.00          {element:>2s}\n"
                    )
                    pdb_lines.append(pdb_line)
                    serial += 1
        pdb_lines.append("END\n")
        pdb_block = "".join(pdb_lines)

        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        return mol
    except Exception:
        return None


def _ad_type_to_element(ad_type: str, atom_name: str) -> str:
    """Map AutoDock atom type to element symbol."""
    t = ad_type.strip().upper()
    mapping = {
        "C": "C", "A": "C", "N": "N", "NA": "N", "OA": "O", "O": "O",
        "S": "S", "SA": "S", "HD": "H", "H": "H", "HS": "H",
        "CL": "Cl", "BR": "Br", "F": "F", "I": "I", "P": "P",
        "FE": "Fe", "ZN": "Zn", "MG": "Mg", "CA": "Ca",
    }
    if t in mapping:
        return mapping[t]
    # Fallback: use first letter of atom name
    for ch in atom_name.strip():
        if ch.isalpha():
            return ch.upper()
    return "C"
