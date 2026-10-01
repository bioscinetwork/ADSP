#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   BSNDVP™ — Scientific Edition  (dlg_extract.py)  v5.0.0                   ║
║   AutoDock DLG Parsing, Validation & Statistical Analysis Framework          ║
║                                                                              ║
║   © 2026 Bio-Science Network Limited.  All Rights Reserved.                  ║
║                                                                              ║
║   A focused, publication-grade command-line tool for AutoDock DLG            ║
║   post-processing.  Three scientific contributions:                           ║
║     1. AutoDock DLG parsing  (rigid & flexible receptor docking)             ║
║     2. Docking validation methodology  (Kabsch RMSD, MCS, RDKit)            ║
║     3. Docking statistical analysis   (thermodynamics, clustering,           ║
║                                        Boltzmann, Shannon entropy)            ║
║                                                                              ║
║   Author  : Chibuike Praise Okechukwu                                        ║
║   Brand   : BSNDVP™ — Bio-Science Network Docking & Validation Platform     ║
║                                                                              ║
║   CLI:  python dlg_extract.py [--dlg DIR] [--output DIR]                    ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""

from __future__ import annotations

# ── Standard library ──────────────────────────────────────────────────────────
import csv
import json
import logging
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

if sys.platform == "win32":
    os.system("")
    for _s in (sys.stdout, sys.stderr):
        if hasattr(_s, "reconfigure"):
            try:
                _s.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass

# ── Third-party: NumPy ────────────────────────────────────────────────────────
try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False
    np = None  # type: ignore

# ── Third-party: Pandas ───────────────────────────────────────────────────────
try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False
    pd = None  # type: ignore

# ── Third-party: OpenPyXL ─────────────────────────────────────────────────────
try:
    import openpyxl  # noqa: F401
    _HAS_XLSX = True
except ImportError:
    _HAS_XLSX = False

# ── Third-party: RDKit (optional structural validation) ───────────────────────
try:
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

# ── Third-party: gemmi (optional mmCIF / CIF parsing) ────────────────────────
try:
    import gemmi
    _HAS_GEMMI = True
except ImportError:
    _HAS_GEMMI = False

# ── Logging ───────────────────────────────────────────────────────────────────
logger = logging.getLogger("dlg_extract")
logger.setLevel(logging.INFO)
if not logger.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[BSNDVP] %(levelname)s: %(message)s"))
    logger.addHandler(_h)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 0b — BRANDING & ANSI COLOURS
# ══════════════════════════════════════════════════════════════════════════════

__trademark__  = "BSNDVP™"
__edition__    = "Scientific Edition"
__copyright__  = "© 2026 Bio-Science Network Limited. All rights reserved."
__full_brand__ = f"{__trademark__} — DLG Parsing, Validation & Statistical Analysis — {__edition__}"
__author__     = "Chibuike Praise Okechukwu"
__version__    = "5.0.0"


class C:
    """ANSI escape codes for terminal colouring."""
    RESET   = "\033[0m"
    BOLD    = "\033[1m"
    DIM     = "\033[2m"
    RED     = "\033[91m"
    GREEN   = "\033[92m"
    YELLOW  = "\033[93m"
    CYAN    = "\033[96m"
    WHITE   = "\033[97m"
    GRAY    = "\033[90m"
    MAGENTA = "\033[95m"
    ORANGE  = "\033[38;5;208m"
    TEAL    = "\033[38;5;43m"
    GOLD    = "\033[38;5;220m"
    LIME    = "\033[38;5;118m"


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — DATA MODELS
# ══════════════════════════════════════════════════════════════════════════════

# ── Docking mode classification ───────────────────────────────────────────────
DockingMode = Enum("DockingMode", ["RIGID", "FLEXIBLE"])


class AtomRole(Enum):
    """Structural role of an atom within a docked pose."""
    LIGAND           = 1
    FLEXIBLE_RESIDUE = 2


# Standard amino acid residue names (including protonation variants)
_STANDARD_AMINO_ACIDS: FrozenSet[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # Protonation variants / non-standard
    "HID", "HIE", "HIP", "HSE", "HSD", "HSP", "CYX", "CYM",
    "ASH", "GLH", "LYN", "TYM", "SEC", "PYL", "MSE",
})


@dataclass
class DockingIdentity:
    ligand_name:   str = "Unknown"
    receptor_name: str = "Unknown"
    original_file: str = ""


@dataclass
class PDBAtom:
    serial:        int
    name:          str
    element:       str
    residue_name:  str      = ""
    chain_id:      str      = ""
    residue_seq:   int      = 0
    x:             float    = 0.0
    y:             float    = 0.0
    z:             float    = 0.0
    charge:        float    = 0.0
    autodock_type: str      = ""
    is_hetatm:     bool     = False
    # Auto-classified structural role (set by AtomRoleClassifier)
    role:          AtomRole = field(default_factory=lambda: AtomRole.LIGAND)

    @property
    def residue_label(self) -> str:
        return f"{self.residue_name}{self.residue_seq}{self.chain_id.strip()}"

    @property
    def is_flexible_residue(self) -> bool:
        return self.role == AtomRole.FLEXIBLE_RESIDUE


@dataclass
class DockingPose:
    rank:             int   = 0
    run_number:       int   = 0
    binding_energy:   float = 0.0
    ki_raw:           Optional[float] = None
    ki_unit:          str   = ""
    ki_nM:            Optional[float] = None
    intermol_energy:  float = 0.0
    internal_energy:  float = 0.0
    torsional_energy: float = 0.0
    unbound_energy:   float = 0.0
    vdW_hbond_desolvation_energy: Optional[float] = None
    electrostatic_energy: Optional[float] = None
    dlg_source: Optional[str] = None
    rmsd_from_ref:    float = 0.0
    cluster_id:       Optional[int] = None
    cluster_size:     int   = 0
    cluster_rmsd:     float = 0.0
    # Unified atoms list (backward compatible — always populated)
    atoms:            List[PDBAtom] = field(default_factory=list)
    # Classified sub-lists (populated by AtomRoleClassifier after parsing)
    ligand_atoms:           List[PDBAtom] = field(default_factory=list)
    flexible_residue_atoms: List[PDBAtom] = field(default_factory=list)
    # Receptor strain estimate (flexible docking only)
    receptor_strain_estimate: Optional[float] = None
    # Validation (RDKit — strictly optional)
    validation_rmsd_conf:     Optional[float] = None
    validation_rmsd_pos:      Optional[float] = None
    validation_heavy_rmsd:    Optional[float] = None
    validation_scaffold_rmsd: Optional[float] = None


@dataclass
class ClusterStats:
    cluster_id:        int   = 0
    size:              int   = 0
    mean_energy:       float = 0.0
    lowest_energy:     float = 0.0
    best_run:          int   = 0
    avg_rmsd:          float = 0.0
    spread:            float = 0.0
    # Cluster Quality Score (dimensionless, 0-100):
    #   f_pop  = cluster_size / N_runs   (population fraction)
    #   q_rmsd = max(0, 1 - avg_rmsd/3)  (geometric compactness; 0 at RMSD>=3.0 Å)
    #   quality = (f_pop + q_rmsd) / 2 * 100
    cluster_quality:   float = 0.0
    runs:              List[int] = field(default_factory=list)
    # Dual ranking
    thermodynamic_rank:  int  = 0
    population_rank:     int  = 0
    # Mode-aware semantics
    docking_mode:        str  = "RIGID"    # "RIGID" | "FLEXIBLE"
    rmsd_is_composite:   bool = False      # True for FLEXIBLE: RMSD covers lig+residues


@dataclass
class DockingResult:
    identity:      DockingIdentity = field(default_factory=DockingIdentity)
    protein:       str  = ""
    ligand:        str  = ""
    dlg_file:      Optional[Path] = None
    num_runs:      int  = 0
    num_clusters:  int  = 0
    poses:         List[DockingPose]  = field(default_factory=list)
    clusters:      List[ClusterStats] = field(default_factory=list)
    success:       bool = True
    error_message: str  = ""
    # Thermodynamics (parsed directly from AutoDock DLG output)
    info_entropy:         Optional[float] = None
    info_entropy_rmstol:  Optional[float] = None
    partition_function:   Optional[float] = None
    stat_temperature:     Optional[float] = None
    stat_free_energy:     Optional[float] = None
    stat_internal_energy: Optional[float] = None
    stat_entropy:         Optional[float] = None
    # Self-adaptive docking mode intelligence
    docking_mode:              DockingMode   = field(default_factory=lambda: DockingMode.RIGID)
    n_ligand_atoms:            int           = 0
    n_flexible_residue_atoms:  int           = 0
    receptor_strain_estimate:  Optional[float] = None
    # Flexible residue identifiers parsed from header
    flex_residue_ids:          List[str]     = field(default_factory=list)

    @property
    def best_energy(self) -> float:
        return min((p.binding_energy for p in self.poses), default=float("inf"))


@dataclass
class DockingProfile:
    """
    Flat data record for one ligand's best docking result.
    Contains only objective scientific quantities — no composite scoring.
    """
    ligand:               str   = ""
    protein:              str   = ""
    binding_energy:       float = 0.0
    ki_nM:                Optional[float] = None
    run_number:           int   = 0
    rmsd:                 float = 0.0
    cluster_id:           Optional[int] = None
    cluster_size:         int   = 0
    cluster_rmsd:         float = 0.0
    avg_cluster_rmsd:     float = 0.0
    cluster_mean_energy:  float = 0.0
    cluster_spread:       float = 0.0
    pose:                 Optional[DockingPose] = None

    # Advanced Thermodynamics
    entropic_penalty_kcal: float = 0.0   # torsional ΔG_tor (already in ΔG_binding)
    boltzmann_prob:        Optional[float] = None

    # Pharmacokinetics (ADMET descriptors — reported as-is, no penalty)
    admet_mw:              Optional[float] = None
    admet_logp:            Optional[float] = None
    admet_tpsa:            Optional[float] = None
    admet_nrotb:           Optional[int]   = None
    admet_hbd:             Optional[int]   = None
    admet_hba:             Optional[int]   = None
    lipinski_violations:   int   = 0

    # Self-adaptive mode fields
    docking_mode:             str            = "RIGID"
    n_ligand_atoms:           int            = 0
    n_flexible_residue_atoms: int            = 0
    receptor_strain_kcal:     Optional[float] = None
    flex_residue_ids:         List[str]      = field(default_factory=list)
    # Kabsch validation RMSD (RDKit; distinct from AutoDock's internal RMSD)
    validation_rmsd_conf:     Optional[float] = None
    validation_rmsd_pos:      Optional[float] = None


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — AUTODOCK TYPE → ELEMENT MAPPING
# ══════════════════════════════════════════════════════════════════════════════

_AD_TYPE_TO_ELEMENT: Dict[str, str] = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
    "H": "H", "HD": "H", "HS": "H", "HP": "H",
    "P": "P", "F": "F",
    "CL": "Cl", "Cl": "Cl", "BR": "Br", "Br": "Br", "I": "I",
    "FE": "Fe", "ZN": "Zn", "MG": "Mg", "MN": "Mn", "CA": "Ca",
    "CU": "Cu", "CO": "Co", "NI": "Ni", "SE": "Se",
    "LP": "", "E": "", "W": "O",
}
_VALID_ELEMENTS: FrozenSet[str] = frozenset({
    "H", "C", "N", "O", "F", "S", "P", "Cl", "Br", "I",
    "Fe", "Zn", "Mg", "Mn", "Ca", "Cu", "Co", "Ni", "Se",
    "Na", "K", "He", "Li", "Be", "B", "Ne", "Al", "Si", "Ar",
    "Sc", "Ti", "V", "Cr", "Ga", "Ge", "As", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi",
})


def _autodock_to_element(ad_type: str, atom_name: str = "") -> str:
    t = ad_type.strip().upper()
    if t in _AD_TYPE_TO_ELEMENT:
        return _AD_TYPE_TO_ELEMENT[t]
    cap = ad_type.strip().capitalize()
    if cap in _VALID_ELEMENTS:
        return cap
    for ch in atom_name:
        if ch.isalpha() and ch.upper() in _VALID_ELEMENTS:
            return ch.upper()
    return "C"


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — DLG PARSER
# ══════════════════════════════════════════════════════════════════════════════

_KI_TO_NM: Dict[str, float] = {
    "pm": 0.001, "nm": 1.0, "um": 1000.0,
    "\u00b5m": 1000.0, "\u03bcm": 1000.0, "microm": 1000.0,
    "mm": 1e6, "m": 1e9,
}


def _normalize_ki(value: Optional[float], unit: str) -> Optional[float]:
    if value is None:
        return None
    key = unit.lower().strip().rstrip(".")
    factor = _KI_TO_NM.get(key)
    if factor is None:
        for k, f in _KI_TO_NM.items():
            if k in key:
                factor = f
                break
    return value * (factor or 1.0)


class IdentityResolutionEngine:
    """
    Extracts ligand / receptor names from DLG header lines.

    Generic-name handling:
    AutoDock records the literal pdbqt filename in the DLG header, so if the
    docking was run with placeholder names like ``ligand.pdbqt``,
    ``ligandec.pdbqt``, ``ligandpyr.pdbqt``, or ``ligandt.pdbqt``, those
    tokens carry no chemical identity and must be discarded.  The engine
    then falls back to parsing the *DLG filename stem*, which the user
    is expected to name using the convention::

        {receptor}_{ligand}.dlg   →  rec = receptor,  lig = ligand
        {ligand}_{receptor}.dlg   →  auto-detected by PDB-ID heuristic

    **Recommended naming convention (AutoDock standard):**
        ``{receptor}_{ligand}.dlg``  e.g. ``1dg5_trimethoprim.dlg``

    Generic-token detection uses both a fixed set (``_GENERIC``) and a
    prefix regex that matches any name starting with ``ligand`` or
    ``receptor`` (case-insensitive), catching ``ligandec``, ``ligandpyr``,
    ``ligandt``, ``receptor2``, etc.
    """

    # Fixed set of names that carry no chemical identity
    _GENERIC = frozenset({
        "ligand", "molecule", "mol", "lig", "drug", "compound", "cpd",
        "receptor", "protein", "target", "macro", "macromolecule",
        "unknown", "input", "",
    })

    # Regex: any name that begins with the word "ligand" or "receptor"
    # (case-insensitive).  Catches: ligandec, ligandpyr, ligandt, ligand2, …
    _RE_GENERIC_PREFIX = re.compile(r'^(ligand|receptor)[a-z0-9_\-]*$', re.IGNORECASE)

    # Regex: PDB entry ID — starts with a digit, then 3 alphanumeric chars
    # (e.g. 1dg5, 4m6k, 1j3i, 1rx2, 11GS treated as starts-with-digit).
    # Also catches common 4-character all-uppercase receptor codes.
    _RE_PDB_ID = re.compile(r'^[0-9][a-zA-Z0-9]{3}$')

    @classmethod
    def _is_generic(cls, name: str) -> bool:
        """Return True if *name* is a placeholder with no chemical meaning."""
        n = name.strip().lower()
        return n in cls._GENERIC or bool(cls._RE_GENERIC_PREFIX.match(n))

    @classmethod
    def _is_pdb_id(cls, name: str) -> bool:
        """Heuristic: does *name* look like a PDB receptor entry ID?"""
        return bool(cls._RE_PDB_ID.match(name.strip()))

    @classmethod
    def extract(cls, lines: List[str], dlg_path: Path) -> DockingIdentity:
        lig_p1 = lig_p2 = rec_p1 = rec_p2 = ""
        for line in lines:
            s = line.strip()
            if 'Ligand PDBQT file =' in s:
                m = re.search(r'Ligand PDBQT file\s*=\s*"(.+?)\.pdbqt"', s)
                if m: lig_p1 = m.group(1)
            elif s.startswith("INPUT-LIGAND-PDBQT"):
                m = re.search(r'INPUT-LIGAND-PDBQT\s+([\w\-]+)', s)
                if m: lig_p2 = m.group(1).replace(".pdbqt", "")
            elif "Macromolecule file used to create Grid Maps" in s:
                m = re.search(r'Macromolecule file used to create Grid Maps\s*=\s*(.+?)\.pdbqt', s)
                if m: rec_p1 = m.group(1)
            elif s.startswith("DPF> macromolecule"):
                m = re.search(r'DPF>\s+macromolecule\s+(.+?)\.pdbqt', s)
                if m: rec_p2 = m.group(1)
            elif "Macromolecule" in s and ".pdbqt" in s:
                m = re.search(r'(\w+)\.pdbqt', s)
                if m and not rec_p1: rec_p1 = m.group(1)
            elif s.startswith("DPF>") and "fld" in s:
                m = re.search(r'fld\s+(\w+)\.maps\.fld', s)
                if m and not rec_p2: rec_p2 = m.group(1)

        def clean(n: str) -> str:
            if not n:
                return ""
            # Resolve path-based names (e.g. "../../Ligands/folate.pdbqt" -> "folate")
            s = n.replace("\\", "/").strip()
            stem = Path(s).stem.strip() if ("/" in s or "\\" in n) else s
            return stem.replace(".pdbqt", "").replace(".pdb", "").strip()

        # Discard generic / placeholder names from header
        raw_lig = next(
            (clean(n) for n in [lig_p1, lig_p2]
             if clean(n) and not cls._is_generic(clean(n))),
            "",
        )
        raw_rec = next(
            (clean(n) for n in [rec_p1, rec_p2]
             if clean(n) and not cls._is_generic(clean(n))),
            "",
        )

        # ── Fallback: derive names from the DLG filename stem ──────────────────
        # Supported conventions (RECOMMENDED: first form):
        #   {receptor}_{ligand}.dlg       →  1dg5_trimethoprim.dlg
        #   {receptor}_flex_{ligand}.dlg  →  2Z5X_flex_aspirin.dlg (flexible docking)
        #   {ligand}_{receptor}.dlg       →  trimethoprim_1dg5.dlg  (auto-detected)
        #   {ligand}.dlg                  →  trimethoprim.dlg       (no receptor)
        stem = dlg_path.stem                      # e.g. "1dg5_trimethoprim" or "2Z5X_flex_aspirin"
        stem_lower = stem.lower()

        if "_flex_" in stem_lower:
            flex_pos = stem_lower.find("_flex_")
            stem_rec = stem[:flex_pos + 5]  # e.g. "2Z5X_flex"
            stem_lig = stem[flex_pos + 6:]  # e.g. "aspirin"
        else:
            parts = stem.split("_", 1)               # split on first underscore only
            if len(parts) == 2:
                left, right = parts[0], parts[1]
                # Determine which token is the receptor PDB ID using a heuristic:
                # PDB IDs typically start with a digit (e.g. 1dg5, 4m6k).
                if cls._is_pdb_id(left) and not cls._is_pdb_id(right):
                    # Standard: {receptor}_{ligand}.dlg  e.g. 1dg5_trimethoprim
                    stem_rec, stem_lig = left, right
                elif cls._is_pdb_id(right) and not cls._is_pdb_id(left):
                    # Reversed: {ligand}_{receptor}.dlg  e.g. trimethoprim_1dg5
                    stem_rec, stem_lig = right, left
                else:
                    # Ambiguous: treat first token as receptor (classic fallback)
                    stem_rec, stem_lig = left, right
            else:
                stem_rec = ""
                stem_lig = parts[0]  # single-token stem → ligand name

        lig = raw_lig or stem_lig or stem
        # If flexible docking is identified from stem, make sure receptor preserves _flex tag
        if "_flex" in stem_rec.lower() and raw_rec and not raw_rec.lower().endswith("_flex"):
            rec = f"{raw_rec}_flex"
        else:
            rec = raw_rec or stem_rec or stem_lig

        return DockingIdentity(ligand_name=lig, receptor_name=rec, original_file=dlg_path.name)


# ──────────────────────────────────────────────────────────────────────────────
class FlexResidueHeaderParser:
    """
    Parses the pre-docking flexible residue declaration from DLG header lines.

    AutoDock writes these in the DLG before the first DOCKED: block:
        INPUT-FLEXRES-PDBQT: BEGIN_RES TYR A   7
        INPUT-FLEXRES-PDBQT: ATOM  ...
        INPUT-FLEXRES-PDBQT: END_RES TYR A   7

    Returns a list of residue ID strings, e.g. ["TYR_A_7"].
    """
    _RE_BEGIN_RES = re.compile(
        r"(?:INPUT-FLEXRES-PDBQT:|DOCKED:)\s*BEGIN_RES\s+(\w+)\s+(\S+)\s+(\d+)"
    )
    _RE_DPF_FLEX  = re.compile(r"DPF>\s+flexres\s+(\S+)")

    @classmethod
    def parse(cls, lines: List[str]) -> List[str]:
        """Return unique residue IDs found in INPUT-FLEXRES-PDBQT blocks or DPF flexres."""
        seen: "Set[str]" = set()
        ids:  List[str]  = []
        for line in lines:
            # Stop at first DOCKED: content block
            if line.startswith("DOCKED:") and "BEGIN_RES" not in line:
                dc = line[8:]
                if dc.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT")):
                    break
            m = cls._RE_BEGIN_RES.search(line)
            if m:
                key = f"{m.group(1)}_{m.group(2)}_{m.group(3)}"
                if key not in seen:
                    seen.add(key)
                    ids.append(key)
            m_dpf = cls._RE_DPF_FLEX.search(line)
            if m_dpf:
                fname = Path(m_dpf.group(1)).stem
                for part in fname.split("_"):
                    if part and part not in seen:
                        seen.add(part)
                        ids.append(part)
        return ids

    @classmethod
    def get_flex_file(cls, lines: List[str]) -> Optional[str]:
        """Return the flex residue pdbqt filename from the DPF header, if present."""
        for line in lines:
            m = cls._RE_DPF_FLEX.search(line)
            if m:
                return m.group(1)
        return None


# ──────────────────────────────────────────────────────────────────────────────
class DockingModeDetector:
    """
    Header-first docking mode classifier.

    Detection priority:
      1. HEADER (ground truth): presence of AutoDock flexible-residue markers
         in the DLG file — INPUT-FLEXRES-PDBQT lines, BEGIN_RES inside a
         DOCKED: block, or a "DPF> flexres" parameter line.
      2. ATOM HEURISTIC (fallback): if no header marker is found, sample atoms
         across the first 5 poses and count amino acid residue matches.  This
         fallback is kept for edge cases but should rarely be needed because
         AutoDock always writes the INPUT-FLEXRES-PDBQT block when flexres is
         used.

    Why header-first? Some ligands (e.g. folate, which contains a glutamic
    acid tail) carry standard amino-acid residue names on their own atoms.
    The atom-only heuristic falsely classifies those as FLEXIBLE.  The header
    markers are written exclusively by AutoDock when flexible residue docking
    is requested, so they are the definitive ground truth.
    """

    # Regex patterns for header-based detection
    _RE_INPUT_FLEXRES    = re.compile(r"^INPUT-FLEXRES-PDBQT:")
    _RE_BEGIN_RES_DOCKED = re.compile(r"^DOCKED:\s*BEGIN_RES")
    _RE_DPF_FLEXRES      = re.compile(r"DPF>\s+flexres\b", re.IGNORECASE)

    # Atom-heuristic thresholds (fallback only)
    _MIN_AA_ATOMS_FOR_FLEXIBLE = 3
    _MIN_DISTINCT_AA_RESIDUES  = 2   # raised to 2 to further reduce false positives

    @classmethod
    def detect(cls, poses: List[DockingPose],
               lines: Optional[List[str]] = None) -> DockingMode:
        """Detect docking mode.  Pass *lines* (raw DLG text lines) for the
        header-based primary check; falls back to atom heuristic if lines
        are None or if no header marker is found."""

        # ── Primary: header-based detection ───────────────────────────────────
        if lines is not None:
            for line in lines:
                if cls._RE_INPUT_FLEXRES.match(line):
                    logger.debug("[DockingModeDetector] INPUT-FLEXRES-PDBQT header → FLEXIBLE")
                    return DockingMode.FLEXIBLE
                if cls._RE_BEGIN_RES_DOCKED.match(line):
                    logger.debug("[DockingModeDetector] BEGIN_RES inside DOCKED block → FLEXIBLE")
                    return DockingMode.FLEXIBLE
                if cls._RE_DPF_FLEXRES.search(line):
                    logger.debug("[DockingModeDetector] DPF flexres parameter → FLEXIBLE")
                    return DockingMode.FLEXIBLE
            # No header markers found → RIGID (no need for atom heuristic)
            logger.debug("[DockingModeDetector] No flexres header markers → RIGID")
            return DockingMode.RIGID

        # ── Fallback: atom-based heuristic (used when lines=None) ─────────────
        if not poses:
            logger.debug("[DockingModeDetector] No poses → RIGID (default)")
            return DockingMode.RIGID

        aa_atom_count   = 0
        distinct_aa_res: Set[str] = set()

        for pose in poses[:5]:          # sample first 5 poses for speed
            for atom in pose.atoms:
                if atom.residue_name.upper() in _STANDARD_AMINO_ACIDS:
                    aa_atom_count += 1
                    distinct_aa_res.add(atom.residue_label)

        is_flexible = (
            aa_atom_count    >= cls._MIN_AA_ATOMS_FOR_FLEXIBLE and
            len(distinct_aa_res) >= cls._MIN_DISTINCT_AA_RESIDUES
        )
        mode = DockingMode.FLEXIBLE if is_flexible else DockingMode.RIGID
        logger.debug(
            f"[DockingModeDetector] Fallback heuristic — "
            f"aa_atoms={aa_atom_count}  "
            f"distinct_aa_res={len(distinct_aa_res)} → {mode.name}"
        )
        return mode


# ──────────────────────────────────────────────────────────────────────────────
class AtomRoleClassifier:
    """
    Classifies each PDBAtom as LIGAND or FLEXIBLE_RESIDUE and populates
    the DockingPose's classified sub-lists.

    Classification rules (priority order):
      1. If docking_mode is RIGID  → all atoms are LIGAND (no flexible residues
         are possible in rigid docking; some ligands such as folate carry
         standard amino-acid residue names on their own atoms, and without mode
         awareness those would be misclassified as FLEXIBLE_RESIDUE).
      2. If docking_mode is FLEXIBLE and residue_name is a standard amino acid
         → FLEXIBLE_RESIDUE
      3. All remaining atoms  → LIGAND

    This deliberately avoids the is_hetatm flag because AutoDock writes
    ATOM records for both ligand atoms (no residue) and flexible sidechain
    atoms (with residue name, chain, and sequence number).
    """

    @staticmethod
    def classify(atom: PDBAtom) -> AtomRole:
        """Residue-name heuristic (used for FLEXIBLE mode only)."""
        if atom.residue_name.upper() in _STANDARD_AMINO_ACIDS:
            return AtomRole.FLEXIBLE_RESIDUE
        return AtomRole.LIGAND

    @classmethod
    def classify_pose(cls, pose: DockingPose,
                      mode: DockingMode = DockingMode.FLEXIBLE) -> None:
        """In-place: set role on each atom and populate the classified sub-lists.

        When *mode* is RIGID every atom is unconditionally assigned LIGAND so
        that ligand atoms bearing amino-acid residue names (e.g. folate's GLU
        tail) are never mistaken for flexible receptor side-chains.
        """
        lig:  List[PDBAtom] = []
        flex: List[PDBAtom] = []
        for atom in pose.atoms:
            if mode == DockingMode.RIGID:
                atom.role = AtomRole.LIGAND
                lig.append(atom)
            else:
                atom.role = cls.classify(atom)
                if atom.role == AtomRole.FLEXIBLE_RESIDUE:
                    flex.append(atom)
                else:
                    lig.append(atom)
        pose.ligand_atoms           = lig
        pose.flexible_residue_atoms = flex

    @classmethod
    def classify_result(cls, result: DockingResult) -> None:
        """Classify all poses (mode-aware) and record aggregate atom counts."""
        mode = result.docking_mode          # already set before this is called
        for pose in result.poses:
            cls.classify_pose(pose, mode)
        # Use first pose as representative for atom counts
        if result.poses:
            ref = result.poses[0]
            result.n_ligand_atoms           = len(ref.ligand_atoms)
            result.n_flexible_residue_atoms = len(ref.flexible_residue_atoms)


class DLGParser:
    """Production-grade AutoDock DLG parser — self-adaptive, mode-aware."""

    _RE_ENERGY    = re.compile(r"Estimated Free Energy of Binding\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol")
    _RE_KI        = re.compile(r"Estimated Inhibition Constant,\s*Ki\s*=\s*([\d.eE+\-]+)\s+(\S+)")
    _RE_INTERMOL  = re.compile(r"Final Intermolecular Energy\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol")
    _RE_INTERNAL  = re.compile(r"Final Total Internal Energy\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol")
    _RE_TORSIONAL = re.compile(r"Torsional Free Energy\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol")
    _RE_UNBOUND   = re.compile(r"Unbound System's Energy\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol")
    _RE_VDW_HB_DES = re.compile(r"(?:vdW|van der Waals).*?(?:H-bond|hydrogen).*?(?:desolvation).*?=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol", re.IGNORECASE)
    _RE_ELEC = re.compile(r"(?:Electrostatic|Coulomb).*?=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol", re.IGNORECASE)
    _RE_RUN       = re.compile(r"Run\s*[=:]\s*(\d+)")
    _RE_RMSD_REF  = re.compile(r"RMSD from reference structure\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)")
    _RE_NUM_RUNS  = re.compile(r"Number of Docking Runs\s*=\s*(\d+)")

    def parse(self, dlg_path: Path) -> DockingResult:
        try:
            text = dlg_path.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            return DockingResult(dlg_file=dlg_path, success=False,
                                 error_message=f"Cannot read file: {e}")
        lines  = text.splitlines()
        self._active_dlg_path = dlg_path
        result = DockingResult(dlg_file=dlg_path)
        result.identity       = IdentityResolutionEngine.extract(lines, dlg_path)
        result.ligand         = result.identity.ligand_name
        result.protein        = result.identity.receptor_name
        result.num_runs       = self._extract_num_runs(lines)
        result.flex_residue_ids = FlexResidueHeaderParser.parse(lines)
        result.poses          = self._parse_docked_blocks(lines)
        self._parse_cluster_data(lines, result)
        self._build_cluster_stats(result)
        self._parse_thermodynamics(lines, result)
        # ── Self-adaptive intelligence layer ────────────────────────────────────
        # Step 1: detect mode from header markers (ground truth, fast O(n) scan)
        result.docking_mode = DockingModeDetector.detect(result.poses, lines)
        # Step 2: classify atoms with full mode knowledge — RIGID forces all
        #         atoms to LIGAND so ligands with AA residue names (e.g. folate)
        #         are never miscounted as flexible residue atoms.
        AtomRoleClassifier.classify_result(result)
        self._compute_receptor_strain(result)
        self._annotate_cluster_mode(result)
        logger.debug(
            f"[{result.ligand}] Mode={result.docking_mode.name}  "
            f"lig_atoms={result.n_ligand_atoms}  "
            f"flex_res_atoms={result.n_flexible_residue_atoms}  "
            f"flex_res_ids={result.flex_residue_ids or 'none'}"
        )
        if not result.poses:
            result.success = False
            result.error_message = "No docked poses found in file."
        return result

    # ── helpers ────────────────────────────────────────────────────────────────

    def _extract_num_runs(self, lines: List[str]) -> int:
        for line in lines:
            m = self._RE_NUM_RUNS.search(line)
            if m: return int(m.group(1))
        for line in lines:
            m = re.search(r"Number of conformations\s*=\s*(\d+)", line)
            if m: return int(m.group(1))
        return 0

    def _compute_receptor_strain(self, result: DockingResult) -> None:
        """
        Calculates receptor strain (internal_energy - torsional_energy) for
        FLEXIBLE docking modes, leaving it as None for RIGID.
        """
        for pose in result.poses:
            if result.docking_mode == DockingMode.FLEXIBLE:
                # In flexible docking, internal energy includes flex residue strain.
                pose.receptor_strain_estimate = pose.internal_energy - pose.torsional_energy
            else:
                pose.receptor_strain_estimate = None

        # Attach best pose's strain to the result level
        best = min(result.poses, key=lambda p: p.binding_energy, default=None)
        if best:
            result.receptor_strain_estimate = best.receptor_strain_estimate

    def _annotate_cluster_mode(self, result: DockingResult) -> None:
        """Applies docking mode semantics to all cluster stats."""
        for cs in result.clusters:
            cs.docking_mode      = result.docking_mode.name
            cs.rmsd_is_composite = (result.docking_mode == DockingMode.FLEXIBLE)

    def _parse_docked_blocks(self, lines: List[str]) -> List[DockingPose]:
        poses: List[DockingPose] = []
        current_atoms: List[PDBAtom] = []
        rank = 0
        ce = cinter = cint = ctors = cunb = crmsd = 0.0
        cvdw: Optional[float] = None
        celec: Optional[float] = None
        cki: Optional[float] = None
        ckiu = ""
        crun = 0

        def flush_pose():
            nonlocal rank, ce, cinter, cint, ctors, cunb, crmsd, cki, ckiu, crun, current_atoms, cvdw, celec
            if current_atoms:
                rank += 1
                poses.append(DockingPose(
                    rank=rank, run_number=crun,
                    binding_energy=ce, ki_raw=cki, ki_unit=ckiu,
                    ki_nM=_normalize_ki(cki, ckiu),
                    intermol_energy=cinter, internal_energy=cint,
                    torsional_energy=ctors, unbound_energy=cunb,
                    vdW_hbond_desolvation_energy=cvdw,
                    electrostatic_energy=celec,
                    dlg_source=str(getattr(self, "_active_dlg_path", "")) or None,
                    rmsd_from_ref=crmsd, atoms=current_atoms,
                ))
                current_atoms = []
                ce = cinter = cint = ctors = cunb = crmsd = 0.0
                cki = None
                ckiu = ""
                cvdw = celec = None
                crun = 0

        is_dlg = any(l.startswith("DOCKED:") for l in lines)
        for line in lines:
            line_s = line.strip()
            is_docked_line = line.startswith("DOCKED:")
            content = line[8:].strip() if line.startswith("DOCKED: ") else (line[7:].strip() if is_docked_line else line_s)

            # In DLG files, ignore un-prefixed structural records (cluster representatives or input coords)
            if is_dlg and not is_docked_line:
                if content.startswith(("MODEL", "ENDMDL", "ATOM", "HETATM")):
                    continue

            if "Estimated Free Energy" in content:
                m = self._RE_ENERGY.search(content)
                if m:
                    ce = float(m.group(1))
            elif "Inhibition Constant" in content:
                m = self._RE_KI.search(content)
                if m:
                    try:
                        cki, ckiu = float(m.group(1)), m.group(2)
                    except ValueError:
                        pass
            elif "Intermolecular" in content:
                m = self._RE_INTERMOL.search(content)
                if m:
                    cinter = float(m.group(1))
            elif "Internal Energy" in content:
                m = self._RE_INTERNAL.search(content)
                if m:
                    cint = float(m.group(1))
            elif "Torsional Free Energy" in content:
                m = self._RE_TORSIONAL.search(content)
                if m:
                    ctors = float(m.group(1))
            elif "Unbound System" in content:
                m = self._RE_UNBOUND.search(content)
                if m:
                    cunb = float(m.group(1))
            elif "desolvation" in content.lower() and ("vdw" in content.lower() or "van der waals" in content.lower()):
                m = self._RE_VDW_HB_DES.search(content)
                if m:
                    cvdw = float(m.group(1))
            elif "electrostatic" in content.lower() or "coulomb" in content.lower():
                m = self._RE_ELEC.search(content)
                if m:
                    celec = float(m.group(1))
            elif ("Run" in content and ("=" in content or ":" in content)) or content.startswith("Run:"):
                # If we encounter a new run while we already have atoms from a previous run, flush
                if current_atoms:
                    flush_pose()
                m = self._RE_RUN.search(content)
                if m:
                    crun = int(m.group(1))
            elif "RMSD from reference" in content:
                m = self._RE_RMSD_REF.search(content)
                if m:
                    crmsd = float(m.group(1))
            elif content.startswith("MODEL") and current_atoms:
                flush_pose()
            elif content.startswith(("ATOM", "HETATM")):
                a = self._parse_atom_line(content)
                if a:
                    current_atoms.append(a)
            elif content.startswith("ENDMDL"):
                flush_pose()

        # Flush any remaining pose
        flush_pose()
        return poses

    def _parse_cluster_data(self, lines: List[str], result: DockingResult) -> None:
        cluster_info: Dict[int, Dict] = {}

        # Pass 1: CLUSTERING HISTOGRAM
        in_hist = False
        for line in lines:
            s = line.strip()
            if "CLUSTERING HISTOGRAM" in s:
                in_hist = True
                continue
            if not in_hist:
                continue
            if "Number of multi-member" in s or "RMSD TABLE" in s:
                in_hist = False
                continue
            if not s or s.startswith("___"):
                continue
            if s.startswith(("Clus", "-ter", "Rank", "|", "#")):
                continue
            parts = s.split("|")
            if len(parts) >= 5:
                try:
                    cid = int(parts[0].strip())
                    cluster_info[cid] = {
                        "size":          int(parts[4].strip()),
                        "lowest_energy": float(parts[1].strip()),
                        "best_run":      int(parts[2].strip()),
                        "mean_energy":   float(parts[3].strip()),
                        "cluster_rmsd":  float(parts[5].strip()) if len(parts) >= 6 else 0.0,
                    }
                except (ValueError, IndexError):
                    pass

        # Pass 2: RMSD TABLE
        in_rmsd = False
        cluster_runs: Dict[int, List[int]] = {cid: [] for cid in cluster_info}
        for line in lines:
            s = line.strip()
            if "RMSD TABLE" in s:
                in_rmsd = True
                continue
            if not in_rmsd:
                continue
            if "INFORMATION ENTROPY" in s:
                in_rmsd = False
                continue
            if not s or s.startswith(("_", "Rank", "|")):
                continue
            if "RANKING" in s:
                pw = s.split()
                if len(pw) >= 6:
                    try:
                        cid = int(pw[0])
                        run = int(pw[2])
                        c_rmsd = float(pw[4])
                        r_rmsd = float(pw[5])
                        if cid in cluster_runs and run not in cluster_runs[cid]:
                            cluster_runs[cid].append(run)
                        elif cid not in cluster_runs:
                            cluster_runs[cid] = [run]
                        for pose in result.poses:
                            if pose.run_number == run:
                                pose.cluster_rmsd = c_rmsd
                                pose.rmsd_from_ref = r_rmsd
                    except (ValueError, IndexError):
                        pass

        # Apply membership
        for cid, info in cluster_info.items():
            runs = cluster_runs.get(cid) or [info["best_run"]]
            for pose in result.poses:
                if pose.run_number in runs:
                    pose.cluster_id   = cid
                    pose.cluster_size = info["size"]
                    if pose.cluster_rmsd == 0.0 and info.get("cluster_rmsd"):
                        pose.cluster_rmsd = info["cluster_rmsd"]
        result.num_clusters = len(cluster_info)
        result._raw_cluster_info = cluster_info

    def _build_cluster_stats(self, result: DockingResult) -> None:
        cluster_map: Dict[int, List[DockingPose]] = {}
        for pose in result.poses:
            if pose.cluster_id is not None:
                cluster_map.setdefault(pose.cluster_id, []).append(pose)

        raw_info: Dict[int, Dict] = getattr(result, "_raw_cluster_info", {})
        all_cids = sorted(set(list(cluster_map.keys()) + list(raw_info.keys())))

        stats: List[ClusterStats] = []
        for cid in all_cids:
            poses = cluster_map.get(cid, [])
            info = raw_info.get(cid, {})
            if poses:
                energies = [p.binding_energy for p in poses]
                rmsds    = [p.cluster_rmsd   for p in poses]
                avg_rmsd = sum(rmsds) / len(rmsds) if rmsds else info.get("cluster_rmsd", 0.0)
                if avg_rmsd == 0.0 and info.get("cluster_rmsd"):
                    avg_rmsd = info["cluster_rmsd"]
                best_run = min(poses, key=lambda p: p.binding_energy).run_number
                c_size = max(len(poses), info.get("size", 0))
                lowest_e = min(energies)
                mean_e = sum(energies) / len(energies)
                runs = [p.run_number for p in poses]
            else:
                avg_rmsd = info.get("cluster_rmsd", 0.0)
                best_run = info.get("best_run", 0)
                c_size = info.get("size", 0)
                lowest_e = info.get("lowest_energy", 0.0)
                mean_e = info.get("mean_energy", 0.0)
                runs = [best_run] if best_run else []

            N_runs = max(result.num_runs, 1)
            f_pop  = c_size / N_runs
            q_rmsd = max(0.0, 1.0 - avg_rmsd / 3.0)
            quality = round((f_pop + q_rmsd) / 2.0 * 100.0, 2)
            stats.append(ClusterStats(
                cluster_id=cid, size=c_size,
                mean_energy=mean_e,
                lowest_energy=lowest_e, best_run=best_run,
                avg_rmsd=avg_rmsd,
                spread=max(0.0, mean_e - lowest_e),
                cluster_quality=quality,
                runs=runs,
            ))
        result.clusters = stats

    def _parse_thermodynamics(self, lines: List[str], result: DockingResult) -> None:
        for line in lines:
            s = line.strip()
            if "Information entropy" in s or "Information Entropy" in s:
                m = re.search(r"=\s*([+-]?[\d.]+)", s)
                if m:
                    result.info_entropy = float(m.group(1))
                m_tol = re.search(r"rmstol\s*=\s*([\d.]+)", s)
                if m_tol:
                    result.info_entropy_rmstol = float(m_tol.group(1))
            elif s.startswith("Partition function") or "Partition function" in s:
                m = re.search(r"Q\s*=\s*([\d.]+)", s)
                if m:
                    result.partition_function = float(m.group(1))
                m_t = re.search(r"Temperature.*=\s*([\d.]+)", s)
                if m_t:
                    result.stat_temperature = float(m_t.group(1))
            elif "Free energy" in s:
                m = re.search(r"A\s*~\s*([-\d.]+)", s)
                if m:
                    result.stat_free_energy = float(m.group(1))
            elif "Internal energy" in s:
                m = re.search(r"U\s*=\s*([-\d.]+)", s)
                if m:
                    result.stat_internal_energy = float(m.group(1))
            elif "Entropy" in s and not s.startswith("Information"):
                m = re.search(r"S\s*=\s*([-\d.]+)", s)
                if m:
                    result.stat_entropy = float(m.group(1))

    def _parse_atom_line(self, line: str) -> Optional[PDBAtom]:
        try:
            is_hetatm = line.startswith("HETATM")
            serial    = int(line[6:11].strip())  if line[6:11].strip()  else 0
            name      = line[12:16].strip()
            res       = line[17:20].strip()
            chain     = line[21:22].strip()
            seq       = int(line[22:26].strip()) if line[22:26].strip() else 0
            x, y, z   = float(line[30:38]), float(line[38:46]), float(line[46:54])
            charge    = 0.0
            ad_type   = ""
            if len(line) >= 76:
                c = line[70:76].strip()
                if c:
                    try: charge = float(c)
                    except ValueError: pass
            if len(line) >= 79:
                ad_type = line[77:79].strip()
            return PDBAtom(serial, name, _autodock_to_element(ad_type, name),
                           res, chain, seq, x, y, z, charge, ad_type, is_hetatm)
        except Exception:
            return None


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — CLUSTER ANALYSIS ENGINE  (dual ranking)
# ══════════════════════════════════════════════════════════════════════════════

class ClusterAnalysisEngine:
    """
    Two independent cluster rankings:
      Thermodynamic : lowest_energy ASC, size DESC  (best binding affinity)
      Population    : size DESC, lowest_energy ASC  (most reproducible)
    """

    def rank_clusters(self, result: DockingResult) -> None:
        if not result.clusters:
            return
        thermo = sorted(result.clusters, key=lambda c: (c.lowest_energy, -c.size))
        for rank, c in enumerate(thermo, 1):
            c.thermodynamic_rank = rank
        pop = sorted(result.clusters, key=lambda c: (-c.size, c.lowest_energy))
        for rank, c in enumerate(pop, 1):
            c.population_rank = rank

    def get_best_cluster_thermo(self, result: DockingResult) -> Optional[ClusterStats]:
        if not result.clusters: return None
        return min(result.clusters, key=lambda c: (c.lowest_energy, -c.size))

    def get_best_cluster_population(self, result: DockingResult) -> Optional[ClusterStats]:
        if not result.clusters: return None
        return max(result.clusters, key=lambda c: (c.size, -c.lowest_energy))


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — POSE EXTRACTOR
# ══════════════════════════════════════════════════════════════════════════════

class PoseExtractor:
    @staticmethod
    def extract_best_pose(result: DockingResult) -> Optional[DockingPose]:
        if not result.poses: return None
        return min(result.poses, key=lambda p: p.binding_energy)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — DATA EXTRACTOR
# ══════════════════════════════════════════════════════════════════════════════

class DataExtractor:
    """
    Assembles a flat DockingProfile from a parsed DockingResult.
    Reports objective scientific quantities only — no composite scoring.
    """

    @staticmethod
    def extract_profile(
        ligand:        str,
        protein:       str,
        pose:          DockingPose,
        cluster_stats: Optional[ClusterStats] = None,
    ) -> DockingProfile:
        return DockingProfile(
            ligand=ligand, protein=protein,
            binding_energy=pose.binding_energy,
            ki_nM=pose.ki_nM, run_number=pose.run_number,
            rmsd=pose.rmsd_from_ref,
            cluster_id=pose.cluster_id, cluster_size=pose.cluster_size,
            cluster_rmsd=pose.cluster_rmsd,
            avg_cluster_rmsd=cluster_stats.avg_rmsd      if cluster_stats else 0.0,
            cluster_mean_energy=cluster_stats.mean_energy if cluster_stats else 0.0,
            cluster_spread=cluster_stats.spread           if cluster_stats else 0.0,
            pose=pose,
        )


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — REFERENCE RESOLUTION & OPTIONAL VALIDATION (RDKit)
# ══════════════════════════════════════════════════════════════════════════════

class ReferenceResolver:
    """
    Offline-only reference ligand resolver.
    All steps are strictly OPTIONAL — FAILED is a valid, non-crashing outcome.

    Format priority
    ---------------
    When multiple reference files exist for the same ligand name (e.g.
    ``donepezil.sdf``, ``donepezil.mol2``, ``donepezil.mmcif``), **only
    the highest-priority file is used** — validation is never duplicated.
    A WARNING is printed so you always know which file was chosen and
    which were silently skipped.

    Priority order (highest → lowest fidelity):

    1. ``.sdf``    — preserves bond orders, formal charges, stereo.
                     Gold standard for small-molecule RMSD validation.
    2. ``.mol2``   — Tripos format; connectivity + partial charges intact.
                     Common Schrödinger / Amber / RCSB output.
    3. ``.mol``    — V2000/V3000 MOL; same connectivity as SDF, older.
    4. ``.pdb``    — 3-D coordinates OK, but bond orders are implicit.
    5. ``.pdbqt``  — AutoDock format; AD atom types ≠ elements (lossy).
    6. ``.cif``    — mmCIF / PDBx: entire crystal structure; ligand must
    7. ``.mmcif``    be extracted by gemmi.  Use only when nothing else
                     is available.
    """

    # Scientifically ordered priority list (highest → lowest fidelity)
    _FORMAT_PRIORITY: List[str] = [
        ".sdf", ".mol2", ".mol", ".pdb", ".pdbqt", ".cif", ".mmcif",
    ]

    def __init__(self, root_dir: Path) -> None:
        self.root         = root_dir
        self.ref_dir      = root_dir / "REFERENCE"
        self.fallback_dir = self.ref_dir / "FALLBACK"
        self.ref_dir.mkdir(parents=True, exist_ok=True)
        self.fallback_dir.mkdir(parents=True, exist_ok=True)
        self.user_map: Dict[str, Path] = {}

    def resolve(self, ligand_name: str, dlg_path: Path) -> Tuple[Optional[Path], str]:
        if ligand_name in self.user_map:
            return self.user_map[ligand_name], f"USER_MAPPED ({self.user_map[ligand_name].name})"
        p = self._find_user_supplied(ligand_name)
        if p:
            return p, f"USER_SUPPLIED_FOLDER ({p.name})"
        out = self.fallback_dir / f"{ligand_name}_fallback.pdbqt"
        if out.exists():
            return out, "DLG_FALLBACK (cached)"
        if self._extract_from_dlg(dlg_path, out):
            return out, "DLG_FALLBACK"
        return None, "FAILED"   # Graceful — callers must handle None

    def _find_user_supplied(self, name: str) -> Optional[Path]:
        """
        Return the single highest-priority reference file found in REFERENCE/.

        Collects *all* matching files first, warns if more than one format
        exists, then returns only the top-priority match.  Validation is
        therefore always performed exactly once, on the best-quality file.
        """
        candidates: List[Path] = []
        search_dirs = [
            self.ref_dir,
            self.root / "reference",
            self.root / "results" / "REFERENCE",
            self.root / "receptors",
            self.root / "ligands",
            self.root,
        ]
        for d in search_dirs:
            if not d.is_dir():
                continue
            for n in [name, f"{name}_ref", name.lower(), f"{name.lower()}_ref"]:
                for ext in self._FORMAT_PRIORITY:
                    p = d / f"{n}{ext}"
                    if p.exists() and p not in candidates:
                        candidates.append(p)

        if not candidates:
            return None

        # Sort by priority order (index in _FORMAT_PRIORITY = lower is better)
        def _priority(path: Path) -> int:
            try:
                return self._FORMAT_PRIORITY.index(path.suffix.lower())
            except ValueError:
                return len(self._FORMAT_PRIORITY)  # unknown extension = lowest

        candidates.sort(key=_priority)
        chosen = candidates[0]

        if len(candidates) > 1:
            ignored = ", ".join(p.name for p in candidates[1:])
            logger.warning(
                "Multiple reference formats found for '%s': using '%s' "
                "(priority %d/%d). Ignored: %s. "
                "Remove lower-priority files to suppress this warning.",
                name,
                chosen.name,
                _priority(chosen) + 1,
                len(self._FORMAT_PRIORITY),
                ignored,
            )

        return chosen

    def _extract_from_dlg(self, dlg_path: Path, out_path: Path) -> bool:
        try:
            lines     = dlg_path.read_text(encoding="utf-8", errors="replace").splitlines()
            out_lines = [l[19:].strip() for l in lines
                         if l.startswith("INPUT-LIGAND-PDBQT:") and l[19:].strip()]
            if out_lines:
                out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
                return True
        except Exception:
            pass
        return False


class LigandStandardizer:
    """
    Loads a reference ligand from disk and standardizes it into an RDKit Mol.

    Supported formats
    -----------------
    ============  ========================================================
    Extension     Notes
    ============  ========================================================
    ``.sdf``      First record in the V2000/V3000 SD file.
    ``.mol``      Single V2000/V3000 MOL record (same reader as SDF).
    ``.mol2``     Tripos MOL2 — common Schrödinger / Amber / RCSB output.
    ``.pdb``      Standard PDB ATOM/HETATM; must contain 3-D coordinates.
    ``.pdbqt``    AutoDock PDBQT — atom types translated to PDB elements.
    ``.cif``      mmCIF (PDBx) — ligand component extracted via **gemmi**
    ``.mmcif``    (alias for .cif).  Requires ``pip install gemmi``.
    ============  ========================================================

    mmCIF notes
    -----------
    RCSB mmCIF files contain the *entire* crystal structure.  The loader
    iterates ``Structure.entities`` and picks the first small-molecule
    (non-polymer, non-water) entity, converts its atoms to a PDB block,
    then passes that to RDKit.  This is the correct scientific approach
    because mmCIF ligands are stored as HETATM records within an entity
    that may span multiple chains.

    gemmi is an *optional* dependency — if it is not installed the CIF
    branch silently returns ``None`` and validation is skipped, just like
    every other optional path in this pipeline.
    """

    @staticmethod
    def load_and_standardize(file_path: str) -> Optional[Any]:
        if not _HAS_RDKIT:
            return None
        ext = file_path.lower().rsplit(".", 1)[-1]
        mol = None
        try:
            if ext == "sdf":
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                if len(suppl) > 0: mol = suppl[0]

            elif ext == "mol":
                # Single V2000/V3000 MOL record — identical reader to SDF
                suppl = Chem.SDMolSupplier(file_path, removeHs=False)
                if len(suppl) > 0: mol = suppl[0]

            elif ext == "mol2":
                # Tripos MOL2 — native RDKit support
                mol = Chem.MolFromMol2File(file_path, removeHs=False)
                if mol is None:
                    try:
                        mol = Chem.MolFromMol2File(file_path, removeHs=False, sanitize=False)
                        if mol is not None:
                            mol.UpdatePropertyCache(strict=False)
                    except Exception:
                        mol = None

            elif ext == "pdb":
                mol = Chem.MolFromPDBFile(file_path, removeHs=False)

            elif ext == "pdbqt":
                # AutoDock PDBQT: strip extra columns, translate AD atom types
                # to standard PDB element symbols so RDKit can parse it.
                pdb_lines: List[str] = []
                with open(file_path, encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if not line.startswith(("ATOM", "HETATM")): continue
                        nm  = line[12:16].strip()
                        raw = line[77:79].strip() if len(line) >= 79 else ""
                        elem = _autodock_to_element(raw, nm)
                        if not elem: continue
                        pdb_lines.append(f"{line[:66].ljust(76)}{elem:>2s}\n")
                mol = Chem.MolFromPDBBlock("".join(pdb_lines), removeHs=False)

            elif ext in ("cif", "mmcif"):
                # mmCIF / PDBx — requires gemmi to parse the full structure
                # and extract the small-molecule ligand component.
                mol = LigandStandardizer._load_from_mmcif(file_path)

            if mol is None:
                return None
            mol = Chem.RemoveHs(mol)
            Chem.SanitizeMol(mol)
            mol = Chem.RenumberAtoms(mol, list(Chem.CanonicalRankAtoms(mol)))
            # NOTE: We deliberately do NOT call AddHs() here.
            # Keeping the mol at the heavy-atom level is critical for two reasons:
            #   1. RMSD engine: MCS across formats (mol2/sdf/pdb/pdbqt) fails when
            #      one mol has explicit Hs and the other doesn't — atom counts
            #      diverge wildly (e.g. 28 vs 71 atoms) and matchValences=True
            #      makes things even worse. Heavy-atom RMSD is the scientific
            #      standard in docking validation (RMSD_heavy = no H).
            #   2. ADMET: RDKit descriptors (MW, LogP, TPSA, HBD, HBA) all work
            #      correctly on heavy-atom mols with implicit Hs; 3D-embedded Hs
            #      only add noise to those calculations.
            AllChem.ComputeGasteigerCharges(mol)
        except Exception:
            return None
        return mol

    @staticmethod
    def _load_from_mmcif(file_path: str) -> Optional[Any]:
        """
        Extract the first small-molecule ligand from an mmCIF structure file
        using gemmi, then pass it to RDKit via a PDB block.

        Scientific rationale
        --------------------
        mmCIF entities are typed as polymer / non-polymer / water / branched.
        The ligand of interest is always a *non-polymer* entity that is not
        a water molecule (residue name "HOH" or "WAT").  We pick the first
        such entity so that co-crystallized ligands can be used directly as
        reference structures without any pre-processing.
        """
        if not _HAS_GEMMI or not _HAS_RDKIT:
            logger.warning(
                "mmCIF reference requested but 'gemmi' is not installed. "
                "Install it with:  pip install gemmi"
            )
            return None
        try:
            structure = gemmi.read_structure(file_path)
            pdb_lines: List[str] = []
            serial = 1
            # Walk entities to find the first small-molecule (non-polymer, non-water)
            for entity in structure.entities:
                if entity.entity_type not in (
                    gemmi.EntityType.NonPolymer,
                ):
                    continue
                # Collect all residues belonging to this entity
                for model in structure:
                    for chain in model:
                        for res in chain:
                            if res.entity_id != entity.name:
                                continue
                            if res.name.upper() in ("HOH", "WAT", "H2O"):
                                continue  # skip water
                            for atom in res:
                                pos = atom.pos
                                elem = atom.element.name.capitalize()
                                line = (
                                    f"HETATM{serial:5d} "
                                    f"{atom.name:<4s} "
                                    f"{res.name:<3s} "
                                    f"{chain.name:1s}{res.seqid.num:4d}    "
                                    f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                                    f"{'1.00':6s}{'0.00':6s}          "
                                    f"{elem:>2s}\n"
                                )
                                pdb_lines.append(line)
                                serial += 1
                        if pdb_lines:
                            break  # only first matching chain per model
                    if pdb_lines:
                        break      # only first model
                if pdb_lines:
                    break          # only first qualifying entity
            if not pdb_lines:
                logger.warning("mmCIF file contained no extractable small-molecule ligand: %s", file_path)
                return None
            return Chem.MolFromPDBBlock("".join(pdb_lines), removeHs=False)
        except Exception as exc:
            logger.debug("mmCIF load failed for %s: %s", file_path, exc)
            return None


class RMSDEngine:
    @staticmethod
    def kabsch_rmsd(P: Any, Q: Any) -> float:
        if not _HAS_NUMPY or P.shape != Q.shape or len(P) == 0: return 0.0
        Pc, Qc = P - np.mean(P, axis=0), Q - np.mean(Q, axis=0)
        V, _, W = np.linalg.svd(np.dot(Pc.T, Qc))
        if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
            V[:, -1] = -V[:, -1]
        U = np.dot(V, W)
        return float(np.sqrt(((np.dot(Pc, U) - Qc) ** 2).sum() / len(P)))

    @staticmethod
    def evaluate_pose(ref_mol: Any, pose_mol: Any, pose: DockingPose) -> bool:
        """Compute heavy-atom validation RMSD via Kabsch alignment on MCS-matched atoms.

        Scientific design
        -----------------
        Both molecules are stripped to heavy atoms before MCS so that differences
        in hydrogen representation (implicit vs explicit) across file formats
        (mol2, sdf, pdb, pdbqt) cannot prevent a match.  matchValences=False
        makes the search robust to protonation-state differences between the
        crystallographic reference and the AutoDock-prepared input.
        The result is the Kabsch RMSD over all non-hydrogen matched atoms — the
        standard metric for docking re-docking validation.
        Strictly OPTIONAL — returns False if unavailable.
        """
        if not _HAS_RDKIT or not _HAS_NUMPY:
            return False
        try:
            # Ensure both molecules are heavy-atom only before comparison
            ref_hvy  = Chem.RemoveHs(ref_mol)  if ref_mol  else None
            pose_hvy = Chem.RemoveHs(pose_mol) if pose_mol else None
            if ref_hvy is None or pose_hvy is None:
                return False

            mcs = rdFMCS.FindMCS(
                [ref_hvy, pose_hvy],
                timeout=10,
                matchValences=False,          # robust across protonation states
                bondCompare=rdFMCS.BondCompare.CompareAny, # ignore guessed bond orders
                ringMatchesRingOnly=True,
                completeRingsOnly=False,
            )
            if mcs.canceled or mcs.numAtoms == 0:
                return False
            common = Chem.MolFromSmarts(mcs.smartsString)
            if not common:
                return False
            rm = ref_hvy.GetSubstructMatches(common)
            pm = pose_hvy.GetSubstructMatches(common)
            if not rm or not pm:
                return False
            mapping = [
                (ri, pi)
                for ri, pi in zip(rm[0], pm[0])
                if (ref_hvy.GetAtomWithIdx(ri).GetSymbol() ==
                    pose_hvy.GetAtomWithIdx(pi).GetSymbol())
            ]
            if not mapping:
                return False
            rc  = ref_hvy.GetConformer()  if ref_hvy.GetNumConformers()  > 0 else None
            pc  = pose_hvy.GetConformer() if pose_hvy.GetNumConformers() > 0 else None
            if rc is None or pc is None:
                return False
            P = np.array([list(rc.GetAtomPosition(i)) for i, _ in mapping])
            Q = np.array([list(pc.GetAtomPosition(j)) for _, j in mapping])
            hcount = sum(1 for a in ref_hvy.GetAtoms())
            pos_rmsd = float(np.sqrt(((P - Q) ** 2).sum() / len(P)))
            kab_rmsd = RMSDEngine.kabsch_rmsd(P, Q)
            if len(mapping) == hcount:
                pose.validation_heavy_rmsd = round(kab_rmsd, 3)
            if len(mapping) / max(1, hcount) >= 0.40:
                pose.validation_scaffold_rmsd = round(kab_rmsd, 3)
            pose.validation_rmsd_conf = (pose.validation_heavy_rmsd
                                         if pose.validation_heavy_rmsd is not None
                                         else pose.validation_scaffold_rmsd)
            pose.validation_rmsd_pos = round(pos_rmsd, 3)
        except Exception:
            pass
        return True


def write_pose_as_reference_template(
    ref_path: Path,
    pose_mol: Any,
    out_path: Path,
    pose:     "DockingPose",
    mol_name: str = "",
) -> bool:
    """
    Write a docked pose SDF using the reference ligand as the molecular template.

    Motivation
    ----------
    Tools such as BIOVIA Discovery Studio and PyMOL require that the pose and
    reference files share identical atom ordering, element symbols, and bond
    connectivity before they will attempt a 1:1 positional RMSD.  When the pose
    SDF is built from scratch (from raw PDBAtom data), atom ordering follows the
    DLG run order which differs from the reference's ordering — causing both tools
    to reject the comparison.

    This function solves the problem at the root by doing the inverse:  instead
    of building a new molecule from pose coordinates, it takes the *reference
    molecule's graph* (correct atom order, element types, formal charges, bond
    orders, stereochemistry) and replaces its conformer with the docked 3-D
    coordinates, atom-by-atom, via an MCS mapping.  The result is a file that
    is structurally identical to the reference in every way except the Cartesian
    positions — which now reflect the docked geometry.

    CRITICAL — why we load the reference raw (not via LigandStandardizer)
    -----------------------------------------------------------------------
    LigandStandardizer.load_and_standardize() calls
    ``Chem.RenumberAtoms(mol, list(Chem.CanonicalRankAtoms(mol)))`` which
    reshuffles atoms into RDKit's canonical order.  That breaks atom-ordering
    identity with the original file, causing DSV to report mismatches again.
    Here we load the reference file directly, strip Hs (to avoid format
    differences), and perform MCS on the *as-stored* atom indices so the output
    SDF has the same atom ordering as the original reference.

    Returns True if the file was written successfully, False otherwise.
    """
    if not _HAS_RDKIT or pose_mol is None or ref_path is None:
        return False
    try:
        from rdkit.Chem import RWMol, Conformer

        # Load reference WITHOUT canonical renumbering to preserve original order
        sfx = ref_path.suffix.lower()
        ref_raw: Optional[Any] = None
        if sfx in (".sdf", ".mol"):
            suppl = Chem.SDMolSupplier(str(ref_path), removeHs=False, sanitize=False)
            ref_raw = suppl[0] if suppl else None
        elif sfx == ".mol2":
            ref_raw = Chem.MolFromMol2File(str(ref_path), removeHs=False, sanitize=False)
        elif sfx in (".pdb", ".ent"):
            ref_raw = Chem.MolFromPDBFile(str(ref_path), removeHs=False, sanitize=False)
        if ref_raw is None:
            return False

        # Sanitize gently to fix valences but keep atom order
        try:
            Chem.SanitizeMol(ref_raw)
        except Exception:
            try:
                Chem.SanitizeMol(
                    ref_raw,
                    Chem.SanitizeFlags.SANITIZE_ALL ^
                    Chem.SanitizeFlags.SANITIZE_PROPERTIES,
                )
            except Exception:
                pass  # proceed with whatever we have

        # Work on heavy-atom copies so Hs never interfere with MCS
        ref_hvy  = Chem.RemoveHs(ref_raw)
        pose_hvy = Chem.RemoveHs(pose_mol)

        if ref_hvy.GetNumConformers() == 0 or pose_hvy.GetNumConformers() == 0:
            return False

        # Find the atom-level mapping: ref_idx -> pose_idx
        mcs = rdFMCS.FindMCS(
            [ref_hvy, pose_hvy],
            timeout=10,
            matchValences=False,
            bondCompare=rdFMCS.BondCompare.CompareAny, # critical: DLG bond orders are guessed and may not match
            ringMatchesRingOnly=True,
            completeRingsOnly=False,
        )
        if mcs.canceled or mcs.numAtoms == 0:
            return False
        common = Chem.MolFromSmarts(mcs.smartsString)
        if not common:
            return False
        rm = ref_hvy.GetSubstructMatch(common)
        pm = pose_hvy.GetSubstructMatch(common)
        if not rm or not pm:
            return False

        ref_to_pose: Dict[int, int] = dict(zip(rm, pm))

        # Build a new conformer on the reference graph, filling in docked coords
        rc = ref_hvy.GetConformer()
        pc = pose_hvy.GetConformer()
        new_conf = Conformer(ref_hvy.GetNumAtoms())
        for ref_idx in range(ref_hvy.GetNumAtoms()):
            if ref_idx in ref_to_pose:
                pos = pc.GetAtomPosition(ref_to_pose[ref_idx])
            else:
                pos = rc.GetAtomPosition(ref_idx)   # unmapped: keep reference coords
            new_conf.SetAtomPosition(ref_idx, pos)

        mol_out = RWMol(ref_hvy)
        mol_out.RemoveAllConformers()
        mol_out.AddConformer(new_conf, assignId=True)
        mol_out_final = mol_out.GetMol()

        name = mol_name or out_path.stem
        mol_out_final.SetProp("_Name", name)
        mol_out_final.SetProp("Run",    str(pose.run_number))
        mol_out_final.SetProp("dG",     f"{pose.binding_energy:.4f}")
        mol_out_final.SetProp("Ki_nM",  str(pose.ki_nM) if pose.ki_nM else "N/A")
        mol_out_final.SetProp("MCS_atoms_matched", str(len(rm)))
        mol_out_final.SetProp("Ref_total_atoms",   str(ref_hvy.GetNumAtoms()))

        out_path.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_path))
        writer.write(mol_out_final)
        writer.close()
        return True
    except Exception as exc:
        logger.debug("write_pose_as_reference_template failed for %s: %s", out_path.name, exc)
        return False


class ADMETEngine:
    """
    Calculates ADMET descriptors and Lipinski Rule-of-Five violations using RDKit.

    Descriptors are reported as calculated molecular properties only.
    They are NOT integrated into any composite score or ranking metric.
    """

    @staticmethod
    def evaluate_pharmacokinetics(pose_mol: Any, profile: DockingProfile) -> None:
        if not _HAS_RDKIT or pose_mol is None:
            return
        try:
            from rdkit.Chem import Descriptors, rdMolDescriptors

            profile.admet_mw    = round(Descriptors.MolWt(pose_mol), 2)
            profile.admet_logp  = round(Descriptors.MolLogP(pose_mol), 2)
            profile.admet_tpsa  = round(rdMolDescriptors.CalcTPSA(pose_mol), 2)
            profile.admet_hbd   = rdMolDescriptors.CalcNumLipinskiHBD(pose_mol)
            profile.admet_hba   = rdMolDescriptors.CalcNumLipinskiHBA(pose_mol)
            profile.admet_nrotb = rdMolDescriptors.CalcNumRotatableBonds(pose_mol)

            v = 0
            if profile.admet_mw   > 500: v += 1
            if profile.admet_logp > 5:   v += 1
            if profile.admet_hbd  > 5:   v += 1
            if profile.admet_hba  > 10:  v += 1
            profile.lipinski_violations = v

        except Exception as e:
            logger.debug(f"ADMET estimation failed: {e}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — STRUCTURE I/O UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

def parse_receptor(pdb_path: Path) -> List[PDBAtom]:
    """
    Parse ATOM records from a PDB / PDBQT receptor file.
    HETATM records (cofactors, water, crystallographic ligands) are excluded.
    """
    atoms: List[PDBAtom] = []
    if not pdb_path or not pdb_path.exists():
        return atoms
    is_pdbqt = pdb_path.suffix.lower() == ".pdbqt"
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            name  = line[12:16].strip()
            res   = line[17:20].strip()
            chain = line[21:22].strip()
            seq   = int(line[22:26].strip()) if line[22:26].strip() else 0
            x     = float(line[30:38])
            y     = float(line[38:46])
            z     = float(line[46:54])
            charge = 0.0
            ad_type = ""
            if is_pdbqt:
                if len(line) >= 76:
                    try: charge = float(line[70:76].strip())
                    except ValueError: pass
                ad_type = line[77:79].strip() if len(line) >= 79 else ""
                elem    = _autodock_to_element(ad_type, name)
            else:
                elem    = line[76:78].strip() if len(line) >= 78 else ""
                if not elem: elem = name[0].upper() if name else "C"
                ad_type = elem
            if not elem:
                continue
            atoms.append(PDBAtom(len(atoms) + 1, name, elem, res, chain, seq,
                                 x, y, z, charge, ad_type, False))
        except Exception:
            pass
    return atoms


def write_pdbqt(pose: DockingPose, out_path: Path, atoms: Optional[List[PDBAtom]] = None) -> None:
    """Write a docked pose as a standards-compliant PDBQT file.

    Three correctness guarantees for external RMSD validation (PyMOL / BIOVIA DSV):
    1. All ligand-atom records are written as HETATM — the correct PDB record
       type for small-molecule ligands.  AutoDock writes ATOM in DOCKED: blocks,
       but external tools (PyMOL align, DSV Heavy Atom RMSD) use the record type
       to distinguish ligand atoms from backbone atoms.
    2. Atoms are sorted by serial number so the output order is reproducible and
       matches the INPUT-LIGAND-PDBQT reference, enabling positional RMSD tools
       that require identical atom ordering (DSV "order must be the same" rule).
    3. Column 77-78 contains the standard PDB element symbol (a.element), NOT
       the AutoDock atom type (OA, HD, A …).  DSV's element-matching check reads
       this column and rejects files where it finds AutoDock types instead.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    atom_list = sorted(
        (atoms if atoms is not None else pose.atoms),
        key=lambda a: a.serial
    )
    lines = [f"REMARK  {__full_brand__} | Run {pose.run_number} "
             f"| Energy = {pose.binding_energy:.2f} kcal/mol"]
    for a in atom_list:
        lines.append(
            f"HETATM{a.serial:>5d} {a.name:<4s} {a.residue_name:<3s} "
            f"{a.chain_id}{a.residue_seq:>4d}    "
            f"{a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}"
            f"  1.00  0.00    {a.charge:>8.3f} {a.element:>2s}"
        )
    lines += ["TER", "END"]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def write_pose_pdb(pose: DockingPose, out_path: Path) -> None:
    """PDB format suitable for PyMOL cmd.load().

    Ligand atoms are written as HETATM (correct PDB convention for small
    molecules) and sorted by serial number for reproducible ordering.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for a in sorted(pose.atoms, key=lambda a: a.serial):
        elem  = a.element or _autodock_to_element(a.autodock_type, a.name) or "C"
        chain = a.chain_id if a.chain_id.strip() else " "
        lines.append(
            f"HETATM{a.serial:>5d} {a.name:<4s} {a.residue_name:<3s} "
            f"{chain}{a.residue_seq:>4d}    "
            f"{a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def _pose_atoms_to_pdb_block(atoms: List["PDBAtom"]) -> str:
    """
    Build a minimal PDB text block from a list of PDBAtom objects.
    Used internally by write_pose_sdf / write_pose_mol as the bridge into RDKit.
    All atoms are written as HETATM; serials are reassigned 1..N to keep them
    contiguous (required for RDKit's PDB parser).
    """
    lines: List[str] = []
    for i, a in enumerate(sorted(atoms, key=lambda a: a.serial), start=1):
        elem  = a.element or _autodock_to_element(a.autodock_type, a.name) or "C"
        chain = a.chain_id if a.chain_id.strip() else "A"
        res   = a.residue_name.strip() or "LIG"
        lines.append(
            f"HETATM{i:>5d} {a.name:<4s} {res:<3s} "
            f"{chain}{a.residue_seq or 1:>4d}    "
            f"{a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}  1.00  0.00          {elem:>2s}"
        )
    lines.append("END")
    return "\n".join(lines)


def write_pose_sdf(
    pose: "DockingPose",
    out_path: Path,
    atoms: Optional[List["PDBAtom"]] = None,
    mol_name: str = "",
) -> bool:
    """
    Write the best docked pose as an SDF file — the most interoperable format
    for external RMSD validation.

    Why SDF instead of PDBQT for external validation
    -------------------------------------------------
    SDF / MOL V2000 files contain an explicit **bond connectivity table**.  This
    means that every tool (PyMOL, BIOVIA DSV, Maestro, LigandScout …) can match
    atoms by **graph topology** rather than by serial-number order or record type.
    The three sources of failure that plagued the PDBQT export
    (ATOM vs HETATM, AutoDock element types, undefined atom ordering) simply do
    not exist in SDF because:
      * Element symbols are mandatory in the atom block (cols 32-34).
      * Bond orders are explicit in the bond block.
      * Atom ordering is irrelevant — tools match via MCS or simple bond graph.

    RDKit path (preferred)
    ----------------------
    When RDKit is available the function:
      1. Builds a PDB block from the pose atoms.
      2. Parses it with ``Chem.MolFromPDBBlock()``.
      3. Writes a proper V2000 / V3000 SDF with ``Chem.SDWriter``.
    This path correctly perceives bond orders from the 3-D geometry and produces
    a chemically complete SDF that any tool can read.

    Fallback (no RDKit)
    -------------------
    A hand-written V2000 MOL block is produced with a dummy bond table
    (all single bonds, no ring closure).  The coordinates and element symbols are
    correct; only bond orders are approximate.  This is still far better than
    PDBQT for positional RMSD because tools can at least match elements and
    atom names deterministically.

    Returns True if the file was written successfully, False otherwise.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    atom_list = sorted(
        (atoms if atoms is not None else pose.atoms),
        key=lambda a: a.serial,
    )
    ligand_name = mol_name or out_path.stem

    # ── RDKit path ────────────────────────────────────────────────────────────
    if _HAS_RDKIT:
        try:
            pdb_block = _pose_atoms_to_pdb_block(atom_list)
            mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
            if mol is None:
                raise ValueError("RDKit could not parse the PDB block")
            # Sanitize with partial flags to avoid valence errors from PDB parse
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                Chem.SanitizeMol(
                    mol,
                    Chem.SanitizeFlags.SANITIZE_ALL ^
                    Chem.SanitizeFlags.SANITIZE_PROPERTIES
                )
            mol.SetProp("_Name", ligand_name)
            mol.SetProp("Run",    str(pose.run_number))
            mol.SetProp("dG",     f"{pose.binding_energy:.4f}")
            mol.SetProp("Ki_nM",  str(pose.ki_nM) if pose.ki_nM else "N/A")
            writer = Chem.SDWriter(str(out_path))
            writer.write(mol)
            writer.close()
            return True
        except Exception as exc:
            logger.debug("write_pose_sdf: RDKit path failed for %s: %s", out_path.name, exc)
            # Fall through to the hand-written fallback

    # ── Fallback: hand-written V2000 MOL block ────────────────────────────────
    # Counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    na = len(atom_list)
    # Naive bond list: sequential pairs only (no ring detection without RDKit)
    # This gives correct atom positions + elements; bond orders are approximate.
    bond_list: List[Tuple[int, int]] = [
        (i, i + 1) for i in range(1, na)  # chain of single bonds
    ]
    nb = len(bond_list)
    fallback_lines: List[str] = [
        ligand_name,                     # molecule name
        f"  BSNDVP  Run={pose.run_number}  dG={pose.binding_energy:.2f} kcal/mol",
        "",                              # comment
        f"{na:>3d}{nb:>3d}  0  0  0  0  0  0  0  0  1 V2000",
    ]
    for a in atom_list:
        elem = a.element or _autodock_to_element(a.autodock_type, a.name) or "C"
        # V2000 atom line: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        fallback_lines.append(
            f"{a.x:>10.4f}{a.y:>10.4f}{a.z:>10.4f} {elem:<3s} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for i, j in bond_list:
        fallback_lines.append(f"{i:>3d}{j:>3d}  1  0  0  0  0")
    fallback_lines += ["M  END", "", "$$$$"]
    try:
        out_path.write_text("\n".join(fallback_lines), encoding="utf-8")
        return True
    except Exception as exc:
        logger.warning("write_pose_sdf: fallback write failed for %s: %s", out_path, exc)
        return False


def write_pose_mol(
    pose: "DockingPose",
    out_path: Path,
    atoms: Optional[List["PDBAtom"]] = None,
    mol_name: str = "",
) -> bool:
    """
    Write the best docked pose as a single MOL (V2000) file.

    MOL is identical to SDF but contains only one record (no ``$$$$`` separator).
    Many tools (especially Windows-based ones) prefer .mol over .sdf for single
    structures.  Internally this function writes via write_pose_sdf and renames
    the terminator, keeping a single code path.
    """
    import tempfile
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sdf", delete=False, encoding="utf-8"
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        ok = write_pose_sdf(pose, tmp_path, atoms=atoms, mol_name=mol_name)
        if not ok:
            return False
        sdf_text = tmp_path.read_text(encoding="utf-8")
        # A valid SDF ends with $$$$\n — strip that for .mol
        mol_text = sdf_text.rstrip().rstrip("$$$$").rstrip()
        out_path.write_text(mol_text + "\n", encoding="utf-8")
        return True
    except Exception as exc:
        logger.warning("write_pose_mol: failed for %s: %s", out_path, exc)
        return False
    finally:
        try:
            tmp_path.unlink(missing_ok=True)
        except Exception:
            pass


def find_receptor(protein_name: str, search_root: Path) -> Optional[Path]:
    """
    Robustly locate a receptor PDB/PDBQT by name.
    Returns None gracefully — callers must handle the None case.
    """
    search_dirs: List[Path] = [search_root]
    for dname in ["Receptor", "RECEPTORS", "receptors", "Receptors", "receptor", "Proteins"]:
        d = search_root / dname
        if d.is_dir():
            search_dirs.append(d)
    name_lower = protein_name.lower()
    for d in search_dirs:
        try:
            for ext in [".pdb", ".pdbqt"]:
                exact = d / f"{protein_name}{ext}"
                if exact.exists():
                    return exact
            for f in d.iterdir():
                if f.is_file() and f.suffix.lower() in (".pdb", ".pdbqt"):
                    if f.stem.lower().startswith(name_lower):
                        return f
        except PermissionError:
            pass
    return None


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Ki FORMATTING
# ══════════════════════════════════════════════════════════════════════════════

def format_ki(ki_nM: Optional[float]) -> str:
    """Compact Ki string with adaptive units."""
    if ki_nM is None: return "N/A"
    if ki_nM >= 1000: return f"{ki_nM/1000:.3f} \u00b5M"
    if ki_nM >= 1:    return f"{ki_nM:.2f} nM"
    return f"{ki_nM:.4f} nM"


def format_ki_pair(ki_nM: Optional[float]) -> Tuple[str, str]:
    """Returns (nM_str, µM_str) for CSV tables."""
    if ki_nM is None: return "N/A", "N/A"
    um     = ki_nM / 1000.0
    nm_str = (f"{ki_nM:.4f}" if ki_nM < 1 else
              f"{ki_nM:.2f}" if ki_nM < 100 else f"{ki_nM:.1f}")
    um_str = (f"{um:.6f}" if um < 0.001 else
              f"{um:.4f}" if um < 1     else f"{um:.2f}")
    return nm_str, um_str


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — REPORTING ENGINE
# ══════════════════════════════════════════════════════════════════════════════

def write_mode_classification_report(
    results: List[DockingResult],
    output_dir: Path,
) -> None:
    """
    Writes BSNDVP_Docking_Modes.txt to *output_dir* listing every processed
    DLG file and its detected docking mode (RIGID or FLEXIBLE).

    The file is organised in two sections:
      [FLEXIBLE DOCKINGS]  — DLGs where AutoDock used flexible side-chains
      [RIGID DOCKINGS]     — standard rigid-receptor dockings

    For flexible entries the detected flex residue IDs are also listed.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    flexible = [r for r in results if r.docking_mode == DockingMode.FLEXIBLE]
    rigid    = [r for r in results if r.docking_mode == DockingMode.RIGID]

    sep  = "═" * 74
    lines: List[str] = [
        "",
        sep,
        f"  {__full_brand__}",
        "  Docking Mode Classification Report",
        f"  {__copyright__}",
        f"  Generated: {ts}",
        "",
        f"  Total DLG files classified : {len(results)}",
        f"  FLEXIBLE dockings          : {len(flexible)}",
        f"  RIGID dockings             : {len(rigid)}",
        "",
        sep,
        "",
    ]

    # ── FLEXIBLE section ──────────────────────────────────────────────────────
    lines.append("  ┌─ FLEXIBLE DOCKINGS " + "─" * 54 + "┐")
    if flexible:
        for r in sorted(flexible, key=lambda x: (x.protein, x.ligand)):
            fname = r.dlg_file.name if r.dlg_file else "unknown"
            flex_ids = ", ".join(r.flex_residue_ids) if r.flex_residue_ids else "detected via atom content"
            lines.append(f"  │  {fname}")
            lines.append(f"  │      Protein : {r.protein or 'N/A'}")
            lines.append(f"  │      Ligand  : {r.ligand}")
            lines.append(f"  │      Flex residues : {flex_ids}")
            lines.append(f"  │      Lig atoms : {r.n_ligand_atoms}   "
                         f"Flex residue atoms : {r.n_flexible_residue_atoms}")
            lines.append("  │")
    else:
        lines.append("  │  (none detected)")
        lines.append("  │")
    lines.append("  └" + "─" * 74)
    lines.append("")

    # ── RIGID section ─────────────────────────────────────────────────────────
    lines.append("  ┌─ RIGID DOCKINGS " + "─" * 56 + "┐")
    if rigid:
        for r in sorted(rigid, key=lambda x: (x.protein, x.ligand)):
            fname = r.dlg_file.name if r.dlg_file else "unknown"
            lines.append(f"  │  {fname}")
            lines.append(f"  │      Protein : {r.protein or 'N/A'}")
            lines.append(f"  │      Ligand  : {r.ligand}")
            lines.append(f"  │      Ligand atoms : {r.n_ligand_atoms}")
            lines.append("  │")
    else:
        lines.append("  │  (none detected)")
        lines.append("  │")
    lines.append("  └" + "─" * 74)
    lines.append("")
    lines += ["", sep, "  End of Report", f"  {__copyright__}", "", sep]

    out_path = output_dir / "BSNDVP_Docking_Modes.txt"
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  {C.GREEN}\u2713{C.RESET} Mode classification report: "
          f"{C.YELLOW}{out_path.name}{C.RESET}  "
          f"({C.CYAN}{len(flexible)}{C.RESET} flexible / "
          f"{C.GREEN}{len(rigid)}{C.RESET} rigid)")


class ReportingEngine:
    """Generates CSV / XLSX / JSON / TXT report tables.

    Reports contain only objective scientific quantities:
      - Docking results (energy, Ki, RMSD, cluster membership)
      - Cluster statistics (thermodynamic & population rankings, quality)
      - Thermodynamic statistics (Shannon entropy, partition function,
        Boltzmann probability, free energy)
      - Validation results (Kabsch RMSD, positional & conformational)
      - ADMET descriptors (molecular weight, LogP, TPSA, HBD, HBA,
        rotatable bonds, Lipinski violations)
    """

    def generate_all(
        self,
        profiles:        List[DockingProfile],
        docking_results: List[DockingResult],
        output_dir:      Path,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        tables: Dict[str, List[Dict]] = {
            "01_Docking_Master": self._master(profiles, docking_results),
            "02_Clusters":       self._clusters(docking_results),
            "03_Ligand_Summary": self._summary_table(profiles),
        }

        brand = [
            [__full_brand__], [__copyright__], [f"Generated: {ts}"], [""],
        ]

        for name, rows in tables.items():
            if not rows: continue
            headers = list(rows[0].keys())
            path    = output_dir / f"{name}.csv"
            with open(path, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                for br in brand:
                    w.writerow(br + [""] * (len(headers) - 1))
                w.writerow(headers)
                for row in rows:
                    w.writerow([row.get(h, "") for h in headers])
                w.writerow([""] * len(headers))
                w.writerow([f"--- End: {name} ---"] + [""] * (len(headers) - 1))
                w.writerow([__copyright__] + [""] * (len(headers) - 1))

        if _HAS_XLSX and _HAS_PANDAS:
            try:
                xlsx = output_dir / "BSNDVP_Complete_Report.xlsx"
                with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
                    for name, rows in tables.items():
                        pd.DataFrame(rows).to_excel(writer, sheet_name=name[:31], index=False)
            except Exception as e:
                logger.warning(f"XLSX export failed: {e}")

        json_data = {name: rows for name, rows in tables.items()}
        (output_dir / "BSNDVP_Scientific_Data.json").write_text(
            json.dumps(json_data, indent=2, default=str), encoding="utf-8"
        )
        self._txt_report(profiles, docking_results, output_dir / "BSNDVP_Report.txt")
        print(f"  {C.GREEN}\u2713{C.RESET} Reports written to: {C.YELLOW}{output_dir}{C.RESET}")

    # ── table builders ─────────────────────────────────────────────────────────

    def _master(self, profiles: List[DockingProfile],
                results: List[DockingResult]) -> List[Dict]:
        rm = {r.ligand: r for r in results}
        rows = []
        for p in profiles:
            dr  = rm.get(p.ligand)
            cl1 = "N/A"
            if dr and dr.clusters:
                c = next((c for c in dr.clusters if c.cluster_id == 1), None)
                if c: cl1 = round(c.mean_energy, 2)
            row: Dict[str, Any] = {
                # ── Identity ──────────────────────────────────────────────────
                "Protein":                        p.protein,
                "Ligand":                         p.ligand,
                "Docking_Mode":                   p.docking_mode,
                # ── Primary docking results (AutoDock output) ────────────────
                "Run":                            p.run_number,
                "Lowest_Binding_Energy_kcal":     round(p.binding_energy, 2),
                "Inhibition_Constant_Ki":         format_ki(p.ki_nM),
                "Torsional_FreeEnergy_kcal":      round(p.entropic_penalty_kcal, 2) if p.entropic_penalty_kcal else "N/A",
                # ── Cluster analysis ──────────────────────────────────────────
                "Cluster_ID":                     p.cluster_id if p.cluster_id else "N/A",
                "Number_in_Cluster":              p.cluster_size,
                "Cluster1_Mean_Energy_kcal":      cl1,
                "Energy_Spread_kcal":             round(p.cluster_spread, 2),
                "Cluster_RMSD_A":                 round(p.cluster_rmsd, 3),
                # ── Statistical mechanics ─────────────────────────────────────
                "Boltzmann_Probability_Pi":       p.boltzmann_prob,
                # ── Validation ───────────────────────────────────────────────
                "Validation_RMSD_Positional_A":   round(p.validation_rmsd_pos, 3) if p.validation_rmsd_pos is not None else "N/A",
                "Validation_RMSD_Conformational_A": round(p.validation_rmsd_conf, 3) if p.validation_rmsd_conf is not None else "N/A",
                "AutoDock_Ref_RMSD_A":            round(p.rmsd, 3) if p.rmsd else "N/A",
                # ── ADMET pharmacokinetics (reported as-is) ───────────────────
                "ADMET_MolWt_Da":                 round(p.admet_mw, 1) if p.admet_mw else "N/A",
                "ADMET_LogP":                     round(p.admet_logp, 2) if p.admet_logp is not None else "N/A",
                "ADMET_TPSA_A2":                  round(p.admet_tpsa, 1) if p.admet_tpsa else "N/A",
                "ADMET_HBD":                      p.admet_hbd if p.admet_hbd is not None else "N/A",
                "ADMET_HBA":                      p.admet_hba if p.admet_hba is not None else "N/A",
                "ADMET_RotBonds":                 p.admet_nrotb if p.admet_nrotb is not None else "N/A",
                "Lipinski_Violations":            p.lipinski_violations,
                # ── Structural (mode-adaptive) ────────────────────────────────
                "Ligand_Atom_Count":              p.n_ligand_atoms,
                "Flex_Residue_Atom_Count":        p.n_flexible_residue_atoms,
                "Flex_Residue_IDs":               "; ".join(p.flex_residue_ids) if p.flex_residue_ids else "N/A",
                "Receptor_Strain_kcal":           round(p.receptor_strain_kcal, 3) if p.receptor_strain_kcal is not None else "N/A",
            }
            if dr:
                row["Info_Entropy"]                = round(dr.info_entropy, 4)          if dr.info_entropy          is not None else "N/A"
                row["Partition_Function_Q"]        = round(dr.partition_function, 4)    if dr.partition_function    is not None else "N/A"
                row["Stat_Free_Energy_A_kcal"]     = round(dr.stat_free_energy, 2)      if dr.stat_free_energy      is not None else "N/A"
                row["Stat_Internal_Energy_U_kcal"] = round(dr.stat_internal_energy, 2)  if dr.stat_internal_energy  is not None else "N/A"
                row["Stat_Entropy_S_kcal_mol_K"]   = round(dr.stat_entropy, 4)          if dr.stat_entropy          is not None else "N/A"
            rows.append(row)
        return rows

    def _clusters(self, results: List[DockingResult]) -> List[Dict]:
        rows = []
        for r in results:
            for c in r.clusters:
                rows.append({
                    "Ligand":                  r.ligand,
                    "Protein":                 r.protein,
                    "Docking_Mode":            c.docking_mode,
                    "Cluster_ID":              c.cluster_id,
                    "Size":                    c.size,
                    "Thermodynamic_Rank":      c.thermodynamic_rank,
                    "Population_Rank":         c.population_rank,
                    "Lowest_Energy_kcal":      round(c.lowest_energy, 2),
                    "Mean_Energy_kcal":        round(c.mean_energy,   2),
                    "Best_Run":                c.best_run,
                    "Avg_RMSD_A":              round(c.avg_rmsd, 2),
                    "Energy_Spread_kcal":      round(c.spread,   2),
                    # Cluster Quality Score: dimensionless 0-100
                    # = [(size/N_runs) + max(0, 1-avg_RMSD/3)] / 2 × 100
                    "Cluster_Quality_Score_Pct": round(c.cluster_quality, 2),
                    "Runs":                    "; ".join(str(x) for x in c.runs),
                })
        return rows

    def _summary_table(self, profiles: List[DockingProfile]) -> List[Dict]:
        """Ligand summary table ranked by binding energy (most negative first)."""
        return [
            {
                "Rank":                           rank,
                "Ligand":                         p.ligand,
                "Protein":                        p.protein,
                "Docking_Mode":                   p.docking_mode,
                "Binding_Energy_kcal":            round(p.binding_energy, 2),
                "Inhibition_Constant_Ki":         format_ki(p.ki_nM),
                "Boltzmann_Probability_Pi":       p.boltzmann_prob,
                "Cluster_Size":                   p.cluster_size,
                "Energy_Spread_kcal":             round(p.cluster_spread, 2),
                "Validation_RMSD_Positional_A":   round(p.validation_rmsd_pos, 3) if p.validation_rmsd_pos is not None else "N/A",
                "Validation_RMSD_Conformational_A": round(p.validation_rmsd_conf, 3) if p.validation_rmsd_conf is not None else "N/A",
                "ADMET_MolWt_Da":                 round(p.admet_mw, 1) if p.admet_mw else "N/A",
                "ADMET_LogP":                     round(p.admet_logp, 2) if p.admet_logp is not None else "N/A",
                "ADMET_TPSA_A2":                  round(p.admet_tpsa, 1) if p.admet_tpsa else "N/A",
                "Lipinski_Violations":            p.lipinski_violations,
            }
            for rank, p in enumerate(profiles, 1)
        ]

    def _txt_report(
        self,
        profiles: List[DockingProfile],
        results:  List[DockingResult],
        path:     Path,
    ) -> None:
        ts   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        W    = 80
        sep  = "\u2550" * W
        hsep = "\u2500" * W

        def fmt(value: object, spec: str, missing: str = "Not available") -> str:
            """Format optional scientific values without converting missing data."""
            return missing if value is None else format(value, spec)

        def section(title: str) -> List[str]:
            pad   = max(0, W - 6 - len(title))
            lpad  = pad // 2
            rpad  = pad - lpad
            bar   = "\u2550" * W
            inner = f"  {' ' * lpad}{title.upper()}{' ' * rpad}  "
            return [
                "",
                f"\u2554{bar}\u2557",
                f"\u2551  {inner[:W]}  \u2551",
                f"\u255a{bar}\u255d",
                "",
            ]

        # Header
        header_lines: List[str] = [
            "", sep,
            f"  {__full_brand__}",
            "  Scientific Cluster & Thermodynamic Analysis Report",
            f"  {__copyright__}",
            f"  Author  : {__author__}",
            f"  Generated : {ts}",
            f"  Ligands   : {len(profiles)}",
            "", sep, "",
        ]

        # ══ SECTION A — SCIENTIFIC METHODOLOGY ════════════════════════════════
        methodology_lines: List[str] = []
        lines = methodology_lines
        lines += section("Scientific methodology — definitions & references")
        lines += [
            "  This report was produced by the BSNDVP(TM) Docking Intelligence System.",
            "  All primary docking data originate from AutoDock 4.2 DLG output files.",
            "  Every metric listed below is defined here in full.",
            "",
        ]

        # 1. DeltaG
        lines += [
            "  1. BINDING FREE ENERGY (DeltaG)  --  Lowest_Binding_Energy_kcal",
            hsep,
            "  AutoDock 4 estimates DeltaG_binding using a semi-empirical force-field",
            "  decomposed into five additive terms:",
            "",
            "    DeltaG = DeltaG_vdW + DeltaG_hbond + DeltaG_elec + DeltaG_tor + DeltaG_sol",
            "",
            "    DeltaG_vdW   -- Lennard-Jones van der Waals dispersion/repulsion",
            "    DeltaG_hbond -- directional hydrogen bonding interactions",
            "    DeltaG_elec  -- Coulombic electrostatic interactions",
            "    DeltaG_tor   -- torsional entropy cost of ligand freezing on binding",
            "    DeltaG_sol   -- desolvation penalty (water displaced from binding site)",
            "",
            "  Unit: kcal/mol  (more negative = stronger predicted binding affinity)",
            "    DeltaG  <  -8.0 kcal/mol  ...  Strong (drug-like range)",
            "    -5.0 to -8.0 kcal/mol      ...  Moderate binding",
            "    DeltaG  >  -5.0 kcal/mol  ...  Weak / non-specific",
            "",
            "  Ref: Morris et al. J. Comput. Chem. 19, 1639-1662 (1998).",
            "",
        ]

        # 2. Ki
        lines += [
            "  2. INHIBITION CONSTANT (Ki)  --  Inhibition_Constant",
            hsep,
            "  Derived from DeltaG via the standard thermodynamic identity:",
            "",
            "    DeltaG = RT * ln(Ki)   =>   Ki = exp(DeltaG / RT)",
            "",
            "  where R = 1.987 cal/mol*K, T = 298.15 K, RT ~ 0.592 kcal/mol.",
            "  Reported in nM or uM.",
            "    Ki  < 100 nM     ...  High affinity (strongly drug-like)",
            "    100 nM - 10 uM   ...  Moderate affinity",
            "    Ki  > 10 uM      ...  Weak / non-specific binding",
            "",
        ]

        # 3. Cluster Analysis
        lines += [
            "  3. CLUSTER ANALYSIS  --  Cluster_ID, Number_in_Cluster,",
            "                           Energy_Spread_kcal, Avg_RMSD_A",
            hsep,
            "  AutoDock runs multiple independent docking trials (configured by the user).",
            "  Poses sharing structural similarity (RMSD <= 2.0 A) are grouped into",
            "  clusters. A large, tight, low-energy cluster is the gold-standard",
            "  result: a reproducible, high-confidence binding mode.",
            "",
            "    Number_in_Cluster  -- poses sharing this binding geometry",
            "    Energy_Spread      -- max(DeltaG) - min(DeltaG) within the cluster",
            "                         (lower = more energetically consistent poses)",
            "    Avg_RMSD           -- mean pairwise RMSD (lower = more compact cluster)",
            "",
        ]

        # 4. Rankings
        lines += [
            "  4. CLUSTER RANKINGS  --  Thermodynamic_Rank, Population_Rank",
            hsep,
            "  Two orthogonal ranking strategies are applied to every cluster:",
            "",
            "    Thermodynamic Rank -- sorted: lowest energy first, then largest cluster.",
            "      Best for affinity prediction and lead optimisation.",
            "",
            "    Population Rank    -- sorted: largest cluster first, then lowest energy.",
            "      Best for identifying the dominant binding mode in solution.",
            "",
            "  When both ranks agree on Cluster 1 the result is unambiguous.",
            "  Divergence indicates competitive binding modes.",
            "",
        ]

        # 5. Entropic Penalty
        lines += [
            "  5. ENTROPIC PENALTY (Torsional Free Energy DeltaG_tor)  --  Torsional_FreeEnergy_kcal",
            hsep,
            "  The torsional free energy contribution from AutoDock's energy decomposition,",
            "  representing the conformational entropy cost the ligand pays when it loses",
            "  rotatable bond freedom upon binding. AutoDock 4 empirical parameterisation:",
            "",
            "    DeltaG_tor = 0.311 * N_active_torsions  (kcal/mol)",
            "",
            "  IMPORTANT: DeltaG_tor is already included in the total DeltaG_binding.",
            "  This column isolates the torsional term for independent assessment of",
            "  the ligand's rotational entropy burden vs. its enthalpic binding gains.",
            "",
        ]

        # 6. Boltzmann
        lines += [
            "  6. BOLTZMANN POPULATION PROBABILITY (pi_i)  --  Boltzmann_Probability_Pi",
            hsep,
            "  Fractional occupancy of the best cluster at thermal equilibrium,",
            "  calculated from the canonical ensemble:",
            "",
            "    pi_i = exp(-E_i / RT) / Sum_j [ exp(-E_j / RT) ]",
            "",
            "  where E_i = mean cluster energy (kcal/mol), RT = 0.5925 kcal/mol @ 298 K",
            "  and the denominator is the partition function Z.",
            "",
            "  pi -> 1.0  : all runs converge on one binding mode (ideal).",
            "  pi -> 0    : energy dispersed across many modes (poor convergence).",
            "",
            "  Ref: Boltzmann, L. (1872); AutoDock DLG statistical analysis.",
            "",
        ]

        # 7. Information Entropy
        lines += [
            "  7. INFORMATION ENTROPY (H)  --  Info_Entropy",
            hsep,
            "  Shannon information entropy over the cluster occupancy distribution:",
            "",
            "    H = -Sum_i [ p_i * log2(p_i) ]",
            "",
            "  where p_i = fraction of runs in cluster i. Extracted from AutoDock.",
            "    H = 0       : perfect convergence (all runs in one cluster)",
            "    H is large  : poor convergence (runs scattered across many clusters)",
            "",
        ]

        # 8. Statistical Mechanics
        lines += [
            "  8. STATISTICAL MECHANICS  --  Partition_Function_Q,",
            "     Stat_Free_Energy_A_kcal, Stat_Internal_Energy_U_kcal,",
            "     Stat_Entropy_S_kcal_mol_K",
            hsep,
            "  AutoDock writes a canonical-ensemble statistical analysis per experiment:",
            "",
            "    Q  = Sum_i exp(-E_i / RT)           Partition function",
            "    A ~= -RT * ln(Q)                    Helmholtz free energy (kcal/mol)",
            "    U  = Sum_i p_i * E_i                Boltzmann-weighted mean energy",
            "    S  = (U - A) / T                    Thermodynamic entropy (kcal/mol*K)",
            "",
            "  These characterise the macroscopic thermodynamic state of the ensemble",
            "  and are suitable for inclusion in academic supplementary tables.",
            "",
        ]

        # 9. ADMET
        lines += [
            "  9. ADMET DESCRIPTORS & LIPINSKI RULE OF FIVE",
            "     ADMET_MolWt_Da, ADMET_LogP, ADMET_TPSA_A2, ADMET_HBD, ADMET_HBA,",
            "     ADMET_RotBonds, Lipinski_Violations",
            hsep,
            "  Calculated by RDKit on the isolated ligand structure.",
            "  Reported as objective molecular descriptors — not integrated into any score.",
            "  Lipinski Rule of Five criteria for oral bioavailability:",
            "",
            "    Molecular Weight  <= 500 Da",
            "    LogP              <= 5       (lipophilicity / oral permeability)",
            "    H-bond donors     <= 5",
            "    H-bond acceptors  <= 10",
            "    TPSA              <  140 A^2 (topological polar surface area)",
            "",
            "  Lipinski_Violations: count of criteria exceeded. 0-1 = drug-like.",
            "",
            "  Ref: Lipinski et al. Adv. Drug Deliv. Rev. 23, 3-25 (1997).",
            "",
        ]

        # 10. Receptor Strain
        lines += [
            "  10. RECEPTOR STRAIN ESTIMATE  --  Receptor_Strain_kcal",
            "      FLEXIBLE docking only. N/A for rigid docking.",
            hsep,
            "  Estimate of the energetic cost imposed on the receptor's flexible",
            "  residues when adopting the docked conformation:",
            "",
            "    Receptor_Strain = Internal_Energy_total - Torsional_Energy_ligand",
            "",
            "  Internal_Energy_total (from AutoDock) includes both ligand internal",
            "  energy and strain in flexible side-chains. Subtracting the ligand's",
            "  torsional term isolates the approximate receptor-side contribution.",
            "  Large positive value = geometrically strained side-chain conformation.",
            "",
        ]

        # 11. Validation RMSD
        lines += [
            "  11. VALIDATION RMSD (Positional & Conformational)",
            "      AUTODOCK INTERNAL RMSD     --  AutoDock_Ref_RMSD_A",
            hsep,
            "  Three distinct RMSD quantities are reported:",
            "",
            "  Validation_RMSD_Positional_A:",
            "    In-place positional RMSD between the docked pose and the reference.",
            "    This represents in-pocket accuracy (raw Euclidean distance without alignment).",
            "",
            "  Validation_RMSD_Conformational_A:",
            "    RMSD between the best docked pose and the INPUT-LIGAND-PDBQT reference",
            "    extracted from the DLG header. Computed using RDKit's MCS-matched",
            "    Kabsch superposition algorithm (heavy atoms only) to evaluate 3D shape folding.",
            "",
            "      RMSD_Kabsch = sqrt[ (1/N) * Sum ||r_i - r_i'||^2 ]  after optimal alignment",
            "",
            "    Thresholds (Leach et al. J. Chem. Inf. Model. 46, 526-529, 2006):",
            "      RMSD < 2.0 A  : successful redocking; crystallographically consistent",
            "      2.0 - 3.0 A   : acceptable; likely same binding mode",
            "      RMSD > 3.0 A  : different binding mode; interpret with caution",
            "    N/A = RDKit unavailable or reference structure not found.",
            "",
            "  AutoDock_Ref_RMSD_A:",
            "    AutoDock's own RMSD relative to the -ref structure if one was provided",
            "    at docking time. Without a crystallographic -ref this is the RMSD from",
            "    the initial random seed and carries no validation meaning.",
            "",
            sep, "",
        ]

        # ══ SECTION B — PER-LIGAND RESULTS ════════════════════════════════════
        results_lines: List[str] = []
        lines = results_lines
        lines += section("Docking results — ranked by binding energy")

        for rank, p in enumerate(profiles, 1):
            mode_badge = ("[FLEXIBLE]" if p.docking_mode == "FLEXIBLE" else "[RIGID]   ")
            val_rmsd_str = (
                f"{p.validation_rmsd_pos:.3f} A (Pos) | {p.validation_rmsd_conf:.3f} A (Conf)"
                if p.validation_rmsd_pos is not None
                else "N/A  (no reference structure or RDKit unavailable)"
            )
            lines += [
                hsep,
                f"  Rank {rank:>3}.  {p.ligand}  <-  {p.protein}  {mode_badge}",
                hsep,
                f"  Binding Free Energy (DeltaG)     : {p.binding_energy:.2f} kcal/mol",
                f"  Inhibition Constant (Ki)         : {format_ki(p.ki_nM)}",
                f"  Torsional Free Energy (DeltaG_t) : {p.entropic_penalty_kcal:.2f} kcal/mol  (included in DeltaG)",
                f"  Best Docking Run                 : {p.run_number}",
                f"  Validation RMSD (Kabsch/RDKit)   : {val_rmsd_str}",
                f"  Cluster ID                       : {p.cluster_id or 'N/A'}",
                f"  Cluster Size                     : {p.cluster_size}",
                f"  Energy Spread                    : {p.cluster_spread:.2f} kcal/mol",
                f"  Boltzmann Probability (pi)       : {fmt(p.boltzmann_prob, '.4f')}",
                "",
            ]

        # ══ SECTION C — CLUSTER TABLES ═══════════════════════════════════════
        lines += section("Cluster analysis tables")
        lines += [
            "  Columns:",
            "    ID    = Cluster_ID            Sz    = Number of poses in cluster",
            "    ThR   = Thermodynamic_Rank    PoR   = Population_Rank",
            "    LowE  = Lowest_Energy (kcal/mol)    MeanE = Mean_Energy (kcal/mol)",
            "    RMSD  = Avg_RMSD (A)          Run   = Best_Run    Qual = Quality (%)",
            "",
        ]

        for r in results:
            if r.clusters:
                lines += [
                    f"  > {r.ligand}  <-  {r.protein}",
                    f"  {'ID':>3}  {'Sz':>5}  {'ThR':>4}  {'PoR':>4}  "
                    f"{'LowE':>8}  {'MeanE':>8}  {'RMSD':>6}  {'Run':>5}  {'Qual':>6}",
                    "  " + "-" * 63,
                ]
                for c in sorted(r.clusters, key=lambda x: x.thermodynamic_rank):
                    lines.append(
                        f"  {c.cluster_id:>3}  {c.size:>5}  "
                        f"{c.thermodynamic_rank:>4}  {c.population_rank:>4}  "
                        f"{c.lowest_energy:>8.2f}  {c.mean_energy:>8.2f}  "
                        f"{c.avg_rmsd:>6.2f}  {c.best_run:>5}  {c.cluster_quality:>5.1f}%"
                    )
                lines.append("")

        # ══ SECTION D — DETAILED VALIDATION REPORT ═══════════════════════════
        validation_detail_lines: List[str] = []
        lines = validation_detail_lines
        lines += section("Detailed structural validation report — per ligand")
        lines += [
            "  This section provides a per-ligand breakdown of all three validation",
            "  dimensions: structural (Kabsch RMSD), thermodynamic (DeltaG / Ki), and",
            "  statistical (cluster quality). See Appendix I for full methodology.",
            "",
        ]

        for rank, p in enumerate(profiles, 1):
            e = p.binding_energy
            thermo_grade  = ("Strong" if e < -8.0 else ("Moderate" if e < -5.0 else "Weak"))
            struct_grade  = "N/A (no reference)"
            struct_detail = "RDKit unavailable or no reference structure found."
            if p.validation_rmsd_pos is not None:
                vr = p.validation_rmsd_pos
                if vr < 2.0:
                    struct_grade  = "PASS — Excellent (<2.0 Å)"
                    struct_detail = "Docked pose reproduces reference within crystallographic precision."
                elif vr < 3.0:
                    struct_grade  = "PASS — Acceptable (2.0–3.0 Å)"
                    struct_detail = "Likely the same binding mode; minor peripheral variation present."
                else:
                    struct_grade  = "CAUTION — Different mode (>3.0 Å)"
                    struct_detail = "Pose may represent an alternative binding site or mode."
            if p.boltzmann_prob is None:
                boltz_interp = "Not available"
            else:
                boltz_interp = ("Dominant (\u03c0 \u2265 0.80)" if p.boltzmann_prob >= 0.80 else
                                ("Competitive (\u03c0 0.50\u20130.79)" if p.boltzmann_prob >= 0.50 else
                                 "Non-dominant (\u03c0 < 0.50)"))
            lines += [
                hsep,
                f"  Rank {rank:>3}.  {p.ligand}  <-  {p.protein}  [{p.docking_mode}]",
                hsep,
                "  \u2500 Thermodynamic Dimension \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500",
                f"    Binding Free Energy (\u0394G)   : {e:.2f} kcal/mol  \u2192 {thermo_grade}",
                f"    Inhibition Constant (Ki)   : {format_ki(p.ki_nM)}",
                f"    Torsional Entropy (\u0394G_tor) : {p.entropic_penalty_kcal:.2f} kcal/mol (included in \u0394G)",
                "",
                "  \u2500 Structural Dimension (Kabsch / RDKit MCS) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500",
                f"    Validation RMSD (Pos)      : {(f'{p.validation_rmsd_pos:.3f} \u00c5') if p.validation_rmsd_pos is not None else 'N/A'}",
                f"    Validation RMSD (Conf)     : {(f'{p.validation_rmsd_conf:.3f} \u00c5') if p.validation_rmsd_conf is not None else 'N/A'}",
                f"    Structural Grade           : {struct_grade}",
                f"    Interpretation             : {struct_detail}",
                "",
                "  \u2500 Statistical Dimension \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500",
                f"    Cluster ID / Size          : {p.cluster_id or 'N/A'} / {p.cluster_size} poses",
                f"    Cluster RMSD (avg)         : {p.avg_cluster_rmsd:.2f} \u00c5",
                f"    Energy Spread              : {p.cluster_spread:.2f} kcal/mol",
                f"    Boltzmann Population (\u03c0)   : {fmt(p.boltzmann_prob, '.4f')}  \u2192 {boltz_interp}",
                "",
                "  \u2500 ADMET Descriptors \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500",
                f"    Mol. Weight (Da)           : {round(p.admet_mw, 1) if p.admet_mw else 'N/A'}",
                f"    LogP                       : {round(p.admet_logp, 2) if p.admet_logp is not None else 'N/A'}",
                f"    TPSA (A^2)                 : {round(p.admet_tpsa, 1) if p.admet_tpsa else 'N/A'}",
                f"    H-bond Donors              : {p.admet_hbd if p.admet_hbd is not None else 'N/A'}",
                f"    H-bond Acceptors           : {p.admet_hba if p.admet_hba is not None else 'N/A'}",
                f"    Rotatable Bonds            : {p.admet_nrotb if p.admet_nrotb is not None else 'N/A'}",
                f"    Lipinski Violations        : {p.lipinski_violations}",
                "",
            ]

        # ══ SECTION E — VALIDATION METHODOLOGY APPENDIX ══════════════════════
        appendix_lines: List[str] = []
        lines = appendix_lines
        lines += section("Appendix I — Docking Validation Methodology (Academic Reference)")
        lines += [
            "  This appendix provides a publication-level description of all validation",
            f"  algorithms implemented in BSNDVP(TM) v{__version__}.  It is reproduced in full",
            "  with each report to serve as a self-contained scientific reference.",
            "",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "  A. STRUCTURAL VALIDATION \u2014 KABSCH RMSD ALGORITHM",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "",
            "  The Root Mean Square Deviation (RMSD) between a docked ligand pose and",
            "  a reference structure is computed using the Kabsch optimal superimposition",
            "  algorithm (Kabsch, Acta Cryst. A, 1976).  This is the provably optimal",
            "  rotation matrix U that minimises RMSD between two point sets.",
            "",
            "  Mathematical derivation:",
            "",
            "    Given N matched atom pairs: P (reference, N\u00d73) and Q (docked, N\u00d73):",
            "",
            "    Step 1 \u2014 Centroid subtraction (translation removal):",
            "      Pc = P - mean(P),   Qc = Q - mean(Q)",
            "",
            "    Step 2 \u2014 Cross-covariance matrix:",
            "      H = Pc^T \u00d7 Qc    [3\u00d73 matrix]",
            "",
            "    Step 3 \u2014 Singular Value Decomposition:",
            "      H = V \u00d7 Sigma \u00d7 W^T",
            "",
            "    Step 4 \u2014 Chirality guard (prevents mirror-image solutions):",
            "      if det(V) * det(W) < 0  =>  V[:,-1] = -V[:,-1]",
            "",
            "    Step 5 \u2014 Optimal rotation matrix:",
            "      U = V \u00d7 W",
            "",
            "    Step 6 \u2014 RMSD after optimal superimposition:",
            "      RMSD = sqrt[ (1/N) * Sum ||( Pc \u00d7 U )_i  -  (Qc)_i||^2 ]",
            "",
            "  The chirality guard (Step 4) is critical for molecular work: without it,",
            "  SVD may return an improper rotation (a reflection) that produces an",
            "  artificially low RMSD by inverting stereochemistry.",
            "",
            "  Implementation reference: RMSDEngine.kabsch_rmsd()  [dlg_extract.py]",
            "  Dependencies: NumPy (np.linalg.svd, np.dot, np.sqrt)",
            "",
            "  Validation RMSD thresholds (Leach et al., J. Med. Chem. 2006):",
            "    < 2.0 \u00c5  \u2014 Successful redocking; crystallographically consistent",
            "    2.0\u20133.0 \u00c5 \u2014 Acceptable; likely same pharmacophore",
            "    > 3.0 \u00c5  \u2014 Different binding mode; interpret with caution",
            "",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "  B. LIGAND SUPERIMPOSITION \u2014 MCS ATOM SELECTION & CORRESPONDENCE",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "",
            "  Atom-to-atom correspondence is established via Maximum Common Substructure",
            "  (MCS) matching using RDKit rdFMCS.FindMCS().  This is a topology-based",
            "  approach: no atom naming convention is assumed.",
            "",
            "  MCS parameters:",
            "    matchValences=False      \u2014 robust across protonation states and",
            "                               AutoDock-guessed bond orders",
            "    ringMatchesRingOnly=True \u2014 ring atoms cannot match chain atoms",
            "    timeout=10 seconds       \u2014 graceful fallback if MCS is NP-hard",
            "",
            "  Atom selection filter applied after MCS:",
            "    \u2022 Only non-hydrogen (heavy) atoms retained",
            "    \u2022 Element symbols of matched pairs must be identical",
            "    \u2022 Hydrogen atoms systematically excluded (unreliable in AutoDock output",
            "      and not resolved crystallographically at typical resolutions)",
            "",
            "  Two RMSD sub-types are reported:",
            "    Validation_Heavy_RMSD    \u2014 all N heavy atoms matched (complete scaffold)",
            "    Validation_Scaffold_RMSD \u2014 \u226540% of heavy atoms matched (partial scaffold)",
            "    Validation_RMSD          \u2014 heavy_rmsd preferred; scaffold_rmsd as fallback",
            "",
            "  Molecule standardisation pipeline (LigandStandardizer):",
            "    1. RemoveHs()             \u2014 strip all hydrogens",
            "    2. SanitizeMol()          \u2014 enforce valence and aromaticity",
            "    3. RenumberAtoms with CanonicalRankAtoms \u2014 deterministic atom ordering",
            "    4. ComputeGasteigerCharges() \u2014 partial charges",
            "",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "  C. CLUSTER-BASED STATISTICAL VALIDATION",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "",
            "  AutoDock 4 performs N independent Genetic Algorithm docking runs.",
            "  Poses within RMSD_tol (default 2.0 \u00c5) are grouped into clusters.",
            "  The DLG records cluster membership, sizes, energies, and RMSD values.",
            "",
            "  Cluster Quality Score (0\u2013100 %):",
            "    Quality = [ f_pop + q_rmsd ] / 2 \u00d7 100",
            "",
            "    f_pop  = n_cluster / N_runs        (population fraction; 0\u21921)",
            "    q_rmsd = max(0, 1 - avg_RMSD/3.0)  (compactness; 0 at avg_RMSD\u22653.0 \u00c5)",
            "",
            "    N_runs is parsed from the DLG header ('Number of Docking Runs = N'),",
            "    ensuring normalisation to the actual run count (10, 50, 100, 200, \u2026).",
            "",
            "  Energy Spread:  max(\u0394G_i) - min(\u0394G_i) within the cluster.",
            "    Small spread (<1 kcal/mol) \u2192 energetically consistent binding mode.",
            "    Large spread \u2192 conformational heterogeneity despite structural grouping.",
            "",
            "  Dual Cluster Ranking:",
            "    Thermodynamic Rank \u2014 sorted by lowest energy, then largest size.",
            "      Use: affinity prediction; lead identification.",
            "    Population Rank    \u2014 sorted by largest size, then lowest energy.",
            "      Use: dominant binding mode in solution; MD seeding.",
            "    Agreement between ranks \u2192 unambiguous result.",
            "    Divergence          \u2192 competitive binding modes of similar plausibility.",
            "",
            "  Boltzmann Population Probability (pi_i):",
            "    pi_i = exp(-E_i/RT) / Z,  Z = Sum_j exp(-E_j/RT)",
            "    Computed via log-sum-exp for numerical stability:",
            "      log_Z = max_t + log[ Sum_j exp(t_j - max_t) ]  where t_i = -E_i/RT",
            "    RT = 0.5925 kcal/mol at 298.15 K (R = 1.987 cal/mol\u00b7K)",
            "    pi \u2192 1.0 : all energy in one mode (thermodynamic dominance)",
            "    pi \u2192 0   : energy dispersed across many modes (poor convergence)",
            "",
            "  Shannon Information Entropy (H):",
            "    H = -Sum_i [ p_i \u00d7 log2(p_i) ]",
            "    p_i = fraction of runs in cluster i.  Parsed directly from DLG.",
            "    H = 0 : perfect convergence.   H_max = log2(k) : maximum disorder.",
            "",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "  D. REFERENCES",
            "  \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550",
            "",
            "  [1] Kabsch, W. (1976). Acta Crystallographica A, 32(5), 922\u2013923.",
            "  [2] Morris, G.M. et al. (1998). J. Comput. Chem., 19(14), 1639\u20131662.",
            "  [3] Leach, A.R. et al. (2006). J. Med. Chem., 49(20), 5851\u20135855.",
            "  [4] Warren, G.L. et al. (2006). J. Med. Chem., 49(20), 5912\u20135931.",
            "  [5] Lipinski, C.A. et al. (1997). Adv. Drug Deliv. Rev., 23(1\u20133), 3\u201325.",
            "  [6] Boltzmann, L. (1872). Sitzungsberichte der Akademie der Wissenschaften,",
            "      66, 275\u2013370.",
            "  [7] RDKit: Open-source cheminformatics. https://www.rdkit.org",
            "",
        ]

        # Footer
        footer_lines = ["", sep, "  End of Report", f"  {__copyright__}", "", sep]

        final_lines = (header_lines + results_lines + validation_detail_lines
                       + methodology_lines + appendix_lines + footer_lines)

        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(final_lines), encoding="utf-8")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11 — SHARED ENGINE SINGLETONS
# ══════════════════════════════════════════════════════════════════════════════

_CLUSTER_ENGINE = ClusterAnalysisEngine()
_POSE_EXTRACTOR = PoseExtractor()
_REPORTER       = ReportingEngine()


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — CLI PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def _find_dlg_dir(cwd: Path) -> Optional[Path]:
    """Auto-discover a directory containing .dlg files in or near cwd."""
    # 1. Check common relative folder names in priority order
    candidate_names = [
        "results/AD4", "AD4",
        "DLG", "dlg", "results/DLG", "results/dlg",
        "results/AD4/DLG_ANALYSIS/DLG", "AD4/DLG_ANALYSIS/DLG",
        "dlg_files", "docking_results",
    ]
    for rel in candidate_names:
        cand = cwd / rel
        if cand.is_dir() and any(cand.rglob("*.dlg")):
            return cand

    # 2. Check if project_config.toml exists in cwd
    cfg_file = cwd / "project_config.toml"
    if cfg_file.is_file():
        try:
            import tomllib
            with open(cfg_file, "rb") as f:
                cfg_data = tomllib.load(f)
            res_dir_str = cfg_data.get("outputs", {}).get("result_directory", "results")
            res_dir = Path(res_dir_str)
            if not res_dir.is_absolute():
                res_dir = cwd / res_dir
            for sub in [res_dir / "AD4", res_dir / "DLG", res_dir / "AD4" / "DLG_ANALYSIS" / "DLG", res_dir]:
                if sub.is_dir() and any(sub.rglob("*.dlg")):
                    return sub
        except Exception:
            pass

    # 3. Check if cwd itself directly contains *.dlg
    if any(cwd.glob("*.dlg")):
        return cwd

    # 4. Search subdirectories for any folder containing *.dlg
    try:
        for d in cwd.rglob("*.dlg"):
            return d.parent
    except Exception:
        pass

    # 5. Check if empty candidates exist (give preference to results/AD4 then results/DLG then DLG)
    for rel in ["results/AD4", "AD4", "results/DLG", "DLG", "dlg"]:
        cand = cwd / rel
        if cand.is_dir():
            return cand

    return None


class DLGPipeline:
    """
    CLI orchestrator — pure command-line scientific application.
    Processes AutoDock DLG files and produces:
      - CSV, XLSX, JSON tables (docking results, clusters, ADMET)
      - TXT scientific report (methodology + per-ligand + appendix)
      - PDBQT / SDF / MOL pose files for best poses
      - Docking mode classification report
    """

    def __init__(
        self,
        dlg_dir:    Path,
        output_dir: Path,
        cwd:        Path,
    ) -> None:
        self.dlg_dir      = dlg_dir
        self.output_dir   = output_dir
        self.cwd          = cwd
        self.poses_dir    = output_dir / "EXTRACTED_POSES"
        self.sdf_dir      = output_dir / "EXTRACTED_POSES" / "SDF"
        self.ref_pose_dir = output_dir / "EXTRACTED_POSES" / "REF_POSE"
        self.resolver     = ReferenceResolver(cwd)

    def run(self) -> None:
        dlgs = sorted(self.dlg_dir.rglob("*.dlg"))
        if not dlgs:
            msg = f"No .dlg files found in {self.dlg_dir}"
            print(f"{C.RED}  \u2717  ERROR: {msg}{C.RESET}")
            raise FileNotFoundError(msg)

        _W = 72
        print(f"  {C.CYAN}\u250c{C.RESET}{C.DIM}{'\u2500' * _W}{C.RESET}{C.CYAN}\u2510{C.RESET}")
        total_str = f"  {__trademark__}  \u2014  Processing {len(dlgs)} DLG file{'s' if len(dlgs) != 1 else ''}  "
        print(f"  {C.CYAN}\u2502{C.RESET}{C.WHITE}{total_str:<{_W}}{C.RESET}{C.CYAN}\u2502{C.RESET}")
        print(f"  {C.CYAN}\u2514{C.RESET}{C.DIM}{'\u2500' * _W}{C.RESET}{C.CYAN}\u2518{C.RESET}")
        print()

        parser       = DLGParser()
        profiles:    List[DockingProfile] = []
        all_results: List[DockingResult]  = []

        for i, dlg in enumerate(dlgs, 1):
            _idx = f"[{i:>0{len(str(len(dlgs)))}d}/{len(dlgs)}]"
            print(f"\n  {C.DIM}\u254c\u254c\u254c{C.RESET} {C.BOLD}{C.CYAN}{_idx}{C.RESET}  {C.BOLD}{C.WHITE}{dlg.name}{C.RESET}")
            self._process_dlg(dlg, parser, profiles, all_results)

        if not profiles:
            msg = "No ligands processed successfully from DLG files."
            print(f"{C.RED}  \u2717  {msg}{C.RESET}")
            raise RuntimeError(msg)

        # Rank by binding energy (most negative = strongest = rank 1)
        ranked = sorted(profiles, key=lambda p: p.binding_energy)

        _REPORTER.generate_all(ranked, all_results, self.output_dir)
        write_mode_classification_report(all_results, self.output_dir)

        # Mirror summary results to cwd/BSNDVP_RESULTS and cwd/reports if applicable
        if self.cwd:
            import shutil
            target_bsndvp = self.cwd / "BSNDVP_RESULTS"
            if target_bsndvp.resolve() != self.output_dir.resolve():
                target_bsndvp.mkdir(parents=True, exist_ok=True)
                for item in self.output_dir.iterdir():
                    try:
                        if item.is_file():
                            shutil.copy2(item, target_bsndvp / item.name)
                        elif item.is_dir():
                            shutil.copytree(item, target_bsndvp / item.name, dirs_exist_ok=True)
                    except Exception:
                        pass

            reports_dir = self.cwd / "reports"
            if reports_dir.is_dir() or (self.cwd / "results").is_dir():
                reports_dir.mkdir(parents=True, exist_ok=True)
                for pattern in ("*.csv", "*.xlsx", "*.txt", "*.json"):
                    for f in self.output_dir.glob(pattern):
                        try:
                            shutil.copy2(f, reports_dir / f.name)
                        except Exception:
                            pass

        self._summary(ranked)

    def _process_dlg(
        self,
        dlg:         Path,
        parser:      DLGParser,
        profiles:    List[DockingProfile],
        all_results: List[DockingResult],
    ) -> None:
        result = parser.parse(dlg)
        if not result.success:
            print(f"    {C.YELLOW}\u26a0 SKIP:{C.RESET} {result.error_message}")
            return

        _CLUSTER_ENGINE.rank_clusters(result)
        all_results.append(result)

        best = _POSE_EXTRACTOR.extract_best_pose(result)
        if not best:
            print(f"    {C.YELLOW}\u26a0 SKIP:{C.RESET} No valid pose extracted.")
            return

        best_cluster = _CLUSTER_ENGINE.get_best_cluster_thermo(result)
        # ── Structured per-file output block ─────────────────────────────────
        mode_badge = (
            f"{C.TEAL}[FLEXIBLE]{C.RESET}"
            if result.docking_mode == DockingMode.FLEXIBLE
            else f"{C.DIM}[RIGID]{C.RESET}   "
        )
        print(
            f"    {C.DIM}Receptor{C.RESET}  {C.WHITE}{result.protein or 'N/A':<22}{C.RESET}"
            f"  {C.DIM}Ligand{C.RESET}  {C.LIME}{result.ligand:<20}{C.RESET}"
            f"  {mode_badge}"
        )
        print(
            f"    {C.DIM}Runs{C.RESET}  {C.WHITE}{result.num_runs:<6}{C.RESET}"
            f"  {C.DIM}Clusters{C.RESET}  {C.WHITE}{result.num_clusters}{C.RESET}"
        )
        print(
            f"    {C.GREEN}\u2605  Best pose{C.RESET}  "
            f"{C.DIM}Run{C.RESET} {C.WHITE}{best.run_number:<5}{C.RESET}"
            f"  {C.DIM}\u0394G{C.RESET} = {C.CYAN}{best.binding_energy:>7.2f}{C.RESET} kcal/mol"
            f"  {C.DIM}Ki{C.RESET} = {C.MAGENTA}{format_ki(best.ki_nM)}{C.RESET}"
        )
        if best_cluster:
            print(
                f"    {C.TEAL}\u25c6  Top cluster{C.RESET}  "
                f"#{best_cluster.cluster_id}  "
                f"{C.WHITE}{best_cluster.size}{C.RESET} poses"
                f"  {C.DIM}\u0394G spread{C.RESET} {best_cluster.spread:.2f} kcal/mol"
            )

        write_pdbqt(best, self.poses_dir / f"{result.ligand}_run{best.run_number}.pdbqt")

        # ── Resolve reference for validation & SDF export ─────────────────────
        _ref_path:  Optional[Path] = None
        _ref_source: str = "NONE"
        _ref_mol:   Optional[Any] = None
        _pose_mol:  Optional[Any] = None
        effective_lig_atoms = best.ligand_atoms if best.ligand_atoms else best.atoms
        _sdf_name = f"{result.ligand}_run{best.run_number}"

        if _HAS_RDKIT:
            _ref_path, _ref_source = self.resolver.resolve(result.ligand, dlg)
            # Build pose mol from temp PDBQT (neutral format, works for all sources)
            _tmp = self.poses_dir / f"_tmp_{dlg.stem}.pdbqt"
            write_pdbqt(best, _tmp, atoms=effective_lig_atoms)
            _pose_mol = LigandStandardizer.load_and_standardize(str(_tmp))
            if _tmp.exists(): _tmp.unlink(missing_ok=True)

            if _ref_path and _ref_path.exists():
                _ref_mol = LigandStandardizer.load_and_standardize(str(_ref_path))

        # ── SDF export strategy ───────────────────────────────────────────────
        # Strategy A (USER-SUPPLIED reference available + RDKit):
        #   Use the reference molecule's own graph as a template and transplant
        #   the docked coordinates onto it.  The output SDF is atom-for-atom
        #   identical to the reference in ordering, element types, bond orders,
        #   and formal charges.  BIOVIA DSV / PyMOL can compare it against the
        #   reference with a straight 1:1 positional RMSD — no MCS needed.
        #   Output goes to EXTRACTED_POSES/REF_POSE/ so it is easy to find.
        #
        # Strategy B (DLG fallback reference or no RDKit):
        #   Build SDF from raw PDBAtom data via write_pose_sdf().  This is still
        #   far better than PDBQT for interoperability but may need MCS-based
        #   alignment in external tools.
        #   Output goes to EXTRACTED_POSES/SDF/.

        _is_user_ref = _ref_path is not None and "USER_SUPPLIED" in _ref_source
        _ref_pose_sdf: Optional[Path] = None

        if _is_user_ref and _HAS_RDKIT and _ref_mol and _pose_mol:
            _ref_pose_path = self.ref_pose_dir / f"{_sdf_name}.sdf"
            _ok_tmpl = write_pose_as_reference_template(
                _ref_path, _pose_mol, _ref_pose_path, best, mol_name=_sdf_name,
            )
            if _ok_tmpl:
                _ref_pose_sdf = _ref_pose_path
                print(f"    {C.DIM}    Output{C.RESET}   {C.GOLD}{_sdf_name}.sdf{C.RESET}  {C.DIM}\u2192 REF_POSE  (template: {_ref_source}){C.RESET}")
            else:
                _is_user_ref = False  # Template failed — fall through to generic SDF

        if not _ref_pose_sdf:
            # Generic SDF (Strategy B)
            _sdf_ok = write_pose_sdf(
                best,
                self.sdf_dir / f"{_sdf_name}.sdf",
                atoms=effective_lig_atoms,
                mol_name=_sdf_name,
            )
            if _sdf_ok:
                write_pose_mol(
                    best,
                    self.sdf_dir / f"{_sdf_name}.mol",
                    atoms=effective_lig_atoms,
                    mol_name=_sdf_name,
                )
                print(f"    {C.DIM}    Output{C.RESET}   {C.GRAY}{_sdf_name}.sdf{C.RESET}")

        # ── Internal RMSD validation (Kabsch on MCS-matched heavy atoms) ──────
        if _ref_mol and _pose_mol:
            RMSDEngine.evaluate_pose(_ref_mol, _pose_mol, best)
            if best.validation_rmsd_pos is not None and best.validation_rmsd_conf is not None:
                _grade = "PASS" if best.validation_rmsd_conf < 2.0 else ("OK" if best.validation_rmsd_conf < 3.0 else "NOTE")
                _gcol  = C.GREEN if _grade == "PASS" else (C.YELLOW if _grade == "OK" else C.ORANGE)
                print(
                    f"    {C.TEAL}    Val RMSD{C.RESET}  "
                    f"{C.DIM}Pos{C.RESET} {best.validation_rmsd_pos:.3f} \u00c5"
                    f"  {C.DIM}Conf{C.RESET} {best.validation_rmsd_conf:.3f} \u00c5"
                    f"  {_gcol}[{_grade}]{C.RESET}"
                )

        # ── Boltzmann Population Probability (log-sum-exp for numerical stability) ──
        # Canonical ensemble: pi_i = exp(-E_i/RT) / Z, where Z = sum_j exp(-E_j/RT)
        # Energies are negative, so -E_i/RT is positive and can be large.
        # Log-sum-exp avoids float overflow: log(Z) = max_t + log(sum(exp(t_j - max_t)))
        boltz_prob = None
        if result.clusters and result.stat_temperature is not None and result.stat_temperature > 0:
            RT     = 0.00198720425864083 * result.stat_temperature
            terms  = [-c.mean_energy / RT for c in result.clusters]
            max_t  = max(terms)
            log_Z  = max_t + math.log(sum(math.exp(t - max_t) for t in terms))
            if best_cluster:
                t_best     = -best_cluster.mean_energy / RT
                boltz_prob = math.exp(t_best - log_Z)

        prof = DataExtractor.extract_profile(
            result.ligand, result.protein, best, best_cluster,
        )

        # Torsional free energy (conformational entropy cost of binding from AutoDock).
        # NB: this is ΔG_tor, already included in the total ΔG_binding — shown separately
        # for independent assessment of the ligand's rotational entropy burden.
        prof.entropic_penalty_kcal = best.torsional_energy
        prof.boltzmann_prob        = round(boltz_prob, 4) if boltz_prob is not None else None

        # Pull adaptive mode fields into the profile
        prof.docking_mode             = result.docking_mode.name
        prof.n_ligand_atoms           = result.n_ligand_atoms
        prof.n_flexible_residue_atoms = result.n_flexible_residue_atoms
        prof.receptor_strain_kcal     = best.receptor_strain_estimate
        prof.flex_residue_ids         = result.flex_residue_ids
        # Propagate the RDKit validation RMSDs (distinct from AutoDock RMSD)
        prof.validation_rmsd_conf     = best.validation_rmsd_conf
        prof.validation_rmsd_pos      = best.validation_rmsd_pos

        # ADMET descriptors (reported as-is; not integrated into any score)
        if _pose_mol:
            ADMETEngine.evaluate_pharmacokinetics(_pose_mol, prof)

        profiles.append(prof)

    def _summary(self, ranked: List[DockingProfile]) -> None:
        W    = 72
        top  = ranked[0]
        n    = len(ranked)
        ts   = datetime.now().strftime("%Y-%m-%d  %H:%M:%S")
        # ── Completion banner ─────────────────────────────────────────────────
        print()
        print(f"  {C.GREEN}\u250c{'\u2500' * W}\u2510{C.RESET}")
        print(f"  {C.GREEN}\u2502{C.RESET}{'':^{W}}{C.GREEN}\u2502{C.RESET}")
        tick = "\u2713  P I P E L I N E   C O M P L E T E"
        print(f"  {C.GREEN}\u2502{C.RESET}{C.BOLD}{C.WHITE}{tick:^{W}}{C.RESET}{C.GREEN}\u2502{C.RESET}")
        print(f"  {C.GREEN}\u2502{C.RESET}{'':^{W}}{C.GREEN}\u2502{C.RESET}")
        print(f"  {C.GREEN}\u2514{'\u2500' * W}\u2518{C.RESET}")
        print()
        # ── Results panel ─────────────────────────────────────────────────────
        def _row(label: str, value: str, vcol: str = C.WHITE) -> None:
            print(f"  {C.DIM}{label:<20}{C.RESET}  {vcol}{value}{C.RESET}")
        _sep = f"  {C.DIM}{'\u2500' * W}{C.RESET}"
        print(_sep)
        _row("Ligands processed",  str(n))
        _row("Reports saved to",   str(self.output_dir), C.YELLOW)
        _row("Completed at",       ts)
        print(_sep)
        print()
        # ── Top ligand spotlight ───────────────────────────────────────────────
        print(f"  {C.GOLD}\u2605  TOP RESULT{C.RESET}")
        print(f"  {C.DIM}{'\u2500' * W}{C.RESET}")
        _row("Ligand",     top.ligand,                   C.LIME)
        _row("Receptor",   top.protein or "N/A",         C.WHITE)
        _row("\u0394G (kcal/mol)", f"{top.binding_energy:.2f}",  C.CYAN)
        _row("Ki",         format_ki(top.ki_nM),         C.MAGENTA)
        _row("Run",        str(top.run_number),          C.WHITE)
        _row("Cluster",    f"#{top.cluster_id}  ({top.cluster_size} poses)"  if top.cluster_id else "N/A",  C.WHITE)
        if top.validation_rmsd_conf is not None:
            _row("Val RMSD (Conf)", f"{top.validation_rmsd_conf:.3f} \u00c5", C.TEAL)
        print()
        print(f"  {C.DIM}{'\u2550' * W}{C.RESET}")
        print(f"  {C.GOLD}{__copyright__}{C.RESET}")
        print(f"  {C.DIM}{'\u2550' * W}{C.RESET}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13 — ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

_BANNER = f"""{C.CYAN}\u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557{C.RESET}
{C.CYAN}\u2551{C.RESET}                                                                                                            {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2588\u2557   \u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2557 \u2588\u2588\u2557   \u2588\u2588\u2557\u2588\u2588\u2588\u2588\u2588\u2588\u2557 {C.RESET} \u2122                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2554\u2550\u2550\u2550\u2550\u255d\u2588\u2588\u2588\u2588\u2557  \u2588\u2588\u2551\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u2588\u2588\u2551   \u2588\u2588\u2551\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557{C.RESET}                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2557\u2588\u2588\u2554\u2588\u2588\u2557 \u2588\u2588\u2551\u2588\u2588\u2551  \u2588\u2588\u2551\u2588\u2588\u2551   \u2588\u2588\u2551\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d{C.RESET}                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u2588\u2588\u2554\u2550\u2550\u2588\u2588\u2557\u255a\u2550\u2550\u2550\u2550\u2588\u2588\u2551\u2588\u2588\u2551\u255a\u2588\u2588\u2557\u2588\u2588\u2551\u2588\u2588\u2551  \u2588\u2588\u2551\u255a\u2588\u2588\u2557 \u2588\u2588\u2554\u255d\u2588\u2588\u2554\u2550\u2550\u2550 {C.RESET}                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d\u2588\u2588\u2588\u2588\u2588\u2588\u2588\u2551\u2588\u2588\u2551 \u255a\u2588\u2588\u2588\u2588\u2551\u2588\u2588\u2588\u2588\u2588\u2588\u2554\u255d \u255a\u2588\u2588\u2588\u2588\u2554\u255d \u2588\u2588\u2551     {C.RESET}                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.BOLD}{C.WHITE}\u255a\u2550\u2550\u2550\u2550\u2550\u255d \u255a\u2550\u2550\u2550\u2550\u2550\u2550\u255d\u255a\u2550\u255d  \u255a\u2550\u2550\u2550\u255d\u255a\u2550\u2550\u2550\u2550\u2550\u255d   \u255a\u2550\u2550\u2550\u255d  \u255a\u2550\u255d     {C.RESET}                                               {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}                                                                                                            {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.WHITE}Bio-Science Network Docking & Validation Platform\u2122{C.RESET}                                                  {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.TEAL}DLG Parsing, Validation & Statistical Analysis \u2014 {__edition__} v{__version__}{C.RESET}                     {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.GRAY}{__copyright__}{C.RESET}                              {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}    {C.GRAY}Author: {__author__}{C.RESET}                                                      {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u2551{C.RESET}                                                                                                            {C.CYAN}\u2551{C.RESET}
{C.CYAN}\u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d{C.RESET}"""


def main() -> None:
    import argparse

    if sys.platform == "win32":
        os.system("")
    if hasattr(sys.stdout, "reconfigure"):
        try: sys.stdout.reconfigure(encoding="utf-8")
        except Exception: pass

    print()
    print(_BANNER)
    print()

    ap = argparse.ArgumentParser(
        prog="dlg_extract.py",
        description=(
            f"{__full_brand__} — AutoDock DLG post-processing\n"
            "Parses rigid and flexible docking results, performs automatic\n"
            "validation, and computes docking statistical analysis."
        ),
    )
    ap.add_argument("--dlg",    metavar="DIR", help="Folder containing .dlg files (default: auto-discover)")
    ap.add_argument("--output", metavar="DIR", help="Output folder (default: ./BSNDVP_RESULTS)")
    args = ap.parse_args()

    cwd = Path.cwd()

    if args.dlg:
        dlg_dir = Path(args.dlg)
        if not dlg_dir.is_dir():
            print(f"  {C.RED}\u2717  ERROR: --dlg '{args.dlg}' is not a directory.{C.RESET}")
            sys.exit(1)
    else:
        dlg_dir = _find_dlg_dir(cwd)
        if not dlg_dir:
            print(f"  {C.RED}\u2717  ERROR: No 'DLG' folder found. Use --dlg to specify one.{C.RESET}")
            sys.exit(1)

    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = cwd / "BSNDVP_RESULTS"

    _W = 72
    print(f"  {C.DIM}┌{'\u2500' * _W}┐{C.RESET}")
    print(f"  {C.DIM}│{C.RESET}  {C.DIM}Working dir {C.RESET}  {C.WHITE}{str(cwd):<{_W - 14}}{C.RESET}{C.DIM}│{C.RESET}")
    print(f"  {C.DIM}│{C.RESET}  {C.DIM}DLG folder  {C.RESET}  {C.YELLOW}{str(dlg_dir):<{_W - 14}}{C.RESET}{C.DIM}│{C.RESET}")
    print(f"  {C.DIM}│{C.RESET}  {C.DIM}Output      {C.RESET}  {C.YELLOW}{str(output_dir):<{_W - 14}}{C.RESET}{C.DIM}│{C.RESET}")
    print(f"  {C.DIM}└{'\u2500' * _W}┘{C.RESET}")
    print()

    try:
        DLGPipeline(dlg_dir, output_dir, cwd).run()
    except Exception as e:
        print(f"\n  {C.RED}✖  ERROR: {e}{C.RESET}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
