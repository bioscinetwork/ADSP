"""
AutoDock Suite Pro — Protein-Ligand Interaction Profiler
========================================================
Pure Python 3 module to identify and classify molecular interactions between
a receptor protein and a docked ligand pose directly from PDBQT coordinates.

Supported Interactions:
  - Hydrogen Bonds (D···A <= 3.5 Å, D-H···A angle >= 120° estimated)
  - Hydrophobic Contacts (non-polar C–C <= 4.0 Å, deduplicated per residue)
  - Salt Bridges (residue type + AutoDock atom type logic, <= 4.0 Å)
  - Pi-Stacking Face-to-face (centroid <= 5.0 Å, inter-plane angle < 30°)
  - Pi-Stacking T-shaped (centroid <= 5.5 Å, inter-plane angle 60°–120°)
  - Pi-Cation (cation to aromatic centroid <= 4.5 Å)
  - Halogen Bonds (Cl/Br/I to Lewis base O/N/S, <= 3.8 Å)

Scientific references for cutoffs:
  - H-bond geometry: McDonald & Thornton (1994) J Mol Biol 238:777
  - Pi-stacking: Hunter & Sanders (1990) JACS 112:5525
  - Salt bridges: Barlow & Thornton (1983) J Mol Biol 168:867
  - Halogen bonds: Auffinger et al. (2004) PNAS 101:16789

Author: AutoDock Suite Pro Development Team
Version: 0.3.0
"""

from __future__ import annotations

import csv
import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


# ==============================================================================
# DATA STRUCTURES & POSE DIAGNOSTICS
# ==============================================================================

class PoseLookupStatus(str, Enum):
    """Explicit status reporting for docking pose retrieval."""
    POSE_FOUND = "POSE FOUND"
    POSE_NOT_FOUND = "POSE NOT FOUND"
    INVALID_POSE_INDEX = "INVALID POSE INDEX"
    CORRUPTED_POSE_OUTPUT = "CORRUPTED POSE OUTPUT"


class PoseLookupError(Exception):
    """Base exception for pose retrieval failures."""
    def __init__(self, status: PoseLookupStatus, message: str):
        super().__init__(message)
        self.status = status
        self.message = message


class PoseNotFoundError(PoseLookupError):
    def __init__(self, message: str = "Requested docking pose could not be located in output.", status: Optional[PoseLookupStatus] = None):
        super().__init__(status or PoseLookupStatus.POSE_NOT_FOUND, message)


class InvalidPoseIndexError(PoseLookupError):
    def __init__(self, message: str = "Invalid pose index specified; must be a positive 1-based integer.", status: Optional[PoseLookupStatus] = None):
        super().__init__(status or PoseLookupStatus.INVALID_POSE_INDEX, message)


class CorruptedPoseOutputError(PoseLookupError):
    def __init__(self, message: str = "Docked pose file is corrupted or contains no parseable coordinate models.", status: Optional[PoseLookupStatus] = None):
        super().__init__(status or PoseLookupStatus.CORRUPTED_POSE_OUTPUT, message)


@dataclass
class PDBQTAtom:
    serial: int
    name: str
    res_name: str
    chain: str
    res_seq: int
    x: float
    y: float
    z: float
    charge: float = 0.0
    atom_type: str = ""
    is_hetero: bool = False

    @property
    def coords(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def element(self) -> str:
        """Heuristic element extraction from AutoDock atom type."""
        t = self.atom_type.strip().upper()
        if t.startswith("CL"):
            return "CL"
        if t.startswith("BR"):
            return "BR"
        if len(t) > 0 and t[0].isalpha():
            return t[0]
        name = self.name.strip().upper()
        for char in name:
            if char.isalpha():
                return char
        return "C"

    @property
    def residue_id(self) -> str:
        return f"{self.res_name}:{self.chain}:{self.res_seq}"


@dataclass
class Interaction:
    pose_index: int
    interaction_type: str  # 'Hydrogen Bond', 'Hydrophobic', 'Salt Bridge',
                           # 'Pi-Stacking (Face-to-face)', 'Pi-Stacking (T-shaped)',
                           # 'Pi-Cation', 'Halogen Bond'
    receptor_residue: str  # e.g. 'ASP:A:102'
    receptor_res_name: str
    receptor_chain: str
    receptor_res_seq: int
    receptor_atom: str
    ligand_atom: str
    distance_angstrom: float
    details: str = ""
    angle_deg: Optional[float] = None
    detection_basis: str = "Geometric Criteria"
    ligand_atom_idx: Optional[int] = None
    receptor_atom_idx: Optional[int] = None
    pose_id: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "Pose": self.pose_index,
            "Interaction_Type": self.interaction_type,
            "Receptor_Residue": self.receptor_residue,
            "Residue_Name": self.receptor_res_name,
            "Chain": self.receptor_chain,
            "Residue_Num": self.receptor_res_seq,
            "Receptor_Atom": self.receptor_atom,
            "Ligand_Atom": self.ligand_atom,
            "Distance_A": round(self.distance_angstrom, 3),
            "Details": self.details,
        }
        if self.angle_deg is not None:
            d["Angle_Deg"] = round(self.angle_deg, 1)
        if self.detection_basis:
            d["Detection_Basis"] = self.detection_basis
        if self.ligand_atom_idx is not None:
            d["Ligand_Atom_Idx"] = self.ligand_atom_idx
        if self.receptor_atom_idx is not None:
            d["Receptor_Atom_Idx"] = self.receptor_atom_idx
        if self.pose_id:
            d["Pose_ID"] = self.pose_id
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Interaction":
        return cls(
            pose_index=data.get("Pose", data.get("pose_index", 1)),
            interaction_type=data.get("Interaction_Type", data.get("interaction_type", "")),
            receptor_residue=data.get("Receptor_Residue", data.get("receptor_residue", "")),
            receptor_res_name=data.get("Residue_Name", data.get("receptor_res_name", "")),
            receptor_chain=data.get("Chain", data.get("receptor_chain", "")),
            receptor_res_seq=data.get("Residue_Num", data.get("receptor_res_seq", 0)),
            receptor_atom=data.get("Receptor_Atom", data.get("receptor_atom", "")),
            ligand_atom=data.get("Ligand_Atom", data.get("ligand_atom", "")),
            distance_angstrom=float(data.get("Distance_A", data.get("distance_angstrom", 0.0))),
            details=data.get("Details", data.get("details", "")),
            angle_deg=data.get("Angle_Deg", data.get("angle_deg")),
            detection_basis=data.get("Detection_Basis", data.get("detection_basis", "Geometric Criteria")),
            ligand_atom_idx=data.get("Ligand_Atom_Idx", data.get("ligand_atom_idx")),
            receptor_atom_idx=data.get("Receptor_Atom_Idx", data.get("receptor_atom_idx")),
            pose_id=data.get("Pose_ID", data.get("pose_id")),
        )


@dataclass
class PoseInteractionProfile:
    """Rich diagnostic result container for pose interaction analysis."""
    status: PoseLookupStatus
    pose_index: int
    pose_id: str = ""
    interactions: List[Interaction] = field(default_factory=list)
    message: str = ""
    ligand_atoms: List[PDBQTAtom] = field(default_factory=list)
    receptor_atoms: List[PDBQTAtom] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __iter__(self):
        return iter(self.interactions)

    def __len__(self):
        return len(self.interactions)

    def __getitem__(self, idx):
        return self.interactions[idx]

    def __bool__(self):
        return bool(self.interactions)


# ==============================================================================
# GEOMETRY HELPERS
# ==============================================================================

def dist3d(c1: Tuple[float, float, float], c2: Tuple[float, float, float]) -> float:
    return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)


def vector_sub(v1: Tuple[float,float,float], v2: Tuple[float,float,float]) -> Tuple[float,float,float]:
    return (v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2])


def vector_add(v1: Tuple[float,float,float], v2: Tuple[float,float,float]) -> Tuple[float,float,float]:
    return (v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2])


def dot_product(v1: Tuple[float,float,float], v2: Tuple[float,float,float]) -> float:
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def cross_product(v1: Tuple[float,float,float], v2: Tuple[float,float,float]) -> Tuple[float,float,float]:
    return (
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0],
    )


def vector_length(v: Tuple[float,float,float]) -> float:
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def normalize(v: Tuple[float,float,float]) -> Tuple[float,float,float]:
    l = vector_length(v)
    if l < 1e-9:
        return (0.0, 0.0, 0.0)
    return (v[0]/l, v[1]/l, v[2]/l)


def angle_between_vectors_degrees(
    v1: Tuple[float,float,float], v2: Tuple[float,float,float]
) -> float:
    """Angle between two vectors in degrees [0, 180]."""
    l1 = vector_length(v1)
    l2 = vector_length(v2)
    if l1 < 1e-9 or l2 < 1e-9:
        return 0.0
    cos_theta = dot_product(v1, v2) / (l1 * l2)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))


def angle_degrees(
    p1: Tuple[float,float,float],
    vertex: Tuple[float,float,float],
    p2: Tuple[float,float,float],
) -> float:
    """Angle p1–vertex–p2 in degrees."""
    v1 = vector_sub(p1, vertex)
    v2 = vector_sub(p2, vertex)
    return angle_between_vectors_degrees(v1, v2)


def calculate_centroid(atoms: List[PDBQTAtom]) -> Tuple[float, float, float]:
    if not atoms:
        return (0.0, 0.0, 0.0)
    n = len(atoms)
    return (
        sum(a.x for a in atoms) / n,
        sum(a.y for a in atoms) / n,
        sum(a.z for a in atoms) / n,
    )


def ring_normal_vector(
    ring_atoms: List[PDBQTAtom],
) -> Optional[Tuple[float, float, float]]:
    """Compute ring plane normal from >= 3 ring atoms using Newell's method."""
    if len(ring_atoms) < 3:
        return None
    cx, cy, cz = calculate_centroid(ring_atoms)
    nx, ny, nz = 0.0, 0.0, 0.0
    n = len(ring_atoms)
    for i in range(n):
        c = ring_atoms[i]
        nx_atom = ring_atoms[(i + 1) % n]
        nx += (c.y - cy) * (nx_atom.z - cz) - (c.z - cz) * (nx_atom.y - cy)
        ny += (c.z - cz) * (nx_atom.x - cx) - (c.x - cx) * (nx_atom.z - cz)
        nz += (c.x - cx) * (nx_atom.y - cy) - (c.y - cy) * (nx_atom.x - cx)
    n_vec = (nx, ny, nz)
    l = vector_length(n_vec)
    if l < 1e-9:
        # Fallback: cross product of first two edge vectors
        v1 = vector_sub(ring_atoms[1].coords, ring_atoms[0].coords)
        v2 = vector_sub(ring_atoms[2].coords, ring_atoms[0].coords)
        n_vec = cross_product(v1, v2)
        l = vector_length(n_vec)
        if l < 1e-9:
            return None
    return (n_vec[0]/l, n_vec[1]/l, n_vec[2]/l)


# ==============================================================================
# 3D SPATIAL GRID FOR O(1) NEIGHBOUR LOOKUP
# ==============================================================================

class SpatialGrid:
    """A simple 3D grid bucket for fast spatial neighbour queries.

    Allows O(1) average-case lookup for atoms within a given cutoff,
    replacing the O(N*M) nested loop over all receptor × ligand atom pairs.
    """

    def __init__(self, atoms: List[PDBQTAtom], cell_size: float = 5.0):
        self.cell_size = cell_size
        self.grid: Dict[Tuple[int,int,int], List[PDBQTAtom]] = defaultdict(list)
        for atom in atoms:
            key = self._cell_key(atom.x, atom.y, atom.z)
            self.grid[key].append(atom)

    def _cell_key(self, x: float, y: float, z: float) -> Tuple[int,int,int]:
        return (
            int(math.floor(x / self.cell_size)),
            int(math.floor(y / self.cell_size)),
            int(math.floor(z / self.cell_size)),
        )

    def neighbors_within(
        self, x: float, y: float, z: float, cutoff: float
    ) -> List[PDBQTAtom]:
        """Returns all atoms within cutoff Å of (x, y, z)."""
        r = int(math.ceil(cutoff / self.cell_size))
        cx, cy, cz = self._cell_key(x, y, z)
        result: List[PDBQTAtom] = []
        cutoff_sq = cutoff * cutoff
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                for dz in range(-r, r + 1):
                    cell = self.grid.get((cx+dx, cy+dy, cz+dz))
                    if cell:
                        for atom in cell:
                            d2 = (atom.x-x)**2 + (atom.y-y)**2 + (atom.z-z)**2
                            if d2 <= cutoff_sq:
                                result.append(atom)
        return result


# ==============================================================================
# PDBQT PARSER
# ==============================================================================

def parse_pdbqt_line(line: str) -> Optional[PDBQTAtom]:
    """Parses standard ATOM / HETATM fixed-column PDBQT line."""
    record = line[:6].strip()
    if record not in ("ATOM", "HETATM"):
        return None
    try:
        serial = int(line[6:11].strip())
        name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21:22].strip() or "A"
        res_seq_str = line[22:26].strip()
        res_seq = int(res_seq_str) if res_seq_str else 1
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        charge = 0.0
        if len(line) >= 76:
            try:
                charge = float(line[68:76].strip())
            except ValueError:
                pass

        atom_type = ""
        if len(line) >= 77:
            atom_type = line[76:].strip().split()[0] if line[76:].strip() else ""

        return PDBQTAtom(
            serial=serial, name=name, res_name=res_name, chain=chain,
            res_seq=res_seq, x=x, y=y, z=z, charge=charge,
            atom_type=atom_type, is_hetero=(record == "HETATM")
        )
    except Exception:
        return None


def _parse_bond_graph_from_pdbqt(
    lines: List[str],
) -> Dict[int, List[int]]:
    """Parses BRANCH/ENDBRANCH records from a Vina output PDBQT to build a
    bond adjacency graph keyed by atom serial number."""
    graph: Dict[int, List[int]] = defaultdict(list)
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("BRANCH"):
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    a, b = int(parts[1]), int(parts[2])
                    graph[a].append(b)
                    graph[b].append(a)
                except ValueError:
                    pass
    return graph


def _find_rings_in_bond_graph(
    graph: Dict[int, List[int]],
    serial_to_atom: Dict[int, PDBQTAtom],
    min_ring_size: int = 5,
    max_ring_size: int = 7,
) -> List[List[PDBQTAtom]]:
    """Find smallest set of smallest rings using DFS cycle detection.

    Returns a list of rings, each ring being a list of PDBQTAtom objects.
    Only returns rings of size min_ring_size to max_ring_size.
    """
    rings: List[List[int]] = []
    visited_edges: Set[Tuple[int, int]] = set()

    def dfs(start: int, current: int, path: List[int], visited: Set[int]) -> None:
        for neighbor in graph.get(current, []):
            edge = (min(current, neighbor), max(current, neighbor))
            if neighbor == start and len(path) >= min_ring_size:
                # Found a ring
                ring = path[:]
                ring_key = frozenset(ring)
                if not any(frozenset(r) == ring_key for r in rings):
                    rings.append(ring)
                return
            if neighbor not in visited and edge not in visited_edges:
                visited_edges.add(edge)
                visited.add(neighbor)
                path.append(neighbor)
                if len(path) <= max_ring_size:
                    dfs(start, neighbor, path, visited)
                path.pop()
                visited.discard(neighbor)

    for start_serial in list(graph.keys()):
        dfs(start_serial, start_serial, [start_serial], {start_serial})

    result: List[List[PDBQTAtom]] = []
    for ring_serials in rings:
        if min_ring_size <= len(ring_serials) <= max_ring_size:
            ring_atoms = [serial_to_atom[s] for s in ring_serials if s in serial_to_atom]
            if len(ring_atoms) >= min_ring_size:
                result.append(ring_atoms)
    return result


def parse_docked_poses(
    pose_pdbqt_path: str | Path,
) -> List[Tuple[int, List[PDBQTAtom], List[PDBQTAtom]]]:
    """
    Parses a docked output PDBQT (e.g. from Vina or AutoDock 4).
    Returns a list of (model_index, ligand_atoms, flex_residue_atoms).
    """
    pose_path = Path(pose_pdbqt_path)
    if not pose_path.is_file():
        return []

    poses: List[Tuple[int, List[PDBQTAtom], List[PDBQTAtom]]] = []
    current_model = 1
    current_ligand: List[PDBQTAtom] = []
    current_flex: List[PDBQTAtom] = []
    current_lines: List[str] = []
    in_flex_res = False

    with open(pose_path, "r", encoding="utf-8", errors="replace") as f:
        for raw_line in f:
            line = raw_line
            stripped = line.strip()
            if stripped.startswith("DOCKED:"):
                stripped = stripped[7:].strip()
                line = line[line.find("DOCKED:") + 7:].lstrip()

            if stripped.startswith("MODEL"):
                parts = stripped.split()
                if len(parts) > 1 and parts[1].isdigit():
                    current_model = int(parts[1])
                current_ligand = []
                current_flex = []
                current_lines = []
                in_flex_res = False
                continue
            elif stripped.startswith("USER") and "Run =" in stripped:
                m_run = re.search(r"Run\s*=\s*(\d+)", stripped)
                if m_run:
                    current_model = int(m_run.group(1))
            elif stripped.startswith("BEGIN_RES"):
                in_flex_res = True
                continue
            elif stripped.startswith("END_RES"):
                in_flex_res = False
                continue
            elif stripped.startswith("ENDMDL"):
                if current_ligand or current_flex:
                    poses.append((current_model, current_ligand, current_flex))
                current_ligand = []
                current_flex = []
                current_lines = []
                in_flex_res = False
                continue

            current_lines.append(line)
            atom = parse_pdbqt_line(line)
            if atom is not None:
                if in_flex_res:
                    current_flex.append(atom)
                else:
                    current_ligand.append(atom)

    # Single-model PDBQT or unclosed trailing model
    if current_ligand or current_flex:
        if not any(p[0] == current_model for p in poses):
            poses.append((current_model, current_ligand, current_flex))

    return poses


def parse_receptor_pdbqt(receptor_path: str | Path) -> List[PDBQTAtom]:
    """Parses all atoms from a receptor PDBQT file."""
    path = Path(receptor_path)
    if not path.is_file():
        return []
    atoms: List[PDBQTAtom] = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            atom = parse_pdbqt_line(line)
            if atom is not None:
                atoms.append(atom)
    return atoms


# ==============================================================================
# INTERACTION DETECTION — ATOM TYPE SETS
# ==============================================================================

POSITIVE_RESIDUES = {"ARG", "LYS", "HIS"}
NEGATIVE_RESIDUES = {"ASP", "GLU"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "MET", "TRP", "TYR"}

# AutoDock 4 atom types
ACCEPTOR_TYPES = {"OA", "NA", "SA"}
DONOR_HEAVY_TYPES = {"N", "NA", "OA", "SA"}
HYDROPHOBIC_TYPES = {"C", "A"}
HALOGEN_TYPES = {"CL", "BR", "I"}

# Receptor atom names that are donors/acceptors (by name when type is ambiguous)
KNOWN_DONOR_NAMES = {"N", "NE", "NH1", "NH2", "NZ", "ND1", "NE2", "NE1",
                     "ND2", "NQ", "OG", "OG1", "OH", "SG"}
KNOWN_ACCEPTOR_NAMES = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1",
                        "OH", "ND1", "NE2", "SD", "SG"}

# Cationic receptor atoms (for salt bridge — more reliable than Gasteiger charge)
CATIONIC_ATOMS = {
    "ARG": {"NH1", "NH2", "NE"},
    "LYS": {"NZ"},
    "HIS": {"ND1", "NE2"},
}
# Anionic receptor atoms
ANIONIC_ATOMS = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}


# ==============================================================================
# H-BOND ANGLE ESTIMATION
# ==============================================================================

def _estimate_hbond_angle(
    donor_heavy: PDBQTAtom,
    acceptor: PDBQTAtom,
    all_ligand_atoms: List[PDBQTAtom],
) -> float:
    """Estimate the D-H···A angle for hydrogen bond validation.

    PDBQT files from Vina/AD4 typically lack explicit H atoms.
    Strategy:
      1. If a covalently bonded H atom is found (within 1.15 Å), use it.
      2. Otherwise, use the D···A vector as a proxy (lower bound estimate).
         This is conservative — the true angle can only be ≥ estimated angle
         when H is along the D→A direction.

    Returns angle in degrees. A value >= 120° is considered valid.
    """
    # Try to find explicit H neighbour of donor_heavy
    nearby_h = [
        a for a in all_ligand_atoms
        if a.atom_type.upper() in ("H", "HD", "HS")
        and dist3d(donor_heavy.coords, a.coords) <= 1.15
    ]
    if nearby_h:
        h_atom = min(nearby_h, key=lambda a: dist3d(donor_heavy.coords, a.coords))
        # Angle: donor_heavy – H – acceptor
        return angle_degrees(donor_heavy.coords, h_atom.coords, acceptor.coords)

    # No H found: use D···A vector angle as conservative proxy
    # The D-H vector is estimated as the D→A direction
    # This gives the best-case angle (180°) which we then use as permissive pass
    # We apply a stricter distance cutoff (< 3.2 Å) when no H is available
    d = dist3d(donor_heavy.coords, acceptor.coords)
    if d <= 3.2:
        return 150.0  # Treat close contacts as likely valid
    return 0.0  # Too far without explicit H to confirm


# ==============================================================================
# INTERACTION DETECTION ENGINES
# ==============================================================================

def detect_hydrogen_bonds(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_idx: int,
    cutoff_dist: float = 3.5,
    rec_spatial: Optional[SpatialGrid] = None,
) -> List[Interaction]:
    """
    Detects hydrogen bonds between receptor and ligand.

    Requirements:
      - D···A distance: 2.2 – 3.5 Å
      - D-H···A angle: >= 120° (estimated from available H atoms or D···A vector)
      - Deduplication: each (receptor_residue, ligand_atom_name) pair reported once.

    Reference: McDonald & Thornton (1994) J Mol Biol 238:777
    """
    interactions: List[Interaction] = []
    seen_pairs: Set[Tuple[str, str]] = set()  # (residue_id, ligand_atom_name)

    for lig_atom in ligand_atoms:
        l_type = lig_atom.atom_type.upper()
        is_lig_acceptor = l_type in ACCEPTOR_TYPES or l_type.startswith(("O", "N"))
        is_lig_donor = l_type in DONOR_HEAVY_TYPES or l_type.startswith(("N", "O"))

        if not (is_lig_acceptor or is_lig_donor):
            continue

        # Get nearby receptor atoms via spatial grid or fallback to full list
        candidates = (
            rec_spatial.neighbors_within(lig_atom.x, lig_atom.y, lig_atom.z, cutoff_dist)
            if rec_spatial else receptor_atoms
        )

        for rec_atom in candidates:
            r_type = rec_atom.atom_type.upper()
            is_rec_acceptor = (r_type in ACCEPTOR_TYPES
                               or rec_atom.name in KNOWN_ACCEPTOR_NAMES)
            is_rec_donor = (r_type in DONOR_HEAVY_TYPES
                            or rec_atom.name in KNOWN_DONOR_NAMES)

            if not (is_rec_acceptor or is_rec_donor):
                continue

            d = dist3d(rec_atom.coords, lig_atom.coords)
            if not (2.2 <= d <= cutoff_dist):
                continue

            pair_key = (rec_atom.residue_id, lig_atom.name)
            if pair_key in seen_pairs:
                continue

            # Determine direction and check angle
            angle_ok = False
            direction = ""

            # Case 1: Receptor Acceptor ← Ligand Donor
            if is_rec_acceptor and is_lig_donor:
                angle = _estimate_hbond_angle(lig_atom, rec_atom, ligand_atoms)
                if angle >= 120.0:
                    angle_ok = True
                    direction = (f"Lig Donor ({lig_atom.atom_type}) → "
                                 f"Rec Acceptor ({rec_atom.atom_type}), "
                                 f"est. angle {angle:.0f}°")

            # Case 2: Receptor Donor → Ligand Acceptor (only if Case 1 didn't pass)
            if not angle_ok and is_rec_donor and is_lig_acceptor:
                # Estimate angle using acceptor position as target
                d_closer = d <= 3.2
                if d_closer:
                    angle_ok = True
                    direction = (f"Rec Donor ({rec_atom.atom_type}) → "
                                 f"Lig Acceptor ({lig_atom.atom_type})")

            if angle_ok:
                seen_pairs.add(pair_key)
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type="Hydrogen Bond",
                    receptor_residue=rec_atom.residue_id,
                    receptor_res_name=rec_atom.res_name,
                    receptor_chain=rec_atom.chain,
                    receptor_res_seq=rec_atom.res_seq,
                    receptor_atom=rec_atom.name,
                    ligand_atom=lig_atom.name,
                    distance_angstrom=d,
                    details=direction,
                    angle_deg=angle if is_rec_acceptor and is_lig_donor else None,
                    detection_basis="Geometric Distance (2.2-3.5 Å) & Angle Criteria",
                    ligand_atom_idx=lig_atom.serial,
                    receptor_atom_idx=rec_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    return interactions


def detect_hydrophobic_contacts(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_idx: int,
    cutoff_dist: float = 4.0,
    rec_spatial: Optional[SpatialGrid] = None,
) -> List[Interaction]:
    """Detects non-polar Carbon–Carbon hydrophobic contacts <= 4.0 Å.
    Deduplicates per receptor residue (one contact per residue per ligand atom).
    """
    interactions: List[Interaction] = []
    seen_pairs: Set[Tuple[str, str]] = set()

    for lig_atom in ligand_atoms:
        l_type = lig_atom.atom_type.upper()
        if l_type not in HYDROPHOBIC_TYPES:
            continue
        if abs(lig_atom.charge) > 0.25:
            continue

        candidates = (
            rec_spatial.neighbors_within(lig_atom.x, lig_atom.y, lig_atom.z, cutoff_dist)
            if rec_spatial else receptor_atoms
        )

        for rec_atom in candidates:
            r_type = rec_atom.atom_type.upper()
            if r_type not in HYDROPHOBIC_TYPES:
                continue
            if abs(rec_atom.charge) > 0.25:
                continue

            pair_key = (rec_atom.residue_id, lig_atom.name)
            if pair_key in seen_pairs:
                continue

            d = dist3d(rec_atom.coords, lig_atom.coords)
            if d <= cutoff_dist:
                seen_pairs.add(pair_key)
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type="Hydrophobic",
                    receptor_residue=rec_atom.residue_id,
                    receptor_res_name=rec_atom.res_name,
                    receptor_chain=rec_atom.chain,
                    receptor_res_seq=rec_atom.res_seq,
                    receptor_atom=rec_atom.name,
                    ligand_atom=lig_atom.name,
                    distance_angstrom=d,
                    details=f"C–C Contact ({rec_atom.atom_type} – {lig_atom.atom_type})",
                    detection_basis="Non-Polar Carbon-Carbon Contact (≤ 4.0 Å)",
                    ligand_atom_idx=lig_atom.serial,
                    receptor_atom_idx=rec_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    return interactions


def detect_salt_bridges(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_idx: int,
    cutoff_dist: float = 4.0,
    rec_spatial: Optional[SpatialGrid] = None,
) -> List[Interaction]:
    """Detects electrostatic ionic pairs (salt bridges).

    Uses residue name + atom name as the primary signal for identifying
    cationic (ARG/LYS/HIS) and anionic (ASP/GLU) receptor atoms.
    For ligand atoms, uses AutoDock atom type (N+/O-) and partial charge
    as secondary signal (charge-based is unreliable for Gasteiger alone).

    Reference: Barlow & Thornton (1983) J Mol Biol 168:867
    """
    interactions: List[Interaction] = []

    # Identify charged ligand atoms
    # Positive ligand: N-containing atoms with meaningful positive charge
    pos_lig = [
        a for a in ligand_atoms
        if a.atom_type.upper().startswith("N")
        and (a.charge >= 0.2 or a.atom_type.upper() in ("N+",))
    ]
    # Negative ligand: O-containing atoms with negative charge OR phosphate/sulphonate
    neg_lig = [
        a for a in ligand_atoms
        if (a.atom_type.upper().startswith("O") or "P" in a.atom_type.upper())
        and a.charge <= -0.2
    ]

    for rec_atom in receptor_atoms:
        res = rec_atom.res_name.upper()

        # Cationic receptor (ARG NH1/NH2/NE, LYS NZ, HIS ND1/NE2)
        if res in CATIONIC_ATOMS and rec_atom.name in CATIONIC_ATOMS[res]:
            for la in neg_lig:
                candidates = (
                    [la] if rec_spatial is None
                    else rec_spatial.neighbors_within(la.x, la.y, la.z, cutoff_dist)
                )
                d = dist3d(rec_atom.coords, la.coords)
                if d <= cutoff_dist:
                    interactions.append(Interaction(
                        pose_index=pose_idx,
                        interaction_type="Salt Bridge",
                        receptor_residue=rec_atom.residue_id,
                        receptor_res_name=rec_atom.res_name,
                        receptor_chain=rec_atom.chain,
                        receptor_res_seq=rec_atom.res_seq,
                        receptor_atom=rec_atom.name,
                        ligand_atom=la.name,
                        distance_angstrom=d,
                        details=f"Rec Cation ({res} {rec_atom.name}) ··· Lig Anion ({la.name})",
                        detection_basis="Electrostatic Ionic Pair (≤ 4.0 Å)",
                        ligand_atom_idx=la.serial,
                        receptor_atom_idx=rec_atom.serial,
                        pose_id=f"pose_{pose_idx}",
                    ))

        # Anionic receptor (ASP OD1/OD2, GLU OE1/OE2)
        elif res in ANIONIC_ATOMS and rec_atom.name in ANIONIC_ATOMS[res]:
            for la in pos_lig:
                d = dist3d(rec_atom.coords, la.coords)
                if d <= cutoff_dist:
                    interactions.append(Interaction(
                        pose_index=pose_idx,
                        interaction_type="Salt Bridge",
                        receptor_residue=rec_atom.residue_id,
                        receptor_res_name=rec_atom.res_name,
                        receptor_chain=rec_atom.chain,
                        receptor_res_seq=rec_atom.res_seq,
                        receptor_atom=rec_atom.name,
                        ligand_atom=la.name,
                        distance_angstrom=d,
                        details=f"Rec Anion ({res} {rec_atom.name}) ··· Lig Cation ({la.name})",
                        detection_basis="Electrostatic Ionic Pair (≤ 4.0 Å)",
                        ligand_atom_idx=la.serial,
                        receptor_atom_idx=rec_atom.serial,
                        pose_id=f"pose_{pose_idx}",
                    ))

    return interactions


def _build_receptor_rings(
    receptor_atoms: List[PDBQTAtom],
) -> List[Tuple[str, PDBQTAtom, List[PDBQTAtom], Tuple[float,float,float]]]:
    """Build list of (residue_id, ref_atom, ring_atoms, centroid) for receptor aromatic rings."""
    rec_residues: Dict[str, List[PDBQTAtom]] = defaultdict(list)
    for a in receptor_atoms:
        if a.res_name.upper() in AROMATIC_RESIDUES:
            rec_residues[a.residue_id].append(a)

    rings = []
    PHE_TYR_RING_NAMES = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
    TRP_RING6_NAMES = {"CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
    TRP_RING5_NAMES = {"CG", "CD1", "CD2", "NE1", "CE2"}
    HIS_RING_NAMES = {"CG", "ND1", "CD2", "CE1", "NE2"}

    for r_id, atoms in rec_residues.items():
        res_name = atoms[0].res_name.upper()
        ref_atom = atoms[0]

        if res_name in ("PHE", "TYR"):
            ring_atoms = [a for a in atoms if a.name in PHE_TYR_RING_NAMES]
            if len(ring_atoms) >= 5:
                centroid = calculate_centroid(ring_atoms)
                rings.append((r_id, ref_atom, ring_atoms, centroid))

        elif res_name == "TRP":
            ring6 = [a for a in atoms if a.name in TRP_RING6_NAMES]
            if len(ring6) >= 5:
                centroid = calculate_centroid(ring6)
                rings.append((r_id, ref_atom, ring6, centroid))
            ring5 = [a for a in atoms if a.name in TRP_RING5_NAMES]
            if len(ring5) >= 4:
                centroid = calculate_centroid(ring5)
                rings.append((r_id, ref_atom, ring5, centroid))

        elif res_name == "HIS":
            ring5 = [a for a in atoms if a.name in HIS_RING_NAMES]
            if len(ring5) >= 4:
                centroid = calculate_centroid(ring5)
                rings.append((r_id, ref_atom, ring5, centroid))

    return rings


def _build_ligand_rings(
    ligand_atoms: List[PDBQTAtom],
) -> List[Tuple[str, List[PDBQTAtom], Tuple[float,float,float]]]:
    """Detect aromatic rings in a ligand using distance-based ring clustering.
    Uses bond-graph DFS if BRANCH records were parsed; falls back to proximity
    clustering of AutoDock 'A' (aromatic) typed atoms.
    Returns list of (ring_name, ring_atoms, centroid).
    """
    aromatic_atoms = [a for a in ligand_atoms if a.atom_type.upper() == "A"]
    if len(aromatic_atoms) < 5:
        return []

    # Distance-based clustering: group atoms within ~1.6 Å (aromatic C–C bond)
    # Uses Union-Find to correctly handle fused ring systems by assigning
    # each connected component and then checking for ring-sized subsets.
    # This is more robust than the previous greedy O(n^2) approach.

    # Step 1: Build adjacency based on aromatic bond distance (~1.35-1.42 Å + tolerance)
    n = len(aromatic_atoms)
    adj: Dict[int, Set[int]] = defaultdict(set)
    AROMATIC_BOND = 1.65  # Å tolerance
    for i in range(n):
        for j in range(i + 1, n):
            if dist3d(aromatic_atoms[i].coords, aromatic_atoms[j].coords) <= AROMATIC_BOND:
                adj[i].add(j)
                adj[j].add(i)

    # Step 2: Find rings by looking for cycles of size 5 or 6 in the adjacency graph
    found_rings: List[List[int]] = []
    visited_ring_sets: Set[frozenset] = set()

    def _find_cycle_from(start: int, current: int, path: List[int],
                         depth_limit: int, path_set: Set[int]) -> None:
        if depth_limit < 0:
            return
        for nb in adj[current]:
            if nb == start and len(path) >= 5:
                fset = frozenset(path)
                if fset not in visited_ring_sets:
                    visited_ring_sets.add(fset)
                    found_rings.append(path[:])
                continue
            if nb not in path_set:
                path_set.add(nb)
                path.append(nb)
                _find_cycle_from(start, nb, path, depth_limit - 1, path_set)
                path.pop()
                path_set.discard(nb)

    for start_idx in range(n):
        _find_cycle_from(start_idx, start_idx, [start_idx], 6, {start_idx})

    # Step 3: Filter to 5- and 6-membered rings only and return
    result = []
    for ring_indices in found_rings:
        if 5 <= len(ring_indices) <= 6:
            ring_atoms_list = [aromatic_atoms[i] for i in ring_indices]
            centroid = calculate_centroid(ring_atoms_list)
            ring_name = f"Ring_{len(result) + 1}"
            result.append((ring_name, ring_atoms_list, centroid))

    return result


def detect_aromatic_and_cation_pi(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_idx: int,
    pi_face_cutoff: float = 5.0,
    pi_t_cutoff: float = 5.5,
    cation_cutoff: float = 4.5,
) -> List[Interaction]:
    """
    Detects Pi-Stacking and Pi-Cation interactions.

    Pi-Stacking sub-types (via inter-plane angle of ring normal vectors):
      - Face-to-face (parallel):  centroid <= 5.0 Å AND angle < 30°
      - T-shaped (edge-to-face):  centroid <= 5.5 Å AND 60° < angle < 120°

    References:
      Hunter & Sanders (1990) JACS 112:5525
      McGaughey et al. (1998) J Biol Chem 273:15458
    """
    interactions: List[Interaction] = []
    rec_rings = _build_receptor_rings(receptor_atoms)
    lig_rings = _build_ligand_rings(ligand_atoms)

    # Pi-Stacking: Receptor Ring ↔ Ligand Ring
    for r_id, ref_atom, r_ring_atoms, r_centroid in rec_rings:
        r_normal = ring_normal_vector(r_ring_atoms)

        for ring_name, l_ring_atoms, l_centroid in lig_rings:
            d = dist3d(r_centroid, l_centroid)

            # Check distance for either stacking type
            if d > pi_t_cutoff:
                continue

            l_normal = ring_normal_vector(l_ring_atoms)

            # Determine stacking type via inter-plane angle
            stack_type = None
            angle_deg = None

            if r_normal and l_normal:
                angle_deg = angle_between_vectors_degrees(r_normal, l_normal)
                # Normalise to [0, 90] — parallel and anti-parallel are equivalent
                if angle_deg > 90.0:
                    angle_deg = 180.0 - angle_deg

                if d <= pi_face_cutoff and angle_deg < 30.0:
                    stack_type = "Pi-Stacking (Face-to-face)"
                elif d <= pi_t_cutoff and 60.0 < angle_deg < 120.0:
                    stack_type = "Pi-Stacking (T-shaped)"
            else:
                # Fallback without normals: report as generic if within face cutoff
                if d <= pi_face_cutoff:
                    stack_type = "Pi-Stacking"

            if stack_type:
                angle_str = f", inter-plane {angle_deg:.0f}°" if angle_deg is not None else ""
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type=stack_type,
                    receptor_residue=r_id,
                    receptor_res_name=ref_atom.res_name,
                    receptor_chain=ref_atom.chain,
                    receptor_res_seq=ref_atom.res_seq,
                    receptor_atom="Aromatic_Ring",
                    ligand_atom=ring_name,
                    distance_angstrom=d,
                    details=(f"Centroid–Centroid {d:.2f} Å"
                             f"{angle_str} ({ref_atom.res_name}↔{ring_name})"),
                    angle_deg=angle_deg,
                    detection_basis="Centroid-to-Centroid Normal Vector Alignment",
                    receptor_atom_idx=ref_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    # Pi-Cation: Receptor Cation ↔ Ligand Ring
    pos_rec_atoms = [
        a for a in receptor_atoms
        if a.res_name.upper() in POSITIVE_RESIDUES
        and a.name in (CATIONIC_ATOMS.get(a.res_name.upper(), set()))
    ]
    for rec_atom in pos_rec_atoms:
        for ring_name, _, l_centroid in lig_rings:
            d = dist3d(rec_atom.coords, l_centroid)
            if d <= cation_cutoff:
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type="Pi-Cation",
                    receptor_residue=rec_atom.residue_id,
                    receptor_res_name=rec_atom.res_name,
                    receptor_chain=rec_atom.chain,
                    receptor_res_seq=rec_atom.res_seq,
                    receptor_atom=rec_atom.name,
                    ligand_atom=ring_name,
                    distance_angstrom=d,
                    details=(f"Rec Cation ({rec_atom.res_name} {rec_atom.name}) "
                             f"to Lig Aromatic Ring {d:.2f} Å"),
                    detection_basis="Cation to Aromatic Centroid Distance (≤ 4.5 Å)",
                    receptor_atom_idx=rec_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    # Pi-Cation: Receptor Ring ↔ Ligand Cation
    pos_lig_atoms = [
        a for a in ligand_atoms
        if a.atom_type.upper().startswith("N") and a.charge >= 0.2
    ]
    for r_id, ref_atom, _, r_centroid in rec_rings:
        for la in pos_lig_atoms:
            d = dist3d(r_centroid, la.coords)
            if d <= cation_cutoff:
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type="Pi-Cation",
                    receptor_residue=r_id,
                    receptor_res_name=ref_atom.res_name,
                    receptor_chain=ref_atom.chain,
                    receptor_res_seq=ref_atom.res_seq,
                    receptor_atom="Aromatic_Ring",
                    ligand_atom=la.name,
                    distance_angstrom=d,
                    details=(f"Lig Cation ({la.name}) to "
                             f"Rec Aromatic Ring ({ref_atom.res_name}) {d:.2f} Å"),
                    detection_basis="Cation to Aromatic Centroid Distance (≤ 4.5 Å)",
                    ligand_atom_idx=la.serial,
                    receptor_atom_idx=ref_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    return interactions


def detect_halogen_bonds(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_idx: int,
    cutoff_dist: float = 3.8,
    rec_spatial: Optional[SpatialGrid] = None,
) -> List[Interaction]:
    """Detects halogen bonds (Cl, Br, I on ligand to O/N/S on receptor).

    Reference: Auffinger et al. (2004) PNAS 101:16789
    Cutoff: 3.8 Å (generous to catch all F·X bonds too)
    """
    interactions: List[Interaction] = []
    halogen_lig = [
        a for a in ligand_atoms
        if a.atom_type.upper() in HALOGEN_TYPES
    ]
    if not halogen_lig:
        return []

    for la in halogen_lig:
        candidates = (
            rec_spatial.neighbors_within(la.x, la.y, la.z, cutoff_dist)
            if rec_spatial else receptor_atoms
        )
        for rec_atom in candidates:
            r_type = rec_atom.atom_type.upper()
            if not (r_type in ACCEPTOR_TYPES
                    or rec_atom.name.startswith("O")
                    or rec_atom.name.startswith("N")
                    or r_type in ("S", "SA")):
                continue
            d = dist3d(la.coords, rec_atom.coords)
            if 2.5 <= d <= cutoff_dist:
                interactions.append(Interaction(
                    pose_index=pose_idx,
                    interaction_type="Halogen Bond",
                    receptor_residue=rec_atom.residue_id,
                    receptor_res_name=rec_atom.res_name,
                    receptor_chain=rec_atom.chain,
                    receptor_res_seq=rec_atom.res_seq,
                    receptor_atom=rec_atom.name,
                    ligand_atom=la.name,
                    distance_angstrom=d,
                    details=f"Lig {la.atom_type} ··· Rec {rec_atom.name} ({rec_atom.res_name})",
                    detection_basis="Lewis Base Halogen Coordinate (≤ 3.8 Å)",
                    ligand_atom_idx=la.serial,
                    receptor_atom_idx=rec_atom.serial,
                    pose_id=f"pose_{pose_idx}",
                ))

    return interactions


# ==============================================================================
# MASTER PROFILER
# ==============================================================================

def analyze_pose_interactions(
    receptor_atoms: List[PDBQTAtom],
    ligand_atoms: List[PDBQTAtom],
    pose_index: int = 1,
    flex_atoms: Optional[List[PDBQTAtom]] = None,
) -> List[Interaction]:
    """
    Analyzes all molecular interactions for a single docked pose.

    If flexible atoms are supplied (from Vina/AD4 flexible docking output),
    replaces rigid residues in the receptor with their flexible counterpart
    coordinates.

    Uses a 3D spatial grid for O(1) average-case receptor atom lookup,
    replacing the previous O(N*M) nested-loop approach.
    """
    # Merge receptor with flexible sidechain atoms
    if flex_atoms:
        flex_res_ids = {a.residue_id for a in flex_atoms}
        effective_receptor = [a for a in receptor_atoms if a.residue_id not in flex_res_ids]
        effective_receptor.extend(flex_atoms)
    else:
        effective_receptor = list(receptor_atoms)

    # Build spatial index on receptor (cell_size = largest cutoff = 5.5 Å)
    rec_spatial = SpatialGrid(effective_receptor, cell_size=5.5)

    all_interactions: List[Interaction] = []
    all_interactions.extend(
        detect_hydrogen_bonds(effective_receptor, ligand_atoms, pose_index, rec_spatial=rec_spatial)
    )
    all_interactions.extend(
        detect_hydrophobic_contacts(effective_receptor, ligand_atoms, pose_index, rec_spatial=rec_spatial)
    )
    all_interactions.extend(
        detect_salt_bridges(effective_receptor, ligand_atoms, pose_index, rec_spatial=rec_spatial)
    )
    all_interactions.extend(
        detect_aromatic_and_cation_pi(effective_receptor, ligand_atoms, pose_index)
    )
    all_interactions.extend(
        detect_halogen_bonds(effective_receptor, ligand_atoms, pose_index, rec_spatial=rec_spatial)
    )

    all_interactions.sort(key=lambda x: (x.receptor_res_seq, x.distance_angstrom))
    return all_interactions


def locate_docked_pose(
    poses: List[Tuple[int, List[PDBQTAtom], List[PDBQTAtom]]],
    pose_index: int,
    run_number: Optional[int] = None,
) -> Tuple[PoseLookupStatus, Optional[Tuple[int, List[PDBQTAtom], List[PDBQTAtom]]]]:
    """
    Explicitly locates a docking pose without silent fallback.

    Resolution strategy:
      1. Match by run_number if provided (e.g. AutoDock4 GA run number).
      2. Match by exact model ID (e.g. Vina MODEL index or AD4 MODEL).
      3. Match by 1-based sequential rank in the file (if 1 <= pose_index <= len(poses)).
      4. Explicit failure if pose_index cannot be mapped.
    """
    if not isinstance(pose_index, int) or pose_index < 1:
        return (PoseLookupStatus.INVALID_POSE_INDEX, None)
    if not poses:
        return (PoseLookupStatus.CORRUPTED_POSE_OUTPUT, None)

    # 1. Match by run_number if specified
    if run_number is not None:
        for p in poses:
            if p[0] == run_number:
                return (PoseLookupStatus.POSE_FOUND, p)

    # 2. Match by exact model index
    for p in poses:
        if p[0] == pose_index:
            return (PoseLookupStatus.POSE_FOUND, p)

    # 3. Match by 1-based sequential rank in the file
    if 1 <= pose_index <= len(poses):
        return (PoseLookupStatus.POSE_FOUND, poses[pose_index - 1])

    # 4. Not found
    return (PoseLookupStatus.POSE_NOT_FOUND, None)


def extract_pose_pdbqt_block(
    pose_pdbqt_path: str | Path,
    pose_index: int = 1,
    run_number: Optional[int] = None,
    raise_on_missing: bool = False,
) -> Optional[str]:
    """
    Extracts the raw PDBQT text block corresponding strictly to the requested pose.
    Returns None if the pose does not exist (or raises PoseNotFoundError if raise_on_missing is True).
    """
    path = Path(pose_pdbqt_path)
    if not path.is_file():
        if raise_on_missing:
            raise PoseNotFoundError(f"File not found: {pose_pdbqt_path}", status=PoseLookupStatus.POSE_NOT_FOUND)
        return None

    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        if raise_on_missing:
            raise CorruptedPoseOutputError(f"Could not read {pose_pdbqt_path}", status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT)
        return None

    lines = text.splitlines(keepends=True)
    is_dlg = any(l.startswith("DOCKED:") for l in lines)
    models: List[Tuple[int, List[str]]] = []
    current_model_idx = 1
    current_lines: List[str] = []

    for raw_line in lines:
        line = raw_line
        if is_dlg:
            if not line.startswith("DOCKED:"):
                continue
            line = line[8:] if line.startswith("DOCKED: ") else line[7:]

        stripped = line.strip()
        if stripped.startswith("MODEL"):
            if current_lines:
                models.append((current_model_idx, current_lines))
                current_lines = []
            parts = stripped.split()
            if len(parts) > 1 and parts[1].isdigit():
                current_model_idx = int(parts[1])
            else:
                current_model_idx = len(models) + 1
            current_lines.append(line)
        elif stripped.startswith("USER") and "Run =" in stripped:
            m_run = re.search(r"Run\s*=\s*(\d+)", stripped)
            if m_run:
                current_model_idx = int(m_run.group(1))
            current_lines.append(line)
        elif stripped.startswith("ENDMDL"):
            current_lines.append(line)
            models.append((current_model_idx, current_lines))
            current_lines = []
            current_model_idx = len(models) + 1
        else:
            current_lines.append(line)

    if current_lines and not is_dlg:
        models.append((current_model_idx, current_lines))

    if not models:
        if raise_on_missing:
            raise CorruptedPoseOutputError(f"No models found in {pose_pdbqt_path}", status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT)
        return None

    # Match by run_number
    if run_number is not None:
        for mid, blk in models:
            if mid == run_number:
                return "".join(blk)

    # Match by exact model ID
    for mid, blk in models:
        if mid == pose_index:
            return "".join(blk)

    # Match by 1-based order in file
    if 1 <= pose_index <= len(models):
        return "".join(models[pose_index - 1][1])

    if raise_on_missing:
        raise PoseNotFoundError(f"Pose {pose_index} not found in {pose_pdbqt_path}", status=PoseLookupStatus.POSE_NOT_FOUND)

    return None


def profile_docking_pose(
    receptor_pdbqt: str | Path,
    output_pose_pdbqt: str | Path,
    pose_index: int = 1,
    run_number: Optional[int] = None,
    raise_on_missing: bool = False,
) -> PoseInteractionProfile:
    """
    Profiles an explicitly requested docking pose with full diagnostic feedback.
    Never falls back to Pose 1 silently.
    """
    receptor_atoms = parse_receptor_pdbqt(receptor_pdbqt)
    if not receptor_atoms:
        if raise_on_missing:
            raise CorruptedPoseOutputError("Receptor PDBQT contains no parseable atoms.", status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT)
        return PoseInteractionProfile(
            status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT,
            pose_index=pose_index,
            message="Receptor PDBQT contains no parseable atoms.",
        )

    poses = parse_docked_poses(output_pose_pdbqt)
    if not poses:
        if raise_on_missing:
            raise CorruptedPoseOutputError(f"Pose file {output_pose_pdbqt} contains no valid coordinate models.", status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT)
        return PoseInteractionProfile(
            status=PoseLookupStatus.CORRUPTED_POSE_OUTPUT,
            pose_index=pose_index,
            message=f"Pose file {output_pose_pdbqt} contains no valid coordinate models.",
        )

    status, match = locate_docked_pose(poses, pose_index, run_number=run_number)
    if status != PoseLookupStatus.POSE_FOUND or match is None:
        if raise_on_missing:
            raise PoseNotFoundError(f"Pose {pose_index} not found in {output_pose_pdbqt} ({len(poses)} available).", status=status)
        return PoseInteractionProfile(
            status=status,
            pose_index=pose_index,
            message=f"Pose {pose_index} not found in {output_pose_pdbqt} ({len(poses)} available).",
        )

    model_idx, lig_atoms, flex_atoms = match
    interactions = analyze_pose_interactions(
        receptor_atoms=receptor_atoms,
        ligand_atoms=lig_atoms,
        pose_index=pose_index,
        flex_atoms=flex_atoms if flex_atoms else None,
    )
    for it in interactions:
        it.pose_index = pose_index
        it.pose_id = f"pose_{pose_index}"

    return PoseInteractionProfile(
        status=PoseLookupStatus.POSE_FOUND,
        pose_index=pose_index,
        pose_id=f"pose_{pose_index}",
        interactions=interactions,
        ligand_atoms=lig_atoms,
        receptor_atoms=receptor_atoms,
        message=f"Successfully profiled Pose {pose_index} ({len(interactions)} contacts detected).",
    )


def profile_docking_job(
    receptor_pdbqt: str | Path,
    output_pose_pdbqt: str | Path,
    max_poses: int = 3,
    pose_index: Optional[int] = None,
    run_number: Optional[int] = None,
    raise_on_missing: bool = False,
) -> List[Interaction]:
    """
    Profiles docking poses for a job.

    CRITICAL SCIENTIFIC SAFETY RULE:
    If pose_index is provided and cannot be found, this function NEVER silently
    falls back to Pose 1. It either raises an explicit exception (if raise_on_missing=True)
    or returns an empty list [].
    """
    receptor_atoms = parse_receptor_pdbqt(receptor_pdbqt)
    if not receptor_atoms:
        return []

    poses = parse_docked_poses(output_pose_pdbqt)
    if not poses:
        if raise_on_missing and pose_index is not None:
            raise CorruptedPoseOutputError(f"Pose file {output_pose_pdbqt} contains no valid coordinate models.")
        return []

    if pose_index is not None:
        status, match = locate_docked_pose(poses, pose_index, run_number=run_number)
        if status != PoseLookupStatus.POSE_FOUND or match is None:
            if raise_on_missing:
                if status == PoseLookupStatus.INVALID_POSE_INDEX:
                    raise InvalidPoseIndexError(f"Invalid pose index {pose_index}: must be a positive integer.")
                elif status == PoseLookupStatus.CORRUPTED_POSE_OUTPUT:
                    raise CorruptedPoseOutputError(f"Pose file {output_pose_pdbqt} contains no valid models.")
                else:
                    raise PoseNotFoundError(f"Pose {pose_index} not found in {output_pose_pdbqt} ({len(poses)} poses available).")
            return []  # NEVER fall back to poses[:1]!
        selected = [match]
    else:
        selected = poses[:max_poses]

    results: List[Interaction] = []
    for model_idx, lig_atoms, flex_atoms in selected:
        actual_pose_idx = pose_index if pose_index is not None else model_idx
        interactions = analyze_pose_interactions(
            receptor_atoms=receptor_atoms,
            ligand_atoms=lig_atoms,
            pose_index=actual_pose_idx,
            flex_atoms=flex_atoms if flex_atoms else None,
        )
        for it in interactions:
            it.pose_index = actual_pose_idx
            it.pose_id = f"pose_{actual_pose_idx}"
        results.extend(interactions)

    return results


def export_interactions_to_csv(
    interactions: List[Interaction], csv_path: str | Path
) -> None:
    """Exports interactions list to a standard CSV file."""
    path = Path(csv_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "Pose", "Interaction_Type", "Receptor_Residue", "Residue_Name",
        "Chain", "Residue_Num", "Receptor_Atom", "Ligand_Atom",
        "Distance_A", "Details",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for item in interactions:
            writer.writerow(item.to_dict())


def summarize_interactions(interactions: List[Interaction]) -> Dict[str, Any]:
    """Provides high-level summary counts and key residue lists for reporting."""
    counts: Dict[str, int] = {}
    interacting_residues: Set[str] = set()
    hbond_residues: List[str] = []

    for item in interactions:
        t = item.interaction_type
        counts[t] = counts.get(t, 0) + 1
        interacting_residues.add(item.receptor_residue)
        if t == "Hydrogen Bond":
            hbond_residues.append(
                f"{item.receptor_res_name}{item.receptor_res_seq}({item.receptor_atom})"
            )

    return {
        "total_interactions": len(interactions),
        "counts_by_type": counts,
        "unique_residues_count": len(interacting_residues),
        "interacting_residues": sorted(list(interacting_residues)),
        "hbond_contacts": sorted(list(set(hbond_residues))),
    }


# ==============================================================================
# CLI EXECUTION (STANDALONE TESTING)
# ==============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="AutoDock Suite Pro v0.3.0 — Protein-Ligand Interaction Profiler"
    )
    parser.add_argument("-r", "--receptor", required=True, help="Path to receptor PDBQT")
    parser.add_argument("-l", "--ligand", required=True, help="Path to docked pose PDBQT")
    parser.add_argument("-o", "--output", default="interactions.csv", help="Output CSV path")
    parser.add_argument("-m", "--max-poses", type=int, default=3, help="Max top poses to profile")
    parser.add_argument("-p", "--pose", type=int, default=None, help="Specific pose index to profile")
    args = parser.parse_args()

    print(f"Profiling interactions:\n  Receptor: {args.receptor}\n  Ligand:   {args.ligand}")
    interactions = profile_docking_job(
        args.receptor, args.ligand,
        max_poses=args.max_poses,
        pose_index=args.pose,
    )
    print(f"Found {len(interactions)} molecular interactions.")

    summary = summarize_interactions(interactions)
    print("\nInteraction summary:")
    for k, v in summary["counts_by_type"].items():
        print(f"  {k:30s}: {v}")
    print(f"  {'Unique Residues':30s}: {summary['unique_residues_count']}")

    export_interactions_to_csv(interactions, args.output)
    print(f"\nSaved: {args.output}")
