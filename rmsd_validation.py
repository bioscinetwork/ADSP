"""Chemically constrained reference-to-pose heavy-atom RMSD utilities.

These functions require complete element/bond graph correspondence.  They do
not guess based on atom order and never return an RMSD for a partial match.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Optional


class MappingStatus(str, Enum):
    SUCCESS = "SUCCESS"
    PARTIAL_MAPPING = "PARTIAL_MAPPING"
    AMBIGUOUS_MAPPING = "AMBIGUOUS_MAPPING"
    MAPPING_FAILED = "MAPPING_FAILED"
    INCOMPATIBLE_STRUCTURES = "INCOMPATIBLE_STRUCTURES"


@dataclass(frozen=True)
class RMSDResult:
    mapping_status: MappingStatus
    mapping_method: str
    rmsd_method: str
    reference_atoms: int
    target_atoms: int
    reference_heavy_atoms: int
    target_heavy_atoms: int
    mapped_heavy_atoms: int
    unmapped_reference_atoms: tuple[int, ...]
    unmapped_target_atoms: tuple[int, ...]
    rmsd_angstrom: Optional[float]
    alternate_mappings: int = 0
    message: str = ""

    def to_dict(self):
        data = asdict(self)
        data["mapping_status"] = self.mapping_status.value
        return data


def heavy_atom_rmsd(reference, target, *, fit: bool = False, max_matches: int = 10000) -> RMSDResult:
    """Map complete heavy-atom molecular graphs and calculate RMSD.

    Atom order is irrelevant. Element, formal charge, aromaticity, bond order,
    and connectivity constrain mapping. Symmetry-equivalent complete graph
    mappings are enumerated and the minimum RMSD is returned with the number
    of alternatives recorded. Multi-component structures are allowed only if
    their complete disconnected graph topology matches.
    """
    import math
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    if reference is None or target is None:
        raise ValueError("Both reference and target molecules are required")
    ref = Chem.RemoveHs(Chem.Mol(reference), sanitize=True)
    tgt = Chem.RemoveHs(Chem.Mol(target), sanitize=True)
    nr, nt = ref.GetNumAtoms(), tgt.GetNumAtoms()
    rh, th = nr, nt

    def fail(status, message):
        return RMSDResult(status, "RDKit complete graph isomorphism", "fitted" if fit else "coordinate",
                          reference.GetNumAtoms(), target.GetNumAtoms(), rh, th, 0,
                          tuple(range(nr)), tuple(range(nt)), None, 0, message)

    if nr != nt:
        return fail(MappingStatus.INCOMPATIBLE_STRUCTURES,
                    f"Heavy-atom counts differ ({nr} reference, {nt} target); partial RMSD is not reported")
    if nr == 0 or ref.GetNumBonds() != tgt.GetNumBonds():
        return fail(MappingStatus.INCOMPATIBLE_STRUCTURES, "Empty or different bond-count molecular graphs")
    if not ref.GetNumConformers() or not tgt.GetNumConformers():
        return fail(MappingStatus.MAPPING_FAILED, "Both molecules need coordinates")

    matches = tgt.GetSubstructMatches(ref, uniquify=False, useChirality=True, maxMatches=max_matches)
    matches = [m for m in matches if len(m) == nr]
    if not matches:
        return fail(MappingStatus.MAPPING_FAILED, "No complete chemically constrained atom mapping exists")
    ref_conf, tgt_conf = ref.GetConformer(), tgt.GetConformer()
    ref_xyz = np.asarray([list(ref_conf.GetAtomPosition(i)) for i in range(nr)], dtype=float)
    results = []
    for match in matches:
        target_xyz = np.asarray([list(tgt_conf.GetAtomPosition(j)) for j in match], dtype=float)
        if fit:
            # RDKit's alignment operation requires a writable target copy and
            # an explicit atom map from reference indices to target indices.
            aligned = Chem.Mol(tgt)
            atom_map = [(int(j), int(i)) for i, j in enumerate(match)]
            value = rdMolAlign.AlignMol(aligned, ref, atomMap=atom_map)
            target_xyz = np.asarray([list(aligned.GetConformer().GetAtomPosition(j)) for j in match], dtype=float)
        else:
            delta = ref_xyz - target_xyz
            value = math.sqrt(float(np.square(delta).sum(axis=1).mean()))
        results.append(float(value))
    best = min(results)
    used_mapping_count = sum(abs(x - best) <= 1e-8 for x in results)
    return RMSDResult(MappingStatus.SUCCESS, "RDKit complete graph isomorphism",
                      "fitted" if fit else "coordinate", reference.GetNumAtoms(), target.GetNumAtoms(),
                      rh, th, nr, (), (), round(best, 6), used_mapping_count,
                      "Symmetry-equivalent mappings were evaluated" if used_mapping_count > 1 else "")


def load_rdkit_molecule(path):
    """Load supported small-molecule files without silently choosing a pose."""
    from pathlib import Path
    from rdkit import Chem
    p = Path(path)
    ext = p.suffix.lower()
    if ext in {".sdf", ".mol"}:
        supplier = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=True)
        return next((m for m in supplier if m is not None), None)
    if ext == ".mol2":
        return Chem.MolFromMol2File(str(p), removeHs=False, sanitize=True)
    if ext in {".pdb", ".cif", ".mmcif"}:
        return Chem.MolFromPDBFile(str(p), removeHs=False, sanitize=True)
    raise ValueError(f"Unsupported reference format: {ext}")
