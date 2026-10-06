"""Reproducible 2NV6 / INH-NAD benchmark inventory and preparation helpers.

The benchmark deliberately keeps crystallographic extraction and independent
SDF input as separate scenarios.  Engine execution is delegated to the same
workflow services used by the CLI; this module does not implement a second
docking path.
"""
from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class BenchmarkInput:
    scenario: str
    path: Path
    sha256: str
    atom_count: int
    residue_name: Optional[str] = None
    chain_id: Optional[str] = None
    residue_seq: Optional[int] = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_pdb_component(source: Path, residue_name: str = "ZID") -> str:
    """Return the exact ATOM/HETATM records for one PDB component.

    The component identity is selected by residue name and the first matching
    chain/residue tuple.  A missing or ambiguous component raises instead of
    silently extracting unrelated HETATM records.
    """
    rows = [line for line in source.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip().upper() == residue_name.upper()]
    keys = {(line[21].strip(), line[22:26].strip()) for line in rows}
    if not rows:
        raise ValueError(f"component {residue_name!r} was not found in {source}")
    if len(keys) != 1:
        raise ValueError(f"component {residue_name!r} is ambiguous: {sorted(keys)}")
    return "\n".join(rows) + "\n"


def discover(root: Path) -> dict[str, BenchmarkInput]:
    pdb = root / "2NV6 (2).pdb"
    cif = root / "2NV6 (1).cif"
    sdf = root / "2nv6_B_ZID.sdf"
    if not pdb.is_file() or not cif.is_file():
        raise FileNotFoundError("2NV6 PDB/CIF benchmark inputs are required")
    component = extract_pdb_component(pdb)
    component_atoms = sum(1 for line in component.splitlines() if line.startswith(("ATOM", "HETATM")))
    result = {
        "crystallographic_pdb": BenchmarkInput("crystallographic", pdb, sha256(pdb), component_atoms, "ZID", "A", 300),
        "crystallographic_cif": BenchmarkInput("crystallographic", cif, sha256(cif), component_atoms, "ZID", "A", 300),
    }
    if sdf.is_file():
        try:
            from rdkit import Chem
            mol = Chem.MolFromMolFile(str(sdf), removeHs=False, sanitize=False)
            atom_count = mol.GetNumAtoms() if mol is not None else 0
        except Exception:
            atom_count = 0
        result["independent_sdf"] = BenchmarkInput("independent_sdf", sdf, sha256(sdf), atom_count)
    return result


def write_inventory(root: Path, output: Path) -> Path:
    import json
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {key: {**asdict(value), "path": str(value.path)} for key, value in discover(root).items()}
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return output


def compare_ligand_representations(root: Path) -> dict:
    """Compare the supplied independent SDF with the ZID component dictionary."""
    from rdkit import Chem
    sdf = root / "2nv6_B_ZID.sdf"
    mol = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True)[0]
    import gemmi
    block = gemmi.cif.read(str(root / "ZID.cif")).sole_block()
    atom_table = block.find_mmcif_category("_chem_comp_atom.")
    bond_table = block.find_mmcif_category("_chem_comp_bond.")
    cif_elements = [row[3] for row in atom_table]
    return {
        "classification": "REPRESENTATION_DIFFERENCE" if len(cif_elements) != mol.GetNumAtoms() else "UNRESOLVED",
        "independent_sdf": {
            "atoms": mol.GetNumAtoms(), "heavy_atoms": mol.GetNumHeavyAtoms(),
            "bonds": mol.GetNumBonds(), "formal_charge": sum(a.GetFormalCharge() for a in mol.GetAtoms()),
            "elements": sorted({a.GetSymbol() for a in mol.GetAtoms()}),
            "connected_components": len(Chem.GetMolFrags(mol)),
        },
        "zid_dictionary": {
            "atoms_including_explicit_hydrogen": len(cif_elements),
            "bonds": len(bond_table), "elements": sorted(set(cif_elements)),
            "formal_charge": 0, "connected_components": "NOT_DECLARED",
        },
        "warning": "ZID.cif contains explicit hydrogen records while the supplied SDF uses implicit hydrogen representation.",
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Inventory the ADSP 2NV6 benchmark inputs")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    out = args.output or args.root / "results" / "2nv6_benchmark_inventory.json"
    print(write_inventory(args.root, out))
