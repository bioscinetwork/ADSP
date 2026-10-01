from pathlib import Path

from benchmark_2nv6 import compare_ligand_representations, discover, extract_pdb_component


ROOT = Path(__file__).resolve().parent


def test_2nv6_component_is_unambiguous_and_complete():
    component = extract_pdb_component(ROOT / "2NV6 (2).pdb")
    atoms = [line for line in component.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert len(atoms) == 52
    assert {(line[17:20].strip(), line[21].strip(), line[22:26].strip()) for line in atoms} == {("ZID", "A", "300")}


def test_2nv6_inventory_keeps_independent_sdf_separate():
    inventory = discover(ROOT)
    assert inventory["crystallographic_pdb"].atom_count == 52
    assert inventory["crystallographic_cif"].residue_name == "ZID"
    assert inventory["independent_sdf"].scenario == "independent_sdf"
    assert inventory["independent_sdf"].atom_count == 52


def test_2nv6_dictionary_comparison_exposes_explicit_hydrogen_difference():
    comparison = compare_ligand_representations(ROOT)
    assert comparison["classification"] == "REPRESENTATION_DIFFERENCE"
    assert comparison["independent_sdf"]["atoms"] == 52
    assert comparison["zid_dictionary"]["atoms_including_explicit_hydrogen"] == 82
