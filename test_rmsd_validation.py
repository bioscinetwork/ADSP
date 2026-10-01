from rdkit import Chem
from rdkit.Chem import AllChem

from rmsd_validation import MappingStatus, heavy_atom_rmsd


def mol(smiles, seed=31):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(m, randomSeed=seed) == 0
    return m


def test_identical_coordinates_have_zero_coordinate_rmsd():
    m = mol("CCO")
    result = heavy_atom_rmsd(m, m)
    assert result.mapping_status is MappingStatus.SUCCESS
    assert result.mapped_heavy_atoms == 3
    assert result.rmsd_method == "coordinate"
    assert result.rmsd_angstrom == 0


def test_translation_differs_without_fit_and_matches_after_fit():
    a = mol("CCO")
    b = Chem.Mol(a)
    conf = b.GetConformer()
    for i in range(b.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 4, p.y - 2, p.z + 7))
    direct = heavy_atom_rmsd(a, b)
    fitted = heavy_atom_rmsd(a, b, fit=True)
    assert direct.rmsd_angstrom > 1
    assert fitted.rmsd_angstrom < 1e-5


def test_atom_order_permutation_is_graph_mapped():
    a = mol("CC(=O)O")
    permutation = list(reversed(range(a.GetNumAtoms())))
    b = Chem.RenumberAtoms(a, permutation)
    result = heavy_atom_rmsd(a, b)
    assert result.mapping_status is MappingStatus.SUCCESS
    assert result.rmsd_angstrom < 1e-5


def test_hydrogens_are_excluded_from_heavy_atom_rmsd():
    a = mol("CCO")
    b = Chem.RemoveHs(a)
    result = heavy_atom_rmsd(a, b)
    assert result.mapping_status is MappingStatus.SUCCESS
    assert result.reference_atoms == a.GetNumAtoms()
    assert result.target_atoms == b.GetNumAtoms()
    assert result.mapped_heavy_atoms == 3


def test_removed_atom_fails_instead_of_partial_rmsd():
    a = mol("CCO")
    b = Chem.RemoveHs(a)
    b = Chem.RWMol(b)
    b.RemoveAtom(2)
    result = heavy_atom_rmsd(a, b.GetMol())
    assert result.mapping_status is MappingStatus.INCOMPATIBLE_STRUCTURES
    assert result.rmsd_angstrom is None


def test_element_substitution_fails_mapping():
    a = mol("CCO")
    b = mol("CCN")
    result = heavy_atom_rmsd(a, b)
    assert result.mapping_status is MappingStatus.MAPPING_FAILED
    assert result.rmsd_angstrom is None


def test_symmetric_mapping_is_examined():
    a = mol("CC(C)C")
    b = Chem.Mol(a)
    result = heavy_atom_rmsd(a, b)
    assert result.mapping_status is MappingStatus.SUCCESS
    assert result.alternate_mappings >= 1
