from pathlib import Path

from prepare import cleanup_provenance


def test_cleanup_provenance_separates_water_and_ligand(tmp_path):
    source = tmp_path / "rec.pdb"
    source.write_text(
        "HETATM    1  O   HOH A 101       0.000   0.000   0.000  1.00 10.00           O  \n"
        "HETATM    2  C1  ZID A 300       1.000   0.000   0.000  1.00 10.00           C  \n"
    )
    result = cleanup_provenance(source)
    assert result["removed_components"][0]["resname"] == "HOH"
    assert result["removed_components"][0]["reason"] == "default_cleanup_rule"
    assert "cleanup_rules" in result
    assert "user_overrides" in result
    assert "warnings" in result
    assert result["retained_components"][0]["resname"] == "ZID"
    assert result["scientific_component_relevance_established"] is False


def test_cleanup_provenance_records_user_removal(tmp_path):
    source = tmp_path / "rec.pdb"
    source.write_text("HETATM    1  C1  ZID A 300       1.000   0.000   0.000  1.00 10.00           C  \n")
    result = cleanup_provenance(source, removed_keys={("ZID", "A", 300)})
    assert result["removed_components"][0]["reason"] == "user_selected_removal"
