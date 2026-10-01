from pathlib import Path

from validators import inspect_pdbqt_charges


def test_charge_completeness_and_total(tmp_path):
    path = tmp_path / "charged.pdbqt"
    path.write_text(
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.250 C  \n"
        "ATOM      2  O   LIG A   1       1.000   0.000   0.000  1.00  0.00    -0.250 OA \n"
    )
    result = inspect_pdbqt_charges(path)
    assert result["completeness_fraction"] == 1.0
    assert result["calculated_total_charge"] == 0.0
    assert result["expected_charge_status"] == "NOT_DECLARED"
    assert result["scientific_correctness_established"] is False


def test_charge_missing_and_malformed_are_reported(tmp_path):
    path = tmp_path / "bad.pdbqt"
    path.write_text(
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00                C  \n"
        "ATOM      2  O   LIG A   1       1.000   0.000   0.000  1.00  0.00       NaN      OA \n"
    )
    result = inspect_pdbqt_charges(path)
    assert result["atoms_without_charge"]
    assert result["atoms_with_invalid_charge"]
    assert result["calculated_total_charge"] is None
