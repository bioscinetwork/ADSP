from pathlib import Path

from validators import validate_pdbqt


def test_pdbqt_fixed_columns_are_validated(tmp_path: Path):
    valid = tmp_path / "valid.pdbqt"
    valid.write_text(
        "ATOM      1  C1  LIG     1       1.000   2.000   3.000  1.00 20.00     0.000 C\nEND\n",
        encoding="utf-8",
    )
    assert validate_pdbqt(valid)[0]

    malformed = tmp_path / "malformed.pdbqt"
    malformed.write_text(
        "ATOM      1  N   UNL     1      50.360 158.754  29.366 +0.46 +2.80    -0.309    165.230\n",
        encoding="utf-8",
    )
    ok, errors = validate_pdbqt(malformed)
    assert not ok
    assert any("AutoDock atom type" in error for error in errors)

