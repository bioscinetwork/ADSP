from pathlib import Path

from dlg_extract import DockingProfile, ReportingEngine


def test_txt_report_accepts_missing_boltzmann_probability(tmp_path: Path):
    profile = DockingProfile(
        ligand="ligand",
        protein="receptor",
        binding_energy=-1.5,
        boltzmann_prob=None,
    )
    output = tmp_path / "report.txt"
    ReportingEngine()._txt_report([profile], [], output)
    text = output.read_text(encoding="utf-8")
    assert "Boltzmann Probability (pi)       : Not available" in text
    assert "Boltzmann Population (π)   : Not available" in text

