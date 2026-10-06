from pathlib import Path

from complex_builder import build_complex_pdbqt
from dlg_extract import DLGParser
from models import DockingJob, Engine, VinaResult


def test_vina_canonical_result_has_no_thermodynamic_ki():
    job = DockingJob(engine=Engine.VINA)
    job.vina_results = [VinaResult(pose=1, binding_affinity=-7.2, rmsd_lower_bound=0.0, rmsd_upper_bound=0.0)]
    result = job.ensure_canonical_result()
    pose = result.poses[0]
    assert pose.estimated_ki is None
    assert "not applicable" in pose.estimated_ki_formatted.lower()
    assert result.thermodynamics is None


def test_ad4_dlg_pose_values_are_pose_specific(tmp_path: Path):
    dlg = tmp_path / "job.dlg"
    dlg.write_text(
        """DOCKED: MODEL 1
DOCKED: Run = 1
DOCKED: Estimated Free Energy of Binding = -5.00 kcal/mol
DOCKED: Estimated Inhibition Constant, Ki = 200.0 nM
DOCKED: Final Intermolecular Energy = -4.00 kcal/mol
DOCKED: Torsional Free Energy = 0.50 kcal/mol
DOCKED: ATOM      1  C1  LIG     1       1.000   2.000   3.000  0.00  C
DOCKED: ENDMDL
DOCKED: MODEL 2
DOCKED: Run = 2
DOCKED: Estimated Free Energy of Binding = -6.00 kcal/mol
DOCKED: Estimated Inhibition Constant, Ki = 50.0 nM
DOCKED: Final Intermolecular Energy = -5.00 kcal/mol
DOCKED: Torsional Free Energy = 0.60 kcal/mol
DOCKED: ATOM      1  C1  LIG     1       4.000   5.000   6.000  0.00  C
DOCKED: ENDMDL
""",
        encoding="utf-8",
    )
    result = DLGParser().parse(dlg)
    assert len(result.poses) == 2
    assert result.poses[0].binding_energy == -5.0
    assert result.poses[1].binding_energy == -6.0
    assert result.poses[0].ki_nM == 200.0
    assert result.poses[1].ki_nM == 50.0
    assert result.poses[0].run_number == 1
    assert result.poses[1].run_number == 2


def test_complex_builder_preserves_selected_pose_coordinates(tmp_path: Path):
    receptor = tmp_path / "rec.pdbqt"
    output = tmp_path / "out.pdbqt"
    complex_file = tmp_path / "complex.pdbqt"
    receptor.write_text(
        "ATOM      1  ZN  ZN  A   1       0.000   0.000   0.000  1.00 20.00     0.000 Zn\nEND\n",
        encoding="utf-8",
    )
    output.write_text(
            "MODEL        1\nATOM      1  C1  LIG     1       1.000   2.000   3.000  1.00 20.00     0.000 C\nENDMDL\n"
            "MODEL        2\nATOM      1  C1  LIG     1       4.000   5.000   6.000  1.00 20.00     0.000 C\nENDMDL\n",
        encoding="utf-8",
    )
    assert build_complex_pdbqt(receptor, output, complex_file, pose_index=2)
    text = complex_file.read_text(encoding="utf-8")
    assert "REMARK Pose     : 2" in text
    assert "  4.000   5.000   6.000" in text
    assert "  1.000   2.000   3.000" not in text
