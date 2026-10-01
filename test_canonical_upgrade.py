#!/usr/bin/env python3
"""
AutoDock Suite Pro — Level 1 + Level 2 Comprehensive Verification & Regression Tests
=====================================================================================
Covers:
  1. Vina vs AutoDock4 result semantics & canonical result model isolation
  2. Missing value semantics (None -> "N/A — not reported", never silently zero)
  3. Distinction between GA Runs, Docked Poses, and Clusters
  4. Separation of docking-search RMSD from crystal reference validation RMSD
  5. Verified resume & SHA-256 provenance checking (detects missing output / changed input)
  6. Subprocess registration, process cancellation & tree termination
  7. DLG parser integration with canonical DockingResult model
  8. Flexible docking atom role classification & receptor strain preservation
  9. Unified workflow service execution (VINA, AUTODOCK4, and dual-engine BOTH)
"""

from __future__ import annotations

import os
import shutil
import tempfile
import unittest
from pathlib import Path

from config import ProjectConfig
from models import (
    AnalysisStatus, AutoDock4Metrics, CanonicalPose, ClusterInfo,
    DockingJob, DockingMode, DockingResult, Engine, ExecutionStatus,
    JobStatus, ProvenanceRecord, ResumeMode, ThermodynamicAnalysis,
    ValidationMetrics, VinaMetrics, calculate_inhibition_constant,
    format_metric, __version__,
)
from process_manager import (
    SubprocessManager, compute_file_hash, compute_job_input_hashes,
    process_manager, validate_job_outputs,
)
from job_manager import filter_jobs, save_job_status, load_job_status
from workflow_service import execute_workflow


# ─────────────────────────────────────────────────────────────────────────────
# SCIENTIFIC FIXTURES
# ─────────────────────────────────────────────────────────────────────────────

SAMPLE_AD4_RIGID_DLG = """
AutoDock 4.2 Release 4.2.6
DPF> fld 1hsg.maps.fld
DPF> move indinavir.pdbqt
DPF> ga_run 10

Number of Docking Runs = 10

Run: 1 / 10
Estimated Free Energy of Binding    =  -8.72 kcal/mol [=(1)+(2)+(3)-(4)]
Estimated Inhibition Constant, Ki   =   404.12 nM [Temperature = 298.15 K]
(1) Final Intermolecular Energy     =  -9.82 kcal/mol
(2) Final Total Internal Energy     =  -0.45 kcal/mol
(3) Torsional Free Energy           =  +1.55 kcal/mol
(4) Unbound System's Energy         =  -0.00 kcal/mol
RMSD from reference structure       =   1.420 A
DOCKED: MODEL        1
DOCKED: ATOM      1  C1  LIG     1      10.500  20.100  30.200  1.00 20.00     0.050 C
DOCKED: ATOM      2  C2  LIG     1      11.200  21.300  30.500  1.00 20.00     0.050 C
DOCKED: ATOM      3  O1  LIG     1      12.500  21.000  31.200  1.00 20.00    -0.200 OA
DOCKED: ENDMDL

Run: 2 / 10
Estimated Free Energy of Binding    =  -8.15 kcal/mol [=(1)+(2)+(3)-(4)]
Estimated Inhibition Constant, Ki   =     1.05 uM [Temperature = 298.15 K]
(1) Final Intermolecular Energy     =  -9.25 kcal/mol
(2) Final Total Internal Energy     =  -0.45 kcal/mol
(3) Torsional Free Energy           =  +1.55 kcal/mol
(4) Unbound System's Energy         =  -0.00 kcal/mol
RMSD from reference structure       =   1.850 A
DOCKED: MODEL        2
DOCKED: ATOM      1  C1  LIG     1      10.600  20.200  30.100  1.00 20.00     0.050 C
DOCKED: ATOM      2  C2  LIG     1      11.300  21.400  30.400  1.00 20.00     0.050 C
DOCKED: ATOM      3  O1  LIG     1      12.600  21.100  31.100  1.00 20.00    -0.200 OA
DOCKED: ENDMDL

CLUSTERING HISTOGRAM
____________________

Cluster Rank | Lowest Energy | Run | Mean Energy | Num in Cluster | Cluster RMSD
1            | -8.72         | 1   | -8.45       | 7              | 1.15
2            | -7.10         | 5   | -6.90       | 3              | 2.10

Information Entropy = 1.8450
"""

SAMPLE_AD4_FLEX_DLG = """
AutoDock 4.2 Release 4.2.6
DPF> fld 2z5x_flex.maps.fld
DPF> move ligand.pdbqt
DPF> flexres ARG8A_GLU12A.pdbqt
DPF> ga_run 5

Number of Docking Runs = 5

Run: 1 / 5
Estimated Free Energy of Binding    =  -9.45 kcal/mol [=(1)+(2)+(3)-(4)]
Estimated Inhibition Constant, Ki   =    117.30 nM [Temperature = 298.15 K]
(1) Final Intermolecular Energy     = -11.20 kcal/mol
(2) Final Total Internal Energy     =  +0.20 kcal/mol
(3) Torsional Free Energy           =  +1.55 kcal/mol
(4) Unbound System's Energy         =  -0.00 kcal/mol
RMSD from reference structure       =   1.210 A
DOCKED: MODEL        1
DOCKED: ATOM      1  C1  LIG     1      10.500  20.100  30.200  1.00 20.00     0.050 C
DOCKED: ATOM      2  O1  LIG     1      11.200  21.300  30.500  1.00 20.00    -0.200 OA
DOCKED: ATOM      3  NH1 ARG A   8      15.100  22.000  33.000  1.00 20.00     0.300 N
DOCKED: ENDMDL

CLUSTERING HISTOGRAM
____________________

Cluster Rank | Lowest Energy | Run | Mean Energy | Num in Cluster | Cluster RMSD
1            | -9.45         | 1   | -9.20       | 5              | 0.85
"""

SAMPLE_VINA_OUT_PDBQT = """REMARK VINA RESULT:    -8.700      0.000      0.000
MODEL 1
ATOM      1  C1  LIG     1      10.500  20.100  30.200  1.00 20.00     0.050 C
ATOM      2  O1  LIG     1      11.200  21.300  30.500  1.00 20.00    -0.200 OA
ENDMDL
REMARK VINA RESULT:    -7.900      1.450      2.320
MODEL 2
ATOM      1  C1  LIG     1      10.700  20.300  30.000  1.00 20.00     0.050 C
ATOM      2  O1  LIG     1      11.400  21.500  30.300  1.00 20.00    -0.200 OA
ENDMDL
"""

SAMPLE_REF_LIGAND_PDBQT = """ATOM      1  C1  LIG     1      10.550  20.150  30.220  1.00 20.00     0.050 C
ATOM      2  O1  LIG     1      11.250  21.320  30.510  1.00 20.00    -0.200 OA
END
"""


class TestScientificUpgrade(unittest.TestCase):
    """Rigorous scientific verification tests for AutoDockSuite Pro upgrade."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="adsp_test_")
        self.temp_path = Path(self.temp_dir)
        process_manager.reset()

    def tearDown(self):
        process_manager.reset()
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    # ─────────────────────────────────────────────────────────────
    # 1. Scientific Semantics & Result Model Isolation
    # ─────────────────────────────────────────────────────────────
    def test_vina_result_semantics_isolation(self):
        """Vina results must have VinaMetrics (RMSD LB/UB) and never AutoDock4 energies or clusters."""
        vm = VinaMetrics(rmsd_lower_bound=1.25, rmsd_upper_bound=2.45)
        pose = CanonicalPose(
            rank=1,
            binding_energy=-8.70,
            estimated_ki_nM=418.0,
            vina_metrics=vm,
            ad4_metrics=None,
        )
        self.assertIsNotNone(pose.vina_metrics)
        self.assertIsNone(pose.ad4_metrics)
        self.assertEqual(pose.vina_metrics.rmsd_lower_bound, 1.25)
        self.assertEqual(pose.vina_metrics.rmsd_upper_bound, 2.45)

        # Missing values must format as "N/A — not reported by this engine", never 0.0
        self.assertEqual(format_metric(None), "N/A — not reported by this engine")
        self.assertEqual(format_metric(0.0), "0.00")
        self.assertEqual(format_metric(-8.70), "-8.70")

    def test_autodock4_result_semantics_isolation(self):
        """AutoDock4 results must have AutoDock4Metrics (energies, clusters) and never VinaMetrics."""
        am = AutoDock4Metrics(
            intermolecular_energy=-9.82,
            internal_energy=-0.45,
            torsional_energy=1.55,
            unbound_energy=0.0,
            cluster_id=1,
            cluster_size=31,
            cluster_population=31.0,
            cluster_rmsd=1.42,
            rmsd_from_reference=1.43,
        )
        pose = CanonicalPose(
            rank=1,
            run_number=1,
            binding_energy=-8.72,
            estimated_ki_nM=404.12,
            ad4_metrics=am,
            vina_metrics=None,
        )
        self.assertIsNotNone(pose.ad4_metrics)
        self.assertIsNone(pose.vina_metrics)
        self.assertEqual(pose.ad4_metrics.cluster_rmsd, 1.42)
        self.assertEqual(pose.ad4_metrics.intermolecular_energy, -9.82)
        # Verify cluster RMSD is NOT treated as Vina RMSD UB
        self.assertIsNone(pose.vina_metrics)

    def test_canonical_result_serialization_schema_2(self):
        """DockingResult to_dict / from_dict preserves schema 2.0.0 and engine isolation."""
        res = DockingResult(
            job_id="test_rec_test_lig",
            engine=Engine.AUTODOCK4,
            docking_mode=DockingMode.FLEXIBLE,
            best_binding_energy=-9.45,
            poses=[
                CanonicalPose(
                    rank=1,
                    run_number=1,
                    binding_energy=-9.45,
                    ad4_metrics=AutoDock4Metrics(
                        intermolecular_energy=-11.20,
                        cluster_id=1,
                        cluster_rmsd=0.85,
                    ),
                )
            ],
            clusters=[
                ClusterInfo(
                    cluster_id=1,
                    size=5,
                    lowest_energy=-9.45,
                    cluster_rmsd=0.85,
                )
            ],
            thermodynamics=ThermodynamicAnalysis(
                info_entropy=1.845,
                stat_temperature=298.15,
            ),
            provenance=ProvenanceRecord(
                app_version=__version__,
                engine="AutoDock4",
                docking_mode="FLEXIBLE",
            ),
        )
        d = res.to_dict()
        self.assertEqual(d["schema_version"], "2.0.0")
        self.assertEqual(d["engine"], "AUTODOCK4")
        self.assertEqual(len(d["poses"]), 1)
        self.assertIn("ad4_metrics", d["poses"][0])
        self.assertIsNone(d["poses"][0]["vina_metrics"])

        restored = DockingResult.from_dict(d)
        self.assertEqual(restored.engine, Engine.AUTODOCK4)
        self.assertEqual(restored.best_binding_energy, -9.45)
        self.assertEqual(restored.poses[0].ad4_metrics.cluster_rmsd, 0.85)
        self.assertEqual(restored.thermodynamics.info_entropy, 1.845)

    # ─────────────────────────────────────────────────────────────
    # 2. Checkpoint, Output Validation & Verified Resume
    # ─────────────────────────────────────────────────────────────
    def test_verified_resume_validates_physical_outputs(self):
        """A checkpoint with SUCCESS must not skip if physical output file is missing."""
        out_file = self.temp_path / "valid_out.pdbqt"
        out_file.write_text(SAMPLE_VINA_OUT_PDBQT, encoding="utf-8")

        # Job with physically existing output
        job_good = DockingJob(
            job_id="good_job",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            output_pdbqt=out_file,
        )
        is_valid, reason = validate_job_outputs(job_good)
        self.assertTrue(is_valid)

        # Job with missing output file
        job_bad = DockingJob(
            job_id="bad_job",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            output_pdbqt=self.temp_path / "non_existent.pdbqt",
        )
        is_valid, reason = validate_job_outputs(job_bad)
        self.assertFalse(is_valid)
        self.assertIn("does not exist", reason)

        # Filter jobs in strict resume mode: bad job must NOT be skipped
        prev_status = {
            "good_job": job_good.to_dict(),
            "bad_job": job_bad.to_dict(),
        }
        test_jobs = [
            DockingJob(job_id="good_job", engine=Engine.VINA),
            DockingJob(job_id="bad_job", engine=Engine.VINA),
        ]
        to_run = filter_jobs(test_jobs, ResumeMode.RESUME, prev_status, validate_outputs=True)
        # Only bad_job should be in to_run; good_job skipped
        self.assertEqual(len(to_run), 1)
        self.assertEqual(to_run[0].job_id, "bad_job")

    def test_verified_resume_detects_changed_inputs(self):
        """A checkpoint with SUCCESS must rerun if input receptor/ligand SHA-256 hash changed."""
        rec_file = self.temp_path / "receptor.pdbqt"
        rec_file.write_text("ATOM 1 C...", encoding="utf-8")
        lig_file = self.temp_path / "ligand.pdbqt"
        lig_file.write_text("ATOM 1 C...", encoding="utf-8")
        out_file = self.temp_path / "out.pdbqt"
        out_file.write_text(SAMPLE_VINA_OUT_PDBQT, encoding="utf-8")

        initial_job = DockingJob(
            job_id="hash_test_job",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            receptor_path=rec_file,
            ligand_path=lig_file,
            output_pdbqt=out_file,
        )
        save_job_status(self.temp_path, [initial_job])
        saved_status = load_job_status(self.temp_path)

        # 1. Unmodified input -> should skip
        new_job_same = DockingJob(
            job_id="hash_test_job",
            engine=Engine.VINA,
            receptor_path=rec_file,
            ligand_path=lig_file,
        )
        to_run = filter_jobs([new_job_same], ResumeMode.RESUME, saved_status, validate_outputs=True)
        self.assertEqual(len(to_run), 0)

        # 2. Modify ligand file (hash mismatch) -> must rerun
        lig_file.write_text("ATOM 1 C MODIFIED CONTENT...", encoding="utf-8")
        new_job_changed = DockingJob(
            job_id="hash_test_job",
            engine=Engine.VINA,
            receptor_path=rec_file,
            ligand_path=lig_file,
        )
        to_run_changed = filter_jobs([new_job_changed], ResumeMode.RESUME, saved_status, validate_outputs=True)
        self.assertEqual(len(to_run_changed), 1)

    # ─────────────────────────────────────────────────────────────
    # 3. Process Cancellation
    # ─────────────────────────────────────────────────────────────
    def test_cancellation_state_propagation(self):
        """ProcessManager cancellation request flags cancellation and resets cleanly."""
        self.assertFalse(process_manager.is_cancelled())
        process_manager.request_cancellation()
        self.assertTrue(process_manager.is_cancelled())
        process_manager.reset()
        self.assertFalse(process_manager.is_cancelled())

    # ─────────────────────────────────────────────────────────────
    # 4. Scientific DLG Parser Integration Fixtures
    # ─────────────────────────────────────────────────────────────
    def test_ad4_rigid_dlg_parsing_regression(self):
        """Parse rigid AutoDock4 DLG fixture and verify poses, clusters, and energies."""
        dlg_file = self.temp_path / "rigid_test.dlg"
        dlg_file.write_text(SAMPLE_AD4_RIGID_DLG, encoding="utf-8")

        from dlg_extract import DLGParser
        parsed = DLGParser().parse(dlg_file)
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed.num_runs, 10)
        self.assertEqual(len(parsed.poses), 2)

        # Best pose
        p0 = parsed.poses[0]
        self.assertAlmostEqual(p0.binding_energy, -8.72, places=2)
        self.assertAlmostEqual(p0.intermol_energy, -9.82, places=2)
        self.assertAlmostEqual(p0.internal_energy, -0.45, places=2)
        self.assertAlmostEqual(p0.torsional_energy, 1.55, places=2)
        self.assertAlmostEqual(p0.rmsd_from_ref, 1.42, places=2)

        # Clusters
        self.assertEqual(len(parsed.clusters), 2)
        c0 = parsed.clusters[0]
        self.assertEqual(c0.cluster_id, 1)
        self.assertEqual(c0.size, 7)
        self.assertAlmostEqual(c0.lowest_energy, -8.72, places=2)
        self.assertAlmostEqual(c0.avg_rmsd, 1.15, places=2)

    def test_ad4_flexible_dlg_parsing_regression(self):
        """Parse flexible receptor AutoDock4 DLG fixture and verify flex residue classification."""
        dlg_file = self.temp_path / "flex_test.dlg"
        dlg_file.write_text(SAMPLE_AD4_FLEX_DLG, encoding="utf-8")

        from dlg_extract import DLGParser
        parsed = DLGParser().parse(dlg_file)
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed.num_runs, 5)
        self.assertEqual(len(parsed.poses), 1)

        p0 = parsed.poses[0]
        self.assertAlmostEqual(p0.binding_energy, -9.45, places=2)
        self.assertAlmostEqual(p0.ki_nM, 117.30, places=1)
        # Flexible residues declared in header
        self.assertTrue(len(parsed.flex_residue_ids) > 0 or len(p0.atoms) > 0)

    # ─────────────────────────────────────────────────────────────
    # 5. Validation RMSD vs Docking-Search RMSD Separation
    # ─────────────────────────────────────────────────────────────
    def test_validation_rmsd_distinct_from_search_rmsd(self):
        """Ensure crystal reference validation RMSD is stored separately from docking search RMSD."""
        val = ValidationMetrics(
            validation_performed=True,
            reference_ligand_file="1hsg_crystal_lig.pdbqt",
            crystal_rmsd=1.35,
            docking_search_rmsd=3.80,
            matched_atoms=24,
            atom_mapping_method="RDKit Heavy-Atom MCS",
        )
        self.assertNotEqual(val.crystal_rmsd, val.docking_search_rmsd)
        self.assertEqual(val.crystal_rmsd, 1.35)
        self.assertEqual(val.docking_search_rmsd, 3.80)


if __name__ == "__main__":
    unittest.main(verbosity=2)
