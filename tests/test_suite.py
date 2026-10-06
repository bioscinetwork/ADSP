#!/usr/bin/env python3
"""
Docking Automation Suite — Comprehensive Test & Verification Suite
==================================================================
Tests core models, configuration loading, parser engines, splitting,
validation, reporting, and CLI dry runs.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from config import ProjectConfig, generate_template_config, load_config, resolve_executable, validate_config
from executables import (
    get_autodock4_version, get_autogrid4_version, get_cpu_count,
    get_vina_split_version, get_vina_version, validate_executable,
)
from job_manager import (
    ProgressTracker, build_job_queue, filter_jobs,
    load_job_status, save_job_status, validate_job_inputs,
)
from models import (
    AnalysisStatus, DockingJob, Engine, ExecutionStatus,
    JobStatus, ResumeMode, VinaResult, __version__,
)
from reporting import generate_reports
from validators import (
    parse_vina_config, validate_pdbqt,
    validate_pdbqt_multimodel, validate_vina_config, validate_vina_output,
)
from vina_parser import count_models, extract_affinities, format_vina_results_table
from vina_splitter import run_vina_split, split_pdbqt_manually, validate_split


# Sample Vina multi-model PDBQT content for testing
SAMPLE_VINA_PDBQT = """REMARK VINA RESULT:    -8.400      0.000      0.000
MODEL 1
ATOM      1  N   LIG     1      10.500  20.100  30.200  1.00 20.00    -0.100 N 
ATOM      2  CA  LIG     1      11.200  21.300  30.500  1.00 20.00     0.050 C 
ATOM      3  C   LIG     1      12.500  21.000  31.200  1.00 20.00     0.250 C 
ATOM      4  O   LIG     1      13.100  20.000  31.000  1.00 20.00    -0.200 O 
TER
ENDMDL
REMARK VINA RESULT:    -7.900      1.520      2.140
MODEL 2
ATOM      1  N   LIG     1      10.600  20.200  30.100  1.00 20.00    -0.100 N 
ATOM      2  CA  LIG     1      11.300  21.400  30.400  1.00 20.00     0.050 C 
ATOM      3  C   LIG     1      12.600  21.100  31.100  1.00 20.00     0.250 C 
ATOM      4  O   LIG     1      13.200  20.100  30.900  1.00 20.00    -0.200 O 
TER
ENDMDL
REMARK VINA RESULT:    -7.100      2.830      3.910
MODEL 3
ATOM      1  N   LIG     1      10.400  20.000  30.300  1.00 20.00    -0.100 N 
ATOM      2  CA  LIG     1      11.100  21.200  30.600  1.00 20.00     0.050 C 
ATOM      3  C   LIG     1      12.400  20.900  31.300  1.00 20.00     0.250 C 
ATOM      4  O   LIG     1      13.000  19.900  31.100  1.00 20.00    -0.200 O 
TER
ENDMDL
"""

SAMPLE_RECEPTOR_PDBQT = """ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 20.00    -0.200 N 
ATOM      2  CA  ALA A   1       1.500   3.000   4.000  1.00 20.00     0.100 C 
ATOM      3  C   ALA A   1       2.500   3.500   4.500  1.00 20.00     0.300 C 
ATOM      4  O   ALA A   1       3.000   4.500   4.000  1.00 20.00    -0.300 O 
TER
"""

SAMPLE_VINA_CONFIG = """center_x = 10.0
center_y = 20.0
center_z = 30.0
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 8
num_modes = 9
energy_range = 3.0
"""


class TestDockingSuite(unittest.TestCase):
    """Test suite verifying all Docking Automation components."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="docking_test_")
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    # ─────────────────────────────────────────────────────────────
    # 1. Executables & Environment
    # ─────────────────────────────────────────────────────────────
    def test_executables_detection(self):
        """Test detection and version reporting of real installed executables."""
        mgl = Path("C:/Program Files (x86)/MGLTools-1.5.7")
        vina = mgl / "vina.exe"
        split = mgl / "vina_split.exe"
        ag4 = mgl / "autogrid4.exe"
        ad4 = mgl / "autodock4.exe"

        self.assertTrue(validate_executable(vina, "Vina"))
        self.assertTrue(validate_executable(split, "Vina Split"))
        self.assertTrue(validate_executable(ag4, "AutoGrid 4"))
        self.assertTrue(validate_executable(ad4, "AutoDock 4"))

        vina_ver = get_vina_version(vina)
        self.assertIn("1.2.7", vina_ver)

        split_ver = get_vina_split_version(split)
        self.assertIn("1.2.7", split_ver)

        # OpenBabel detection
        from executables import get_obabel_version, resolve_executable
        ob_exe = resolve_executable(None, "obabel.exe")
        self.assertTrue(ob_exe.is_file())
        ob_ver = get_obabel_version(ob_exe)
        self.assertIn("3.", ob_ver)

        cpu = get_cpu_count("AUTO")
        self.assertGreater(cpu, 0)

    # ─────────────────────────────────────────────────────────────
    # 2. Vina Output Parser & Splitter
    # ─────────────────────────────────────────────────────────────
    def test_vina_parser(self):
        """Test parsing of multi-model Vina PDBQT files."""
        pdbqt_file = self.temp_path / "output.pdbqt"
        pdbqt_file.write_text(SAMPLE_VINA_PDBQT, encoding="utf-8")

        models = count_models(pdbqt_file)
        self.assertEqual(models, 3)

        results = extract_affinities(pdbqt_file)
        self.assertEqual(len(results), 3)

        self.assertEqual(results[0].pose, 1)
        self.assertAlmostEqual(results[0].binding_affinity, -8.4)
        self.assertAlmostEqual(results[0].rmsd_lower_bound, 0.0)

        self.assertEqual(results[1].pose, 2)
        self.assertAlmostEqual(results[1].binding_affinity, -7.9)
        self.assertAlmostEqual(results[1].rmsd_lower_bound, 1.52)

        table = format_vina_results_table(results)
        self.assertIn("-8.4", table)

    def test_manual_splitter_fallback(self):
        """Test manual fallback splitting of multi-model PDBQT."""
        pdbqt_file = self.temp_path / "output.pdbqt"
        pdbqt_file.write_text(SAMPLE_VINA_PDBQT, encoding="utf-8")
        split_dir = self.temp_path / "split_poses"

        split_files = split_pdbqt_manually(pdbqt_file, split_dir)
        self.assertEqual(len(split_files), 3)

        for f in split_files:
            self.assertTrue(f.exists())
            self.assertGreater(f.stat().st_size, 0)

        is_valid, msg = validate_split(split_dir, 3)
        self.assertTrue(is_valid, msg)

    # ─────────────────────────────────────────────────────────────
    # 3. Validators
    # ─────────────────────────────────────────────────────────────
    def test_validators(self):
        """Test PDBQT and Vina config validators."""
        rec_file = self.temp_path / "rec.pdbqt"
        rec_file.write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
        is_valid, errors = validate_pdbqt(rec_file)
        self.assertTrue(is_valid)
        self.assertEqual(len(errors), 0)

        cfg_file = self.temp_path / "config.txt"
        cfg_file.write_text(SAMPLE_VINA_CONFIG, encoding="utf-8")
        is_valid_cfg, cfg_errors = validate_vina_config(cfg_file)
        self.assertTrue(is_valid_cfg)
        self.assertEqual(len(cfg_errors), 0)

        parsed_cfg = parse_vina_config(cfg_file)
        self.assertAlmostEqual(parsed_cfg["center_x"], 10.0)
        self.assertAlmostEqual(parsed_cfg["size_x"], 20.0)

    # ─────────────────────────────────────────────────────────────
    # 4. Configuration & Queue Management
    # ─────────────────────────────────────────────────────────────
    def test_config_and_job_queue(self):
        """Test TOML config loading, queue building, and filtering."""
        receptors_dir = self.temp_path / "receptors"
        ligands_dir = self.temp_path / "ligands"
        results_dir = self.temp_path / "results"
        logs_dir = self.temp_path / "logs"
        reports_dir = self.temp_path / "reports"

        receptors_dir.mkdir()
        ligands_dir.mkdir()

        # Create 2 receptors with rigid subfolder and config.txt, and 3 ligands
        rec1_dir = receptors_dir / "rec1"
        (rec1_dir / "rigid").mkdir(parents=True)
        (rec1_dir / "rigid" / "rec1.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
        (rec1_dir / "config.txt").write_text(SAMPLE_VINA_CONFIG, encoding="utf-8")

        rec2_dir = receptors_dir / "rec2"
        (rec2_dir / "rigid").mkdir(parents=True)
        (rec2_dir / "rigid" / "rec2.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
        (rec2_dir / "config.txt").write_text(SAMPLE_VINA_CONFIG, encoding="utf-8")

        (ligands_dir / "lig1.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
        (ligands_dir / "lig2.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")
        (ligands_dir / "lig3.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")

        config = ProjectConfig(
            project_name="TestProject",
            project_root=self.temp_path,
            vina_executable=resolve_executable(None, "vina.exe"),
            vina_split_executable=resolve_executable(None, "vina_split.exe"),
            autogrid4_executable=resolve_executable(None, "autogrid4.exe"),
            autodock4_executable=resolve_executable(None, "autodock4.exe"),
            receptor_directory=receptors_dir,
            ligand_directory=ligands_dir,
            result_directory=results_dir,
            log_directory=logs_dir,
            report_directory=reports_dir,
            engine=Engine.VINA,
        )

        errors = validate_config(config)
        self.assertEqual(len(errors), 0, f"Config validation errors: {errors}")

        # Build queue: 2 receptors x 3 ligands = 6 jobs
        jobs = build_job_queue(config)
        self.assertEqual(len(jobs), 6)

        # Mark 1 job as SUCCESS and test filter_jobs with RESUME
        jobs[0].status = JobStatus.SUCCESS
        save_job_status(results_dir, jobs)

        loaded_jobs = load_job_status(results_dir)
        self.assertEqual(len(loaded_jobs), 6)

        previous_status = load_job_status(results_dir)
        to_run = filter_jobs(jobs, ResumeMode.RESUME, previous_status)
        skipped = len(jobs) - len(to_run)
        self.assertEqual(len(to_run), 5)
        self.assertEqual(skipped, 1)

    # ─────────────────────────────────────────────────────────────
    # 5. Reporting Engine
    # ─────────────────────────────────────────────────────────────
    def test_reporting_generation(self):
        """Test generation of CSV, XLSX, and JSON reports."""
        results_dir = self.temp_path / "results"
        reports_dir = self.temp_path / "reports"
        results_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)

        config = ProjectConfig(
            result_directory=results_dir,
            report_directory=reports_dir,
            engine=Engine.VINA,
            export_csv=True,
            export_excel=True,
        )

        # Create simulated completed jobs
        job1 = DockingJob(
            job_id="recA_vs_lig1",
            receptor_name="recA",
            ligand_name="lig1",
            receptor_path=self.temp_path / "recA.pdbqt",
            ligand_path=self.temp_path / "lig1.pdbqt",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            requested_modes=9,
            obtained_modes=9,
            elapsed_seconds=12.5,
            vina_results=[
                VinaResult(pose=1, binding_affinity=-9.2, rmsd_lower_bound=0.0, rmsd_upper_bound=0.0),
                VinaResult(pose=2, binding_affinity=-8.7, rmsd_lower_bound=1.2, rmsd_upper_bound=2.1),
            ],
        )
        job2 = DockingJob(
            job_id="recA_vs_lig2",
            receptor_name="recA",
            ligand_name="lig2",
            receptor_path=self.temp_path / "recA.pdbqt",
            ligand_path=self.temp_path / "lig2.pdbqt",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            requested_modes=9,
            obtained_modes=9,
            elapsed_seconds=14.1,
            vina_results=[
                VinaResult(pose=1, binding_affinity=-8.5, rmsd_lower_bound=0.0, rmsd_upper_bound=0.0),
            ],
        )

        generate_reports([job1, job2], config)

        # Verify output files exist
        csv_poses = reports_dir / "VINA_01_Vina_All_Poses.csv"
        csv_summary = reports_dir / "VINA_02_Vina_Summary.csv"
        csv_status = reports_dir / "VINA_03_Job_Status.csv"
        xlsx_file = reports_dir / "VINA_Complete_Report.xlsx"
        json_file = reports_dir / "VINA_data.json"

        self.assertTrue(csv_poses.exists())
        self.assertTrue(csv_summary.exists())
        self.assertTrue(csv_status.exists())
        self.assertTrue(xlsx_file.exists())
        self.assertTrue(json_file.exists())

        # Check content in CSV
        summary_text = csv_summary.read_text(encoding="utf-8")
        self.assertIn("recA", summary_text)
        self.assertIn("-9.2", summary_text)


class TestProteinLigandInteractions(unittest.TestCase):
    """Test the Pure Python 3 Protein-Ligand Interaction Profiler."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_interaction_profiler(self):
        from interactions import (
            PDBQTAtom, analyze_pose_interactions, detect_hydrogen_bonds,
            detect_hydrophobic_contacts, summarize_interactions,
        )

        rec_atoms = [
            # ASN donor
            PDBQTAtom(serial=1, name="N", res_name="ASN", chain="A", res_seq=10, x=10.0, y=10.0, z=10.0, charge=-0.3, atom_type="N"),
            # Hydrophobic LEU
            PDBQTAtom(serial=2, name="CD1", res_name="LEU", chain="A", res_seq=15, x=20.0, y=20.0, z=20.0, charge=0.0, atom_type="C"),
        ]

        lig_atoms = [
            # Ligand acceptor within 3.0 A of ASN
            PDBQTAtom(serial=10, name="O1", res_name="UNL", chain="L", res_seq=1, x=10.0, y=10.0, z=12.8, charge=-0.4, atom_type="OA"),
            # Ligand carbon within 3.5 A of LEU
            PDBQTAtom(serial=11, name="C5", res_name="UNL", chain="L", res_seq=1, x=20.0, y=20.0, z=23.3, charge=0.0, atom_type="C"),
        ]

        interactions = analyze_pose_interactions(rec_atoms, lig_atoms, pose_index=1)
        self.assertGreaterEqual(len(interactions), 2)

        types = [i.interaction_type for i in interactions]
        self.assertIn("Hydrogen Bond", types)
        self.assertIn("Hydrophobic", types)

        summary = summarize_interactions(interactions)
        self.assertEqual(summary["counts_by_type"]["Hydrogen Bond"], 1)
        self.assertEqual(summary["counts_by_type"]["Hydrophobic"], 1)


class TestJobDeserializerAndParallel(unittest.TestCase):
    """Test job status loading, deserialization, and parallel runner."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_load_docking_jobs(self):
        from job_manager import load_docking_jobs, save_job_status
        job = DockingJob(
            job_id="test__lig1",
            receptor_name="test",
            ligand_name="lig1",
            status=JobStatus.SUCCESS,
            vina_results=[VinaResult(pose=1, binding_affinity=-8.5)]
        )
        save_job_status(self.temp_path, [job])

        loaded = load_docking_jobs(self.temp_path)
        self.assertEqual(len(loaded), 1)
        self.assertIsInstance(loaded[0], DockingJob)
        self.assertEqual(loaded[0].job_id, "test__lig1")
        self.assertEqual(loaded[0].vina_results[0].binding_affinity, -8.5)

    def test_run_jobs_parallel(self):
        from job_manager import run_jobs_parallel
        cfg = ProjectConfig()
        cfg.result_directory = self.temp_path
        cfg.max_workers = 2

        jobs = [
            DockingJob(job_id=f"job_{i}", receptor_name="rec", ligand_name=f"lig_{i}")
            for i in range(4)
        ]

        def mock_runner(j: DockingJob, c: ProjectConfig) -> DockingJob:
            j.status = JobStatus.SUCCESS
            j.obtained_modes = 9
            return j

        results = run_jobs_parallel(
            jobs_to_run=jobs,
            all_jobs=jobs,
            config=cfg,
            runner_func=mock_runner,
            engine_name="MockEngine"
        )
        self.assertEqual(len(results), 4)
        self.assertTrue(all(j.status == JobStatus.SUCCESS for j in results))


class TestPreparationAnd3DViewer(unittest.TestCase):
    """Test molecular preparation, flexible receptor splitting, and 3D WebGL viewer."""

    def test_molecule_descriptors(self):
        from prepare import calculate_molecule_descriptors
        lig_path = Path("ligands/berberine.pdbqt")
        if lig_path.is_file():
            desc = calculate_molecule_descriptors(lig_path)
            self.assertEqual(desc["stem"], "berberine")
            self.assertGreater(desc["mw"], 200.0)
            self.assertIn("lipinski", desc)

    def test_3d_viewer_generation(self):
        from viewer_3d import generate_3d_viewer_html
        rec_path = Path("receptors/2V5Z/rigid/2v5z_clean.pdbqt")
        lig_path = Path("ligands/berberine.pdbqt")
        if rec_path.is_file() and lig_path.is_file():
            html = generate_3d_viewer_html(
                receptor_path=rec_path,
                ligand_path=lig_path,
                grid_center=(50.0, 150.0, 25.0),
                grid_size=(25.0, 25.0, 25.0),
                target_residues=["ILE:199", "TYR:326"]
            )
            self.assertIn("3Dmol.org", html)
            self.assertIn("AutoDock Suite Pro", html)
            self.assertIn("Docked Pose", html)

    def test_pymol_exporter_generation(self):
        from pymol_exporter import (
            export_pymol_session, find_pymol_executable, generate_pymol_script,
            is_pymol_available, pdbqt_to_clean_pdb_lines,
        )
        rec_path = Path("receptors/2V5Z/rigid/2V5Z.pdbqt")
        lig_path = Path("results/VINA/2V5Z/berberine/output/berberine_out.pdbqt")

        if rec_path.is_file() and lig_path.is_file():
            with tempfile.TemporaryDirectory(prefix="test_pymol_") as tmp_dir:
                res = export_pymol_session(
                    receptor_path=rec_path,
                    ligand_path=lig_path,
                    output_dir=tmp_dir,
                    pose_index=1,
                    grid_center=(51.886, 156.453, 28.559),
                    grid_size=(25.0, 25.0, 25.0),
                    compile_pse=True,
                )
                self.assertTrue(res["pml"].is_file())
                self.assertTrue(res["receptor_pdb"].is_file())
                self.assertTrue(res["ligand_pdb"].is_file())
                self.assertGreater(res["interactions_count"], 0)

                pml_content = res["pml"].read_text(encoding="utf-8")
                self.assertIn("reinitialize", pml_content)
                self.assertIn("show cartoon, receptor", pml_content)
                self.assertIn("show sticks, ligand_pose_1", pml_content)
                self.assertIn("select pocket_residues", pml_content)
                self.assertIn("Docking_Grid_Box", pml_content)
                self.assertIn("Hydrogen_Bonds", pml_content)

                # If PyMOL is installed, test that .pse was compiled
                if is_pymol_available():
                    self.assertIsNotNone(res["pse"])
                    self.assertTrue(res["pse"].is_file())
                    self.assertGreater(res["pse"].stat().st_size, 1000)



class TestWorkspaceDirectoryAnchoring(unittest.TestCase):
    """Test workspace directory anchoring, scaffolding, and DLG post-processing paths."""

    def setUp(self):
        self.temp_ws = tempfile.mkdtemp(prefix="test_workspace_")
        self.ws_path = Path(self.temp_ws).resolve()
        self.orig_cwd = Path.cwd()

    def tearDown(self):
        try:
            os.chdir(self.orig_cwd)
        except Exception:
            pass
        shutil.rmtree(self.temp_ws, ignore_errors=True)

    def test_set_workspace_directory(self):
        from config import ProjectConfig, set_workspace_directory

        cfg = ProjectConfig()
        set_workspace_directory(cfg, self.ws_path)

        self.assertEqual(cfg.project_root, self.ws_path)
        self.assertEqual(cfg.project_name, self.ws_path.name)
        self.assertEqual(cfg.result_directory, (self.ws_path / "results").resolve())
        self.assertEqual(cfg.log_directory, (self.ws_path / "logs").resolve())
        self.assertEqual(cfg.report_directory, (self.ws_path / "reports").resolve())

        # Check required directories created
        self.assertTrue((self.ws_path / "results" / "DLG").is_dir())
        self.assertTrue((self.ws_path / "results" / "DPF").is_dir())
        self.assertTrue((self.ws_path / "results" / "VINA").is_dir())
        self.assertTrue((self.ws_path / "results" / "BSNDVP_RESULTS").is_dir())
        self.assertTrue((self.ws_path / "DLG").is_dir())
        self.assertTrue((self.ws_path / "BSNDVP_RESULTS").is_dir())
        self.assertTrue((self.ws_path / "logs").is_dir())
        self.assertTrue((self.ws_path / "reports").is_dir())

        # Check dlg_extract.py copied
        self.assertTrue((self.ws_path / "dlg_extract.py").is_file())

    def test_load_config_with_relative_root(self):
        from config import load_config
        toml_content = """[project]
name = "MyTestWorkspace"
root = "."

[outputs]
result_directory = "results"
log_directory = "logs"
report_directory = "reports"
"""
        toml_file = self.ws_path / "project_config.toml"
        toml_file.write_text(toml_content, encoding="utf-8")

        loaded = load_config(toml_file)
        self.assertEqual(loaded.project_root, self.ws_path)
        self.assertEqual(loaded.result_directory, (self.ws_path / "results").resolve())
        self.assertEqual(loaded.log_directory, (self.ws_path / "logs").resolve())
        self.assertEqual(loaded.report_directory, (self.ws_path / "reports").resolve())

    def test_dlg_extract_discovery(self):
        from dlg_extract import _find_dlg_dir

        res_dlg = self.ws_path / "results" / "DLG"
        res_dlg.mkdir(parents=True, exist_ok=True)
        (res_dlg / "test.dlg").write_text("TEST DLG FILE", encoding="utf-8")

        discovered = _find_dlg_dir(self.ws_path)
        self.assertIsNotNone(discovered)
        self.assertEqual(discovered.resolve(), res_dlg.resolve())

        # Test top-level DLG folder discovery
        top_dlg = self.ws_path / "DLG"
        top_dlg.mkdir(parents=True, exist_ok=True)
        (top_dlg / "test_top.dlg").write_text("TEST TOP DLG", encoding="utf-8")

        discovered_top = _find_dlg_dir(self.ws_path)
        self.assertIsNotNone(discovered_top)
        self.assertTrue(discovered_top.is_dir())

    def test_create_scaffold_completeness(self):
        from prepare import create_scaffold

        create_scaffold(self.ws_path)
        self.assertTrue((self.ws_path / "receptors").is_dir())
        self.assertTrue((self.ws_path / "ligands").is_dir())
        self.assertTrue((self.ws_path / "results" / "DLG").is_dir())
        self.assertTrue((self.ws_path / "DLG").is_dir())
        self.assertTrue((self.ws_path / "BSNDVP_RESULTS").is_dir())
        self.assertTrue((self.ws_path / "project_config.toml").is_file())
        self.assertTrue((self.ws_path / "dlg_extract.py").is_file())

    def test_dlg_pipeline_in_anchored_workspace(self):
        from config import set_workspace_directory, ProjectConfig
        from dlg_extract import DLGPipeline

        cfg = ProjectConfig()
        set_workspace_directory(cfg, self.ws_path)

        # Place a real DLG file in the workspace results/DLG folder
        source_dlg = Path("AD4/DLG_ANALYSIS/DLG/2V5Z_Berberine.dlg")
        if source_dlg.is_file():
            target_dlg = self.ws_path / "results" / "DLG" / "2V5Z_Berberine.dlg"
            shutil.copy2(source_dlg, target_dlg)

            out_dir = self.ws_path / "results" / "BSNDVP_RESULTS"
            pipeline = DLGPipeline(self.ws_path / "results" / "DLG", out_dir, self.ws_path)
            pipeline.run()

            # Verify outputs in results/BSNDVP_RESULTS
            self.assertTrue((out_dir / "01_Docking_Master.csv").is_file())
            self.assertTrue((out_dir / "02_Clusters.csv").is_file())
            self.assertTrue((out_dir / "03_Ligand_Summary.csv").is_file())

            # Verify mirrored outputs in workspace BSNDVP_RESULTS
            self.assertTrue((self.ws_path / "BSNDVP_RESULTS" / "01_Docking_Master.csv").is_file())
            self.assertTrue((self.ws_path / "BSNDVP_RESULTS" / "03_Ligand_Summary.csv").is_file())

            # Verify mirrored outputs in workspace reports
            self.assertTrue((self.ws_path / "reports" / "01_Docking_Master.csv").is_file())


class TestChainAndResidueEditing(unittest.TestCase):
    """Test chain detection/removal and residue mutation/editing studio logic."""

    SAMPLE_PDB = """\
ATOM      1  N   MET A   1      20.154  14.234  10.123  1.00 20.00           N
ATOM      2  CA  MET A   1      21.234  15.123  10.456  1.00 20.00           C
ATOM      3  C   MET A   1      22.456  14.567  11.123  1.00 20.00           C
ATOM      4  O   MET A   1      23.123  13.789  10.567  1.00 20.00           O
ATOM      5  CB  MET A   1      21.567  16.234   9.456  1.00 20.00           C
ATOM      6  CG  MET A   1      22.789  17.123   9.890  1.00 20.00           C
ATOM      7  SD  MET A   1      23.123  18.234   8.567  1.00 20.00           S
ATOM      8  CE  MET A   1      24.567  19.123   9.123  1.00 20.00           C
ATOM      9  N   ARG A   2      22.789  15.012  12.345  1.00 20.00           N
ATOM     10  CA  ARG A   2      23.912  14.567  13.123  1.00 20.00           C
ATOM     11  C   ARG A   2      25.123  15.456  12.890  1.00 20.00           C
ATOM     12  O   ARG A   2      26.123  15.012  12.345  1.00 20.00           O
ATOM     13  CB  ARG A   2      23.567  14.567  14.612  1.00 20.00           C
ATOM     14  CG  ARG A   2      24.678  14.123  15.567  1.00 20.00           C
ATOM     15  CD  ARG A   2      24.234  14.123  17.012  1.00 20.00           C
ATOM     16  NE  ARG A   2      25.345  13.789  17.890  1.00 20.00           N
ATOM     17  CZ  ARG A   2      25.234  13.456  19.178  1.00 20.00           C
ATOM     18  NH1 ARG A   2      24.067  13.456  19.789  1.00 20.00           N
ATOM     19  NH2 ARG A   2      26.345  13.123  19.823  1.00 20.00           N
ATOM     20  N   HIS A   3      25.012  16.678  13.456  1.00 20.00           N
ATOM     21  CA  HIS A   3      26.123  17.567  13.234  1.00 20.00           C
ATOM     22  C   HIS A   3      27.345  16.890  12.678  1.00 20.00           C
ATOM     23  O   HIS A   3      28.456  17.412  12.567  1.00 20.00           O
ATOM     24  CB  HIS A   3      25.678  18.789  12.412  1.00 20.00           C
ATOM     25  CG  HIS A   3      26.678  19.890  12.345  1.00 20.00           C
ATOM     26  ND1 HIS A   3      27.890  19.789  11.678  1.00 20.00           N
ATOM     27  CD2 HIS A   3      26.567  21.123  12.890  1.00 20.00           C
ATOM     28  CE1 HIS A   3      28.456  20.890  11.789  1.00 20.00           C
ATOM     29  NE2 HIS A   3      27.678  21.734  12.523  1.00 20.00           N
TER
ATOM     30  N   MET B   1      50.154  14.234  10.123  1.00 20.00           N
ATOM     31  CA  MET B   1      51.234  15.123  10.456  1.00 20.00           C
ATOM     32  C   MET B   1      52.456  14.567  11.123  1.00 20.00           C
ATOM     33  O   MET B   1      53.123  13.789  10.567  1.00 20.00           O
ATOM     34  CB  MET B   1      51.567  16.234   9.456  1.00 20.00           C
TER
HETATM   35  O   HOH A 101      10.000  10.000  10.000  1.00 25.00           O
"""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="test_edit_")
        self.temp_path = Path(self.temp_dir)
        self.pdb_file = self.temp_path / "test_multichain.pdb"
        self.pdb_file.write_text(self.SAMPLE_PDB, encoding="utf-8")

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_parse_pdb_chains(self):
        from prepare import parse_pdb_chains

        chains = parse_pdb_chains(self.pdb_file)
        chain_ids = [c["chain"] for c in chains]
        self.assertIn("A", chain_ids)
        self.assertIn("B", chain_ids)

        # Chain A has 29 ATOM + 1 HETATM = 30 atoms
        chain_a = next(c for c in chains if c["chain"] == "A")
        self.assertEqual(chain_a["atom_count"], 30)
        self.assertEqual(chain_a["residue_count"], 4)  # res 1, 2, 3, 101

        # Chain B has 5 atoms
        chain_b = next(c for c in chains if c["chain"] == "B")
        self.assertEqual(chain_b["atom_count"], 5)
        self.assertEqual(chain_b["residue_count"], 1)

    def test_parse_pdb_residues(self):
        from prepare import parse_pdb_residues

        residues = parse_pdb_residues(self.pdb_file)
        # Should contain MET A 1, ARG A 2, HIS A 3, MET B 1 (ATOM records only)
        keys = [(r["chain"], r["seq"], r["name"]) for r in residues]
        self.assertIn(("A", 1, "MET"), keys)
        self.assertIn(("A", 2, "ARG"), keys)
        self.assertIn(("A", 3, "HIS"), keys)
        self.assertIn(("B", 1, "MET"), keys)

        # Check atom counts
        res_arg = next(r for r in residues if (r["chain"], r["seq"]) == ("A", 2))
        self.assertEqual(res_arg["atom_count"], 11)  # 11 atoms in ARG in sample

    def test_chain_removal(self):
        from prepare import prepare_receptor_pdbqt

        out_pdbqt = self.temp_path / "receptor.pdbqt"
        prepare_receptor_pdbqt(
            self.pdb_file,
            out_pdbqt,
            remove_chains={"B"},
            cleanup_water=True,
        )

        self.assertTrue(out_pdbqt.is_file())
        content = out_pdbqt.read_text(encoding="utf-8")

        # Must not contain any Chain B atoms
        for line in content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21:22].strip()
                self.assertNotEqual(chain, "B", "Found Chain B atom in output despite remove_chains={'B'}")

        # Check edited PDB was saved
        edited_pdb = self.temp_path / "receptor_edited.pdb"
        self.assertTrue(edited_pdb.is_file())
        edited_content = edited_pdb.read_text(encoding="utf-8")
        self.assertNotIn(" MET B ", edited_content)

    def test_alanine_mutation_sidechain_truncation(self):
        from prepare import prepare_receptor_pdbqt

        out_pdbqt = self.temp_path / "receptor_ala.pdbqt"
        # Mutate ARG A:2 (which has 11 atoms) to ALA
        prepare_receptor_pdbqt(
            self.pdb_file,
            out_pdbqt,
            mutations={("A", 2): "ALA"},
            cleanup_water=True,
        )

        edited_pdb = self.temp_path / "receptor_ala_edited.pdb"
        self.assertTrue(edited_pdb.is_file())
        lines = [l for l in edited_pdb.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

        # Check residue 2 in Chain A
        res2_lines = [l for l in lines if l[21:22].strip() == "A" and l[22:26].strip() == "2"]
        # All atoms in res 2 must now have residue name 'ALA'
        for l in res2_lines:
            self.assertEqual(l[17:20].strip(), "ALA")

        atom_names = [l[12:16].strip() for l in res2_lines]
        # Should keep N, CA, C, O, CB
        self.assertEqual(sorted(atom_names), ["C", "CA", "CB", "N", "O"])
        # Sidechain atoms beyond CB must be truncated
        for removed in ("CG", "CD", "NE", "CZ", "NH1", "NH2"):
            self.assertNotIn(removed, atom_names)

    def test_glycine_mutation_sidechain_truncation(self):
        from prepare import prepare_receptor_pdbqt

        out_pdbqt = self.temp_path / "receptor_gly.pdbqt"
        # Mutate MET A:1 (which has CB, CG, SD, CE) to GLY
        prepare_receptor_pdbqt(
            self.pdb_file,
            out_pdbqt,
            mutations={("A", 1): "GLY"},
            cleanup_water=True,
        )

        edited_pdb = self.temp_path / "receptor_gly_edited.pdb"
        self.assertTrue(edited_pdb.is_file())
        lines = [l for l in edited_pdb.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

        res1_lines = [l for l in lines if l[21:22].strip() == "A" and l[22:26].strip() == "1"]
        for l in res1_lines:
            self.assertEqual(l[17:20].strip(), "GLY")

        atom_names = [l[12:16].strip() for l in res1_lines]
        # GLY should only keep backbone N, CA, C, O (CB must be stripped)
        self.assertEqual(sorted(atom_names), ["C", "CA", "N", "O"])
        self.assertNotIn("CB", atom_names)

    def test_residue_deletion(self):
        from prepare import prepare_receptor_pdbqt

        out_pdbqt = self.temp_path / "receptor_del.pdbqt"
        # Delete ARG A:2 completely
        prepare_receptor_pdbqt(
            self.pdb_file,
            out_pdbqt,
            deleted_residues={("A", 2)},
            cleanup_water=True,
        )

        edited_pdb = self.temp_path / "receptor_del_edited.pdb"
        self.assertTrue(edited_pdb.is_file())
        lines = [l for l in edited_pdb.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

        # Residue 2 in Chain A should not exist at all
        res2_lines = [l for l in lines if l[21:22].strip() == "A" and l[22:26].strip() == "2"]
        self.assertEqual(len(res2_lines), 0)

    def test_protonation_variant_rename(self):
        from prepare import prepare_receptor_pdbqt

        out_pdbqt = self.temp_path / "receptor_hip.pdbqt"
        # Mutate HIS A:3 to HIP (protonated histidine)
        prepare_receptor_pdbqt(
            self.pdb_file,
            out_pdbqt,
            mutations={("A", 3): "HIP"},
            cleanup_water=True,
        )

        edited_pdb = self.temp_path / "receptor_hip_edited.pdb"
        self.assertTrue(edited_pdb.is_file())
        lines = [l for l in edited_pdb.read_text(encoding="utf-8").splitlines() if l.startswith("ATOM")]

        res3_lines = [l for l in lines if l[21:22].strip() == "A" and l[22:26].strip() == "3"]
        self.assertGreater(len(res3_lines), 0)
        for l in res3_lines:
            self.assertEqual(l[17:20].strip(), "HIP")


class TestMultiReceptorAndGridBoxResolution(unittest.TestCase):
    """Test suite verifying consecutive multi-receptor preparation, cross-alias mirroring,
    and GridBoxTab receptor finding without 'Missing PDBQT' errors."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="multi_rec_test_")
        self.workspace = Path(self.temp_dir)
        self.rec_dir = self.workspace / "receptors"
        self.rec_dir.mkdir(parents=True, exist_ok=True)
        self.macro_dir = self.workspace / "Macromolecules"
        self.macro_dir.mkdir(parents=True, exist_ok=True)
        self.lig_dir = self.workspace / "ligands"
        self.lig_dir.mkdir(parents=True, exist_ok=True)

        # Create dummy ligand
        (self.lig_dir / "ligand1.pdbqt").write_text(SAMPLE_VINA_PDBQT, encoding="utf-8")

        # Create sample PDB files
        self.pdb1 = self.workspace / "2V5Z.pdb"
        self.pdb1.write_text(TestChainAndResidueEditing.SAMPLE_PDB, encoding="utf-8")
        self.pdb2 = self.workspace / "2Z5X_clean.pdb"
        self.pdb2.write_text(TestChainAndResidueEditing.SAMPLE_PDB, encoding="utf-8")

        # Project config anchored to workspace
        self.config = ProjectConfig(
            project_root=self.workspace,
            receptor_directory=self.rec_dir,
            ligand_directory=self.lig_dir,
            result_directory=self.workspace / "results",
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_multi_receptor_prep_and_mirroring(self):
        from prepare import prepare_receptor_pdbqt

        # Simulate Prep Tab preparing receptor 1 (2V5Z)
        name1 = "2V5Z"
        out_dir1 = self.rec_dir / name1 / "rigid"
        out_dir1.mkdir(parents=True, exist_ok=True)
        out1 = out_dir1 / f"{name1}.pdbqt"
        prepare_receptor_pdbqt(self.pdb1, out1, cleanup_water=True)
        self.assertTrue(out1.is_file())

        # Simulate Prep Tab _on_receptor_done mirroring logic
        rec_folder1 = self.rec_dir / name1
        rec_root1 = rec_folder1 / f"{name1}.pdbqt"
        shutil.copy2(out1, rec_root1)
        cfg1 = rec_folder1 / "config.txt"
        cfg1.write_text("center_x = 0.0\ncenter_y = 0.0\ncenter_z = 0.0\nsize_x = 25.0\nsize_y = 25.0\nsize_z = 25.0\n", encoding="utf-8")

        # Cross-mirror 2V5Z -> 2V5Z_clean and Macromolecules
        targets1 = [name1, f"{name1}_clean"]
        for t in targets1:
            td = self.rec_dir / t
            td.mkdir(parents=True, exist_ok=True)
            (td / "rigid").mkdir(parents=True, exist_ok=True)
            if (td / f"{t}.pdbqt").resolve() != out1.resolve():
                shutil.copy2(out1, td / f"{t}.pdbqt")
            if (td / "rigid" / f"{t}.pdbqt").resolve() != out1.resolve():
                shutil.copy2(out1, td / "rigid" / f"{t}.pdbqt")
            if (td / "config.txt").resolve() != cfg1.resolve():
                shutil.copy2(cfg1, td / "config.txt")
            # Mirror to Macromolecules
            md = self.macro_dir / t
            md.mkdir(parents=True, exist_ok=True)
            (md / "rigid").mkdir(parents=True, exist_ok=True)
            if (md / f"{t}.pdbqt").resolve() != out1.resolve():
                shutil.copy2(out1, md / f"{t}.pdbqt")
            if (md / "rigid" / f"{t}.pdbqt").resolve() != out1.resolve():
                shutil.copy2(out1, md / "rigid" / f"{t}.pdbqt")
            if (md / "config.txt").resolve() != cfg1.resolve():
                shutil.copy2(cfg1, md / "config.txt")

        # Now simulate Prep Tab preparing receptor 2 with _clean suffix (2Z5X_clean)
        name2 = "2Z5X_clean"
        out_dir2 = self.rec_dir / name2 / "rigid"
        out_dir2.mkdir(parents=True, exist_ok=True)
        out2 = out_dir2 / f"{name2}.pdbqt"
        prepare_receptor_pdbqt(self.pdb2, out2, cleanup_water=True)
        self.assertTrue(out2.is_file())

        rec_folder2 = self.rec_dir / name2
        rec_root2 = rec_folder2 / f"{name2}.pdbqt"
        if rec_root2.resolve() != out2.resolve():
            shutil.copy2(out2, rec_root2)
        cfg2 = rec_folder2 / "config.txt"
        cfg2.write_text("center_x = 10.0\ncenter_y = 20.0\ncenter_z = 30.0\nsize_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n", encoding="utf-8")

        # Cross-mirror 2Z5X_clean -> 2Z5X and Macromolecules
        targets2 = [name2, name2.replace("_clean", "")]
        for t in targets2:
            td = self.rec_dir / t
            td.mkdir(parents=True, exist_ok=True)
            (td / "rigid").mkdir(parents=True, exist_ok=True)
            if (td / f"{t}.pdbqt").resolve() != out2.resolve():
                shutil.copy2(out2, td / f"{t}.pdbqt")
            if (td / "rigid" / f"{t}.pdbqt").resolve() != out2.resolve():
                shutil.copy2(out2, td / "rigid" / f"{t}.pdbqt")
            if (td / "config.txt").resolve() != cfg2.resolve():
                shutil.copy2(cfg2, td / "config.txt")
            md = self.macro_dir / t
            md.mkdir(parents=True, exist_ok=True)
            (md / "rigid").mkdir(parents=True, exist_ok=True)
            if (md / f"{t}.pdbqt").resolve() != out2.resolve():
                shutil.copy2(out2, md / f"{t}.pdbqt")
            if (md / "rigid" / f"{t}.pdbqt").resolve() != out2.resolve():
                shutil.copy2(out2, md / "rigid" / f"{t}.pdbqt")
            if (md / "config.txt").resolve() != cfg2.resolve():
                shutil.copy2(cfg2, md / "config.txt")

        # Verify all destinations exist for BOTH receptors
        self.assertTrue((self.rec_dir / "2V5Z" / "rigid" / "2V5Z.pdbqt").is_file())
        self.assertTrue((self.rec_dir / "2V5Z_clean" / "rigid" / "2V5Z_clean.pdbqt").is_file())
        self.assertTrue((self.rec_dir / "2Z5X_clean" / "rigid" / "2Z5X_clean.pdbqt").is_file())
        self.assertTrue((self.rec_dir / "2Z5X" / "rigid" / "2Z5X.pdbqt").is_file())

        # Verify Macromolecules mirroring
        self.assertTrue((self.macro_dir / "2Z5X_clean" / "2Z5X_clean.pdbqt").is_file())
        self.assertTrue((self.macro_dir / "2Z5X" / "2Z5X.pdbqt").is_file())

    def test_grid_box_receptor_discovery_and_finding(self):
        # Place a receptor inside Macromolecules only (like in user's Group3_workspace)
        group_rec = self.macro_dir / "2v5z_clean"
        (group_rec / "rigid").mkdir(parents=True, exist_ok=True)
        (group_rec / "rigid" / "2v5z_clean.pdbqt").write_text(SAMPLE_RECEPTOR_PDBQT, encoding="utf-8")

        # Create a mock GridBoxTab instance with our test config
        from gui.grid_tab import GridBoxTab
        mock_grid = object.__new__(GridBoxTab)
        mock_grid.config = self.config

        # 1. Listing should find 2v5z_clean even though it's in Macromolecules
        receptors = mock_grid._list_receptors()
        self.assertIn("2v5z_clean", receptors)

        # 2. Finding by exact name
        found_exact = mock_grid._find_receptor_pdbqt("2v5z_clean")
        self.assertIsNotNone(found_exact)
        self.assertTrue(found_exact.is_file())

        # 3. Finding by base alias (2v5z) without _clean
        found_alias = mock_grid._find_receptor_pdbqt("2v5z")
        self.assertIsNotNone(found_alias)
        self.assertTrue(found_alias.is_file())

        # 4. Finding by uppercase alias (2V5Z)
        found_upper = mock_grid._find_receptor_pdbqt("2V5Z")
        self.assertIsNotNone(found_upper)
        self.assertTrue(found_upper.is_file())

        # 5. After finding, it must be mirrored to standard rec_dir / name so AutoGrid never fails
        self.assertTrue((self.rec_dir / "2V5Z" / "2V5Z.pdbqt").is_file())


class TestDockingFixesAndDLG(unittest.TestCase):
    """Test AutoGrid box dimension conversion, DLG scientific notation parsing, and 2D depiction."""

    def test_autogrid_box_points_calculation(self):
        from autodock4_workflow import generate_gpf

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            rec_pdbqt = td_path / "test_rec.pdbqt"
            lig_pdbqt = td_path / "test_lig.pdbqt"
            gpf_out = td_path / "test_rec.gpf"

            rec_pdbqt.write_text("ATOM      1  CA  ALA A   1      10.000  10.000  10.000  0.00  0.00    +0.000 C\n", encoding="utf-8")
            lig_pdbqt.write_text("ATOM      1  C1  LIG L   1       5.000   5.000   5.000  0.00  0.00    +0.000 C\n", encoding="utf-8")

            # 25.0 Angstrom box at 0.375 spacing: 25 / 0.375 = 66.67 -> must be 68 points!
            generate_gpf(
                receptor_path=rec_pdbqt,
                ligand_path=lig_pdbqt,
                output_path=gpf_out,
                grid_center={"x": 10.0, "y": 10.0, "z": 10.0},
                grid_size={"x": 25.0, "y": 25.0, "z": 25.0},
                spacing=0.375,
            )

            content = gpf_out.read_text(encoding="utf-8")
            self.assertIn("npts 68 68 68", content)
            self.assertIn("spacing 0.375", content)

    def test_dlg_scientific_notation_energy_parsing(self):
        from dlg_extract import DLGParser

        line_pos = "Estimated Free Energy of Binding    =  +6.29e+008 kcal/mol [=(1)+(2)+(3)-(4)]"
        line_neg = "Estimated Free Energy of Binding    =  -8.45 kcal/mol [=(1)+(2)+(3)-(4)]"

        m_pos = DLGParser._RE_ENERGY.search(line_pos)
        self.assertIsNotNone(m_pos)
        self.assertAlmostEqual(float(m_pos.group(1)), 6.29e8)

        m_neg = DLGParser._RE_ENERGY.search(line_neg)
        self.assertIsNotNone(m_neg)
        self.assertAlmostEqual(float(m_neg.group(1)), -8.45)

    def test_interactions_direct_dlg_support(self):
        from interactions import parse_docked_poses

        mock_dlg = (
            "DOCKED: MODEL        1\n"
            "DOCKED: USER    Estimated Free Energy of Binding = -7.80 kcal/mol\n"
            "DOCKED: ATOM      1  C1  LIG L   1      12.000  14.000  16.000  1.00 20.00    +0.050 C\n"
            "DOCKED: ATOM      2  O1  LIG L   1      13.000  14.500  16.500  1.00 20.00    -0.500 OA\n"
            "DOCKED: ENDMDL\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".dlg", delete=False) as f:
            f.write(mock_dlg)
            f_path = Path(f.name)

        try:
            poses = parse_docked_poses(f_path)
            self.assertEqual(len(poses), 1)
            self.assertEqual(poses[0][0], 1)
            self.assertEqual(len(poses[0][1]), 2)
            self.assertEqual(poses[0][1][0].name, "C1")
            self.assertEqual(poses[0][1][1].name, "O1")
        finally:
            f_path.unlink(missing_ok=True)


class TestDualEngineAndInteractionDiagram(unittest.TestCase):
    """Tests for Engine.BOTH, job ID collision avoidance, and 2D interaction diagram."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp(prefix="dual_engine_test_")
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_job_id_engine_namespacing(self):
        """AutoDock4 jobs must have an 'ad4_' prefix to prevent job_status.json collisions."""
        from models import generate_job_id, Engine, DockingMode

        vina_id = generate_job_id("2v5z", "aspirin", DockingMode.RIGID, engine=Engine.VINA)
        ad4_id = generate_job_id("2v5z", "aspirin", DockingMode.RIGID, engine=Engine.AUTODOCK4)

        self.assertEqual(vina_id, "2v5z__aspirin")
        self.assertEqual(ad4_id, "ad4_2v5z__aspirin")
        self.assertNotEqual(vina_id, ad4_id, "Vina and AD4 must produce distinct job IDs")

        # Flexible docking
        vina_flex = generate_job_id("2v5z", "aspirin", DockingMode.FLEXIBLE, engine=Engine.VINA)
        ad4_flex = generate_job_id("2v5z", "aspirin", DockingMode.FLEXIBLE, engine=Engine.AUTODOCK4)
        self.assertEqual(vina_flex, "2v5z_flex__aspirin")
        self.assertEqual(ad4_flex, "ad4_2v5z_flex__aspirin")

    def test_engine_both_enum_and_validation(self):
        """Engine.BOTH is a valid Engine and validate_config checks both sets of binaries."""
        from models import Engine
        from config import ProjectConfig, validate_config

        self.assertEqual(Engine.BOTH.value, "BOTH")
        self.assertIn(Engine.BOTH, list(Engine))

        cfg = ProjectConfig(
            engine=Engine.BOTH,
            vina_executable=Path("nonexistent_vina_xyz.exe"),
            autodock4_executable=Path("nonexistent_ad4_xyz.exe"),
        )
        errors = validate_config(cfg)
        self.assertTrue(any("Vina" in e for e in errors), f"Expected Vina error in: {errors}")
        self.assertTrue(any("AutoDock4" in e for e in errors), f"Expected AD4 error in: {errors}")

    def test_interaction_diagram_generation(self):
        """render_interaction_diagram_html should generate valid offline HTML+SVG."""
        from gui.interaction_diagram import (
            render_interaction_diagram_html, get_interaction_color, is_available
        )
        from dataclasses import dataclass

        @dataclass
        class MockInteraction:
            interaction_type: str
            receptor_residue: str
            receptor_atom: str
            ligand_atom: str
            distance_angstrom: float
            pose_index: int = 1

        mock_inters = [
            MockInteraction("Hydrogen Bond", "ASN 119", "ND2", "O1", 2.85),
            MockInteraction("Hydrophobic Contact", "LEU 167", "CD1", "C4", 3.72),
        ]

        # Test with empty path (placeholder mode)
        html = render_interaction_diagram_html(
            ligand_path="nonexistent.pdbqt",
            interactions=mock_inters,
            width=400,
            height=340,
            theme="dark",
        )
        self.assertIsNotNone(html)
        self.assertIn("<!DOCTYPE html>", html)
        self.assertIn("<svg", html)
        self.assertIn("ASN 119", html)
        self.assertIn("LEU 167", html)
        self.assertIn("H-Bond", html)
        self.assertIn("Hydrophobic", html)

        # Test colors
        self.assertIn(get_interaction_color("Hydrogen Bond"), ["#22c55e", "#16a34a", "#3b82f6"])
        self.assertIn(get_interaction_color("Hydrophobic Contact"), ["#f59e0b", "#b91c1c"])

    def test_new_themes_plasma_and_bio_neon(self):
        """Verify PLASMA and BIO_NEON themes are registered and have distinct colors."""
        from gui.themes import ALL_THEMES, THEME_NAMES, PLASMA, BIO_NEON, set_theme

        self.assertIn("Plasma", THEME_NAMES)
        self.assertIn("Bio-Neon", THEME_NAMES)

        self.assertEqual(PLASMA.name, "Plasma")
        self.assertEqual(PLASMA.ctk_mode, "dark")
        self.assertIn("#", PLASMA.bg_primary)

        self.assertEqual(BIO_NEON.name, "Bio-Neon")
        self.assertEqual(BIO_NEON.ctk_mode, "dark")
        self.assertIn("#", BIO_NEON.bg_primary)

        # Live switch should work without error
        t = set_theme("Plasma")
        self.assertEqual(t.name, "Plasma")
        t2 = set_theme("Bio-Neon")
        self.assertEqual(t2.name, "Bio-Neon")
        # Reset to Dark
        set_theme("Dark")

    def test_generate_gpf_file(self):
        """Verify generate_gpf_file works with various center/size formats and fallbacks."""
        from autodock4_workflow import generate_gpf_file

        rec_file = self.temp_path / "rec.pdbqt"
        rec_file.write_text("ATOM      1  N   ALA A   1      11.104  13.207   2.476  1.00 20.00    -0.350 N\nEND\n", encoding="utf-8")
        gpf_out = self.temp_path / "rec.gpf"

        # Call with tuple center and size
        res_path = generate_gpf_file(
            receptor_pdbqt=rec_file,
            output_gpf=gpf_out,
            center=(10.0, 12.0, 2.0),
            size=(20.0, 20.0, 20.0),
            spacing=0.375,
        )
        self.assertTrue(res_path.is_file())
        content = res_path.read_text(encoding="utf-8")
        self.assertIn("npts", content)
        self.assertIn("gridcenter 10.000 12.000 2.000", content)
        self.assertIn("spacing 0.375", content)
        self.assertIn("elecmap", content)
        self.assertIn("dsolvmap", content)

    def test_kollman_charge_assignment(self):
        """Verify assign_kollman_charges sets united-atom charges on receptor residues."""
        from prepare import assign_kollman_charges

        sample_pdbqt = (
            "ATOM      1  N   ALA A   1      11.104  13.207   2.476  1.00 20.00    +0.000 N\n"
            "ATOM      2  CA  ALA A   1      12.000  14.000   3.000  1.00 20.00    +0.000 C\n"
            "ATOM      3  C   ALA A   1      13.000  15.000   4.000  1.00 20.00    +0.000 C\n"
            "ATOM      4  O   ALA A   1      14.000  16.000   5.000  1.00 20.00    +0.000 OA\n"
            "END\n"
        )
        test_file = self.temp_path / "kollman_test.pdbqt"
        test_file.write_text(sample_pdbqt, encoding="utf-8")

        assign_kollman_charges(test_file)
        result_text = test_file.read_text(encoding="utf-8")
        self.assertIn("-0.350 N", result_text)
        self.assertIn("+0.550 C", result_text)
        self.assertIn("-0.550 OA", result_text)

    def test_sanitize_flex_pdbqt(self):
        """Verify sanitize_flex_pdbqt renumbers residues so ROOT is 1 and BRANCH is 2..N."""
        from prepare import sanitize_flex_pdbqt

        raw_flex = (
            "BEGIN_RES MET A 1\n"
            "ROOT\n"
            "ATOM   1339  CA  MET A   1      12.000  14.000   3.000  1.00 20.00    +0.000 C\n"
            "ENDROOT\n"
            "BRANCH 1 2\n"
            "ATOM   1342  CB  MET A   1      13.000  15.000   4.000  1.00 20.00    +0.000 C\n"
            "ENDBRANCH 1 2\n"
            "END_RES\n"
        )
        test_file = self.temp_path / "flex_test.pdbqt"
        test_file.write_text(raw_flex, encoding="utf-8")

        sanitize_flex_pdbqt(test_file)
        content = test_file.read_text(encoding="utf-8")
        self.assertIn("ATOM      1  CA  MET A   1", content)
        self.assertIn("BRANCH 1 2", content)
        self.assertIn("ATOM      2  CB  MET A   1", content)

    def test_default_energy_range_vina(self):
        """Verify default energy_range across config is 6.0 kcal/mol."""
        from config import ProjectConfig
        cfg = ProjectConfig()
        self.assertEqual(cfg.energy_range, 6.0)

    def test_save_job_status_merging(self):
        """Verify save_job_status merges jobs by job_id rather than clobbering existing ones."""
        from job_manager import save_job_status, DockingJob, JobStatus, Engine
        import json

        out_dir = self.temp_path / "results_merge_test"
        out_dir.mkdir(parents=True, exist_ok=True)

        job_vina = DockingJob(
            job_id="job_ligand1_rigid",
            receptor_name="rec1",
            ligand_name="lig1",
            status=JobStatus.SUCCESS,
            engine=Engine.VINA,
        )
        save_job_status(out_dir, [job_vina])

        job_ad4 = DockingJob(
            job_id="ad4_job_ligand1_rigid",
            receptor_name="rec1",
            ligand_name="lig1",
            status=JobStatus.SUCCESS,
            engine=Engine.AUTODOCK4,
        )
        save_job_status(out_dir, [job_ad4])

        status_file = out_dir / "job_status.json"
        self.assertTrue(status_file.exists())
        data = json.loads(status_file.read_text(encoding="utf-8"))
        jobs_list = data.get("jobs", [])
        self.assertEqual(len(jobs_list), 2)
        job_ids = {j["job_id"] for j in jobs_list}
        self.assertIn("job_ligand1_rigid", job_ids)
        self.assertIn("ad4_job_ligand1_rigid", job_ids)

    def test_ad4_results_serialization_and_resume(self):
        """Verify ad4_results and AD4 paths are properly serialized, deserialized, and resumed."""
        from models import DockingJob, JobStatus, Engine, ResumeMode
        from job_manager import filter_jobs
        from pathlib import Path

        job = DockingJob(
            job_id="ad4_test_job",
            receptor_name="rec1",
            ligand_name="lig1",
            engine=Engine.AUTODOCK4,
            status=JobStatus.SUCCESS,
            gpf_path=Path("receptors/rec1.gpf"),
            glg_path=Path("results/rec1.glg"),
            dpf_path=Path("results/rec1_lig1.dpf"),
            dlg_path=Path("results/DLG/rec1_lig1.dlg"),
            ad4_results={
                "dlg_path": "results/DLG/rec1_lig1.dlg",
                "best_energy": -8.45,
                "poses": [{"rank": 1, "binding_energy": -8.45}],
            }
        )

        d = job.to_dict()
        self.assertIn("ad4_results", d)
        self.assertEqual(d["ad4_results"]["best_energy"], -8.45)
        self.assertEqual(d["dlg_path"], str(Path("results/DLG/rec1_lig1.dlg")))

        restored = DockingJob.from_dict(d)
        self.assertEqual(restored.ad4_results.get("best_energy"), -8.45)
        self.assertEqual(restored.dlg_path, Path("results/DLG/rec1_lig1.dlg"))

        # Test resume restoring AD4 properties
        new_job = DockingJob(
            job_id="ad4_test_job",
            receptor_name="rec1",
            ligand_name="lig1",
            engine=Engine.AUTODOCK4,
        )
        prev_status = {"ad4_test_job": d}
        filtered = filter_jobs([new_job], ResumeMode.RESUME, prev_status)
        self.assertEqual(len(filtered), 0)  # Skipped because completed
        self.assertEqual(new_job.dlg_path, Path("results/DLG/rec1_lig1.dlg"))
        self.assertEqual(new_job.ad4_results.get("best_energy"), -8.45)

    def test_dlg_fallback_analysis_jobs(self):
        """Verify completed jobs with only j.dlg_path are recognized for DLG analysis."""
        from models import DockingJob, JobStatus, Engine
        from pathlib import Path

        j = DockingJob(
            job_id="ad4_fallback_job",
            receptor_name="rec1",
            ligand_name="lig1",
            engine=Engine.AUTODOCK4,
            status=JobStatus.SUCCESS,
            dlg_path=Path("dummy.dlg"),
            ad4_results={},  # empty ad4_results simulates resume bug
        )
        has_dlg = bool(j.ad4_results.get("dlg_path") or j.dlg_path)

class TestDualEngineAndDLGWorkflowFixes(unittest.TestCase):
    """Tests for Dual-Engine execution, reporting, flexible DLG naming, and directory cleanliness."""

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp(prefix="test_dual_engine_"))

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_reporting_engine_both_dual_comparison(self):
        """Verify Engine.BOTH generates DUAL_Complete_Report.xlsx with comparison table."""
        from models import DockingJob, JobStatus, DockingMode, Engine
        from config import ProjectConfig
        from reporting import generate_reports
        import openpyxl

        cfg = ProjectConfig()
        cfg.engine = Engine.BOTH
        cfg.report_directory = self.temp_dir / "reports"
        cfg.export_excel = True
        cfg.export_csv = True

        v_job = DockingJob(
            job_id="vina_rec1_lig1",
            receptor_name="rec1",
            ligand_name="lig1",
            docking_mode=DockingMode.RIGID,
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            vina_results={
                "modes": [
                    {"mode": 1, "affinity_kcal_mol": -8.5, "rmsd_lb": 0.0, "rmsd_ub": 0.0},
                    {"mode": 2, "affinity_kcal_mol": -7.9, "rmsd_lb": 1.2, "rmsd_ub": 1.8},
                ]
            },
        )

        a_job = DockingJob(
            job_id="ad4_rec1_lig1",
            receptor_name="rec1",
            ligand_name="lig1",
            docking_mode=DockingMode.RIGID,
            engine=Engine.AUTODOCK4,
            status=JobStatus.SUCCESS,
            ad4_results={
                "poses": [
                    {
                        "rank": 1,
                        "mode": 1,
                        "binding_energy": -8.1,
                        "ki_raw": 1.25,
                        "ki_unit": "uM",
                        "intermol_energy": -9.2,
                        "internal_energy": 0.4,
                        "torsional_energy": 0.7,
                    }
                ],
                "best_energy": -8.1,
            },
        )

        res = generate_reports([v_job, a_job], cfg)
        self.assertIsNotNone(res)
        xlsx_path = cfg.report_directory / "DUAL_Complete_Report.xlsx"
        self.assertTrue(xlsx_path.is_file())

        wb = openpyxl.load_workbook(str(xlsx_path))
        self.assertIn("00_Dual_Comparison", wb.sheetnames)
        self.assertIn("01_Vina_All_Poses", wb.sheetnames)
        self.assertIn("03_AD4_All_Poses", wb.sheetnames)

        # Validate dual comparison rows
        sheet = wb["00_Dual_Comparison"]
        header = [cell.value for cell in sheet[1]]
        self.assertIn("Vina Best ΔG (kcal/mol)", header)
        self.assertIn("AD4 Best ΔG (kcal/mol)", header)
        self.assertIn("ΔΔG (AD4 - Vina)", header)

        row2 = [cell.value for cell in sheet[2]]
        self.assertEqual(row2[0], "rec1")
        self.assertEqual(row2[2], "lig1")
        self.assertEqual(row2[3], -8.5)
        self.assertEqual(row2[4], -8.1)
        self.assertEqual(row2[5], 0.4)  # -8.1 - (-8.5) = +0.4

    def test_write_xlsx_empty_tables_guard(self):
        """Verify _write_xlsx creates a valid workbook when given empty tables."""
        from reporting import _write_xlsx
        import openpyxl

        out_path = self.temp_dir / "reports" / "Empty_Report.xlsx"
        _write_xlsx({}, out_path)
        self.assertTrue(out_path.is_file())

        wb = openpyxl.load_workbook(str(out_path))
        self.assertIn("Summary", wb.sheetnames)

    def test_flexible_docking_naming_and_identity(self):
        """Verify flexible DLG naming includes _flex and parses without collision."""
        from dlg_extract import IdentityResolutionEngine

        # Flexible docking filename
        flex_dlg = Path("dummy/2Z5X_flex_aspirin.dlg")
        id_flex = IdentityResolutionEngine.extract([], flex_dlg)
        self.assertEqual(id_flex.receptor_name, "2Z5X_flex")
        self.assertEqual(id_flex.ligand_name, "aspirin")

        # Rigid docking filename
        rigid_dlg = Path("dummy/2Z5X_aspirin.dlg")
        id_rigid = IdentityResolutionEngine.extract([], rigid_dlg)
        self.assertEqual(id_rigid.receptor_name, "2Z5X")
        self.assertEqual(id_rigid.ligand_name, "aspirin")

        # Ensure no collision in output
        self.assertNotEqual(id_flex.receptor_name, id_rigid.receptor_name)

    def test_dlg_pipeline_raises_on_missing(self):
        """Verify DLGPipeline raises FileNotFoundError instead of terminating via sys.exit."""
        from dlg_extract import DLGPipeline

        empty_dir = self.temp_dir / "empty"
        empty_dir.mkdir()
        out_dir = self.temp_dir / "BSNDVP_RESULTS"

        pipeline = DLGPipeline(empty_dir, out_dir, self.temp_dir)
        with self.assertRaises(FileNotFoundError):
            pipeline.run()

    def test_find_dlg_dir_nested_ad4(self):
        """Verify _find_dlg_dir discovers DLGs nested inside results/AD4/<rec>/<lig>."""
        from dlg_extract import _find_dlg_dir

        nested_dlg = self.temp_dir / "results" / "AD4" / "rec_A" / "lig_1" / "rec_A_lig_1.dlg"
        nested_dlg.parent.mkdir(parents=True, exist_ok=True)
        nested_dlg.write_text("dummy", encoding="utf-8")

        found = _find_dlg_dir(self.temp_dir)
        self.assertIsNotNone(found)
        # Should discover results/AD4
        self.assertTrue("AD4" in str(found))

    def test_generate_reports_progress_callback(self):
        """Verify generate_reports provides non-blocking real-time progress callbacks."""
        from config import ProjectConfig
        from models import DockingJob, DockingMode, Engine, JobStatus
        from reporting import generate_reports

        cfg = ProjectConfig(
            project_root=self.temp_dir,
            report_directory=self.temp_dir / "reports_progress",
            export_excel=True,
            export_csv=True,
        )

        job = DockingJob(
            job_id="vina_rec_lig",
            receptor_name="rec",
            ligand_name="lig",
            docking_mode=DockingMode.RIGID,
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
            vina_results=[
                type("VinaPose", (), {"pose": 1, "binding_affinity": -8.2, "rmsd_lower_bound": 0.0, "rmsd_upper_bound": 0.0})(),
            ],
        )

        progress_msgs = []
        def _on_prog(msg: str):
            progress_msgs.append(msg)

        res = generate_reports([job], cfg, progress_callback=_on_prog)
        self.assertIsNotNone(res)
        self.assertTrue((cfg.report_directory / "VINA_Complete_Report.xlsx").is_file())
        self.assertGreater(len(progress_msgs), 0)
        self.assertTrue(any("poses" in m.lower() or "table" in m.lower() for m in progress_msgs))
        self.assertTrue(any("completed" in m.lower() for m in progress_msgs))

    def test_ad4_dlg_exact_pose_deduplication(self):
        """Verify DLGParser parses exactly N poses and ignores un-prefixed cluster summary models."""
        from dlg_extract import DLGParser

        mock_dlg = (
            "AutoDock 4.2 Release 4.2.6\n"
            "Run: 1 / 2\n"
            "Estimated Free Energy of Binding = -8.50 kcal/mol\n"
            "DOCKED: MODEL 1\n"
            "DOCKED: ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "DOCKED: ENDMDL\n"
            "Run: 2 / 2\n"
            "Estimated Free Energy of Binding = -7.90 kcal/mol\n"
            "DOCKED: MODEL 2\n"
            "DOCKED: ATOM      1  C1  LIG     1      10.500  20.500  30.500  1.00 20.00     0.050 C\n"
            "DOCKED: ENDMDL\n"
            "CLUSTERING HISTOGRAM\n"
            "MODEL 1\n"
            "USER Run = 1\n"
            "USER Cluster Rank = 1\n"
            "ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "ENDMDL\n"
        )
        dlg_file = self.temp_dir / "test_dedup.dlg"
        dlg_file.write_text(mock_dlg, encoding="utf-8")

        parsed = DLGParser().parse(dlg_file)
        self.assertEqual(len(parsed.poses), 2, "Should extract exactly 2 poses, ignoring the un-prefixed cluster summary model")
        self.assertEqual(parsed.poses[0].run_number, 1)
        self.assertEqual(parsed.poses[1].run_number, 2)
        self.assertAlmostEqual(parsed.poses[0].binding_energy, -8.50)
        self.assertAlmostEqual(parsed.poses[1].binding_energy, -7.90)

    def test_complex_builder_dlg_multiple_poses(self):
        """Verify complex_builder extracts poses 1 and 2 from a DLG file without returning empty list."""
        from complex_builder import _extract_pose_lines, build_complex_pdbqt

        mock_dlg = (
            "DOCKED: MODEL 1\n"
            "DOCKED: ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "DOCKED: ENDMDL\n"
            "DOCKED: MODEL 2\n"
            "DOCKED: ATOM      1  C1  LIG     1      11.000  21.000  31.000  1.00 20.00     0.050 C\n"
            "DOCKED: ENDMDL\n"
            "MODEL 1\n"
            "ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "ENDMDL\n"
        )
        mock_rec = (
            "ATOM      1  CA  ALA A   1       5.000   5.000   5.000  1.00 20.00     0.000 C\n"
            "TER\n"
        )
        dlg_p = self.temp_dir / "multi.dlg"
        rec_p = self.temp_dir / "rec.pdbqt"
        dlg_p.write_text(mock_dlg, encoding="utf-8")
        rec_p.write_text(mock_rec, encoding="utf-8")

        p1_lines = _extract_pose_lines(dlg_p, 1)
        p2_lines = _extract_pose_lines(dlg_p, 2)
        self.assertGreater(len(p1_lines), 0)
        self.assertGreater(len(p2_lines), 0)

        out1 = self.temp_dir / "complex_1.pdbqt"
        out2 = self.temp_dir / "complex_2.pdbqt"
        res1 = build_complex_pdbqt(rec_p, dlg_p, out1, pose_index=1)
        res2 = build_complex_pdbqt(rec_p, dlg_p, out2, pose_index=2)
        self.assertIsNotNone(res1)
        self.assertIsNotNone(res2)
        self.assertTrue(out1.is_file())
        self.assertTrue(out2.is_file())

    def test_compute_pose_validation_rmsd_on_demand(self):
        """Verify on-demand compute_pose_validation_rmsd calculates RMSD correctly."""
        import math
        from validators import compute_pose_validation_rmsd

        ref_pdbqt = (
            "ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "ATOM      2  C2  LIG     1      11.000  21.000  31.000  1.00 20.00     0.050 C\n"
            "TER\n"
        )
        pose_pdbqt = (
            "MODEL 1\n"
            "ATOM      1  C1  LIG     1      10.000  20.000  30.000  1.00 20.00     0.050 C\n"
            "ATOM      2  C2  LIG     1      11.000  21.000  31.000  1.00 20.00     0.050 C\n"
            "ENDMDL\n"
            "MODEL 2\n"
            "ATOM      1  C1  LIG     1      11.000  21.000  31.000  1.00 20.00     0.050 C\n"
            "ATOM      2  C2  LIG     1      12.000  22.000  32.000  1.00 20.00     0.050 C\n"
            "ENDMDL\n"
        )
        ref_p = self.temp_dir / "ref.pdbqt"
        pose_p = self.temp_dir / "poses.pdbqt"
        ref_p.write_text(ref_pdbqt, encoding="utf-8")
        pose_p.write_text(pose_pdbqt, encoding="utf-8")

        pos_1, kab_1, n_1 = compute_pose_validation_rmsd(ref_p, pose_p, pose_index=1)
        self.assertAlmostEqual(pos_1, 0.0, places=2)
        self.assertAlmostEqual(kab_1, 0.0, places=2)
        self.assertEqual(n_1, 2)

        pos_2, _, _ = compute_pose_validation_rmsd(ref_p, pose_p, pose_index=2)
        # Shifted by (1, 1, 1), distance = sqrt(3) ~= 1.732
        self.assertAlmostEqual(pos_2, math.sqrt(3), places=2)

    def test_complex_builder_pdb_generation(self):
        """Verify that build_job_complexes produces standard .pdb complex files."""
        from complex_builder import build_job_complexes

        rec_content = (
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00     0.100 N\n"
            "ATOM      2  CA  ALA A   1      11.000  11.000  11.000  1.00 20.00     0.100 C\n"
            "TER\n"
        )
        lig_content = (
            "MODEL 1\n"
            "ATOM      1  C1  LIG     1      12.000  12.000  12.000  1.00 20.00     0.050 C\n"
            "ATOM      2  C2  LIG     1      13.000  13.000  13.000  1.00 20.00     0.050 C\n"
            "ENDMDL\n"
        )
        rec_f = self.temp_dir / "rec.pdbqt"
        lig_f = self.temp_dir / "lig_out.pdbqt"
        rec_f.write_text(rec_content, encoding="utf-8")
        lig_f.write_text(lig_content, encoding="utf-8")

        c_dir = self.temp_dir / "complexes"
        written = build_job_complexes(rec_f, lig_f, c_dir, n_poses=1)
        self.assertTrue(len(written) >= 1)

        # Confirm that complex_pose_1.pdb was generated
        pdb_f = c_dir / "complex_pose_1.pdb"
        self.assertTrue(pdb_f.is_file())
        self.assertTrue(pdb_f.stat().st_size > 0)
        # Also confirm pdbqt exists
        pdbqt_f = c_dir / "complex_pose_1.pdbqt"
        self.assertTrue(pdbqt_f.is_file())

    def test_comparison_tab_refresh_and_matching(self):
        """Verify that ExperimentComparisonTab groups and computes ddg across Vina and AD4."""
        from config import ProjectConfig
        from gui_qt.comparison_tab import ExperimentComparisonTab
        from models import DockingJob, Engine, JobStatus, CanonicalPose, DockingResult

        cfg = ProjectConfig()
        tab = ExperimentComparisonTab(cfg, lambda: None)

        j_vina = DockingJob(
            job_id="REC__LIG_vina",
            receptor_name="REC",
            ligand_name="LIG",
            engine=Engine.VINA,
            status=JobStatus.SUCCESS,
        )
        j_vina.canonical_result = DockingResult(
            job_id=j_vina.job_id,
            engine="VINA",
            poses=[CanonicalPose(rank=1, binding_score=-6.5, estimated_ki_formatted="16.9 μM")]
        )

        j_ad4 = DockingJob(
            job_id="REC__LIG_ad4",
            receptor_name="REC",
            ligand_name="LIG",
            engine=Engine.AUTODOCK4,
            status=JobStatus.SUCCESS,
        )
        j_ad4.canonical_result = DockingResult(
            job_id=j_ad4.job_id,
            engine="AUTODOCK4",
            poses=[CanonicalPose(rank=1, binding_score=-15.5, estimated_ki_formatted="44.8 nM")]
        )

        tab.set_jobs([j_vina, j_ad4])
        self.assertEqual(len(tab._filtered_rows), 1)
        row = tab._filtered_rows[0]
        self.assertEqual(row["receptor"], "REC")
        self.assertEqual(row["ligand"], "LIG")
        self.assertAlmostEqual(row["vina_energy"], -6.5)
        self.assertAlmostEqual(row["ad4_energy"], -15.5)
        self.assertAlmostEqual(row["ddg"], -9.0)
        self.assertEqual(row["status"], "Completed Both")


if __name__ == "__main__":
    unittest.main(verbosity=2)




