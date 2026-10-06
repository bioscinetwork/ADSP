"""
AutoDockSuite Pro — Level 3 Scientific Platform & UI Upgrade Test Suite
========================================================================
Comprehensive regression and functional test suite validating:
  1. Design system, color tokens, and 5 research themes.
  2. Publication-quality HTML report generation with inline SVG charts.
  3. Experiment Comparison studio logic and ΔΔG scoring differences.
  4. Platform grouped settings round-trip serialization.
  5. Modern MainWindow sidebar navigation and stacked workspace indexing.
  6. Laboratory KPI cards and scientific workflow progression steppers.
  7. Runner execution monitor, timeline stages, and dual-engine badges.
"""

import sys
import tempfile
from pathlib import Path
import pytest

from PySide6.QtWidgets import QApplication

# Ensure QApplication singleton for GUI widget tests
app = QApplication.instance()
if not app:
    app = QApplication(sys.argv)

from config import ProjectConfig
from models import (
    DockingJob,
    DockingMode,
    DockingResult,
    Engine,
    JobStatus,
    CanonicalPose,
    ClusterInfo,
    ThermodynamicAnalysis,
    ValidationMetrics,
    ProvenanceRecord,
    VinaMetrics,
    AutoDock4Metrics,
)
from gui_qt.styles import (
    THEME_NAMES,
    THEME_DEFINITIONS,
    get_stylesheet,
    resolve_theme_name,
    ThemePreviewWidget,
)
from reporting import generate_html_report, _generate_energy_chart_svg, _generate_cluster_chart_svg
from gui_qt.comparison_tab import ExperimentComparisonTab
from gui_qt.settings_tab import SettingsTab
from gui_qt.app import MainWindow
from gui_qt.workspace_tab import WorkspaceTab
from gui_qt.preparation_tab import PreparationTab
from gui_qt.runner_tab import RunnerTab
from gui_qt.results_tab import ResultsTab


class TestThemeSystem:
    """Validate Level 3 Design Tokens and Theme Architecture."""

    def test_all_five_themes_registered(self):
        assert len(THEME_NAMES) == 5
        assert "Dark Studio" in THEME_NAMES
        assert "Scientific Light" in THEME_NAMES
        assert "Molecular Plasma" in THEME_NAMES
        assert "Bio-Neutral" in THEME_NAMES
        assert "Noir" in THEME_NAMES

    def test_legacy_theme_name_resolution(self):
        assert resolve_theme_name("Dark") == "Dark Studio"
        assert resolve_theme_name("Scientific") == "Scientific Light"
        assert resolve_theme_name("Plasma") == "Molecular Plasma"
        assert resolve_theme_name("Bio-Neon") == "Bio-Neutral"
        assert resolve_theme_name("Noir") == "Noir"
        assert resolve_theme_name("UnknownCustom") == "Dark Studio"

    def test_all_themes_generate_valid_stylesheets(self):
        for name in THEME_NAMES:
            qss = get_stylesheet(name)
            assert isinstance(qss, str)
            assert len(qss) > 2000
            assert "QMainWindow" in qss
            assert "#cardFrame" in qss
            assert "#navButton" in qss
            assert "#sidebarFrame" in qss
            assert "QTableWidget" in qss

    def test_theme_preview_widget(self):
        widget = ThemePreviewWidget()
        for name in THEME_NAMES:
            widget.set_theme(name)
            assert widget._theme.name == name


class TestPublicationReporting:
    """Validate Level 3 standalone research reports and inline vector charts."""

    def test_svg_chart_generators(self):
        # Energy chart
        sample_rows = [
            {"Ligand": "Compound_A", "Best_Affinity_kcal_mol": -9.42},
            {"Ligand": "Compound_B", "Best_Affinity_kcal_mol": -8.15},
            {"Ligand": "Compound_C", "Best_Affinity_kcal_mol": -6.50},
        ]
        svg_e = _generate_energy_chart_svg(sample_rows)
        assert "<svg" in svg_e
        assert "</svg>" in svg_e
        assert "Compound_A" in svg_e
        assert "-9.42" in svg_e

        # Cluster chart
        cluster_rows = [
            {"Cluster_ID": 1, "Population_%": 45.0, "Size": 45, "Lowest_Energy": -8.90},
            {"Cluster_ID": 2, "Population_%": 25.0, "Size": 25, "Lowest_Energy": -7.60},
        ]
        svg_c = _generate_cluster_chart_svg(cluster_rows)
        assert "<svg" in svg_c
        assert "</svg>" in svg_c
        assert "Cluster 1" in svg_c
        assert "45.0%" in svg_c

    def test_generate_html_report_complete(self, tmp_path):
        cfg = ProjectConfig()
        cfg.project_name = "Benchmark_Project"

        job = DockingJob(
            job_id="job_001",
            receptor_name="Target_Kinase",
            ligand_name="Inhibitor_X",
            engine=Engine.AUTODOCK4,
            docking_mode=DockingMode.RIGID,
            status=JobStatus.SUCCESS,
            elapsed_seconds=12.5,
        )
        res = DockingResult(
            job_id="job_001",
            engine=Engine.AUTODOCK4,
            docking_mode=DockingMode.RIGID,
            best_binding_energy=-9.15,
            poses=[
                CanonicalPose(
                    rank=1,
                    binding_energy=-9.15,
                    estimated_ki_nM=195.0,
                    estimated_ki_formatted="195.00 nM",
                    ad4_metrics=AutoDock4Metrics(
                        run_number=14,
                        cluster_id=1,
                        cluster_size=42,
                        cluster_rmsd=0.65,
                        intermolecular_energy=-10.20,
                        internal_energy=-0.40,
                        torsional_energy=1.45,
                        unbound_energy=-0.40,
                    ),
                )
            ],
            clusters=[
                ClusterInfo(
                    cluster_id=1,
                    size=42,
                    lowest_energy=-9.15,
                    mean_energy=-8.85,
                    cluster_rmsd=0.65,
                    representative_run=14,
                    population_percent=42.0,
                )
            ],
            validation=ValidationMetrics(
                validation_performed=True,
                reference_ligand="1ABC_ligand.pdbqt",
                crystal_rmsd=1.24,
                matched_atoms=24,
            ),
            thermodynamics=ThermodynamicAnalysis(
                info_entropy=2.14,
                partition_function=1.54,
                stat_temperature=298.15,
            ),
            provenance=ProvenanceRecord(
                engine="AUTODOCK4",
                engine_version="4.2.6",
                receptor_hash="abc123def456",
                ligand_hash="fed654cba321",
            ),
        )
        job.canonical_result = res

        out_html = tmp_path / "docking_analysis_report.html"
        generated = generate_html_report([job], cfg, out_html)

        assert generated.is_file()
        content = generated.read_text(encoding="utf-8")
        assert "AutoDock Suite Pro" in content
        assert "Target_Kinase" in content
        assert "Inhibitor_X" in content
        assert "-9.15" in content
        assert "195.00 nM" in content
        assert "Cluster Analysis" in content
        assert "Crystallographic Reference Ligand Validation" in content
        assert "Thermodynamic &amp; Statistical Descriptors" in content
        assert "Reproducibility &amp; Provenance" in content
        assert "Scientific Disclaimers" in content
        assert "<svg" in content


class TestExperimentComparison:
    """Validate Cross-Job Comparison and scoring delta logic."""

    def test_comparison_delta_calculation(self):
        v_job = DockingJob(
            job_id="v_001",
            receptor_name="Receptor_A",
            ligand_name="Ligand_1",
            engine=Engine.VINA,
            docking_mode=DockingMode.RIGID,
            status=JobStatus.SUCCESS,
        )
        v_res = DockingResult(
            job_id="v_001",
            engine=Engine.VINA,
            docking_mode=DockingMode.RIGID,
            best_binding_energy=-7.50,
            poses=[CanonicalPose(rank=1, binding_energy=-7.50, vina_metrics=VinaMetrics(rmsd_lower_bound=0.0))],
        )
        v_job.canonical_result = v_res

        a_job = DockingJob(
            job_id="a_001",
            receptor_name="Receptor_A",
            ligand_name="Ligand_1",
            engine=Engine.AUTODOCK4,
            docking_mode=DockingMode.RIGID,
            status=JobStatus.SUCCESS,
        )
        a_res = DockingResult(
            job_id="a_001",
            engine=Engine.AUTODOCK4,
            docking_mode=DockingMode.RIGID,
            best_binding_energy=-8.70,
            poses=[CanonicalPose(rank=1, binding_energy=-8.70, ad4_metrics=AutoDock4Metrics(run_number=1))],
        )
        a_job.canonical_result = a_res

        cfg = ProjectConfig()
        tab = ExperimentComparisonTab(cfg, lambda: None)
        tab._jobs = [v_job, a_job]
        tab._apply_filters()

        assert tab.table_comparison.rowCount() == 1
        item_diff = tab.table_comparison.item(0, 5)  # ΔΔG column
        assert item_diff is not None
        # ΔΔG = AD4 (-8.70) - Vina (-7.50) = -1.20 kcal/mol
        assert "-1.20" in item_diff.text()


class TestMainWindowAndNavigation:
    """Validate MainWindow sidebar architecture and stacked workspaces."""

    def test_mainwindow_sidebar_and_workspaces(self):
        cfg = ProjectConfig()
        win = MainWindow(config=cfg)

        # 8 Workspaces in stack
        assert win.stack.count() == 8
        assert len(win.nav_buttons) == 8

        # Switch to each workspace and verify sync
        for idx in range(8):
            win.set_active_tab(idx)
            assert win.stack.currentIndex() == idx
            assert win.nav_buttons[idx].isChecked()

        # Check backward compatibility alias
        assert win.tab_widget is win.stack

        # Check global status updates
        win.set_global_status("READY", level="ready")
        assert "READY" in win.lbl_global_status.text()
        win.set_global_status("RUNNING", level="run")
        assert "RUNNING" in win.lbl_global_status.text()
        win.set_global_status("COMPLETED", level="success")
        assert "COMPLETED" in win.lbl_global_status.text()

        # Check notification banner
        win.notify("Test Notification", level="success", duration_ms=1000)
        assert not win.notification_banner.isHidden()
        assert "Test Notification" in win.notification_banner.lbl_msg.text()


class TestWorkspaceAndRunnerEnhancements:
    """Validate KPI summary cards, workflow progression stepper, and runner timeline."""

    def test_workspace_kpi_cards(self, tmp_path):
        cfg = ProjectConfig()
        cfg.project_root = tmp_path
        tab = WorkspaceTab(cfg, lambda: None)

        assert hasattr(tab, "lbl_kpi_receptors")
        assert hasattr(tab, "lbl_kpi_ligands")
        assert hasattr(tab, "lbl_kpi_completed")
        assert hasattr(tab, "lbl_kpi_status")

    def test_preparation_stepper(self):
        cfg = ProjectConfig()
        tab = PreparationTab(cfg, lambda: None)

        assert hasattr(tab, "prep_stepper")
        assert hasattr(tab, "stepper_labels")
        assert len(tab.stepper_labels) == 7

        # Test stage transitions
        tab._update_stepper_stage(0)
        tab._update_stepper_stage(2)
        assert "✓" in tab.stepper_labels[0].text()
        assert "✓" in tab.stepper_labels[1].text()
        tab._update_stepper_stage(6)
        assert "✓" in tab.stepper_labels[5].text()

    def test_runner_timeline_and_badges(self):
        cfg = ProjectConfig()
        tab = RunnerTab(cfg, lambda: None, lambda jobs: None)

        assert hasattr(tab, "lbl_mon_engine_mode")
        assert hasattr(tab, "lbl_mon_vina_badge")
        assert hasattr(tab, "lbl_mon_ad4_badge")
        assert hasattr(tab, "timeline_labels")
        assert len(tab.timeline_labels) == 6

        # Test timeline updates
        tab._update_timeline_stage(2)
        assert "✓" in tab.timeline_labels[0][0].text()
        assert "●" in tab.timeline_labels[2][0].text()
        assert "○" in tab.timeline_labels[3][0].text()

    def test_results_vector_charts_integrated(self):
        cfg = ProjectConfig()
        tab = ResultsTab(cfg, lambda: None)

        assert hasattr(tab, "ov_chart_card")
        assert hasattr(tab, "cluster_chart_card")

        poses = [
            CanonicalPose(rank=1, binding_energy=-8.9),
            CanonicalPose(rank=2, binding_energy=-7.8),
        ]
        svg_poses = tab._generate_poses_chart_svg(poses, "VINA")
        assert "<svg" in svg_poses
        assert "-8.9" in svg_poses

        clusters = [
            ClusterInfo(cluster_id=1, size=30, lowest_energy=-8.9, mean_energy=-8.5, population_percent=60.0),
        ]
        svg_clusters = tab._generate_clusters_chart_svg(clusters)
        assert "<svg" in svg_clusters
        assert "Cluster 1" in svg_clusters
