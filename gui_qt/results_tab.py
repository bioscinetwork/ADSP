"""
AutoDock Suite Pro — PySide6 Results & Analysis Tab (gui_qt/results_tab.py)
===========================================================================
Engine-aware, publication-quality molecular docking analysis suite:
  - VIEW 1: Overview (concise scientific summary with engine-aware metrics)
  - VIEW 2: Poses (engine-aware columns + advanced energetics toggle)
  - VIEW 3: Clusters (dedicated AutoDock4 LGA clustering interface)
  - VIEW 4: Validation (docking-search RMSD separated from crystal reference validation)
  - VIEW 5: Interactions (thread-safe 2D diagram + non-covalent contacts table)
  - VIEW 6: Thermodynamics (statistical mechanics descriptors from DLG analysis)
  - VIEW 7: ADMET (compound physicochemical properties & rule-of-5 compliance)
  - VIEW 8: Provenance (reproducible execution manifest with SHA-256 hashes)
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import threading
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List, Optional, Tuple, Any

from PySide6.QtCore import QByteArray, QObject, Qt, Signal
from PySide6.QtGui import QColor, QFont, QGuiApplication
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

# Optional WebEngine / SVG imports
try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    HAS_WEBENGINE = True
except ImportError:
    HAS_WEBENGINE = False

try:
    from PySide6.QtSvgWidgets import QSvgWidget
    HAS_SVG = True
except ImportError:
    HAS_SVG = False

if TYPE_CHECKING:
    from config import ProjectConfig
    from models import DockingJob, DockingResult, CanonicalPose, ClusterInfo


class _NumericTableItem(QTableWidgetItem):
    """QTableWidgetItem supporting numeric sorting."""
    def __init__(self, text: str, sort_val: float):
        super().__init__(text)
        self.sort_val = sort_val

    def __lt__(self, other):
        if isinstance(other, _NumericTableItem):
            return self.sort_val < other.sort_val
        return super().__lt__(other)


class _ResultsSignals(QObject):
    diagram_ready = Signal(object, str, int, int)  # (InteractionDiagramResult, job_id, pose_idx, generation)
    interactions_ready = Signal(list, str, int, int)  # (interactions, job_id, pose_idx, generation)
    admet_ready = Signal(dict)
    export_progress = Signal(str)
    export_complete = Signal(dict, str)
    dlg_complete = Signal(str, str)


class ResultsTab(QWidget):
    """Engine-Aware Scientific Results & 2D/3D Interaction Analysis in PySide6."""

    JOB_COLS = ("Job ID", "Receptor", "Ligand", "Engine", "Mode", "Status", "Best Engine Score (kcal/mol)")
    INTER_COLS = ("Pose", "Type", "Residue", "Rec.Atom", "Lig.Atom", "Dist (Å)", "Detection Basis")

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self._jobs: List["DockingJob"] = []
        self._filtered_jobs: List["DockingJob"] = []
        self._current_job: Optional["DockingJob"] = None
        self._current_pose: int = 1
        self._current_interactions: List = []
        self._current_diagram_res: Optional[Any] = None
        self._current_svg: str = ""
        self._current_html: str = ""
        self._render_generation = 0
        self._interactions_cache: Dict[Tuple[str, int, Optional[int], str], List] = {}

        self._signals = _ResultsSignals()
        self._signals.diagram_ready.connect(self._on_diagram_ready)
        self._signals.interactions_ready.connect(self._on_interactions_ready)
        self._signals.admet_ready.connect(self._on_admet_ready)
        self._signals.export_progress.connect(self._on_export_progress)
        self._signals.export_complete.connect(self._on_export_complete)
        self._signals.dlg_complete.connect(self._on_dlg_complete)

        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        # ── Header & Action Toolbar ────────────────────────────
        hdr_layout = QHBoxLayout()
        lbl_title = QLabel("Virtual Screening Results & Scientific Profiler")
        lbl_title.setObjectName("sectionHeader")
        hdr_layout.addWidget(lbl_title)
        hdr_layout.addStretch()

        self.btn_refresh = QPushButton("⟳ Refresh Results")
        self.btn_refresh.setObjectName("accentButton")
        self.btn_refresh.clicked.connect(self.refresh)
        hdr_layout.addWidget(self.btn_refresh)

        self.btn_run_dlg = QPushButton("🔬 Run DLG Extract")
        self.btn_run_dlg.setToolTip("Run BSNDVP™ DLG analysis pipeline to extract clusters and thermodynamic descriptors")
        self.btn_run_dlg.clicked.connect(self._run_dlg_extract)
        hdr_layout.addWidget(self.btn_run_dlg)

        self.btn_export = QPushButton("📤 Export Reports")
        self.btn_export.clicked.connect(self._export_reports)
        hdr_layout.addWidget(self.btn_export)

        self.btn_3d_webgl = QPushButton("🌐 3D WebGL Inspector")
        self.btn_3d_webgl.setObjectName("accentButton")
        self.btn_3d_webgl.clicked.connect(self._launch_3d_viewer)
        hdr_layout.addWidget(self.btn_3d_webgl)

        self.btn_pymol = QPushButton("🔮 PyMOL Inspector")
        self.btn_pymol.clicked.connect(self._launch_pymol)
        hdr_layout.addWidget(self.btn_pymol)

        self.lbl_status_action = QLabel("")
        self.lbl_status_action.setStyleSheet("color: #38bdf8; font-size: 11px; font-weight: bold;")
        hdr_layout.addWidget(self.lbl_status_action)

        main_layout.addLayout(hdr_layout)

        # ── KPI Card Bar ───────────────────────────────────────
        kpi_card = QFrame()
        kpi_card.setObjectName("cardFrame")
        kpi_l = QHBoxLayout(kpi_card)
        kpi_l.setContentsMargins(16, 10, 16, 10)

        self.lbl_kpi_total = self._create_kpi_cell(kpi_l, "Total Jobs", "—")
        self.lbl_kpi_success = self._create_kpi_cell(kpi_l, "Completed", "—", "#10b981")
        self.lbl_kpi_failed = self._create_kpi_cell(kpi_l, "Failed", "—", "#ef4444")
        self.lbl_kpi_lead = self._create_kpi_cell(kpi_l, "Top Lead", "Run docking to see results", "#38bdf8")

        main_layout.addWidget(kpi_card)

        # ── Main Splitter: Left (Jobs Table), Right (8 Scientific Views) ──
        splitter = QSplitter(Qt.Horizontal)

        # ── Left Column: Jobs Table ──
        left_widget = QWidget()
        left_l = QVBoxLayout(left_widget)
        left_l.setContentsMargins(0, 0, 0, 0)
        left_l.setSpacing(6)

        job_bar = QHBoxLayout()
        lbl_jobs_title = QLabel("DOCKING JOBS:")
        lbl_jobs_title.setObjectName("subHeader")
        job_bar.addWidget(lbl_jobs_title)
        self.ent_search_jobs = QLineEdit()
        self.ent_search_jobs.setPlaceholderText("Filter by ligand, receptor, engine...")
        self.ent_search_jobs.textChanged.connect(self._apply_jobs_filter)
        job_bar.addWidget(self.ent_search_jobs, 1)
        left_l.addLayout(job_bar)

        self.table_jobs = QTableWidget()
        self.table_jobs.setColumnCount(len(self.JOB_COLS))
        self.table_jobs.setHorizontalHeaderLabels(self.JOB_COLS)
        self.table_jobs.setAlternatingRowColors(True)
        self.table_jobs.setSelectionBehavior(QTableWidget.SelectRows)
        self.table_jobs.horizontalHeader().setStretchLastSection(True)
        self.table_jobs.setSortingEnabled(True)
        self.table_jobs.itemSelectionChanged.connect(self._on_job_selected)
        left_l.addWidget(self.table_jobs, 1)

        splitter.addWidget(left_widget)

        # ── Right Column: Scientific Views TabWidget ──
        right_widget = QWidget()
        right_l = QVBoxLayout(right_widget)
        right_l.setContentsMargins(0, 0, 0, 0)
        right_l.setSpacing(6)

        self.tab_views = QTabWidget()
        self.tab_views.setObjectName("scientificTabs")

        # 8 Scientific Views
        self._build_view_overview()
        self._build_view_poses()
        self._build_view_clusters()
        self._build_view_validation()
        self._build_view_interactions()
        self._build_view_thermo()
        self._build_view_admet()
        self._build_view_provenance()

        right_l.addWidget(self.tab_views)
        splitter.addWidget(right_widget)

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        main_layout.addWidget(splitter, 1)

    def _create_kpi_cell(self, layout, title: str, default: str, val_color: str = "") -> QLabel:
        cell = QVBoxLayout()
        lbl_t = QLabel(title.upper())
        lbl_t.setObjectName("kpiTitle")
        lbl_v = QLabel(default)
        lbl_v.setObjectName("kpiValue")
        if val_color:
            lbl_v.setStyleSheet(f"color: {val_color};")
        cell.addWidget(lbl_t)
        cell.addWidget(lbl_v)
        layout.addLayout(cell)
        layout.addStretch()
        return lbl_v

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 1 — OVERVIEW
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_overview(self) -> None:
        self.view_overview = QWidget()
        layout = QVBoxLayout(self.view_overview)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        # Engine & Status Banner Card
        self.ov_banner_card = QFrame()
        self.ov_banner_card.setObjectName("cardFrame")
        b_l = QHBoxLayout(self.ov_banner_card)
        b_l.setContentsMargins(14, 10, 14, 10)

        self.lbl_ov_job_title = QLabel("Select a docking job to inspect scientific results")
        self.lbl_ov_job_title.setStyleSheet("font-size: 14px; font-weight: bold; color: #f8fafc;")
        b_l.addWidget(self.lbl_ov_job_title)
        b_l.addStretch()

        self.lbl_ov_engine_badge = QLabel("Engine: —")
        self.lbl_ov_engine_badge.setStyleSheet("background: #0284c7; color: white; padding: 4px 10px; border-radius: 4px; font-weight: bold; font-size: 11px;")
        b_l.addWidget(self.lbl_ov_engine_badge)

        self.lbl_ov_mode_badge = QLabel("Mode: —")
        self.lbl_ov_mode_badge.setStyleSheet("background: #334155; color: #94a3b8; padding: 4px 10px; border-radius: 4px; font-weight: bold; font-size: 11px;")
        b_l.addWidget(self.lbl_ov_mode_badge)

        layout.addWidget(self.ov_banner_card)

        # Metrics Card Grid
        self.ov_grid_card = QFrame()
        self.ov_grid_card.setObjectName("cardFrame")
        self.ov_grid = QGridLayout(self.ov_grid_card)
        self.ov_grid.setContentsMargins(16, 14, 16, 14)
        self.ov_grid.setHorizontalSpacing(24)
        self.ov_grid.setVerticalSpacing(14)

        self.ov_cells: Dict[str, Tuple[QLabel, QLabel]] = {}
        metric_keys = [
            ("Runs / Modes", "Total sampling steps performed by the docking engine"),
            ("Docked Poses", "Number of distinct coordinate conformations recorded"),
            ("Clusters", "Conformational clusters identified by RMSD tolerance"),
            ("Engine Score / Affinity", "Engine-reported score; semantics depend on the selected engine"),
            ("AutoDock4 Ki", "Estimated inhibition constant reported by the AutoDock4 DLG; not applicable to Vina"),
            ("Dominant Cluster", "Largest cluster population fraction and run size"),
            ("Best Cluster RMSD", "Average internal RMSD of the lowest-energy cluster"),
            ("Validation RMSD", "Structural RMSD from reference crystal ligand (if provided)"),
        ]

        for idx, (title, tooltip) in enumerate(metric_keys):
            row = idx // 2
            col = (idx % 2) * 2

            t_lbl = QLabel(title.upper())
            t_lbl.setStyleSheet("color: #94a3b8; font-size: 10px; font-weight: bold; letter-spacing: 0.5px;")
            t_lbl.setToolTip(tooltip)

            v_lbl = QLabel("—")
            v_lbl.setStyleSheet("font-size: 13px; font-weight: bold; color: #f8fafc;")

            self.ov_grid.addWidget(t_lbl, row, col)
            self.ov_grid.addWidget(v_lbl, row, col + 1)
            self.ov_cells[title] = (t_lbl, v_lbl)

        layout.addWidget(self.ov_grid_card)

        # Scientific Notes / Description Card
        self.ov_notes_card = QFrame()
        self.ov_notes_card.setObjectName("cardFrame")
        n_l = QVBoxLayout(self.ov_notes_card)
        n_l.setContentsMargins(14, 10, 14, 10)
        self.lbl_ov_notes = QLabel("No active job selected.")
        self.lbl_ov_notes.setWordWrap(True)
        self.lbl_ov_notes.setStyleSheet("color: #cbd5e1; font-size: 11px; line-height: 1.4;")
        n_l.addWidget(self.lbl_ov_notes)
        layout.addWidget(self.ov_notes_card)

        # Conformation Energy Spectrum Chart Card
        self.ov_chart_card = QFrame()
        self.ov_chart_card.setObjectName("cardFrame")
        c_l = QVBoxLayout(self.ov_chart_card)
        c_l.setContentsMargins(12, 10, 12, 10)
        c_l.setSpacing(6)
        lbl_c_title = QLabel("CONFORMATION ENERGETIC SPECTRUM:")
        lbl_c_title.setObjectName("subHeader")
        c_l.addWidget(lbl_c_title)

        if HAS_SVG:
            self.ov_chart_svg = QSvgWidget()
            self.ov_chart_svg.setFixedHeight(170)
            c_l.addWidget(self.ov_chart_svg)
        else:
            self.lbl_ov_chart_fallback = QLabel("SVG visualization requires QtSvgWidgets.")
            c_l.addWidget(self.lbl_ov_chart_fallback)

        layout.addWidget(self.ov_chart_card)

        layout.addStretch()
        self.tab_views.addTab(self.view_overview, "📊 Overview")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 2 — POSES
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_poses(self) -> None:
        self.view_poses = QWidget()
        layout = QVBoxLayout(self.view_poses)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        tb = QHBoxLayout()
        self.lbl_poses_header = QLabel("DOCKED CONFORMATIONS:")
        self.lbl_poses_header.setObjectName("subHeader")
        tb.addWidget(self.lbl_poses_header)
        tb.addStretch()

        self.chk_show_advanced = QCheckBox("Show Advanced Energetics & Diagnostics")
        self.chk_show_advanced.setToolTip("Toggle display of intermolecular, internal, torsional, and unbound energy components")
        self.chk_show_advanced.toggled.connect(self._on_toggle_advanced_poses)
        tb.addWidget(self.chk_show_advanced)

        layout.addLayout(tb)

        self.table_poses = QTableWidget()
        self.table_poses.setAlternatingRowColors(True)
        self.table_poses.setSelectionBehavior(QTableWidget.SelectRows)
        self.table_poses.horizontalHeader().setStretchLastSection(True)
        self.table_poses.setSortingEnabled(True)
        self.table_poses.itemSelectionChanged.connect(self._on_pose_selected)
        layout.addWidget(self.table_poses)

        self.tab_views.addTab(self.view_poses, "🎯 Poses")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 3 — CLUSTERS (AutoDock4 LGA Dedicated)
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_clusters(self) -> None:
        self.view_clusters = QWidget()
        layout = QVBoxLayout(self.view_clusters)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        self.lbl_cluster_banner = QLabel("CONFORMATIONAL CLUSTERING ANALYSIS (AutoDock 4 / LGA):")
        self.lbl_cluster_banner.setObjectName("subHeader")
        layout.addWidget(self.lbl_cluster_banner)

        # Cluster Table
        self.CLUSTER_COLS = (
            "Cluster ID", "Members", "Population %", "Lowest Energy (kcal/mol)",
            "Mean Energy (kcal/mol)", "Cluster RMSD (Å)", "Rep. Run", "Rep. Pose"
        )
        self.table_clusters = QTableWidget()
        self.table_clusters.setColumnCount(len(self.CLUSTER_COLS))
        self.table_clusters.setHorizontalHeaderLabels(self.CLUSTER_COLS)
        self.table_clusters.setAlternatingRowColors(True)
        self.table_clusters.setSelectionBehavior(QTableWidget.SelectRows)
        self.table_clusters.horizontalHeader().setStretchLastSection(True)
        self.table_clusters.setSortingEnabled(True)
        self.table_clusters.itemSelectionChanged.connect(self._on_cluster_selected)
        layout.addWidget(self.table_clusters, 1)

        # Cluster details card
        self.cluster_detail_card = QFrame()
        self.cluster_detail_card.setObjectName("cardFrame")
        cd_l = QVBoxLayout(self.cluster_detail_card)
        cd_l.setContentsMargins(12, 10, 12, 10)
        self.lbl_cluster_detail = QLabel("Select a cluster to view run membership and energetic distribution.")
        self.lbl_cluster_detail.setStyleSheet("color: #cbd5e1; font-size: 11px;")
        self.lbl_cluster_detail.setWordWrap(True)
        cd_l.addWidget(self.lbl_cluster_detail)

        btn_bar = QHBoxLayout()
        self.btn_inspect_cluster_pose = QPushButton("🎯 Inspect Cluster Representative Pose")
        self.btn_inspect_cluster_pose.setEnabled(False)
        self.btn_inspect_cluster_pose.clicked.connect(self._inspect_cluster_rep_pose)
        btn_bar.addWidget(self.btn_inspect_cluster_pose)
        btn_bar.addStretch()
        cd_l.addLayout(btn_bar)

        layout.addWidget(self.cluster_detail_card)

        # Cluster Population Chart Card
        self.cluster_chart_card = QFrame()
        self.cluster_chart_card.setObjectName("cardFrame")
        cc_l = QVBoxLayout(self.cluster_chart_card)
        cc_l.setContentsMargins(12, 10, 12, 10)
        cc_l.setSpacing(6)
        lbl_cc = QLabel("CONFORMATIONAL CLUSTER POPULATION DISTRIBUTION:")
        lbl_cc.setObjectName("subHeader")
        cc_l.addWidget(lbl_cc)

        if HAS_SVG:
            self.cluster_chart_svg = QSvgWidget()
            self.cluster_chart_svg.setFixedHeight(160)
            cc_l.addWidget(self.cluster_chart_svg)
        else:
            self.lbl_cluster_chart_fallback = QLabel("SVG visualization requires QtSvgWidgets.")
            cc_l.addWidget(self.lbl_cluster_chart_fallback)

        layout.addWidget(self.cluster_chart_card)
        self.tab_views.addTab(self.view_clusters, "🧬 Clusters")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 4 — VALIDATION
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_validation(self) -> None:
        self.view_validation = QWidget()
        layout = QVBoxLayout(self.view_validation)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        hdr = QLabel("CRYSTAL STRUCTURE RE-DOCKING & RMSD VALIDATION (OPTIONAL):")
        hdr.setObjectName("subHeader")
        layout.addWidget(hdr)

        # ── Interactive Reference Input Card ──────────────────────
        ref_card = QFrame()
        ref_card.setObjectName("cardFrame")
        r_layout = QVBoxLayout(ref_card)
        r_layout.setContentsMargins(16, 12, 16, 12)
        r_layout.setSpacing(10)

        r_top = QHBoxLayout()
        self.lbl_val_target = QLabel("Target: No docking job selected")
        self.lbl_val_target.setStyleSheet("font-weight: bold; font-size: 13px; color: #38bdf8;")
        r_top.addWidget(self.lbl_val_target)
        r_top.addStretch()

        self.lbl_val_opt_status = QLabel("Validation Status: Optional (Select reference ligand to evaluate)")
        self.lbl_val_opt_status.setStyleSheet("color: #94a3b8; font-size: 11px;")
        r_top.addWidget(self.lbl_val_opt_status)
        r_layout.addLayout(r_top)

        r_row = QHBoxLayout()
        r_row.setSpacing(8)
        r_row.addWidget(QLabel("Reference Ligand:"))
        self.txt_val_ref = QLineEdit()
        self.txt_val_ref.setPlaceholderText("Path to crystallographic native ligand (SDF, MOL2, PDB, PDBQT, CIF)...")
        r_row.addWidget(self.txt_val_ref)

        self.btn_val_browse = QPushButton("📁 Browse Reference...")
        self.btn_val_browse.clicked.connect(self._browse_validation_reference)
        r_row.addWidget(self.btn_val_browse)

        self.btn_val_run = QPushButton("⚡ Validate Selected Ligand Poses")
        self.btn_val_run.setStyleSheet("background: #0284c7; color: white; font-weight: bold; padding: 6px 14px; border-radius: 4px;")
        self.btn_val_run.clicked.connect(self._run_on_demand_validation)
        r_row.addWidget(self.btn_val_run)

        self.btn_val_clear = QPushButton("Reset")
        self.btn_val_clear.clicked.connect(self._clear_validation)
        r_row.addWidget(self.btn_val_clear)

        r_layout.addLayout(r_row)
        layout.addWidget(ref_card)

        # ── Summary Metrics Card ──────────────────────────────────
        val_card = QFrame()
        val_card.setObjectName("cardFrame")
        v_grid = QGridLayout(val_card)
        v_grid.setContentsMargins(16, 14, 16, 14)
        v_grid.setHorizontalSpacing(24)
        v_grid.setVerticalSpacing(14)

        self.val_cells: Dict[str, Tuple[QLabel, QLabel]] = {}
        val_fields = [
            ("Reference Ligand Source", "Co-crystallized native ligand structure used for RMSD validation"),
            ("Validation Status", "Whether experimental reference structure was supplied and matched"),
            ("Validation RMSD", "All-heavy-atom RMSD relative to co-crystallized reference ligand pose"),
            ("Search RMSD (Docking)", "Intra-run conformational coordinate divergence from search origin"),
            ("Matched Heavy Atoms", "Count of structurally mapped atoms used in RMSD computation"),
            ("Alignment / Mapping Method", "Isomorphic heavy-atom mapping (RDKit / Graph Isomorphism)"),
        ]

        for idx, (title, tooltip) in enumerate(val_fields):
            row = idx // 2
            col = (idx % 2) * 2

            t_lbl = QLabel(title.upper())
            t_lbl.setStyleSheet("color: #94a3b8; font-size: 10px; font-weight: bold;")
            t_lbl.setToolTip(tooltip)

            v_lbl = QLabel("—")
            v_lbl.setStyleSheet("font-size: 13px; font-weight: bold; color: #f8fafc;")

            v_grid.addWidget(t_lbl, row, col)
            v_grid.addWidget(v_lbl, row, col + 1)
            self.val_cells[title] = (t_lbl, v_lbl)

        layout.addWidget(val_card)

        # ── Pose-by-Pose Validation Table ─────────────────────────
        tbl_hdr = QLabel("POSE-BY-POSE RE-DOCKING RMSD VALIDATION:")
        tbl_hdr.setObjectName("subHeader")
        layout.addWidget(tbl_hdr)

        self.table_val_poses = QTableWidget()
        self.table_val_poses.setColumnCount(6)
        self.table_val_poses.setHorizontalHeaderLabels([
            "Pose Rank", "Binding ΔG (kcal/mol)", "Validation RMSD (Å)",
            "Shape / Kabsch RMSD (Å)", "Matched Atoms", "Quality Interpretation"
        ])
        self.table_val_poses.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table_val_poses.setAlternatingRowColors(True)
        self.table_val_poses.setEditTriggers(QTableWidget.NoEditTriggers)
        layout.addWidget(self.table_val_poses)

        # ── Validation Standards Note Card ────────────────────────
        std_card = QFrame()
        std_card.setObjectName("cardFrame")
        s_l = QVBoxLayout(std_card)
        s_l.setContentsMargins(14, 10, 14, 10)
        lbl_std = QLabel(
            "<b>Scientific Benchmark Guidelines:</b><br>"
            "• <b>RMSD ≤ 2.0 Å</b>: High-accuracy reproduction of crystallographic binding mode (standard docking benchmark threshold).<br>"
            "• <b>2.0 Å &lt; RMSD ≤ 3.0 Å</b>: Moderate accuracy; general pose topology captured but local contact orientation diverges.<br>"
            "• <b>RMSD &gt; 3.0 Å</b>: Inverted or alternate site binding mode.<br>"
            "<i>Note: Validation is strictly optional. Docking-search RMSD is always reported independently in Poses view.</i>"
        )
        lbl_std.setStyleSheet("color: #cbd5e1; font-size: 11px; line-height: 1.4;")
        lbl_std.setWordWrap(True)
        s_l.addWidget(lbl_std)
        layout.addWidget(std_card)

        self.tab_views.addTab(self.view_validation, "🔬 Validation")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 5 — INTERACTIONS
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_interactions(self) -> None:
        self.view_interactions = QWidget()
        layout = QVBoxLayout(self.view_interactions)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Pose Awareness Info Card & Navigation (Sections 8 & 9)
        self.card_inter_pose = QFrame()
        self.card_inter_pose.setObjectName("cardFrame")
        pose_info_l = QHBoxLayout(self.card_inter_pose)
        pose_info_l.setContentsMargins(12, 8, 12, 8)

        self.lbl_inter_header_info = QLabel("<b>Interaction Analysis:</b> Select a docking job and pose to view details.")
        self.lbl_inter_header_info.setObjectName("subHeader")
        pose_info_l.addWidget(self.lbl_inter_header_info)
        pose_info_l.addStretch()

        self.btn_prev_pose = QPushButton("◀ Prev Pose")
        self.btn_prev_pose.setEnabled(False)
        self.btn_prev_pose.setToolTip("Switch to previous ranked docking pose")
        self.btn_prev_pose.clicked.connect(lambda: self._select_relative_pose(-1))
        pose_info_l.addWidget(self.btn_prev_pose)

        self.lbl_pose_nav_indicator = QLabel("Pose 1 of 1")
        self.lbl_pose_nav_indicator.setStyleSheet("font-weight: bold; color: #38bdf8; padding: 0 6px;")
        pose_info_l.addWidget(self.lbl_pose_nav_indicator)

        self.btn_next_pose = QPushButton("Next Pose ▶")
        self.btn_next_pose.setEnabled(False)
        self.btn_next_pose.setToolTip("Switch to next ranked docking pose")
        self.btn_next_pose.clicked.connect(lambda: self._select_relative_pose(1))
        pose_info_l.addWidget(self.btn_next_pose)

        layout.addWidget(self.card_inter_pose)

        # 2D Interaction Diagram Box
        diag_box = QGroupBox("2D PROTEIN–LIGAND INTERACTION MAP")
        diag_l = QVBoxLayout(diag_box)
        diag_l.setContentsMargins(6, 6, 6, 6)

        diag_tb = QHBoxLayout()
        diag_tb.addWidget(QLabel("Theme:"))
        self.cmb_diag_theme = QComboBox()
        self.cmb_diag_theme.addItems(["🌙 Dark Studio", "📄 Publication (White)"])
        self.cmb_diag_theme.currentIndexChanged.connect(self._on_diag_theme_changed)
        diag_tb.addWidget(self.cmb_diag_theme)
        diag_tb.addStretch()

        self.btn_export_diagram = QPushButton("💾 Save Diagram (SVG/HTML/PNG)...")
        self.btn_export_diagram.clicked.connect(self._save_diagram_file)
        diag_tb.addWidget(self.btn_export_diagram)
        diag_l.addLayout(diag_tb)

        if HAS_WEBENGINE:
            self.web_view = QWebEngineView()
            self.web_view.setMinimumHeight(240)
            diag_l.addWidget(self.web_view)
        elif HAS_SVG:
            self.svg_view = QSvgWidget()
            self.svg_view.setMinimumHeight(240)
            diag_l.addWidget(self.svg_view)
        else:
            self.lbl_diag_fallback = QLabel("Select a docking pose to render 2D interaction diagram")
            self.lbl_diag_fallback.setAlignment(Qt.AlignCenter)
            self.lbl_diag_fallback.setMinimumHeight(200)
            diag_l.addWidget(self.lbl_diag_fallback)

        layout.addWidget(diag_box, 1)

        # Contacts Table
        inter_bar = QHBoxLayout()
        lbl_inter_title = QLabel("NON-COVALENT RESIDUE CONTACTS:")
        lbl_inter_title.setObjectName("subHeader")
        inter_bar.addWidget(lbl_inter_title)
        inter_bar.addStretch()

        self.cmb_filter_inter = QComboBox()
        self.cmb_filter_inter.addItems([
            "All Types", "Hydrogen Bond", "Hydrophobic Contact",
            "Salt Bridge", "Pi-Stacking", "Halogen Bond",
        ])
        self.cmb_filter_inter.currentTextChanged.connect(self._filter_interactions_table)
        inter_bar.addWidget(self.cmb_filter_inter)

        btn_copy_inter = QPushButton("📋 Copy Contacts")
        btn_copy_inter.clicked.connect(self._copy_contacts_to_clipboard)
        inter_bar.addWidget(btn_copy_inter)
        layout.addLayout(inter_bar)

        self.table_inter = QTableWidget()
        self.table_inter.setColumnCount(len(self.INTER_COLS))
        self.table_inter.setHorizontalHeaderLabels(self.INTER_COLS)
        self.table_inter.setAlternatingRowColors(True)
        self.table_inter.horizontalHeader().setStretchLastSection(True)
        self.table_inter.setSortingEnabled(True)
        layout.addWidget(self.table_inter, 1)

        self.tab_views.addTab(self.view_interactions, "⚡ Interactions")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 6 — THERMODYNAMICS
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_thermo(self) -> None:
        self.view_thermo = QWidget()
        layout = QVBoxLayout(self.view_thermo)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)

        hdr = QLabel("THERMODYNAMIC / STATISTICAL ANALYSIS")
        self.thermo_header = hdr
        hdr.setObjectName("subHeader")
        layout.addWidget(hdr)

        th_card = QFrame()
        th_card.setObjectName("cardFrame")
        th_grid = QGridLayout(th_card)
        th_grid.setContentsMargins(16, 14, 16, 14)
        th_grid.setHorizontalSpacing(24)
        th_grid.setVerticalSpacing(14)

        self.thermo_cells: Dict[str, Tuple[QLabel, QLabel]] = {}
        thermo_fields = [
            ("Information Entropy (S_inf)", "Conformational dispersion metric across GA run ensemble (Shannon entropy in nats/bits)"),
            ("Partition Function (Q)", "Canonical partition function summing Boltzmann weights: Q = Σ exp(-ΔG_i / k_B T)"),
            ("Statistical Temperature", "Effective ensemble temperature parameter for population weighting"),
            ("Statistical Free Energy (A)", "Helmholtz free energy of the ensemble: A = -k_B T ln(Q)"),
            ("Statistical Internal Energy (U)", "Thermal internal energy expectation value of docked ensemble"),
            ("Statistical Entropy (S_stat)", "Statistical mechanical entropy: S = (U - A) / T"),
            ("Boltzmann Population Fraction", "AutoDock4 DLG-derived ensemble descriptor; only shown when temperature is reported"),
            ("Analysis Engine", "BSNDVP™ AutoDock4 DLG statistical mechanics pipeline"),
        ]

        for idx, (title, tooltip) in enumerate(thermo_fields):
            row = idx // 2
            col = (idx % 2) * 2

            t_lbl = QLabel(title.upper())
            t_lbl.setStyleSheet("color: #94a3b8; font-size: 10px; font-weight: bold;")
            t_lbl.setToolTip(tooltip)

            v_lbl = QLabel("—")
            v_lbl.setStyleSheet("font-size: 13px; font-weight: bold; color: #f8fafc;")

            th_grid.addWidget(t_lbl, row, col)
            th_grid.addWidget(v_lbl, row, col + 1)
            self.thermo_cells[title] = (t_lbl, v_lbl)

        layout.addWidget(th_card)

        # Informative note
        n_card = QFrame()
        n_card.setObjectName("cardFrame")
        n_l = QVBoxLayout(n_card)
        n_l.setContentsMargins(14, 10, 14, 10)
        lbl_info = QLabel(
            "<b>Descriptor Attribution:</b><br>"
            "Thermodynamic and information entropy values are derived by AutoDockSuite Pro (BSNDVP™) "
            "from the complete AutoDock4 DLG conformational search history. For AutoDock Vina runs, "
            "partition function data is not applicable because Vina operates via Monte Carlo iterated local search."
        )
        self.thermo_info_label = lbl_info
        lbl_info.setStyleSheet("color: #cbd5e1; font-size: 11px; line-height: 1.4;")
        lbl_info.setWordWrap(True)
        n_l.addWidget(lbl_info)
        layout.addWidget(n_card)

        layout.addStretch()
        self.tab_views.addTab(self.view_thermo, "📈 Thermodynamics")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 7 — ADMET
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_admet(self) -> None:
        self.view_admet = QWidget()
        layout = QVBoxLayout(self.view_admet)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        lbl_hdr = QLabel("PHYSICOCHEMICAL & DRUG-LIKENESS PROFILE (ADMET / Lipinski Rule of 5):")
        lbl_hdr.setObjectName("subHeader")
        layout.addWidget(lbl_hdr)

        self.ADMET_COLS = ("Property", "Calculated Value", "Standard Threshold", "Interpretation", "Model / Source")
        self.table_admet = QTableWidget()
        self.table_admet.setColumnCount(len(self.ADMET_COLS))
        self.table_admet.setHorizontalHeaderLabels(self.ADMET_COLS)
        self.table_admet.setAlternatingRowColors(True)
        self.table_admet.horizontalHeader().setStretchLastSection(True)
        self.table_admet.setSelectionBehavior(QTableWidget.SelectRows)
        layout.addWidget(self.table_admet)

        self.tab_views.addTab(self.view_admet, "🧪 ADMET")

    # ═════════════════════════════════════════════════════════════════════════
    # VIEW 8 — PROVENANCE
    # ═════════════════════════════════════════════════════════════════════════
    def _build_view_provenance(self) -> None:
        self.view_provenance = QWidget()
        layout = QVBoxLayout(self.view_provenance)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        tb = QHBoxLayout()
        lbl_hdr = QLabel("EXECUTION PROVENANCE & REPRODUCIBILITY MANIFEST:")
        lbl_hdr.setObjectName("subHeader")
        tb.addWidget(lbl_hdr)
        tb.addStretch()

        self.btn_copy_provenance = QPushButton("📋 Copy Provenance Manifest (JSON)")
        self.btn_copy_provenance.clicked.connect(self._copy_provenance_json)
        tb.addWidget(self.btn_copy_provenance)
        layout.addLayout(tb)

        self.PROV_COLS = ("Attribute", "Value")
        self.table_provenance = QTableWidget()
        self.table_provenance.setColumnCount(len(self.PROV_COLS))
        self.table_provenance.setHorizontalHeaderLabels(self.PROV_COLS)
        self.table_provenance.setAlternatingRowColors(True)
        self.table_provenance.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table_provenance)

        self.tab_views.addTab(self.view_provenance, "📜 Provenance")

    # ─────────────────────────────────────────────────────────
    # Data Population & Selection
    # ─────────────────────────────────────────────────────────

    def _get_job_best_energy(self, job: "DockingJob") -> Optional[float]:
        if job.canonical_result and job.canonical_result.best_binding_energy is not None:
            return float(job.canonical_result.best_binding_energy)
        if job.vina_results:
            return float(job.vina_results[0].binding_affinity)
        ad4 = getattr(job, "ad4_results", {})
        if isinstance(ad4, dict):
            if "best_energy" in ad4:
                return float(ad4["best_energy"])
            if "poses" in ad4 and ad4["poses"]:
                p0 = ad4["poses"][0]
                return float(p0.get("binding_affinity", p0.get("binding_energy", p0.get("energy", 0.0))))
        dlg_p = getattr(job, "dlg_path", None) or (ad4.get("dlg_path") if isinstance(ad4, dict) else None)
        if dlg_p and Path(dlg_p).is_file():
            try:
                from dlg_extract import DLGParser
                res = DLGParser().parse(Path(dlg_p))
                if res and res.poses:
                    return float(res.poses[0].binding_energy)
            except Exception:
                pass
        return None

    def refresh(self) -> None:
        """Scan results directory and load jobs."""
        from job_manager import load_docking_jobs
        res_dir = getattr(self.config, "result_directory", Path("results"))
        try:
            self._jobs = load_docking_jobs(res_dir)
        except Exception:
            self._jobs = []

        success = 0
        failed = 0
        best_affinity = None
        best_lead = ""

        for job in self._jobs:
            st = job.status.value if hasattr(job.status, "value") else str(job.status)
            dg = self._get_job_best_energy(job)
            if dg is not None:
                if best_affinity is None or dg < best_affinity:
                    best_affinity = dg
                    best_lead = f"{job.ligand_name} ({dg:.2f} kcal/mol)"

            if "SUCCESS" in st.upper():
                success += 1
            elif "FAIL" in st.upper() or "ERROR" in st.upper():
                failed += 1

        self.lbl_kpi_total.setText(str(len(self._jobs)))
        self.lbl_kpi_success.setText(str(success))
        self.lbl_kpi_failed.setText(str(failed))
        self.lbl_kpi_lead.setText(best_lead if best_lead else "No completed jobs yet")

        self._apply_jobs_filter()

    def _apply_jobs_filter(self) -> None:
        q = self.ent_search_jobs.text().strip().lower()
        filtered = []
        for job in self._jobs:
            st = job.status.value if hasattr(job.status, "value") else str(job.status)
            eng = job.engine.value if hasattr(job.engine, "value") else str(job.engine)
            txt = f"{job.job_id} {job.receptor_name} {job.ligand_name} {eng} {st}".lower()
            if not q or q in txt:
                filtered.append(job)

        self._filtered_jobs = filtered
        self.table_jobs.setSortingEnabled(False)
        self.table_jobs.setRowCount(len(filtered))

        for r, job in enumerate(filtered):
            st = job.status.value if hasattr(job.status, "value") else str(job.status)
            dg_val = self._get_job_best_energy(job)
            if dg_val is not None:
                best_dg_val = dg_val
                best_dg_str = f"{dg_val:.2f}"
            else:
                best_dg_val = 999.0
                best_dg_str = "—"

            mode_str = job.docking_mode.value if hasattr(job.docking_mode, "value") else str(job.docking_mode)
            eng_str = job.engine.value if hasattr(job.engine, "value") else str(job.engine)

            it_id = QTableWidgetItem(job.job_id)
            it_id.setData(Qt.UserRole, job)
            it_rec = QTableWidgetItem(job.receptor_name)
            it_lig = QTableWidgetItem(job.ligand_name)
            it_eng = QTableWidgetItem(eng_str)
            it_eng.setTextAlignment(Qt.AlignCenter)
            it_mode = QTableWidgetItem(mode_str)
            it_mode.setTextAlignment(Qt.AlignCenter)

            it_st = QTableWidgetItem(st)
            it_st.setTextAlignment(Qt.AlignCenter)
            it_st.setForeground(QColor("#34d399") if "SUCCESS" in st.upper() else QColor("#f87171"))

            it_dg = _NumericTableItem(best_dg_str, best_dg_val)
            it_dg.setTextAlignment(Qt.AlignCenter)
            it_dg.setForeground(QColor("#38bdf8"))

            row_items = [it_id, it_rec, it_lig, it_eng, it_mode, it_st, it_dg]
            for c, it in enumerate(row_items):
                self.table_jobs.setItem(r, c, it)

        self.table_jobs.setSortingEnabled(True)
        if filtered and not self.table_jobs.selectedItems():
            self.table_jobs.selectRow(0)

    def _on_job_selected(self) -> None:
        sel = self.table_jobs.selectedItems()
        if not sel:
            return
        row = sel[0].row()
        item_0 = self.table_jobs.item(row, 0)
        job: Optional["DockingJob"] = item_0.data(Qt.UserRole) if item_0 else None
        if not job:
            return

        self._current_job = job
        # Ensure canonical scientific result model is populated without redundant parsing
        job.ensure_canonical_result()

        # Update all 8 Scientific Views
        self._update_overview_view(job)
        self._populate_poses_view(job)
        self._update_clusters_view(job)
        self._update_validation_view(job)
        self._update_thermo_view(job)
        self._update_admet_view(job)
        self._update_provenance_view(job)

        # Trigger interaction analysis for top pose
        self._current_pose = 1
        self._load_pose_interactions(job, 1)

    # ─────────────────────────────────────────────────────────
    # VIEW UPDATERS & SCIENTIFIC CHARTS
    # ─────────────────────────────────────────────────────────

    def _generate_poses_chart_svg(self, poses: List["CanonicalPose"], engine: str) -> str:
        """Generate a clean vector SVG bar chart of ranked pose binding affinities."""
        if not poses:
            return ""
        top_poses = poses[:10]
        w = 640
        h = 160
        pad_l = 60
        pad_r = 30
        pad_t = 30
        pad_b = 30

        energies = [p.binding_energy for p in top_poses if p.binding_energy is not None]
        if not energies:
            return ""

        min_e = min(energies)
        span = abs(min_e) if abs(min_e) > 0.01 else 1.0

        bar_w = max(18, int((w - pad_l - pad_r) / (len(top_poses) * 1.5)))
        spacing = int(bar_w * 1.5)

        svg = [
            f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" style="background:transparent; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">',
            f'<text x="{pad_l}" y="18" font-size="11" font-weight="700" fill="#94a3b8">Conformation Energy Spectrum (kcal/mol) — {engine}</text>',
            f'<line x1="{pad_l}" y1="{h - pad_b}" x2="{w - pad_r}" y2="{h - pad_b}" stroke="#475569" stroke-width="1" />',
        ]

        max_bar_h = h - pad_t - pad_b
        for i, p in enumerate(top_poses):
            e = p.binding_energy if p.binding_energy is not None else 0.0
            x = pad_l + i * spacing + 10
            bar_h = max(6, int((abs(e) / span) * max_bar_h))
            y = (h - pad_b) - bar_h

            fill = "#0284c7" if i == 0 else "#38bdf8"
            svg.append(f'<rect x="{x}" y="{y}" width="{bar_w}" height="{bar_h}" rx="3" fill="{fill}" opacity="0.88" />')
            svg.append(f'<text x="{x + bar_w // 2}" y="{h - pad_b + 14}" font-size="9" fill="#94a3b8" text-anchor="middle">P{p.rank}</text>')
            svg.append(f'<text x="{x + bar_w // 2}" y="{y - 4}" font-size="9" font-weight="600" fill="#f8fafc" text-anchor="middle">{e:.1f}</text>')

        svg.append('</svg>')
        return "\n".join(svg)

    def _generate_clusters_chart_svg(self, clusters: List["ClusterInfo"]) -> str:
        """Generate a clean vector SVG bar chart of cluster population percentages."""
        if not clusters:
            return ""
        top_c = clusters[:8]
        w = 640
        h = max(130, 30 + len(top_c) * 18)
        bar_x = 90
        max_w = 460

        svg = [
            f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" style="background:transparent; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">',
            f'<text x="20" y="18" font-size="11" font-weight="700" fill="#94a3b8">Conformational Cluster Population Distribution (%)</text>',
        ]

        for i, c in enumerate(top_c):
            y = 30 + i * 18
            pop = c.population_percent if c.population_percent is not None else 0.0
            bw = max(6, int((pop / 100.0) * max_w))
            fill = "#10b981" if i == 0 else "#0ea5e9"

            svg.append(f'<text x="{bar_x - 10}" y="{y + 11}" font-size="10" fill="#94a3b8" text-anchor="end">Cluster {c.cluster_id}</text>')
            svg.append(f'<rect x="{bar_x}" y="{y}" width="{bw}" height="13" rx="2" fill="{fill}" opacity="0.85" />')
            svg.append(f'<text x="{bar_x + bw + 8}" y="{y + 11}" font-size="10" font-weight="600" fill="#f8fafc">{pop:.1f}% ({c.size} runs, {c.lowest_energy:.2f} kcal/mol)</text>')

        svg.append('</svg>')
        return "\n".join(svg)

    def _update_overview_view(self, job: "DockingJob") -> None:
        from models import Engine
        res = job.canonical_result
        is_ad4 = (job.engine == Engine.AUTODOCK4)

        self.lbl_ov_job_title.setText(f"{job.receptor_name}  ×  {job.ligand_name}")
        self.lbl_ov_engine_badge.setText(f"Engine: {job.engine.value}")
        self.lbl_ov_engine_badge.setStyleSheet(
            "background: #7c3aed; color: white; padding: 4px 10px; border-radius: 4px; font-weight: bold; font-size: 11px;"
            if is_ad4 else
            "background: #0284c7; color: white; padding: 4px 10px; border-radius: 4px; font-weight: bold; font-size: 11px;"
        )
        self.lbl_ov_mode_badge.setText(f"Mode: {job.docking_mode.value}")

        # Engine-specific score semantics; Vina affinity is not ΔG and has no Ki.
        best_dg = res.best_binding_energy if res else self._get_job_best_energy(job)
        dg_str = f"{best_dg:.2f} kcal/mol" if best_dg is not None else "N/A — not calculated"

        ki_str = "DLG-reported only" if is_ad4 else "Not applicable — Vina does not report Ki"

        # Update Grid Cells
        if is_ad4:
            total_runs = res.summary.get("total_runs", getattr(self.config, "ga_run", 100)) if res else getattr(self.config, "ga_run", 100)
            poses_count = len(res.poses) if (res and res.poses) else job.obtained_modes or 0
            clusters_count = len(res.clusters) if (res and res.clusters) else 0

            self.ov_cells["Runs / Modes"][0].setText("GA RUNS:")
            self.ov_cells["Runs / Modes"][1].setText(f"{total_runs} runs")

            self.ov_cells["Docked Poses"][0].setText("DOCKED CONFORMATIONS:")
            self.ov_cells["Docked Poses"][1].setText(f"{poses_count} poses")

            self.ov_cells["Clusters"][0].setText("RMSD CLUSTERS:")
            self.ov_cells["Clusters"][1].setText(f"{clusters_count} clusters" if clusters_count > 0 else "N/A — single run / unclustered")

            # Dominant cluster
            if res and res.clusters:
                c0 = res.clusters[0]
                pop_pct = c0.population_percent if c0.population_percent is not None else (c0.size / total_runs * 100.0 if total_runs else 0)
                self.ov_cells["Dominant Cluster"][1].setText(f"{pop_pct:.1f}% ({c0.size} runs)")
                self.ov_cells["Best Cluster RMSD"][1].setText(f"{c0.cluster_rmsd:.2f} Å" if c0.cluster_rmsd is not None else "N/A — not reported")
            else:
                self.ov_cells["Dominant Cluster"][1].setText("N/A — unclustered")
                self.ov_cells["Best Cluster RMSD"][1].setText("N/A — unclustered")

            self.lbl_ov_notes.setText(
                "<b>AutoDock 4.2.6 Analysis:</b> Genetic Algorithm (LGA) conformational sampling with local search. "
                "Docking poses are clustered by root-mean-square coordinate tolerance. "
                "AutoDock4 computes explicit Coulomb electrostatic interactions from Gasteiger partial charges "
                "(which explains differences from Vina's purely empirical non-electrostatic scoring)."
            )
        else:
            # Vina
            modes_count = len(res.poses) if (res and res.poses) else job.obtained_modes or (len(job.vina_results) if job.vina_results else 0)
            self.ov_cells["Runs / Modes"][0].setText("SAMPLING MODES:")
            self.ov_cells["Runs / Modes"][1].setText(f"{modes_count} modes")

            self.ov_cells["Docked Poses"][0].setText("DOCKED POSES:")
            self.ov_cells["Docked Poses"][1].setText(f"{modes_count} poses")

            self.ov_cells["Clusters"][0].setText("CLUSTERS:")
            self.ov_cells["Clusters"][1].setText("N/A — not reported by this engine")

            self.ov_cells["Dominant Cluster"][1].setText("N/A — not applicable (Vina)")
            self.ov_cells["Best Cluster RMSD"][1].setText("N/A — not applicable (Vina)")

            self.lbl_ov_notes.setText(
                "<b>AutoDock Vina Analysis:</b> Iterated local search with Broyden-Fletcher-Goldfarb-Shanno (BFGS) optimization. "
                "Conformations are ranked strictly by empirical binding affinity. "
                "RMSD lower/upper bounds represent geometric coordinate distance between predicted modes."
            )

        self.ov_cells["Engine Score / Affinity"][1].setText(dg_str)
        # A Vina result has no Ki field.  Hide the AutoDock4-only metric
        # rather than showing an engine-specific field in a Vina view.
        ki_title, ki_value = self.ov_cells["AutoDock4 Ki"]
        ki_title.setVisible(is_ad4)
        ki_value.setVisible(is_ad4)
        if is_ad4:
            ki_value.setText(ki_str)

        # Validation RMSD
        val_rmsd = res.validation.crystal_rmsd if (res and res.validation and res.validation.validation_performed) else None
        if val_rmsd is not None:
            self.ov_cells["Validation RMSD"][1].setText(f"{val_rmsd:.2f} Å")
            self.ov_cells["Validation RMSD"][1].setStyleSheet("font-size: 13px; font-weight: bold; color: #10b981;" if val_rmsd <= 2.0 else "font-size: 13px; font-weight: bold; color: #f59e0b;")
        else:
            self.ov_cells["Validation RMSD"][1].setText("N/A — validation not performed")
            self.ov_cells["Validation RMSD"][1].setStyleSheet("font-size: 13px; font-weight: bold; color: #94a3b8;")

        # Update SVG energy chart
        if HAS_SVG and hasattr(self, "ov_chart_svg"):
            poses = res.poses if res else []
            if poses:
                chart_svg = self._generate_poses_chart_svg(poses, job.engine.value)
                self.ov_chart_svg.load(QByteArray(chart_svg.encode("utf-8")))
                self.ov_chart_card.show()
            else:
                self.ov_chart_card.hide()

    def _populate_poses_view(self, job: "DockingJob") -> None:
        from models import Engine
        res = job.canonical_result
        is_ad4 = (job.engine == Engine.AUTODOCK4)
        is_advanced = self.chk_show_advanced.isChecked()

        if is_ad4:
            if is_advanced:
                headers = [
                    "Rank", "GA Run", "Estimated Free Energy of Binding (kcal/mol)",
                    "Estimated Inhibition Constant, Ki", "Cluster ID",
                    "Cluster Size", "Cluster RMSD (Å)", "Intermol ΔG", "Internal ΔG",
                    "Torsional ΔG", "Unbound ΔG", "Ref RMSD (Å)", "Validation RMSD (Å)"
                ]
            else:
                headers = [
                    "Rank", "GA Run", "Estimated Free Energy of Binding (kcal/mol)",
                    "Estimated Inhibition Constant, Ki", "Cluster ID",
                    "Cluster Size", "Cluster RMSD (Å)", "Validation RMSD (Å)"
                ]
        else:
            headers = [
                "Rank", "Affinity (kcal/mol)", "Ki (not applicable)",
                "RMSD l.b. (Å)", "RMSD u.b. (Å)", "Validation RMSD (Å)"
            ]

        self.table_poses.setSortingEnabled(False)
        self.table_poses.setColumnCount(len(headers))
        self.table_poses.setHorizontalHeaderLabels(headers)

        poses: List[CanonicalPose] = res.poses if res else []
        self.table_poses.setRowCount(len(poses))

        for row, p in enumerate(poses):
            ki_str = p.estimated_ki_formatted if is_ad4 else "Not applicable — Vina does not report Ki"

            if is_ad4:
                ad4 = p.ad4_metrics
                it_rank = _NumericTableItem(str(p.rank), float(p.rank))
                it_rank.setData(Qt.UserRole, p)
                it_run = _NumericTableItem(str(p.run_number or p.rank), float(p.run_number or p.rank))
                it_dg = _NumericTableItem(f"{p.binding_energy:.2f}" if p.binding_energy is not None else "N/A", p.binding_energy or 0.0)
                it_dg.setForeground(QColor("#38bdf8"))
                it_ki = QTableWidgetItem(ki_str)

                c_id_str = str(ad4.cluster_id) if (ad4 and ad4.cluster_id is not None) else "N/A"
                c_id_val = float(ad4.cluster_id) if (ad4 and ad4.cluster_id is not None) else 999.0
                it_cid = _NumericTableItem(c_id_str, c_id_val)

                c_sz_str = str(ad4.cluster_size) if (ad4 and ad4.cluster_size is not None) else "N/A"
                c_sz_val = float(ad4.cluster_size) if (ad4 and ad4.cluster_size is not None) else 0.0
                it_csz = _NumericTableItem(c_sz_str, c_sz_val)

                c_rmsd_str = f"{ad4.cluster_rmsd:.2f}" if (ad4 and ad4.cluster_rmsd is not None) else "N/A"
                c_rmsd_val = float(ad4.cluster_rmsd) if (ad4 and ad4.cluster_rmsd is not None) else 999.0
                it_crmsd = _NumericTableItem(c_rmsd_str, c_rmsd_val)

                val_str = f"{p.validation_rmsd:.2f}" if p.validation_rmsd is not None else "N/A — not performed"
                val_num = p.validation_rmsd if p.validation_rmsd is not None else 999.0
                it_val = _NumericTableItem(val_str, val_num)

                if is_advanced:
                    inter_str = f"{ad4.intermolecular_energy:.2f}" if (ad4 and ad4.intermolecular_energy is not None) else "N/A"
                    it_inter = _NumericTableItem(inter_str, ad4.intermolecular_energy if ad4 and ad4.intermolecular_energy is not None else 0.0)

                    int_str = f"{ad4.internal_energy:.2f}" if (ad4 and ad4.internal_energy is not None) else "N/A"
                    it_int = _NumericTableItem(int_str, ad4.internal_energy if ad4 and ad4.internal_energy is not None else 0.0)

                    tor_str = f"{ad4.torsional_energy:.2f}" if (ad4 and ad4.torsional_energy is not None) else "N/A"
                    it_tor = _NumericTableItem(tor_str, ad4.torsional_energy if ad4 and ad4.torsional_energy is not None else 0.0)

                    unb_str = f"{ad4.unbound_energy:.2f}" if (ad4 and ad4.unbound_energy is not None) else "N/A"
                    it_unb = _NumericTableItem(unb_str, ad4.unbound_energy if ad4 and ad4.unbound_energy is not None else 0.0)

                    ref_rmsd_str = f"{ad4.rmsd_from_reference:.2f}" if (ad4 and ad4.rmsd_from_reference is not None) else "N/A — not reported"
                    it_ref = _NumericTableItem(ref_rmsd_str, ad4.rmsd_from_reference if (ad4 and ad4.rmsd_from_reference is not None) else 999.0)

                    row_items = [it_rank, it_run, it_dg, it_ki, it_cid, it_csz, it_crmsd, it_inter, it_int, it_tor, it_unb, it_ref, it_val]
                else:
                    row_items = [it_rank, it_run, it_dg, it_ki, it_cid, it_csz, it_crmsd, it_val]

            else:
                # Vina
                vm = p.vina_metrics
                it_rank = _NumericTableItem(str(p.rank), float(p.rank))
                it_rank.setData(Qt.UserRole, p)
                it_dg = _NumericTableItem(f"{p.binding_energy:.2f}" if p.binding_energy is not None else "N/A", p.binding_energy or 0.0)
                it_dg.setForeground(QColor("#38bdf8"))
                it_ki = QTableWidgetItem(ki_str)

                lb_str = f"{vm.rmsd_lower_bound:.2f}" if (vm and vm.rmsd_lower_bound is not None) else "0.00"
                lb_val = vm.rmsd_lower_bound if (vm and vm.rmsd_lower_bound is not None) else 0.0
                it_lb = _NumericTableItem(lb_str, lb_val)

                ub_str = f"{vm.rmsd_upper_bound:.2f}" if (vm and vm.rmsd_upper_bound is not None) else "0.00"
                ub_val = vm.rmsd_upper_bound if (vm and vm.rmsd_upper_bound is not None) else 0.0
                it_ub = _NumericTableItem(ub_str, ub_val)

                val_str = f"{p.validation_rmsd:.2f}" if p.validation_rmsd is not None else "N/A — not performed"
                it_val = _NumericTableItem(val_str, p.validation_rmsd if p.validation_rmsd is not None else 999.0)

                row_items = [it_rank, it_dg, it_ki, it_lb, it_ub, it_val]

            for c, it in enumerate(row_items):
                it.setTextAlignment(Qt.AlignCenter)
                self.table_poses.setItem(row, c, it)

        self.table_poses.setSortingEnabled(True)
        if poses:
            self.table_poses.selectRow(0)

    def _on_toggle_advanced_poses(self, checked: bool) -> None:
        if self._current_job:
            self._populate_poses_view(self._current_job)

    def _on_pose_selected(self) -> None:
        sel = self.table_poses.selectedItems()
        if not sel or not self._current_job:
            return
        row = sel[0].row()
        item_rank = self.table_poses.item(row, 0)
        pose_obj: Optional["CanonicalPose"] = item_rank.data(Qt.UserRole) if item_rank else None
        if pose_obj is not None:
            pose_idx = pose_obj.rank
            run_num = pose_obj.ad4_metrics.run_number if (pose_obj.ad4_metrics and pose_obj.ad4_metrics.run_number is not None) else None
        else:
            try:
                pose_idx = int(item_rank.text()) if item_rank else row + 1
            except ValueError:
                pose_idx = row + 1
            run_num = None

        self._current_pose = pose_idx
        self._load_pose_interactions(self._current_job, pose_idx, run_number=run_num)

    # ─────────────────────────────────────────────────────────
    # CLUSTERS VIEW
    # ─────────────────────────────────────────────────────────

    def _update_clusters_view(self, job: "DockingJob") -> None:
        from models import Engine
        res = job.canonical_result
        is_ad4 = (job.engine == Engine.AUTODOCK4)

        self.table_clusters.setSortingEnabled(False)
        self.table_clusters.setRowCount(0)
        self.btn_inspect_cluster_pose.setEnabled(False)

        if not is_ad4:
            if hasattr(self, "cluster_chart_card"):
                self.cluster_chart_card.hide()
            self.lbl_cluster_detail.setText(
                "ℹ <b>Clustering is an AutoDock 4 feature:</b> AutoDock 4 groups conformational states across "
                "Lamarckian Genetic Algorithm runs based on RMS tolerance.<br>"
                "AutoDock Vina conducts Monte Carlo search and produces independent binding modes without cluster partitioning."
            )
            return

        clusters: List[ClusterInfo] = res.clusters if (res and res.clusters) else []
        if not clusters:
            if hasattr(self, "cluster_chart_card"):
                self.cluster_chart_card.hide()
            self.lbl_cluster_detail.setText("No clustering data found for this AutoDock 4 run (single run or DLG clustering section missing).")
            return

        # Update SVG cluster distribution chart
        if HAS_SVG and hasattr(self, "cluster_chart_svg"):
            c_svg = self._generate_clusters_chart_svg(clusters)
            self.cluster_chart_svg.load(QByteArray(c_svg.encode("utf-8")))
            self.cluster_chart_card.show()

        self.table_clusters.setRowCount(len(clusters))
        total_runs = res.summary.get("total_runs", 100) if res else 100

        for r, c in enumerate(clusters):
            pop_pct = c.population_percent if c.population_percent is not None else (c.size / total_runs * 100.0 if total_runs else 0)

            it_cid = _NumericTableItem(str(c.cluster_id), float(c.cluster_id))
            it_cid.setData(Qt.UserRole, c)
            it_sz = _NumericTableItem(str(c.size), float(c.size))
            it_pop = _NumericTableItem(f"{pop_pct:.1f}%", float(pop_pct))

            low_e_str = f"{c.lowest_energy:.2f}" if c.lowest_energy is not None else "N/A"
            it_low = _NumericTableItem(low_e_str, c.lowest_energy or 0.0)
            it_low.setForeground(QColor("#38bdf8"))

            mean_e_str = f"{c.mean_energy:.2f}" if c.mean_energy is not None else "N/A"
            it_mean = _NumericTableItem(mean_e_str, c.mean_energy or 0.0)

            rmsd_str = f"{c.cluster_rmsd:.2f}" if c.cluster_rmsd is not None else "N/A"
            it_rmsd = _NumericTableItem(rmsd_str, c.cluster_rmsd or 0.0)

            rep_run_str = str(c.representative_run) if c.representative_run is not None else "N/A"
            it_reprun = _NumericTableItem(rep_run_str, float(c.representative_run or 0))

            rep_pose_str = str(c.representative_pose) if c.representative_pose is not None else str(r + 1)
            it_reppose = _NumericTableItem(rep_pose_str, float(c.representative_pose or (r + 1)))

            items = [it_cid, it_sz, it_pop, it_low, it_mean, it_rmsd, it_reprun, it_reppose]
            for col, item in enumerate(items):
                item.setTextAlignment(Qt.AlignCenter)
                self.table_clusters.setItem(r, col, item)

        self.table_clusters.setSortingEnabled(True)
        if clusters:
            self.table_clusters.selectRow(0)

    def _on_cluster_selected(self) -> None:
        sel = self.table_clusters.selectedItems()
        if not sel:
            self.btn_inspect_cluster_pose.setEnabled(False)
            return
        row = sel[0].row()
        item_0 = self.table_clusters.item(row, 0)
        c: Optional[ClusterInfo] = item_0.data(Qt.UserRole) if item_0 else None
        if not c:
            return

        runs_str = ", ".join(str(r) for r in c.runs[:25]) if c.runs else "N/A"
        if c.runs and len(c.runs) > 25:
            runs_str += f" ... (+{len(c.runs) - 25} more)"

        detail = (
            f"<b>Cluster ID {c.cluster_id}:</b> {c.size} member conformations ({c.population_percent or 0:.1f}%).<br>"
            f"<b>Lowest Binding Energy:</b> {c.lowest_energy or 0:.2f} kcal/mol  |  "
            f"<b>Mean Energy:</b> {c.mean_energy or 0:.2f} kcal/mol  |  "
            f"<b>Cluster RMSD:</b> {c.cluster_rmsd or 0:.2f} Å<br>"
            f"<b>Representative Run:</b> #{c.representative_run}  |  <b>Member Runs:</b> [{runs_str}]"
        )
        self.lbl_cluster_detail.setText(detail)
        self.btn_inspect_cluster_pose.setEnabled(True)

    def _inspect_cluster_rep_pose(self) -> None:
        sel = self.table_clusters.selectedItems()
        if not sel or not self._current_job:
            return
        row = sel[0].row()
        item_0 = self.table_clusters.item(row, 0)
        c: Optional[ClusterInfo] = item_0.data(Qt.UserRole) if item_0 else None
        if not c:
            return

        target_pose = c.representative_pose or (row + 1)
        # Switch to Poses tab
        self.tab_views.setCurrentWidget(self.view_poses)
        # Select row matching target_pose
        for r in range(self.table_poses.rowCount()):
            it = self.table_poses.item(r, 0)
            if it and it.text() == str(target_pose):
                self.table_poses.selectRow(r)
                break

    # ─────────────────────────────────────────────────────────
    # VALIDATION VIEW (ON-DEMAND & OPTIONAL)
    # ─────────────────────────────────────────────────────────

    def _browse_validation_reference(self) -> None:
        init_dir = str(getattr(self.config, "ligand_directory", Path.cwd()))
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Select Crystallographic Reference Ligand Structure",
            init_dir,
            "Ligand Structures (*.sdf *.mol2 *.pdb *.pdbqt *.mol *.cif *.mmcif);;All Files (*.*)"
        )
        if f:
            self.txt_val_ref.setText(f)

    def _clear_validation(self) -> None:
        self.txt_val_ref.clear()
        if self._current_job:
            self._update_validation_view(self._current_job)

    def _run_on_demand_validation(self) -> None:
        job = self._current_job
        if not job:
            QMessageBox.warning(self, "No Job Selected", "Please select a docking job from the table first.")
            return

        ref_str = self.txt_val_ref.text().strip()
        if not ref_str:
            QMessageBox.warning(
                self, "No Reference Ligand",
                "Please click 'Browse Reference...' to select an experimental crystallographic reference ligand file (SDF, MOL2, PDB, PDBQT, CIF)."
            )
            return

        ref_path = Path(ref_str)
        if not ref_path.is_file():
            QMessageBox.warning(self, "File Not Found", f"The reference ligand file was not found:\n{ref_path}")
            return

        pose_file = self._resolve_ligand_path(job)
        if not pose_file or not Path(pose_file).is_file():
            QMessageBox.warning(self, "Pose File Missing", f"Could not locate the docked pose output file for job {job.job_id}.")
            return

        from validators import compute_pose_validation_rmsd
        from models import ValidationResult

        poses = job.canonical_result.poses if (job.canonical_result and job.canonical_result.poses) else []
        n_poses = len(poses) or job.obtained_modes or 1

        self.table_val_poses.setRowCount(0)
        best_rmsd = 999.0
        best_pose_idx = 1
        matched_heavy = 0

        self.table_val_poses.setSortingEnabled(False)
        self.table_val_poses.setRowCount(n_poses)

        for i in range(1, n_poses + 1):
            p = poses[i - 1] if i - 1 < len(poses) else None
            dg_val = p.binding_energy if (p and p.binding_energy is not None) else None
            dg_str = f"{dg_val:.2f}" if dg_val is not None else "N/A"

            try:
                pos_rmsd, kab_rmsd, matched = compute_pose_validation_rmsd(ref_path, pose_file, pose_index=i)
                matched_heavy = max(matched_heavy, matched)
                if pos_rmsd < best_rmsd:
                    best_rmsd = pos_rmsd
                    best_pose_idx = i

                if p:
                    p.validation_rmsd = pos_rmsd

                if pos_rmsd <= 2.0:
                    badge_str = "✔ High Accuracy (≤ 2.0 Å Standard)"
                    color = "#10b981"
                elif pos_rmsd <= 3.0:
                    badge_str = "⚠ Moderate Accuracy (2.0–3.0 Å)"
                    color = "#f59e0b"
                else:
                    badge_str = "✖ Divergent / Alternate Pose (> 3.0 Å)"
                    color = "#ef4444"

                it_rank = _NumericTableItem(f"Pose {i}", float(i))
                it_rank.setTextAlignment(Qt.AlignCenter)
                it_dg = _NumericTableItem(dg_str, dg_val if dg_val is not None else 0.0)
                it_dg.setTextAlignment(Qt.AlignCenter)
                it_pos = _NumericTableItem(f"{pos_rmsd:.2f} Å", pos_rmsd)
                it_pos.setTextAlignment(Qt.AlignCenter)
                it_pos.setForeground(QColor(color))
                it_kab = _NumericTableItem(f"{kab_rmsd:.2f} Å", kab_rmsd)
                it_kab.setTextAlignment(Qt.AlignCenter)
                it_m = _NumericTableItem(f"{matched} atoms", float(matched))
                it_m.setTextAlignment(Qt.AlignCenter)
                it_q = QTableWidgetItem(badge_str)
                it_q.setForeground(QColor(color))

                for col, item in enumerate([it_rank, it_dg, it_pos, it_kab, it_m, it_q]):
                    self.table_val_poses.setItem(i - 1, col, item)

            except Exception as err:
                it_rank = _NumericTableItem(f"Pose {i}", float(i))
                it_err = QTableWidgetItem(f"Failed: {err}")
                it_err.setForeground(QColor("#ef4444"))
                self.table_val_poses.setItem(i - 1, 0, it_rank)
                self.table_val_poses.setItem(i - 1, 5, it_err)

        self.table_val_poses.setSortingEnabled(True)

        if best_rmsd < 900.0:
            if job.canonical_result:
                job.canonical_result.validation = ValidationResult(
                    validation_performed=True,
                    reference_ligand_file=str(ref_path.name),
                    crystal_rmsd=best_rmsd,
                    matched_atoms=matched_heavy,
                    atom_mapping_method="RDKit MCS Substructure Isomorphism",
                    docking_search_rmsd=getattr(poses[0], "validation_rmsd", None) if poses else None,
                )

            self.val_cells["Reference Ligand Source"][1].setText(ref_path.name)
            self.val_cells["Validation Status"][1].setText(f"✔ Validated against {ref_path.name}")
            self.val_cells["Validation Status"][1].setStyleSheet("color: #10b981; font-weight: bold;")
            self.val_cells["Validation RMSD"][1].setText(f"{best_rmsd:.2f} Å (Best: Pose #{best_pose_idx})")
            self.val_cells["Search RMSD (Docking)"][1].setText("See Poses view for intra-conformation RMSD")
            self.val_cells["Matched Heavy Atoms"][1].setText(f"{matched_heavy} atoms")
            self.val_cells["Alignment / Mapping Method"][1].setText("MCS Substructure Graph Isomorphism")

            self._populate_poses_view(job)
            self._update_overview_view(job)

            QMessageBox.information(
                self,
                "RMSD Validation Complete",
                f"Successfully computed re-docking RMSD against crystallographic reference:\n\n"
                f"Reference: {ref_path.name}\n"
                f"Best Pose RMSD: {best_rmsd:.2f} Å (Pose #{best_pose_idx})\n"
                f"Matched Heavy Atoms: {matched_heavy}\n\n"
                f"Pose table and Overview have been updated."
            )
        else:
            QMessageBox.warning(self, "Validation Incomplete", "Could not compute RMSD for the poses.")

    def _update_validation_view(self, job: "DockingJob") -> None:
        res = job.canonical_result
        val = res.validation if res else None

        if hasattr(self, "lbl_val_target"):
            self.lbl_val_target.setText(f"Target Complex: {job.receptor_name}  ×  {job.ligand_name} ({job.engine.value})")

        if val and val.validation_performed:
            ref_name = val.reference_ligand_file or "Supplied Reference"
            if hasattr(self, "lbl_val_opt_status"):
                self.lbl_val_opt_status.setText(f"Validated against: {ref_name}")
                self.lbl_val_opt_status.setStyleSheet("color: #10b981; font-weight: bold; font-size: 11px;")

            self.val_cells["Reference Ligand Source"][1].setText(ref_name)
            self.val_cells["Validation Status"][1].setText("✔ Validated against Crystal Structure")
            self.val_cells["Validation Status"][1].setStyleSheet("color: #10b981; font-weight: bold;")
            self.val_cells["Validation RMSD"][1].setText(f"{val.crystal_rmsd:.2f} Å" if val.crystal_rmsd is not None else "N/A — not calculated")
            self.val_cells["Search RMSD (Docking)"][1].setText(f"{val.docking_search_rmsd:.2f} Å" if val.docking_search_rmsd is not None else "N/A — not reported")
            self.val_cells["Matched Heavy Atoms"][1].setText(str(val.matched_atoms) if val.matched_atoms else "All Heavy Atoms")
            self.val_cells["Alignment / Mapping Method"][1].setText(val.atom_mapping_method or "RDKit Graph Invariant Isomorphism")

            if hasattr(self, "table_val_poses") and res and res.poses:
                self.table_val_poses.setRowCount(len(res.poses))
                self.table_val_poses.setSortingEnabled(False)
                for i, p in enumerate(res.poses):
                    it_r = _NumericTableItem(f"Pose {p.rank}", float(p.rank))
                    it_r.setTextAlignment(Qt.AlignCenter)
                    dg_val = p.binding_energy if p.binding_energy is not None else 0.0
                    it_dg = _NumericTableItem(f"{dg_val:.2f}", dg_val)
                    it_dg.setTextAlignment(Qt.AlignCenter)
                    v_val = p.validation_rmsd if p.validation_rmsd is not None else (val.crystal_rmsd if i == 0 else None)
                    v_str = f"{v_val:.2f} Å" if v_val is not None else "—"
                    it_v = _NumericTableItem(v_str, v_val if v_val is not None else 999.0)
                    it_v.setTextAlignment(Qt.AlignCenter)

                    if v_val is not None:
                        if v_val <= 2.0:
                            badge, col = "✔ High Accuracy (≤ 2.0 Å)", "#10b981"
                        elif v_val <= 3.0:
                            badge, col = "⚠ Moderate (2.0–3.0 Å)", "#f59e0b"
                        else:
                            badge, col = "✖ Divergent (> 3.0 Å)", "#ef4444"
                    else:
                        badge, col = "Not evaluated", "#94a3b8"
                    it_v.setForeground(QColor(col))
                    it_k = _NumericTableItem("—", 0.0)
                    it_k.setTextAlignment(Qt.AlignCenter)
                    it_m = _NumericTableItem(str(val.matched_atoms or "—"), float(val.matched_atoms or 0))
                    it_m.setTextAlignment(Qt.AlignCenter)
                    it_q = QTableWidgetItem(badge)
                    it_q.setForeground(QColor(col))
                    for col_idx, item in enumerate([it_r, it_dg, it_v, it_k, it_m, it_q]):
                        self.table_val_poses.setItem(i, col_idx, item)
                self.table_val_poses.setSortingEnabled(True)
        else:
            if hasattr(self, "lbl_val_opt_status"):
                self.lbl_val_opt_status.setText("Validation Status: Optional (Select reference ligand to evaluate)")
                self.lbl_val_opt_status.setStyleSheet("color: #94a3b8; font-size: 11px;")

            self.val_cells["Reference Ligand Source"][1].setText("None provided")
            self.val_cells["Validation Status"][1].setText("Optional — not performed")
            self.val_cells["Validation Status"][1].setStyleSheet("color: #94a3b8;")
            self.val_cells["Validation RMSD"][1].setText("N/A — validation optional")
            self.val_cells["Search RMSD (Docking)"][1].setText("See Poses view for intra-conformation RMSD")
            self.val_cells["Matched Heavy Atoms"][1].setText("N/A")
            self.val_cells["Alignment / Mapping Method"][1].setText("N/A")
            if hasattr(self, "table_val_poses"):
                self.table_val_poses.setRowCount(0)

    # ─────────────────────────────────────────────────────────
    # THERMODYNAMICS VIEW
    # ─────────────────────────────────────────────────────────

    def _update_thermo_view(self, job: "DockingJob") -> None:
        from models import Engine
        res = job.canonical_result
        th = res.thermodynamics if res else None
        is_ad4 = (job.engine == Engine.AUTODOCK4)

        if is_ad4 and th:
            self.thermo_header.setText("STATISTICAL MECHANICS & THERMODYNAMIC ANALYSIS (AUTODOCK4 DLG)")
            self.thermo_info_label.setText(
                "<b>Descriptor Attribution:</b><br>"
                "Values shown here are reported by the selected AutoDock4 DLG. "
                "Missing descriptors are shown as <i>Not reported in this DLG</i>; ADSP does not infer them."
            )
            for key, (title, value) in self.thermo_cells.items():
                title.setVisible(True)
                value.setVisible(True)
            self.thermo_cells["Information Entropy (S_inf)"][1].setText(f"{th.info_entropy:.4f} nats" if th.info_entropy is not None else "Not reported in this DLG")
            self.thermo_cells["Partition Function (Q)"][1].setText(f"{th.partition_function:.4e}" if th.partition_function is not None else "Not reported in this DLG")
            self.thermo_cells["Statistical Temperature"][1].setText(f"{th.stat_temperature:.2f} K" if th.stat_temperature is not None else "Not reported in this DLG")
            self.thermo_cells["Statistical Free Energy (A)"][1].setText(f"{th.stat_free_energy:.2f} kcal/mol" if th.stat_free_energy is not None else "Not reported in this DLG")
            self.thermo_cells["Statistical Internal Energy (U)"][1].setText(f"{th.stat_internal_energy:.2f} kcal/mol" if th.stat_internal_energy is not None else "Not reported in this DLG")
            self.thermo_cells["Statistical Entropy (S_stat)"][1].setText(f"{th.stat_entropy:.4f} kcal/mol·K" if th.stat_entropy is not None else "Not reported in this DLG")

            # This is an ensemble-level DLG descriptor, not a selected-pose value.
            # ADSP does not recompute it from the canonical pose list.
            self.thermo_cells["Boltzmann Population Fraction"][1].setText(
                "Not reported in this DLG"
            )

            self.thermo_cells["Analysis Engine"][1].setText("BSNDVP™ AutoDock4 DLG Trajectory Parser")
        else:
            self.thermo_header.setText("THERMODYNAMIC / STATISTICAL ANALYSIS")
            self.thermo_info_label.setText(
                "<b>Not applicable for AutoDock Vina.</b><br>"
                "Vina affinity is reported directly from the Vina scoring output. "
                "ADSP derives no Ki or thermodynamic quantities from Vina."
            )
            # Do not present AD4-only descriptors under a Vina result.
            for key, (title, value) in self.thermo_cells.items():
                visible = key == "Analysis Engine"
                title.setVisible(visible)
                value.setVisible(visible)
            self.thermo_cells["Analysis Engine"][0].setText("ANALYSIS SCOPE")
            self.thermo_cells["Analysis Engine"][1].setText("Not applicable for AutoDock Vina")

    # ─────────────────────────────────────────────────────────
    # ADMET VIEW
    # ─────────────────────────────────────────────────────────

    def _update_admet_view(self, job: "DockingJob") -> None:
        self.table_admet.setRowCount(0)
        lig_path = self._resolve_ligand_path(job)
        if not lig_path:
            return

        def _worker():
            from prepare import calculate_molecule_descriptors
            try:
                d = calculate_molecule_descriptors(lig_path, obabel_exe=self.config.obabel_executable)
                if d:
                    self._signals.admet_ready.emit(d)
            except Exception:
                pass

        threading.Thread(target=_worker, daemon=True).start()

    def _on_admet_ready(self, d: Dict) -> None:
        mw = d.get("mw", 0.0)
        logp = d.get("logp", 0.0)
        psa = d.get("psa", 0.0)
        hbd = d.get("hbd", 0)
        hba = d.get("hba", 0)
        rotb = d.get("rotb", 0)
        ro5 = d.get("ro5_pass", True)

        rules = [
            ("Molecular Weight (MW)", f"{mw:.1f} Da", "≤ 500 Da (Lipinski)", "Pass" if mw <= 500 else "High (Violation)", "OpenBabel / Descriptor Engine"),
            ("Octanol/Water Partition (LogP)", f"{logp:.2f}", "≤ 5.0 (Lipinski)", "Pass" if logp <= 5.0 else "Lipophilic (Violation)", "OpenBabel / Descriptor Engine"),
            ("Polar Surface Area (PSA)", f"{psa:.1f} Å²", "≤ 140 Å² (Veber)", "Good Permeability" if psa <= 140 else "Poor Permeability", "OpenBabel / Descriptor Engine"),
            ("H-Bond Donors (HBD)", str(hbd), "≤ 5 (Lipinski)", "Pass" if hbd <= 5 else "Excess Donors", "OpenBabel / Descriptor Engine"),
            ("H-Bond Acceptors (HBA)", str(hba), "≤ 10 (Lipinski)", "Pass" if hba <= 10 else "Excess Acceptors", "OpenBabel / Descriptor Engine"),
            ("Rotatable Bonds", str(rotb), "≤ 10 (Veber)", "Conformationally Rigid" if rotb <= 10 else "High Flexibility", "OpenBabel / Descriptor Engine"),
            ("Lipinski Rule of 5 Compliance", "✔ PASS" if ro5 else "✖ VIOLATION", "Max 1 Violation Allowed", "Oral Drug-Like" if ro5 else "Sub-Optimal Bioavailability", "Lipinski et al. Consensus"),
        ]

        self.table_admet.setRowCount(len(rules))
        for r, (prop, val, thresh, interp, src) in enumerate(rules):
            it_prop = QTableWidgetItem(prop)
            it_val = QTableWidgetItem(val)
            it_val.setTextAlignment(Qt.AlignCenter)
            it_thresh = QTableWidgetItem(thresh)
            it_thresh.setTextAlignment(Qt.AlignCenter)

            it_interp = QTableWidgetItem(interp)
            it_interp.setTextAlignment(Qt.AlignCenter)
            if "Pass" in interp or "Good" in interp or "Rigid" in interp or "Drug-Like" in interp:
                it_interp.setForeground(QColor("#34d399"))
            else:
                it_interp.setForeground(QColor("#f87171"))

            it_src = QTableWidgetItem(src)

            for col, item in enumerate([it_prop, it_val, it_thresh, it_interp, it_src]):
                self.table_admet.setItem(r, col, item)

    # ─────────────────────────────────────────────────────────
    # PROVENANCE VIEW
    # ─────────────────────────────────────────────────────────

    def _update_provenance_view(self, job: "DockingJob") -> None:
        from models import __version__
        res = job.canonical_result
        prov = res.provenance if res else None

        records = [
            ("Application Version", prov.app_version if prov else __version__),
            ("Docking Engine", prov.engine if prov else job.engine.value),
            ("Engine Version", prov.engine_version if prov else "Auto-detected"),
            ("Docking Mode", prov.docking_mode if prov else job.docking_mode.value),
            ("Receptor Target File", str(job.receptor_path) if job.receptor_path else "N/A"),
            ("Receptor SHA-256 Hash", prov.receptor_hash if prov else "Computed at runtime"),
            ("Ligand File", str(job.ligand_path) if job.ligand_path else "N/A"),
            ("Ligand SHA-256 Hash", prov.ligand_hash if prov else "Computed at runtime"),
            ("Grid / Parameter File", str(job.config_path or job.gpf_path or "N/A")),
            ("Parameter SHA-256 Hash", prov.config_hash if prov else "N/A"),
            ("Executable Invoked", prov.executable_path if prov else "Resolved binary"),
            ("Run Completion Timestamp", (prov.completion_time or prov.start_time or "Recorded in run log") if prov else "Recorded in run log"),
            ("Output DLG / PDBQT Path", str(job.dlg_path or job.output_pdbqt or "N/A")),
            ("Elapsed Calculation Time", f"{job.elapsed_seconds:.1f} s"),
        ]

        self.table_provenance.setRowCount(len(records))
        for r, (attr, val) in enumerate(records):
            it_a = QTableWidgetItem(str(attr))
            _font_a = it_a.font()
            _font_a.setBold(True)
            it_a.setFont(_font_a)
            it_v = QTableWidgetItem(str(val) if val is not None else "N/A")
            self.table_provenance.setItem(r, 0, it_a)
            self.table_provenance.setItem(r, 1, it_v)

    def _copy_provenance_json(self) -> None:
        if not self._current_job:
            QMessageBox.information(self, "No Job Selected", "Select a docking job first.")
            return

        res = self._current_job.canonical_result
        if res and res.provenance:
            data = res.provenance.to_dict()
        else:
            data = {
                "job_id": self._current_job.job_id,
                "engine": self._current_job.engine.value,
                "receptor": self._current_job.receptor_name,
                "ligand": self._current_job.ligand_name,
                "status": self._current_job.status.value,
            }

        text = json.dumps(data, indent=2)
        QGuiApplication.clipboard().setText(text)
        QMessageBox.information(self, "Copied", "Provenance manifest copied to clipboard as JSON.")

    # ─────────────────────────────────────────────────────────
    # INTERACTIONS VIEW (Thread-Safe 2D Diagram & Contacts)
    # ─────────────────────────────────────────────────────────

    def _load_pose_interactions(self, job: "DockingJob", pose_idx: int, run_number: Optional[int] = None) -> None:
        if not job:
            return

        self._current_pose = pose_idx
        self._render_generation += 1
        gen = self._render_generation

        # Auto-resolve run_number from canonical poses if not provided
        if run_number is None and job.canonical_result and job.canonical_result.poses:
            for p in job.canonical_result.poses:
                p_pose_id = getattr(p, "pose_id", f"pose_{p.rank}")
                if p.rank == pose_idx or p_pose_id == f"pose_{pose_idx}":
                    if p.ad4_metrics and p.ad4_metrics.run_number is not None:
                        run_number = p.ad4_metrics.run_number
                    break

        self._update_interaction_header()

        rec_path = self._resolve_receptor_path(job)
        lig_path = self._resolve_ligand_path(job)

        if not rec_path or not lig_path:
            return

        cache_key = (job.job_id, pose_idx, run_number, "v2.0")
        if cache_key in self._interactions_cache:
            cached_inters = self._interactions_cache[cache_key]
            self._signals.interactions_ready.emit(cached_inters, job.job_id, pose_idx, gen)
            return

        job_id = job.job_id

        def _worker():
            try:
                from interactions import profile_docking_pose
                inters = profile_docking_pose(rec_path, lig_path, pose_index=pose_idx, run_number=run_number)
                self._interactions_cache[cache_key] = inters or []
                self._signals.interactions_ready.emit(inters or [], job_id, pose_idx, gen)
            except Exception:
                self._signals.interactions_ready.emit([], job_id, pose_idx, gen)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_interactions_ready(self, inters: List, job_id: str, pose_idx: int, gen: int) -> None:
        if gen != self._render_generation:
            return
        if not self._current_job or self._current_job.job_id != job_id:
            return
        if pose_idx != self._current_pose:
            return

        self._current_interactions = inters
        self._filter_interactions_table()
        if self._current_job:
            self._update_2d_diagram(self._current_job, interactions=inters, pose_idx=pose_idx, gen=gen)

    def _filter_interactions_table(self) -> None:
        filter_type = self.cmb_filter_inter.currentText()
        matching = []
        for it in self._current_interactions:
            itype = it.interaction_type
            if filter_type != "All Types" and filter_type.lower() not in itype.lower():
                continue
            matching.append(it)

        self.table_inter.setSortingEnabled(False)
        self.table_inter.setRowCount(len(matching))
        for r, it in enumerate(matching):
            it_pose = _NumericTableItem(str(it.pose_index), float(it.pose_index))
            it_pose.setTextAlignment(Qt.AlignCenter)

            it_type = QTableWidgetItem(it.interaction_type)
            if "Hydrogen" in it.interaction_type:
                it_type.setForeground(QColor("#38bdf8"))
            elif "Salt" in it.interaction_type:
                it_type.setForeground(QColor("#ec4899"))
            elif "Hydrophobic" in it.interaction_type:
                it_type.setForeground(QColor("#f59e0b"))
            elif "Pi" in it.interaction_type:
                it_type.setForeground(QColor("#a855f7"))
            elif "Halogen" in it.interaction_type:
                it_type.setForeground(QColor("#10b981"))
            else:
                it_type.setForeground(QColor("#94a3b8"))

            it_res = QTableWidgetItem(it.receptor_residue)
            it_rec_atom = QTableWidgetItem(it.receptor_atom)
            it_lig_atom = QTableWidgetItem(it.ligand_atom)

            it_dist = _NumericTableItem(f"{it.distance_angstrom:.2f}", float(it.distance_angstrom))
            it_dist.setTextAlignment(Qt.AlignCenter)

            basis_text = getattr(it, "detection_basis", "Computational Prediction (Geometric Criteria)")
            if getattr(it, "angle_deg", None) is not None:
                basis_text += f" [Angle {it.angle_deg:.1f}°]"
            it_basis = QTableWidgetItem(basis_text)
            it_basis.setForeground(QColor("#94a3b8"))

            row_items = [it_pose, it_type, it_res, it_rec_atom, it_lig_atom, it_dist, it_basis]
            for c, item in enumerate(row_items):
                self.table_inter.setItem(r, c, item)

        self.table_inter.setSortingEnabled(True)

    def _on_diag_theme_changed(self, idx: int) -> None:
        if self._current_job:
            self._update_2d_diagram(self._current_job, interactions=self._current_interactions, pose_idx=self._current_pose)

    def _update_2d_diagram(self, job: "DockingJob", interactions: List, pose_idx: int = 1, run_number: Optional[int] = None, gen: Optional[int] = None) -> None:
        if gen is None:
            self._render_generation += 1
            gen = self._render_generation
        lig_path = self._resolve_ligand_path(job)
        orig_lig_path = getattr(job, "ligand_path", None)
        theme_idx = self.cmb_diag_theme.currentIndex() if hasattr(self, "cmb_diag_theme") else 0
        is_pub = (theme_idx == 1)
        theme_str = "light" if is_pub else "dark"

        if run_number is None and job.canonical_result and job.canonical_result.poses:
            for p in job.canonical_result.poses:
                if p.rank == pose_idx or p.pose_id == f"pose_{pose_idx}":
                    if p.ad4_metrics and p.ad4_metrics.run_number is not None:
                        run_number = p.ad4_metrics.run_number
                    break

        job_id = job.job_id
        pose_meta = {
            "engine": job.engine.value if hasattr(job.engine, "value") else str(job.engine),
            "docking_mode": job.docking_mode.value if hasattr(job.docking_mode, "value") else str(job.docking_mode),
            "job_id": job.job_id,
        }

        def _worker():
            try:
                from gui.interaction_diagram import render_interaction_diagram
                res = render_interaction_diagram(
                    ligand_path=lig_path,
                    interactions=interactions or [],
                    pose_index=pose_idx,
                    run_number=run_number,
                    original_ligand_path=orig_lig_path,
                    pose_metadata=pose_meta,
                    obabel_exe=self.config.obabel_executable,
                    width=480, height=360,
                    theme=theme_str,
                    publication_mode=is_pub,
                )
                if res:
                    self._signals.diagram_ready.emit(res, job_id, pose_idx, gen)
            except Exception:
                pass

        threading.Thread(target=_worker, daemon=True).start()

    def _on_diagram_ready(self, diag_res: Any, job_id: str, pose_idx: int, gen: int) -> None:
        if gen != self._render_generation:
            return
        if not self._current_job or self._current_job.job_id != job_id:
            return
        if pose_idx != self._current_pose:
            return

        self._current_diagram_res = diag_res
        self._current_html = getattr(diag_res, "html", str(diag_res))
        self._current_svg = getattr(diag_res, "svg", "")

        if HAS_WEBENGINE:
            self.web_view.setHtml(self._current_html)
        elif HAS_SVG:
            svg_content = self._current_svg or self._current_html
            self.svg_view.load(QByteArray(svg_content.encode("utf-8")))

    def _save_diagram_file(self) -> None:
        if not self._current_svg and not self._current_html:
            QMessageBox.information(self, "No Diagram", "Select a docking pose to generate a diagram first.")
            return

        job = self._current_job
        lig_name = Path(job.ligand_name).stem if job else "ligand"
        eng_name = job.engine.value if (job and hasattr(job.engine, "value")) else "docking"
        mode_str = job.docking_mode.value if (job and hasattr(job.docking_mode, "value")) else "rigid"
        pose_str = f"pose{self._current_pose:02d}"
        default_base = f"{lig_name}_{eng_name}_{mode_str}_{pose_str}_interactions"
        default_base = re.sub(r'[\\/*?:"<>| ]', "_", default_base)

        initial_path = str(getattr(self.config, "project_root", Path.cwd()) / f"{default_base}.svg")
        f, sel_filter = QFileDialog.getSaveFileName(
            self, "Save 2D Interaction Diagram",
            initial_path,
            "SVG Vector Image (*.svg);;HTML Document (*.html);;PNG Image (*.png);;All Files (*.*)"
        )
        if not f:
            return

        p = Path(f)
        ext = p.suffix.lower()
        try:
            if ext == ".svg":
                if not self._current_svg:
                    raise ValueError("SVG representation not available for this diagram.")
                p.write_text(self._current_svg, encoding="utf-8")
            elif ext == ".html":
                if not self._current_html:
                    raise ValueError("HTML representation not available for this diagram.")
                p.write_text(self._current_html, encoding="utf-8")
            elif ext == ".png":
                if not self._current_svg:
                    raise ValueError("Cannot generate PNG: SVG representation missing.")
                try:
                    from PySide6.QtSvg import QSvgRenderer
                    from PySide6.QtGui import QImage, QPainter
                    renderer = QSvgRenderer(QByteArray(self._current_svg.encode("utf-8")))
                    if not renderer.isValid():
                        raise ValueError("SVG data is invalid for PNG rendering.")
                    img = QImage(1200, 900, QImage.Format_ARGB32)
                    theme_idx = self.cmb_diag_theme.currentIndex() if hasattr(self, "cmb_diag_theme") else 0
                    bg = Qt.white if theme_idx == 1 else QColor("#0b0f1a")
                    img.fill(bg)
                    painter = QPainter(img)
                    renderer.render(painter)
                    painter.end()
                    if not img.save(str(p), "PNG"):
                        raise RuntimeError("Failed to write PNG file.")
                except Exception as raster_err:
                    raise RuntimeError(f"PNG rasterization failed: {raster_err}")
            else:
                p.write_text(self._current_html, encoding="utf-8")

            if not p.exists() or p.stat().st_size == 0:
                raise IOError("Diagram export failed because the output file is missing or empty.")

            QMessageBox.information(
                self,
                "Export Complete",
                f"Successfully exported Pose {self._current_pose} diagram to:\n{p.name}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Export Failed", f"Diagram export failed:\n{e}")

    def _copy_contacts_to_clipboard(self) -> None:
        if not self._current_interactions:
            QMessageBox.information(self, "No Contacts", "No contacts available to copy.")
            return

        lines = ["\t".join(self.INTER_COLS)]
        for it in self._current_interactions:
            dist_val = f"{it.distance_angstrom:.2f}" if getattr(it, "distance_angstrom", None) is not None else "N/A"
            basis_val = getattr(it, "detection_basis", "Computational Prediction")
            if getattr(it, "angle_deg", None) is not None:
                basis_val += f" [Angle {it.angle_deg:.1f}°]"
            lines.append(f"{it.pose_index}\t{it.interaction_type}\t{it.receptor_residue}\t{it.receptor_atom}\t{it.ligand_atom}\t{dist_val}\t{basis_val}")

        QGuiApplication.clipboard().setText("\n".join(lines))
        QMessageBox.information(self, "Copied", f"Copied {len(self._current_interactions)} contact(s) to clipboard as TSV.")

    def _select_relative_pose(self, delta: int) -> None:
        if not self._current_job or not self._current_job.canonical_result:
            return
        poses = self._current_job.canonical_result.poses
        if not poses:
            return
        current_idx = 0
        for i, p in enumerate(poses):
            if p.rank == self._current_pose:
                current_idx = i
                break
        new_idx = max(0, min(len(poses) - 1, current_idx + delta))
        if new_idx == current_idx:
            return
        target_pose = poses[new_idx]
        self._current_pose = target_pose.rank
        run_num = target_pose.ad4_metrics.run_number if target_pose.ad4_metrics else None

        self.table_poses.blockSignals(True)
        self.table_poses.selectRow(new_idx)
        self.table_poses.blockSignals(False)

        self._load_pose_interactions(self._current_job, target_pose.rank, run_number=run_num)

    def _update_interaction_header(self) -> None:
        if not hasattr(self, "lbl_inter_header_info"):
            return
        if not self._current_job:
            self.lbl_inter_header_info.setText("<b>Interaction Analysis:</b> Select a docking job.")
            self.lbl_pose_nav_indicator.setText("Pose - of -")
            self.btn_prev_pose.setEnabled(False)
            self.btn_next_pose.setEnabled(False)
            return

        job = self._current_job
        lig_name = Path(job.ligand_name).stem
        eng_name = job.engine.value.upper() if hasattr(job.engine, "value") else str(job.engine).upper()
        mode_name = job.docking_mode.value.capitalize() if hasattr(job.docking_mode, "value") else str(job.docking_mode).capitalize()

        poses = job.canonical_result.poses if (job.canonical_result and job.canonical_result.poses) else []
        total_poses = len(poses) if poses else 1

        curr_idx = 0
        for i, p in enumerate(poses):
            if p.rank == self._current_pose:
                curr_idx = i
                break

        self.lbl_inter_header_info.setText(
            f"<b>Interaction Analysis</b> &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Ligand:</b> {lig_name} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Pose:</b> <span style='color:#38bdf8; font-weight:bold;'>{self._current_pose}</span> (of {total_poses}) &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Engine:</b> {eng_name} &nbsp;&nbsp;|&nbsp;&nbsp; "
            f"<b>Docking Mode:</b> {mode_name}"
        )
        self.lbl_pose_nav_indicator.setText(f"Pose {self._current_pose} of {total_poses}")
        self.btn_prev_pose.setEnabled(curr_idx > 0)
        self.btn_next_pose.setEnabled(curr_idx < total_poses - 1)

    # ─────────────────────────────────────────────────────────
    # Path Resolvers & Integrations
    # ─────────────────────────────────────────────────────────

    def _resolve_receptor_path(self, job: "DockingJob") -> Optional[Path]:
        p = getattr(job, "receptor_path", None)
        if p and Path(p).is_file():
            return Path(p)
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        candidates = [
            rec_dir / f"{job.receptor_name}.pdbqt",
            rec_dir / job.receptor_name / "rigid" / f"{job.receptor_name}.pdbqt",
            rec_dir / job.receptor_name / "rigid" / "receptor.pdbqt",
            rec_dir / job.receptor_name / f"{job.receptor_name}.pdbqt",
        ]
        if getattr(self.config, "project_root", None):
            ws_rec = self.config.project_root / "receptors"
            candidates.extend([
                ws_rec / f"{job.receptor_name}.pdbqt",
                ws_rec / job.receptor_name / "rigid" / f"{job.receptor_name}.pdbqt",
                ws_rec / job.receptor_name / "rigid" / "receptor.pdbqt",
                ws_rec / job.receptor_name / f"{job.receptor_name}.pdbqt",
            ])
        for c in candidates:
            if c.is_file():
                return c
        return None

    def _resolve_ligand_path(self, job: "DockingJob") -> Optional[Path]:
        from models import Engine
        # For AD4 jobs: docked coordinates live in the DLG file, not a PDBQT
        if getattr(job, "engine", None) == Engine.AUTODOCK4:
            dlg = getattr(job, "dlg_path", None)
            if dlg and Path(dlg).is_file():
                return Path(dlg)
            # Also try ad4_results dict
            ad4 = getattr(job, "ad4_results", None)
            if isinstance(ad4, dict):
                dlg2 = ad4.get("dlg_path")
                if dlg2 and Path(dlg2).is_file():
                    return Path(dlg2)
            # Fallback: search output_dir for .dlg
            if getattr(job, "output_dir", None) and Path(job.output_dir).is_dir():
                dlg_files = list(Path(job.output_dir).glob("*.dlg"))
                if dlg_files:
                    return dlg_files[0]

        # For Vina jobs (or AD4 fallback): output_pdbqt
        p = getattr(job, "output_pdbqt", None)
        if p and Path(p).is_file():
            return Path(p)
        split_poses = getattr(job, "split_poses", None)
        if split_poses and isinstance(split_poses, list):
            for sp in split_poses:
                if Path(sp).is_file():
                    return Path(sp)
        if getattr(job, "output_dir", None) and Path(job.output_dir).is_dir():
            for c in Path(job.output_dir).glob(f"*{job.ligand_name}*.pdbqt"):
                if c.is_file():
                    return c
        res_dir = Path(getattr(self.config, "result_directory", "results"))
        candidates = list(res_dir.rglob(f"*{job.ligand_name}*.pdbqt"))
        if candidates:
            return candidates[0]
        if getattr(self.config, "project_root", None):
            ws_res = self.config.project_root / "results"
            ws_candidates = list(ws_res.rglob(f"*{job.ligand_name}*.pdbqt"))
            if ws_candidates:
                return ws_candidates[0]
        if getattr(job, "ligand_path", None) and Path(job.ligand_path).is_file():
            return Path(job.ligand_path)
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        orig = lig_dir / f"{job.ligand_name}.pdbqt"
        if orig.is_file():
            return orig
        return None

    def _launch_3d_viewer(self) -> None:
        if not self._current_job:
            QMessageBox.information(self, "Selection Required", "Select a docking job first.")
            return
        rec_path = self._resolve_receptor_path(self._current_job)
        lig_path = self._resolve_ligand_path(self._current_job)
        if not rec_path:
            QMessageBox.warning(self, "Receptor Missing", "Receptor file could not be located.")
            return

        from viewer_3d import launch_3d_viewer
        try:
            launch_3d_viewer(receptor_path=rec_path, ligand_path=lig_path)
            QMessageBox.information(self, "3D Inspector", "Launched interactive 3D WebGL viewer in your web browser!")
        except Exception as e:
            QMessageBox.critical(self, "3D Inspector Error", str(e))

    def _launch_pymol(self) -> None:
        if not self._current_job:
            QMessageBox.information(self, "Selection Required", "Select a docking job first.")
            return

        job = self._current_job
        rec_path = self._resolve_receptor_path(job)
        lig_path = self._resolve_ligand_path(job)

        # Check for pre-built complex PDB / PDBQT files first
        complex_dir = job.output_dir / "complexes" if job.output_dir else None
        has_complexes = complex_dir and complex_dir.exists() and (
            any(complex_dir.glob("complex_pose_*.pdb")) or any(complex_dir.glob("complex_pose_*.pdbqt"))
        )

        # Try PyMol if available
        try:
            from pymol_exporter import export_pymol_session, find_pymol_executable, launch_pymol
            exe = find_pymol_executable()
            if exe and rec_path and lig_path:
                res_dir = Path(getattr(self.config, "result_directory", "results"))
                out_dir = res_dir / "pymol_sessions" / f"{job.job_id}_pose{self._current_pose}"
                try:
                    bundle = export_pymol_session(rec_path, lig_path, output_dir=out_dir, pose_index=self._current_pose)
                    target = bundle.get("pse") or bundle.get("pml")
                    if target and target.is_file():
                        launch_pymol(target, custom_pymol_path=exe)
                        QMessageBox.information(self, "PyMOL Launched", f"Launched PyMOL session:\n{target.name}")
                        return
                except Exception as e:
                    QMessageBox.warning(self, "PyMOL Error", f"PyMOL launch failed:\n{e}\n\nOpening complex folder instead.")
        except ImportError:
            pass  # PyMol not installed — fall through to complex folder

        # Fall back: open the complexes folder in the file explorer
        if has_complexes:
            import os
            try:
                os.startfile(str(complex_dir))
                QMessageBox.information(
                    self, "Complex Files (PDB / PDBQT)",
                    f"Opened the complex folder:\n{complex_dir}\n\n"
                    f"Standard PDB complexes (via OpenBabel) and PDBQT files are ready.\n"
                    f"Open any 'complex_pose_N.pdb' in your preferred viewer:\n"
                    f"  • UCSF Chimera / ChimeraX\n"
                    f"  • PyMol (File > Open)\n"
                    f"  • VMD\n"
                    f"  • Discovery Studio\n"
                    f"  • Maestro / Swiss-PdbViewer"
                )
            except Exception as e:
                QMessageBox.information(
                    self, "Complex Files (PDB / PDBQT)",
                    f"Complex PDB and PDBQT files are in:\n{complex_dir}\n\n"
                    f"Open any 'complex_pose_N.pdb' in your preferred viewer."
                )
        else:
            QMessageBox.information(
                self, "No Visualization Files",
                "No complex PDB files found for this job yet.\n\n"
                "Complex files are generated automatically after each successful docking run.\n"
                "If docking succeeded, check the job's output folder for a 'complexes/' subfolder."
            )

    def _run_dlg_extract(self) -> None:
        if not self.btn_run_dlg.isEnabled():
            return
        self.btn_run_dlg.setEnabled(False)
        self.btn_run_dlg.setText("⏳ Extracting...")
        if hasattr(self, "lbl_status_action"):
            self.lbl_status_action.setText("Extracting DLG results...")

        def _worker():
            try:
                from autodock4_workflow import _run_dlg_analysis
                _run_dlg_analysis(self.config, self._jobs)
                self._signals.dlg_complete.emit("", "")
            except Exception as e:
                self._signals.dlg_complete.emit("", str(e))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_dlg_complete(self, _msg: str, err: str) -> None:
        self.btn_run_dlg.setEnabled(True)
        self.btn_run_dlg.setText("🔬 Run DLG Extract")
        if hasattr(self, "lbl_status_action"):
            self.lbl_status_action.setText("")
        if err:
            QMessageBox.critical(self, "Extract Error", err)
        else:
            QMessageBox.information(
                self, "DLG Extract Complete",
                "BSNDVP™ DLG analysis executed successfully.\nPose metrics and reports have been updated."
            )
            self.refresh()

    def _export_reports(self) -> None:
        if not self.btn_export.isEnabled():
            return
        self.btn_export.setEnabled(False)
        self.btn_export.setText("⏳ Exporting...")
        if hasattr(self, "lbl_status_action"):
            self.lbl_status_action.setText("⏳ Preparing reports...")

        def _worker():
            try:
                jobs = list(self._jobs) if self._jobs else []
                if not jobs and getattr(self.config, "result_directory", None):
                    self._signals.export_progress.emit("Scanning workspace for completed jobs...")
                    from job_manager import load_job_status
                    saved_status = load_job_status(self.config.result_directory)
                    if saved_status:
                        from job_manager import build_job_queue
                        from models import JobStatus
                        all_jobs = build_job_queue(self.config)
                        for j in all_jobs:
                            if j.job_id in saved_status:
                                st = saved_status[j.job_id]
                                try:
                                    j.status = JobStatus(st.get("status", "SUCCESS"))
                                except Exception:
                                    j.status = JobStatus.SUCCESS
                                j.vina_results = st.get("vina_results", {})
                                j.ad4_results = st.get("ad4_results", {})
                        jobs = all_jobs

                if not jobs:
                    self._signals.export_complete.emit({}, "NO_JOBS")
                    return

                from reporting import generate_reports
                def _prog(msg: str):
                    self._signals.export_progress.emit(msg)

                res = generate_reports(jobs, self.config, progress_callback=_prog)
                self._signals.export_complete.emit(res, "")
            except Exception as e:
                self._signals.export_complete.emit({}, str(e))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_export_progress(self, msg: str) -> None:
        if hasattr(self, "lbl_status_action"):
            self.lbl_status_action.setText(msg)

    def _on_export_complete(self, res: dict, err: str) -> None:
        self.btn_export.setEnabled(True)
        self.btn_export.setText("📤 Export Reports")
        if hasattr(self, "lbl_status_action"):
            self.lbl_status_action.setText("")

        if err == "NO_JOBS":
            QMessageBox.information(self, "No Jobs", "No docking jobs available to export.")
            return
        if err:
            QMessageBox.critical(self, "Export Error", err)
            return

        xlsx_path = res.get("xlsx_path") if isinstance(res, dict) else None
        rep_dir = getattr(self.config, "report_directory", Path("reports"))
        msg = f"Generated CSV, Excel, and JSON reports in:\n{rep_dir}"
        if xlsx_path:
            msg += f"\n\nExcel Report:\n{xlsx_path}"
        QMessageBox.information(self, "Reports Exported", msg)
