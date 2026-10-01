"""
AutoDock Suite Pro — PySide6 Runner Tab (gui_qt/runner_tab.py)
==============================================================
Docking console with engine selection (Vina, AutoDock4, BOTH),
dual parameter controls, live progress tracking, and color-coded streaming log.
"""

from __future__ import annotations

import copy
import os
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Optional

from PySide6.QtCore import QObject, QThread, QTimer, Qt, Signal
from PySide6.QtGui import QColor, QFont, QTextCursor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig


class _StreamTee:
    """Tees stdout/stderr to both original console and Qt signal."""
    def __init__(self, orig_stream, emit_fn):
        self.orig_stream = orig_stream
        self.emit_fn = emit_fn
        self._buffer = ""

    def write(self, s: str):
        if self.orig_stream:
            try:
                self.orig_stream.write(s)
                self.orig_stream.flush()
            except Exception:
                pass
        self._buffer += str(s)
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            line = line.rstrip("\r")
            if line:
                self.emit_fn(line)

    def flush(self):
        if self._buffer:
            line = self._buffer.rstrip("\r\n")
            self._buffer = ""
            if line:
                self.emit_fn(line)
        if self.orig_stream:
            try:
                self.orig_stream.flush()
            except Exception:
                pass


class _PipelineWorker(QThread):
    """Background worker thread running docking workflows."""
    progress_signal = Signal(int, int, str)
    log_signal = Signal(str)
    done_signal = Signal(list, object)

    def __init__(self, config: "ProjectConfig"):
        super().__init__()
        self.config = config
        self._is_stopped = False

    def run(self):
        from models import Engine, ResumeMode
        orig_stdout = sys.stdout
        orig_stderr = sys.stderr
        tee_out = _StreamTee(orig_stdout, self.log_signal.emit)
        tee_err = _StreamTee(orig_stderr, self.log_signal.emit)
        sys.stdout = tee_out
        sys.stderr = tee_err

        try:
            from workflow_service import execute_workflow
            from models import ResumeMode

            jobs = execute_workflow(
                self.config,
                resume_mode=ResumeMode.RESUME,
                progress_callback=self._on_progress,
            )
            self.done_signal.emit(jobs or [], None)
        except Exception as e:
            self.done_signal.emit([], e)
        finally:
            tee_out.flush()
            tee_err.flush()
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr

    def _on_progress(self, completed: int, total: int, label: str = ""):
        self.progress_signal.emit(completed, total, label)


class RunnerTab(QWidget):
    """Docking execution console in PySide6 with real-time feedback."""

    def __init__(self, config: "ProjectConfig", on_config_change: Callable,
                 on_run_complete: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self.on_run_complete = on_run_complete
        self._worker: Optional[_PipelineWorker] = None
        self._console_font_size = 12
        self._elapsed_seconds = 0
        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.timeout.connect(self._on_timer_tick)
        self._current_progress_label = "Ready"
        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        # ── Control Bar ────────────────────────────────────────
        ctrl_card = QFrame()
        ctrl_card.setObjectName("cardFrame")
        ctrl_l = QHBoxLayout(ctrl_card)
        ctrl_l.setContentsMargins(12, 10, 12, 10)
        ctrl_l.setSpacing(10)

        ctrl_l.addWidget(QLabel("Docking Engine:"))
        self.cmb_engine = QComboBox()
        self.cmb_engine.addItems(["VINA", "AUTODOCK4", "BOTH"])
        cur_eng = self.config.engine.value if hasattr(self.config.engine, "value") else str(self.config.engine)
        if cur_eng in ["VINA", "AUTODOCK4", "BOTH"]:
            self.cmb_engine.setCurrentText(cur_eng)
        self.cmb_engine.setToolTip(
            "VINA: AutoDock Vina (fast Monte Carlo conformational search)\n"
            "AUTODOCK4: AutoDock 4.2.6 (Lamarckian Genetic Algorithm using AutoGrid4 maps)\n"
            "BOTH: Dual-Engine (runs AutoDock Vina, then AutoDock 4 sequentially)"
        )
        self.cmb_engine.currentTextChanged.connect(self._on_engine_change)
        ctrl_l.addWidget(self.cmb_engine)

        ctrl_l.addWidget(QLabel("Mode:"))
        self.cmb_mode = QComboBox()
        self.cmb_mode.addItems(["AUTO", "RIGID", "FLEXIBLE", "BOTH"])
        cur_mode = self.config.docking_mode.value if hasattr(self.config.docking_mode, "value") else str(self.config.docking_mode)
        if cur_mode in ["AUTO", "RIGID", "FLEXIBLE", "BOTH"]:
            self.cmb_mode.setCurrentText(cur_mode)
        ctrl_l.addWidget(self.cmb_mode)

        self.btn_toggle_params = QPushButton("⚙ Parameters ▲")
        self.btn_toggle_params.clicked.connect(self._toggle_params)
        ctrl_l.addWidget(self.btn_toggle_params)

        ctrl_l.addStretch()

        self.btn_start = QPushButton("▶  Start Docking")
        self.btn_start.setObjectName("accentButton")
        self.btn_start.setMinimumHeight(32)
        self.btn_start.setMinimumWidth(150)
        self.btn_start.clicked.connect(self._start_docking)
        ctrl_l.addWidget(self.btn_start)

        self.btn_stop = QPushButton("⏹  Stop")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self._stop_docking)
        ctrl_l.addWidget(self.btn_stop)

        main_layout.addWidget(ctrl_card)

        # ── Collapsible Parameters Container ───────────────────
        self.params_container = QFrame()
        self.params_container.setObjectName("cardFrame")
        self.params_layout = QVBoxLayout(self.params_container)
        self.params_layout.setContentsMargins(12, 10, 12, 10)
        self.params_layout.setSpacing(10)

        # Vina parameters group
        self.vina_group = QGroupBox("AutoDock Vina Parameters")
        vg_l = QHBoxLayout(self.vina_group)

        vg_l.addWidget(QLabel("Scoring Weights:"))
        self.cmb_vina_scoring = QComboBox()
        self.cmb_vina_scoring.addItems(["vina", "vinardo", "ad4"])
        self.cmb_vina_scoring.setCurrentText(getattr(self.config, "vina_scoring", "vina"))
        self.cmb_vina_scoring.setToolTip(
            "Scoring weights used *inside* AutoDock Vina search:\n"
            "• vina: Standard empirical scoring (default)\n"
            "• vinardo: Recalibrated scoring function (Quiroga & Villarreal 2016)\n"
            "• ad4: Vina Monte Carlo search using AutoDock4 scoring weights.\n"
            "(To run the full AutoDock 4.2.6 engine with Lamarckian GA, select AUTODOCK4 in the Docking Engine dropdown above)."
        )
        vg_l.addWidget(self.cmb_vina_scoring)

        vg_l.addWidget(QLabel("Exhaustiveness:"))
        self.spn_exhaust = QSpinBox()
        self.spn_exhaust.setRange(1, 128)
        self.spn_exhaust.setValue(getattr(self.config, "exhaustiveness", 8))
        vg_l.addWidget(self.spn_exhaust)

        vg_l.addWidget(QLabel("Num Modes:"))
        self.spn_modes = QSpinBox()
        self.spn_modes.setRange(1, 100)
        self.spn_modes.setValue(getattr(self.config, "num_modes", 9))
        vg_l.addWidget(self.spn_modes)

        vg_l.addWidget(QLabel("Energy Range (kcal):"))
        self.spn_energy = QDoubleSpinBox()
        self.spn_energy.setRange(0.1, 20.0)
        self.spn_energy.setValue(getattr(self.config, "energy_range", 6.0))
        vg_l.addWidget(self.spn_energy)

        vg_l.addWidget(QLabel("Min RMSD (Å):"))
        self.spn_min_rmsd = QDoubleSpinBox()
        self.spn_min_rmsd.setRange(0.1, 10.0)
        self.spn_min_rmsd.setValue(getattr(self.config, "vina_min_rmsd", 1.0))
        vg_l.addWidget(self.spn_min_rmsd)

        vg_l.addStretch()
        self.params_layout.addWidget(self.vina_group)

        # AutoDock 4 parameters group
        self.ad4_group = QGroupBox("AutoDock 4.2.6 (ADT) Full Parameter Suite")
        ad4_vbox = QVBoxLayout(self.ad4_group)
        ad4_vbox.setSpacing(6)

        # Row 1: Core Genetic Algorithm & AutoGrid settings
        ad4_row1 = QHBoxLayout()
        ad4_row1.addWidget(QLabel("Algorithm:"))
        self.cmb_ad4_alg = QComboBox()
        self.cmb_ad4_alg.addItems(["LGA", "GA", "LS"])
        self.cmb_ad4_alg.setCurrentText(getattr(self.config, "ad4_algorithm", "LGA"))
        self.cmb_ad4_alg.setToolTip(
            "LGA: Lamarckian Genetic Algorithm (GA + Solis-Wets Local Search; recommended)\n"
            "GA: Standard Genetic Algorithm without local search\n"
            "LS: Pure Solis-Wets Local Search"
        )
        ad4_row1.addWidget(self.cmb_ad4_alg)

        ad4_row1.addWidget(QLabel("ga_run (Runs):"))
        self.spn_ad4_runs = QSpinBox()
        self.spn_ad4_runs.setRange(1, 500)
        self.spn_ad4_runs.setValue(getattr(self.config, "ga_run", 100))
        ad4_row1.addWidget(self.spn_ad4_runs)

        ad4_row1.addWidget(QLabel("ga_pop_size:"))
        self.spn_ad4_pop = QSpinBox()
        self.spn_ad4_pop.setRange(10, 1000)
        self.spn_ad4_pop.setValue(getattr(self.config, "ga_pop_size", 150))
        ad4_row1.addWidget(self.spn_ad4_pop)

        ad4_row1.addWidget(QLabel("ga_num_evals:"))
        self.spn_ad4_evals = QSpinBox()
        self.spn_ad4_evals.setRange(10000, 50000000)
        self.spn_ad4_evals.setSingleStep(250000)
        self.spn_ad4_evals.setValue(getattr(self.config, "ga_num_evals", 2500000))
        ad4_row1.addWidget(self.spn_ad4_evals)

        ad4_row1.addWidget(QLabel("Generations:"))
        self.spn_ad4_gens = QSpinBox()
        self.spn_ad4_gens.setRange(1000, 100000)
        self.spn_ad4_gens.setSingleStep(1000)
        self.spn_ad4_gens.setValue(getattr(self.config, "ga_num_generations", 27000))
        ad4_row1.addWidget(self.spn_ad4_gens)

        ad4_row1.addWidget(QLabel("RMSD Tol (Å):"))
        self.spn_ad4_rmstol = QDoubleSpinBox()
        self.spn_ad4_rmstol.setRange(0.2, 10.0)
        self.spn_ad4_rmstol.setSingleStep(0.5)
        self.spn_ad4_rmstol.setValue(getattr(self.config, "rmstol", 2.0))
        ad4_row1.addWidget(self.spn_ad4_rmstol)

        self.chk_reuse_maps = QCheckBox("Reuse Maps")
        self.chk_reuse_maps.setChecked(getattr(self.config, "reuse_existing_maps", False))
        self.chk_reuse_maps.setToolTip("Skip AutoGrid4 grid generation if .map files already exist.")
        ad4_row1.addWidget(self.chk_reuse_maps)

        ad4_row1.addStretch()
        ad4_vbox.addLayout(ad4_row1)

        # Row 2: Advanced GA & Solis-Wets Local Search Fine-Tuning
        ad4_row2 = QHBoxLayout()
        ad4_row2.addWidget(QLabel("Elitism:"))
        self.spn_ad4_elitism = QSpinBox()
        self.spn_ad4_elitism.setRange(0, 10)
        self.spn_ad4_elitism.setValue(getattr(self.config, "ga_elitism", 1))
        self.spn_ad4_elitism.setToolTip("Number of top individuals that survive unchanged into the next generation (default 1).")
        ad4_row2.addWidget(self.spn_ad4_elitism)

        ad4_row2.addWidget(QLabel("Mutation Rate:"))
        self.spn_ad4_mutation = QDoubleSpinBox()
        self.spn_ad4_mutation.setRange(0.00, 1.00)
        self.spn_ad4_mutation.setSingleStep(0.01)
        self.spn_ad4_mutation.setDecimals(3)
        self.spn_ad4_mutation.setValue(getattr(self.config, "ga_mutation_rate", 0.02))
        self.spn_ad4_mutation.setToolTip("Rate of gene mutation per gene (default 0.02).")
        ad4_row2.addWidget(self.spn_ad4_mutation)

        ad4_row2.addWidget(QLabel("Crossover Rate:"))
        self.spn_ad4_crossover = QDoubleSpinBox()
        self.spn_ad4_crossover.setRange(0.00, 1.00)
        self.spn_ad4_crossover.setSingleStep(0.05)
        self.spn_ad4_crossover.setDecimals(2)
        self.spn_ad4_crossover.setValue(getattr(self.config, "ga_crossover_rate", 0.80))
        self.spn_ad4_crossover.setToolTip("Crossover rate between individuals (default 0.80).")
        ad4_row2.addWidget(self.spn_ad4_crossover)

        ad4_row2.addWidget(QLabel("LS Freq:"))
        self.spn_ad4_ls_freq = QDoubleSpinBox()
        self.spn_ad4_ls_freq.setRange(0.00, 1.00)
        self.spn_ad4_ls_freq.setSingleStep(0.01)
        self.spn_ad4_ls_freq.setDecimals(3)
        self.spn_ad4_ls_freq.setValue(getattr(self.config, "ls_search_freq", 0.06))
        self.spn_ad4_ls_freq.setToolTip("Probability of performing Solis-Wets local search on an individual (default 0.06).")
        ad4_row2.addWidget(self.spn_ad4_ls_freq)

        ad4_row2.addWidget(QLabel("SW Max Its:"))
        self.spn_ad4_sw_its = QSpinBox()
        self.spn_ad4_sw_its.setRange(10, 2000)
        self.spn_ad4_sw_its.setSingleStep(50)
        self.spn_ad4_sw_its.setValue(getattr(self.config, "sw_max_its", 300))
        self.spn_ad4_sw_its.setToolTip("Maximum iterations per Solis-Wets local search (default 300).")
        ad4_row2.addWidget(self.spn_ad4_sw_its)

        ad4_row2.addWidget(QLabel("SW Rho:"))
        self.spn_ad4_sw_rho = QDoubleSpinBox()
        self.spn_ad4_sw_rho.setRange(0.1, 10.0)
        self.spn_ad4_sw_rho.setSingleStep(0.1)
        self.spn_ad4_sw_rho.setValue(getattr(self.config, "sw_rho", 1.0))
        self.spn_ad4_sw_rho.setToolTip("Initial variance / step size for Solis-Wets local search (default 1.0).")
        ad4_row2.addWidget(self.spn_ad4_sw_rho)

        ad4_row2.addWidget(QLabel("Unbound:"))
        self.cmb_ad4_unbound = QComboBox()
        self.cmb_ad4_unbound.addItems(["bound", "extended", "compact"])
        self.cmb_ad4_unbound.setCurrentText(getattr(self.config, "unbound_model", "bound"))
        self.cmb_ad4_unbound.setToolTip("Model for unbound conformation free energy calculation (default 'bound').")
        ad4_row2.addWidget(self.cmb_ad4_unbound)

        ad4_row2.addStretch()
        ad4_vbox.addLayout(ad4_row2)

        self.params_layout.addWidget(self.ad4_group)

        main_layout.addWidget(self.params_container)
        self._update_param_view()

        # ── Execution Monitor & Timeline Card ──────────────────
        mon_card = QFrame()
        mon_card.setObjectName("cardFrame")
        mon_l = QVBoxLayout(mon_card)
        mon_l.setContentsMargins(14, 10, 14, 10)
        mon_l.setSpacing(8)

        # Header metrics row
        hdr_row = QHBoxLayout()
        hdr_row.setSpacing(14)

        self.lbl_mon_engine_mode = QLabel("ENGINE: VINA  ·  MODE: AUTO")
        self.lbl_mon_engine_mode.setStyleSheet("font-size: 11px; font-weight: 700; color: #0284c7;")
        hdr_row.addWidget(self.lbl_mon_engine_mode)

        hdr_row.addStretch()

        self.lbl_mon_vina_badge = QLabel("Vina: ○ Pending")
        self.lbl_mon_vina_badge.setStyleSheet(
            "font-size: 11px; font-weight: 600; color: #94a3b8; background: #1e293b; "
            "border: 1px solid #334155; border-radius: 4px; padding: 2px 8px;"
        )
        hdr_row.addWidget(self.lbl_mon_vina_badge)

        self.lbl_mon_ad4_badge = QLabel("AutoDock4: ○ Pending")
        self.lbl_mon_ad4_badge.setStyleSheet(
            "font-size: 11px; font-weight: 600; color: #94a3b8; background: #1e293b; "
            "border: 1px solid #334155; border-radius: 4px; padding: 2px 8px;"
        )
        hdr_row.addWidget(self.lbl_mon_ad4_badge)

        mon_l.addLayout(hdr_row)

        # Scientific Execution Timeline row
        self.timeline_labels = []
        timeline_row = QHBoxLayout()
        timeline_row.setSpacing(6)
        timeline_stages = [
            "Preparation",
            "Grid Verification",
            "Docking Search",
            "DLG Parsing",
            "Biophysical Analysis",
            "Report Generation",
        ]
        for idx, stage in enumerate(timeline_stages):
            lbl = QLabel(f"○ {stage}")
            lbl.setStyleSheet(
                "font-size: 10px; font-weight: 700; color: #64748b; background: #1e293b; "
                "border: 1px solid #334155; border-radius: 4px; padding: 3px 8px;"
            )
            timeline_row.addWidget(lbl)
            self.timeline_labels.append((lbl, stage))
            if idx < len(timeline_stages) - 1:
                arr = QLabel("→")
                arr.setStyleSheet("color: #475569; font-weight: bold;")
                timeline_row.addWidget(arr)
        timeline_row.addStretch()
        mon_l.addLayout(timeline_row)

        # Progress bar & label
        self.lbl_progress = QLabel("Ready — configure workspace & parameters, then click Start Docking")
        self.lbl_progress.setObjectName("subHeader")
        mon_l.addWidget(self.lbl_progress)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        mon_l.addWidget(self.progress_bar)

        # Status footer row
        stat_row = QHBoxLayout()
        self.lbl_mon_current = QLabel("Current Job: None")
        self.lbl_mon_current.setStyleSheet("font-size: 11px; color: #94a3b8;")
        stat_row.addWidget(self.lbl_mon_current)

        stat_row.addStretch()

        self.lbl_mon_counts = QLabel("Jobs: 0 / 0")
        self.lbl_mon_counts.setStyleSheet("font-size: 11px; font-weight: 600; color: #94a3b8;")
        stat_row.addWidget(self.lbl_mon_counts)

        self.lbl_mon_eta = QLabel("ETA: —")
        self.lbl_mon_eta.setStyleSheet("font-size: 11px; font-weight: 600; color: #94a3b8;")
        stat_row.addWidget(self.lbl_mon_eta)

        mon_l.addLayout(stat_row)

        main_layout.addWidget(mon_card)

        # ── Live Console Output ────────────────────────────────
        self.console = QTextEdit()
        self.console.setObjectName("consoleBox")
        self.console.setReadOnly(True)
        self._apply_console_font()
        main_layout.addWidget(self.console, 1)

        # Console Toolbar
        console_tb = QHBoxLayout()
        console_tb.addWidget(QLabel("Font Size:"))

        btn_zoom_out = QPushButton("A -")
        btn_zoom_out.setFixedWidth(40)
        btn_zoom_out.clicked.connect(self._zoom_out)
        console_tb.addWidget(btn_zoom_out)

        self.lbl_font_size = QLabel(f"{self._console_font_size} pt")
        console_tb.addWidget(self.lbl_font_size)

        btn_zoom_in = QPushButton("A +")
        btn_zoom_in.setFixedWidth(40)
        btn_zoom_in.clicked.connect(self._zoom_in)
        console_tb.addWidget(btn_zoom_in)

        console_tb.addStretch()

        btn_open_dlg = QPushButton("📁 Open DLG Folder")
        btn_open_dlg.clicked.connect(self._open_dlg_folder)
        console_tb.addWidget(btn_open_dlg)

        btn_open_res = QPushButton("📁 Open Results Folder")
        btn_open_res.clicked.connect(self._open_results_folder)
        console_tb.addWidget(btn_open_res)

        btn_save_log = QPushButton("💾 Save Log...")
        btn_save_log.clicked.connect(self._save_log_file)
        console_tb.addWidget(btn_save_log)

        btn_clear = QPushButton("🗑 Clear Log")
        btn_clear.clicked.connect(self.console.clear)
        console_tb.addWidget(btn_clear)

        main_layout.addLayout(console_tb)

    # ─────────────────────────────────────────────────────────
    # Param View & Engine Toggling
    # ─────────────────────────────────────────────────────────

    def _update_param_view(self) -> None:
        eng = self.cmb_engine.currentText().upper()
        if eng == "VINA":
            self.vina_group.setVisible(True)
            self.ad4_group.setVisible(False)
        elif eng == "AUTODOCK4":
            self.vina_group.setVisible(False)
            self.ad4_group.setVisible(True)
        elif eng == "BOTH":
            self.vina_group.setVisible(True)
            self.ad4_group.setVisible(True)

    def _on_engine_change(self, val: str) -> None:
        from models import Engine
        try:
            self.config.engine = Engine(val.upper())
        except ValueError:
            pass
        self._update_param_view()
        self.on_config_change()

    def _toggle_params(self) -> None:
        vis = not self.params_container.isVisible()
        self.params_container.setVisible(vis)
        self.btn_toggle_params.setText("⚙ Parameters ▲" if vis else "⚙ Parameters ▼")

    def refresh(self) -> None:
        """Sync UI controls with current configuration."""
        cur_eng = self.config.engine.value if hasattr(self.config.engine, "value") else str(self.config.engine)
        if cur_eng in ["VINA", "AUTODOCK4", "BOTH"]:
            self.cmb_engine.blockSignals(True)
            self.cmb_engine.setCurrentText(cur_eng)
            self.cmb_engine.blockSignals(False)

        cur_mode = self.config.docking_mode.value if hasattr(self.config.docking_mode, "value") else str(self.config.docking_mode)
        if cur_mode in ["AUTO", "RIGID", "FLEXIBLE", "BOTH"]:
            self.cmb_mode.blockSignals(True)
            self.cmb_mode.setCurrentText(cur_mode)
            self.cmb_mode.blockSignals(False)

        self.spn_energy.setValue(getattr(self.config, "energy_range", 6.0))
        self.spn_exhaust.setValue(getattr(self.config, "exhaustiveness", 8))
        self.spn_modes.setValue(getattr(self.config, "num_modes", 9))
        self.spn_min_rmsd.setValue(getattr(self.config, "vina_min_rmsd", 1.0))
        self._update_param_view()

    def _sync_params_to_config(self) -> None:
        from models import DockingMode, Engine
        try:
            self.config.engine = Engine(self.cmb_engine.currentText().upper())
            self.config.docking_mode = DockingMode(self.cmb_mode.currentText().upper())
        except Exception:
            pass

        self.config.vina_scoring = self.cmb_vina_scoring.currentText()
        self.config.exhaustiveness = self.spn_exhaust.value()
        self.config.num_modes = self.spn_modes.value()
        self.config.energy_range = self.spn_energy.value()
        self.config.vina_min_rmsd = self.spn_min_rmsd.value()

        self.config.ad4_algorithm = self.cmb_ad4_alg.currentText()
        self.config.ga_run = self.spn_ad4_runs.value()
        self.config.ga_pop_size = self.spn_ad4_pop.value()
        self.config.ga_num_evals = self.spn_ad4_evals.value()
        self.config.ga_num_generations = self.spn_ad4_gens.value()
        self.config.rmstol = self.spn_ad4_rmstol.value()
        self.config.reuse_existing_maps = self.chk_reuse_maps.isChecked()

        self.config.ga_elitism = self.spn_ad4_elitism.value()
        self.config.ga_mutation_rate = self.spn_ad4_mutation.value()
        self.config.ga_crossover_rate = self.spn_ad4_crossover.value()
        self.config.ls_search_freq = self.spn_ad4_ls_freq.value()
        self.config.sw_max_its = self.spn_ad4_sw_its.value()
        self.config.sw_rho = self.spn_ad4_sw_rho.value()
        self.config.unbound_model = self.cmb_ad4_unbound.currentText()

        self.on_config_change()

    def _update_timeline_stage(self, stage_idx: int) -> None:
        """Update scientific execution timeline stages: 0=Prep, 1=Grid, 2=Docking, 3=Parsing, 4=Analysis, 5=Reporting."""
        if not hasattr(self, "timeline_labels"):
            return
        for i, (lbl, name) in enumerate(self.timeline_labels):
            if i < stage_idx:
                lbl.setText(f"✓ {name}")
                lbl.setStyleSheet(
                    "font-size: 10px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 3px 8px;"
                )
            elif i == stage_idx:
                lbl.setText(f"● {name}")
                lbl.setStyleSheet(
                    "font-size: 10px; font-weight: 700; color: #0284c7; background: rgba(2, 132, 199, 0.15); "
                    "border: 1px solid rgba(2, 132, 199, 0.4); border-radius: 4px; padding: 3px 8px;"
                )
            else:
                lbl.setText(f"○ {name}")
                lbl.setStyleSheet(
                    "font-size: 10px; font-weight: 700; color: #64748b; background: #1e293b; "
                    "border: 1px solid #334155; border-radius: 4px; padding: 3px 8px;"
                )

    # ─────────────────────────────────────────────────────────
    # Docking Pipeline
    # ─────────────────────────────────────────────────────────

    def _start_docking(self) -> None:
        if self._worker and self._worker.isRunning():
            return
        self._sync_params_to_config()

        eng = self.cmb_engine.currentText().upper()
        mode = self.cmb_mode.currentText().upper()
        self.lbl_mon_engine_mode.setText(f"ENGINE: {eng}  ·  MODE: {mode}")

        if eng == "BOTH":
            self.lbl_mon_vina_badge.setText("Vina: ● Running")
            self.lbl_mon_vina_badge.setStyleSheet(
                "font-size: 11px; font-weight: 700; color: #0284c7; background: rgba(2, 132, 199, 0.15); "
                "border: 1px solid rgba(2, 132, 199, 0.4); border-radius: 4px; padding: 2px 8px;"
            )
            self.lbl_mon_ad4_badge.setText("AutoDock4: ○ Pending")
            self.lbl_mon_ad4_badge.setStyleSheet(
                "font-size: 11px; font-weight: 600; color: #94a3b8; background: #1e293b; "
                "border: 1px solid #334155; border-radius: 4px; padding: 2px 8px;"
            )
        elif eng == "AUTODOCK4":
            self.lbl_mon_vina_badge.setText("Vina: —")
            self.lbl_mon_ad4_badge.setText("AutoDock4: ● Running")
            self.lbl_mon_ad4_badge.setStyleSheet(
                "font-size: 11px; font-weight: 700; color: #7c3aed; background: rgba(124, 58, 237, 0.15); "
                "border: 1px solid rgba(124, 58, 237, 0.4); border-radius: 4px; padding: 2px 8px;"
            )
        else:
            self.lbl_mon_vina_badge.setText("Vina: ● Running")
            self.lbl_mon_vina_badge.setStyleSheet(
                "font-size: 11px; font-weight: 700; color: #0284c7; background: rgba(2, 132, 199, 0.15); "
                "border: 1px solid rgba(2, 132, 199, 0.4); border-radius: 4px; padding: 2px 8px;"
            )
            self.lbl_mon_ad4_badge.setText("AutoDock4: —")

        self._update_timeline_stage(0)

        self.btn_start.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.progress_bar.setValue(0)
        self.lbl_progress.setStyleSheet("color: #94a3b8;")
        self.lbl_progress.setText("Initialising docking pipeline...")
        self._elapsed_seconds = 0
        self._current_progress_label = "Initialising..."
        self._elapsed_timer.start(1000)

        # Setup worker thread
        self._worker = _PipelineWorker(self.config)
        self._worker.progress_signal.connect(self._on_worker_progress)
        self._worker.log_signal.connect(self._append_log)
        self._worker.done_signal.connect(self._on_worker_done)
        self._worker.start()

    def _stop_docking(self) -> None:
        from process_manager import process_manager
        if self._worker and self._worker.isRunning():
            self._worker._is_stopped = True
            process_manager.request_cancellation()
            self.lbl_progress.setText("Cancellation requested — terminating active process...")
            self.lbl_progress.setStyleSheet("color: #f87171;")
            self.btn_stop.setEnabled(False)

    def _on_timer_tick(self) -> None:
        self._elapsed_seconds += 1
        m, s = divmod(self._elapsed_seconds, 60)
        h, m = divmod(m, 60)
        time_str = f"{h:02d}:{m:02d}:{s:02d}" if h > 0 else f"{m:02d}:{s:02d}"
        if self._worker and self._worker.isRunning():
            self.lbl_progress.setText(f"{self._current_progress_label}  ·  Elapsed: {time_str}")

    def _on_worker_progress(self, completed: int, total: int, label: str):
        lbl_lower = label.lower()
        if "prep" in lbl_lower or "scaffold" in lbl_lower:
            self._update_timeline_stage(0)
        elif "grid" in lbl_lower or "map" in lbl_lower or "autogrid" in lbl_lower:
            self._update_timeline_stage(1)
        elif "dock" in lbl_lower or "vina" in lbl_lower or "ad4" in lbl_lower or "run" in lbl_lower:
            self._update_timeline_stage(2)
        elif "parse" in lbl_lower or "dlg" in lbl_lower or "extract" in lbl_lower:
            self._update_timeline_stage(3)
        elif "interact" in lbl_lower or "thermo" in lbl_lower or "admet" in lbl_lower or "analys" in lbl_lower:
            self._update_timeline_stage(4)
        elif "report" in lbl_lower or "export" in lbl_lower:
            self._update_timeline_stage(5)

        # Update engine badges in BOTH mode
        if "autodock4" in lbl_lower or "ad4" in lbl_lower:
            self.lbl_mon_vina_badge.setText("Vina: ✓ Complete")
            self.lbl_mon_vina_badge.setStyleSheet(
                "font-size: 11px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 2px 8px;"
            )
            self.lbl_mon_ad4_badge.setText("AutoDock4: ● Running")
            self.lbl_mon_ad4_badge.setStyleSheet(
                "font-size: 11px; font-weight: 700; color: #7c3aed; background: rgba(124, 58, 237, 0.15); "
                "border: 1px solid rgba(124, 58, 237, 0.4); border-radius: 4px; padding: 2px 8px;"
            )

        self.lbl_mon_current.setText(f"Current: {label}")
        self.lbl_mon_counts.setText(f"Jobs: {completed} / {total} ({max(0, total - completed)} left)")
        if completed > 0 and total > completed:
            eta_s = int((self._elapsed_seconds / completed) * (total - completed))
            em, es = divmod(eta_s, 60)
            eh, em = divmod(em, 60)
            self.lbl_mon_eta.setText(f"ETA: {eh:02d}:{em:02d}:{es:02d}" if eh > 0 else f"ETA: {em:02d}:{es:02d}")
        elif completed == total and total > 0:
            self.lbl_mon_eta.setText("ETA: Complete")

        if total > 0:
            pct = int(100 * completed / total)
            self.progress_bar.setValue(pct)
            m, s = divmod(self._elapsed_seconds, 60)
            h, m = divmod(m, 60)
            time_str = f"{h:02d}:{m:02d}:{s:02d}" if h > 0 else f"{m:02d}:{s:02d}"
            self._current_progress_label = f"Job {completed} / {total} — {label} [{pct}%]"
            self.lbl_progress.setText(f"{self._current_progress_label}  ·  Elapsed: {time_str}")

    def _append_log(self, text: str):
        # Color-coded append
        upper = text.upper()
        if any(m in upper for m in ("[OK]", "✔", "[SUCCESS]")):
            color = "#34d399"
        elif any(m in upper for m in ("[WARN]", "⚠", "[WARNING]")):
            color = "#fbbf24"
        elif any(m in upper for m in ("[FAIL]", "[ERROR]", "✖", "ERROR")):
            color = "#f87171"
        elif any(m in upper for m in ("DUAL-ENGINE", "PHASE", "ENGINE:")):
            color = "#38bdf8"
        else:
            color = "#e2e8f0"

        self.console.setTextColor(QColor(color))
        self.console.append(text)
        self.console.moveCursor(QTextCursor.End)

    def _on_worker_done(self, jobs: List, error: Optional[Exception]):
        self._elapsed_timer.stop()
        self.btn_start.setEnabled(True)
        self.btn_stop.setEnabled(False)
        m, s = divmod(self._elapsed_seconds, 60)
        h, m = divmod(m, 60)
        time_str = f"{h:02d}:{m:02d}:{s:02d}" if h > 0 else f"{m:02d}:{s:02d}"

        if error:
            self.progress_bar.setValue(0)
            self.lbl_progress.setText(f"Pipeline error after {time_str}: {error}")
            self.lbl_progress.setStyleSheet("color: #f87171; font-weight: bold;")
        else:
            self._update_timeline_stage(6)
            eng_val = self.cmb_engine.currentText().upper()
            if eng_val == "BOTH":
                self.lbl_mon_vina_badge.setText("Vina: ✓ Complete")
                self.lbl_mon_vina_badge.setStyleSheet(
                    "font-size: 11px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 2px 8px;"
                )
                self.lbl_mon_ad4_badge.setText("AutoDock4: ✓ Complete")
                self.lbl_mon_ad4_badge.setStyleSheet(
                    "font-size: 11px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 2px 8px;"
                )
            elif eng_val == "AUTODOCK4":
                self.lbl_mon_ad4_badge.setText("AutoDock4: ✓ Complete")
                self.lbl_mon_ad4_badge.setStyleSheet(
                    "font-size: 11px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 2px 8px;"
                )
            else:
                self.lbl_mon_vina_badge.setText("Vina: ✓ Complete")
                self.lbl_mon_vina_badge.setStyleSheet(
                    "font-size: 11px; font-weight: 700; color: #10b981; background: rgba(16, 185, 129, 0.15); "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 2px 8px;"
                )
            self.lbl_mon_eta.setText("ETA: Complete")
            self.progress_bar.setValue(100)
            self.lbl_progress.setText(f"✔ Completed in {time_str} — {len(jobs)} job(s) processed")
            self.lbl_progress.setStyleSheet("color: #34d399; font-weight: bold;")

            # Automatic DLG analysis for AD4 and BOTH engines
            from models import Engine
            eng = getattr(self.config, "engine", Engine.VINA)
            if eng in (Engine.AUTODOCK4, Engine.BOTH):
                try:
                    self._append_log("\n[INFO] Running automatic DLG extraction...")
                    from autodock4_workflow import _run_dlg_analysis
                    _run_dlg_analysis(self.config, jobs)
                    self._append_log("[OK] Automatic DLG extraction completed successfully.")
                except Exception as ex:
                    self._append_log(f"[WARNING] Automatic DLG extraction failed: {ex}")

        self.on_run_complete(jobs)

    def _save_log_file(self) -> None:
        f, _ = QFileDialog.getSaveFileName(
            self, "Save Docking Execution Log",
            str(getattr(self.config, "project_root", Path.cwd()) / "docking_execution.log"),
            "Log Files (*.log *.txt);;All Files (*.*)"
        )
        if f:
            try:
                Path(f).write_text(self.console.toPlainText(), encoding="utf-8")
                QMessageBox.information(self, "Log Saved", f"Saved console log to:\n{f}")
            except Exception as e:
                QMessageBox.critical(self, "Save Error", str(e))

    # ─────────────────────────────────────────────────────────
    # Helpers
    # ─────────────────────────────────────────────────────────

    def _apply_console_font(self) -> None:
        f = QFont("Consolas", max(8, min(24, self._console_font_size)))
        f.setStyleHint(QFont.Monospace)
        self.console.setFont(f)

    def _zoom_in(self) -> None:
        if self._console_font_size < 24:
            self._console_font_size += 1
            self._apply_console_font()
            self.lbl_font_size.setText(f"{self._console_font_size} pt")

    def _zoom_out(self) -> None:
        if self._console_font_size > 8:
            self._console_font_size -= 1
            self._apply_console_font()
            self.lbl_font_size.setText(f"{self._console_font_size} pt")

    def _open_dlg_folder(self) -> None:
        res = Path(getattr(self.config, "result_directory", "results")) / "DLG"
        res.mkdir(parents=True, exist_ok=True)
        if os.name == "nt":
            os.startfile(str(res))
        else:
            subprocess.run(["xdg-open", str(res)])

    def _open_results_folder(self) -> None:
        res = Path(getattr(self.config, "result_directory", "results"))
        res.mkdir(parents=True, exist_ok=True)
        if os.name == "nt":
            os.startfile(str(res))
        else:
            subprocess.run(["xdg-open", str(res)])
