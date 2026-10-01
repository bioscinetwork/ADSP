"""
AutoDock Suite Pro — PySide6 Grid Box Tab (gui_qt/grid_tab.py)
==============================================================
Grid center / size inputs with real-time grid point and volume calculation,
preset pocket dimensions, bounding-box calculation from co-crystallized ligand,
Vina config.txt and AutoDock4 GPF generation, and integrated AutoGrid4 runner.
"""

from __future__ import annotations

import math
import os
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Optional, Tuple

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig


class _AutoGridWorkerSignals(QObject):
    log_signal = Signal(str)
    done_signal = Signal(int)


class GridBoxTab(QWidget):
    """Grid Box configuration, file generation, and AutoGrid4 runner in PySide6."""

    PRESETS = {
        "Custom Dimensions": None,
        "Small Pocket (16 × 16 × 16 Å)": (16.0, 16.0, 16.0),
        "Standard Active Site (20 × 20 × 20 Å)": (20.0, 20.0, 20.0),
        "Extended Binding Site (26 × 26 × 26 Å)": (26.0, 26.0, 26.0),
        "Multi-Domain / Macro (32 × 32 × 32 Å)": (32.0, 32.0, 32.0),
        "Blind Docking / Whole Protein (45 × 45 × 45 Å)": (45.0, 45.0, 45.0),
    }

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change

        self._signals = _AutoGridWorkerSignals()
        self._signals.log_signal.connect(self._on_autogrid_log)
        self._signals.done_signal.connect(self._on_autogrid_done)

        self._build_ui()
        self.refresh()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        # Top Bar: Receptor selector
        top_card = QFrame()
        top_card.setObjectName("cardFrame")
        top_bar = QHBoxLayout(top_card)
        top_bar.setContentsMargins(10, 8, 10, 8)
        top_bar.setSpacing(10)

        lbl_rec_tag = QLabel("Target Receptor:")
        lbl_rec_tag.setObjectName("subHeader")
        top_bar.addWidget(lbl_rec_tag)

        self.cmb_receptors = QComboBox()
        self.cmb_receptors.setMinimumWidth(260)
        self.cmb_receptors.currentIndexChanged.connect(self._on_receptor_changed)
        top_bar.addWidget(self.cmb_receptors)

        self.btn_refresh_rec = QPushButton("⟳ Refresh Receptors")
        self.btn_refresh_rec.clicked.connect(self.refresh)
        top_bar.addWidget(self.btn_refresh_rec)

        top_bar.addStretch()

        top_bar.addWidget(QLabel("Preset Size:"))
        self.cmb_presets = QComboBox()
        self.cmb_presets.addItems(list(self.PRESETS.keys()))
        self.cmb_presets.setCurrentText("Standard Active Site (20 × 20 × 20 Å)")
        self.cmb_presets.currentTextChanged.connect(self._on_preset_changed)
        top_bar.addWidget(self.cmb_presets)

        main_layout.addWidget(top_card)

        # Split Layout: Left for Grid Parameters & Generators, Right for AutoGrid Runner
        body_layout = QHBoxLayout()
        body_layout.setSpacing(14)

        # ── Left Column ──
        left_col = QVBoxLayout()

        # Grid Box Parameters Box
        grid_box = QGroupBox("Grid Box Dimensions & Center")
        g_layout = QVBoxLayout(grid_box)
        g_layout.setSpacing(8)

        # Center Coordinates (X, Y, Z)
        lbl_center = QLabel("Center Coordinates (Å):")
        lbl_center.setObjectName("subHeader")
        g_layout.addWidget(lbl_center)

        row_c = QHBoxLayout()
        row_c.addWidget(QLabel("X:"))
        self.spn_cx = QDoubleSpinBox()
        self.spn_cx.setRange(-9999.0, 9999.0)
        self.spn_cx.setDecimals(3)
        self.spn_cx.setSingleStep(1.0)
        self.spn_cx.valueChanged.connect(self._update_grid_metrics)
        row_c.addWidget(self.spn_cx)

        row_c.addWidget(QLabel("Y:"))
        self.spn_cy = QDoubleSpinBox()
        self.spn_cy.setRange(-9999.0, 9999.0)
        self.spn_cy.setDecimals(3)
        self.spn_cy.setSingleStep(1.0)
        self.spn_cy.valueChanged.connect(self._update_grid_metrics)
        row_c.addWidget(self.spn_cy)

        row_c.addWidget(QLabel("Z:"))
        self.spn_cz = QDoubleSpinBox()
        self.spn_cz.setRange(-9999.0, 9999.0)
        self.spn_cz.setDecimals(3)
        self.spn_cz.setSingleStep(1.0)
        self.spn_cz.valueChanged.connect(self._update_grid_metrics)
        row_c.addWidget(self.spn_cz)
        g_layout.addLayout(row_c)

        # Size Dimensions (X, Y, Z)
        lbl_size = QLabel("Box Size (Å):")
        lbl_size.setObjectName("subHeader")
        g_layout.addWidget(lbl_size)

        row_s = QHBoxLayout()
        row_s.addWidget(QLabel("Size X:"))
        self.spn_sx = QDoubleSpinBox()
        self.spn_sx.setRange(1.0, 999.0)
        self.spn_sx.setValue(20.0)
        self.spn_sx.setSingleStep(1.0)
        self.spn_sx.valueChanged.connect(self._update_grid_metrics)
        row_s.addWidget(self.spn_sx)

        row_s.addWidget(QLabel("Size Y:"))
        self.spn_sy = QDoubleSpinBox()
        self.spn_sy.setRange(1.0, 999.0)
        self.spn_sy.setValue(20.0)
        self.spn_sy.setSingleStep(1.0)
        self.spn_sy.valueChanged.connect(self._update_grid_metrics)
        row_s.addWidget(self.spn_sy)

        row_s.addWidget(QLabel("Size Z:"))
        self.spn_sz = QDoubleSpinBox()
        self.spn_sz.setRange(1.0, 999.0)
        self.spn_sz.setValue(20.0)
        self.spn_sz.setSingleStep(1.0)
        self.spn_sz.valueChanged.connect(self._update_grid_metrics)
        row_s.addWidget(self.spn_sz)
        g_layout.addLayout(row_s)

        # Spacing
        row_sp = QHBoxLayout()
        row_sp.addWidget(QLabel("Grid Spacing (Å):"))
        self.spn_spacing = QDoubleSpinBox()
        self.spn_spacing.setRange(0.1, 5.0)
        self.spn_spacing.setValue(0.375)
        self.spn_spacing.setSingleStep(0.025)
        self.spn_spacing.valueChanged.connect(self._update_grid_metrics)
        row_sp.addWidget(self.spn_spacing)
        row_sp.addStretch()
        g_layout.addLayout(row_sp)

        # Live Real-time Grid Metric Card
        metric_card = QFrame()
        metric_card.setObjectName("cardFrame")
        mc_l = QVBoxLayout(metric_card)
        mc_l.setContentsMargins(10, 8, 10, 8)
        mc_l.setSpacing(4)

        self.lbl_points_metric = QLabel("Total Points: 0 (0 × 0 × 0)")
        self.lbl_points_metric.setStyleSheet("font-weight: bold; color: #38bdf8;")
        mc_l.addWidget(self.lbl_points_metric)

        self.lbl_volume_metric = QLabel("Volume: 0.0 Å³  ·  Est. Map Memory: 0.0 MB")
        self.lbl_volume_metric.setStyleSheet("color: #94a3b8; font-size: 11px;")
        mc_l.addWidget(self.lbl_volume_metric)

        self.lbl_points_warning = QLabel("")
        self.lbl_points_warning.setStyleSheet("color: #f59e0b; font-weight: bold; font-size: 11px;")
        mc_l.addWidget(self.lbl_points_warning)

        g_layout.addWidget(metric_card)
        left_col.addWidget(grid_box)

        # Co-crystallized Ligand Bounding Box Auto-Calculation
        colig_box = QGroupBox("Auto-Detect from Co-crystallized Ligand / Active Site")
        cl_layout = QVBoxLayout(colig_box)

        cl_row = QHBoxLayout()
        self.ent_coligand = QLineEdit()
        self.ent_coligand.setPlaceholderText("Select bound ligand PDB / PDBQT / MOL2...")
        cl_row.addWidget(self.ent_coligand, 1)

        btn_browse_colig = QPushButton("📂 Browse")
        btn_browse_colig.clicked.connect(self._browse_coligand)
        cl_row.addWidget(btn_browse_colig)
        cl_layout.addLayout(cl_row)

        btn_autofill = QPushButton("🎯 Auto-Calculate Grid Box (Center + 5Å Padding)")
        btn_autofill.clicked.connect(self._autofill_from_coligand)
        cl_layout.addWidget(btn_autofill)

        left_col.addWidget(colig_box)

        # Generation Buttons Box
        gen_box = QGroupBox("Generate Grid Configuration Files")
        gen_layout = QHBoxLayout(gen_box)

        self.btn_gen_vina = QPushButton("📄 Generate Vina config.txt")
        self.btn_gen_vina.setObjectName("accentButton")
        self.btn_gen_vina.clicked.connect(self._generate_vina_config)
        gen_layout.addWidget(self.btn_gen_vina)

        self.btn_gen_gpf = QPushButton("🧬 Generate AutoDock4 GPF")
        self.btn_gen_gpf.setObjectName("accentButton")
        self.btn_gen_gpf.clicked.connect(self._generate_ad4_gpf)
        gen_layout.addWidget(self.btn_gen_gpf)

        left_col.addWidget(gen_box)
        left_col.addStretch()
        body_layout.addLayout(left_col, 1)

        # ── Right Column: AutoGrid Runner & Preview ──
        right_col = QVBoxLayout()
        ag_box = QGroupBox("AutoGrid 4.2 Grid Map Generation (Pre-docking)")
        ag_layout = QVBoxLayout(ag_box)

        self.btn_run_autogrid = QPushButton("▶ Run AutoGrid4 Map Precalculation")
        self.btn_run_autogrid.setObjectName("accentButton")
        self.btn_run_autogrid.setMinimumHeight(34)
        self.btn_run_autogrid.clicked.connect(self._run_autogrid)
        ag_layout.addWidget(self.btn_run_autogrid)

        self.txt_autogrid_log = QTextEdit()
        self.txt_autogrid_log.setObjectName("consoleBox")
        self.txt_autogrid_log.setReadOnly(True)
        self.txt_autogrid_log.setPlaceholderText("AutoGrid4 execution log will stream here...")
        ag_layout.addWidget(self.txt_autogrid_log)

        right_col.addWidget(ag_box)
        body_layout.addLayout(right_col, 1)

        main_layout.addLayout(body_layout)
        self._update_grid_metrics()

    # ─────────────────────────────────────────────────────────
    # Logic & Calculations
    # ─────────────────────────────────────────────────────────

    def _update_grid_metrics(self) -> None:
        sx = self.spn_sx.value()
        sy = self.spn_sy.value()
        sz = self.spn_sz.value()
        sp = max(self.spn_spacing.value(), 0.05)

        # Points in AutoGrid4 must be even integers
        n_x = int(math.ceil(sx / sp))
        if n_x % 2 != 0:
            n_x += 1
        n_y = int(math.ceil(sy / sp))
        if n_y % 2 != 0:
            n_y += 1
        n_z = int(math.ceil(sz / sp))
        if n_z % 2 != 0:
            n_z += 1

        total_pts = n_x * n_y * n_z
        vol = sx * sy * sz
        # Est memory: 4 bytes per float * ~6 map types
        mem_mb = (total_pts * 4 * 6) / (1024 * 1024)

        self.lbl_points_metric.setText(f"Total Points: {total_pts:,} ({n_x} × {n_y} × {n_z})")
        self.lbl_volume_metric.setText(f"Volume: {vol:,.1f} Å³  ·  Est. Map Memory: ~{mem_mb:.2f} MB")

        if total_pts > 1_500_000:
            self.lbl_points_warning.setText("⚠ Large search space — consider reducing box size for faster docking.")
        elif total_pts < 8_000:
            self.lbl_points_warning.setText("ℹ Box size is quite compact — verify ligand fits with adequate clearance.")
        else:
            self.lbl_points_warning.setText("✔ Optimal search space dimensions for accurate screening.")

    def _on_preset_changed(self, preset_name: str) -> None:
        dims = self.PRESETS.get(preset_name)
        if dims:
            self.spn_sx.setValue(dims[0])
            self.spn_sy.setValue(dims[1])
            self.spn_sz.setValue(dims[2])
            self._update_grid_metrics()

    def refresh(self) -> None:
        """Scan receptor directory and populate dropdown while preserving user selection."""
        current_sel = self.cmb_receptors.currentText().strip()
        self.cmb_receptors.blockSignals(True)
        self.cmb_receptors.clear()
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        if not rec_dir.is_dir():
            self.cmb_receptors.blockSignals(False)
            return
        recs = []
        for entry in sorted(rec_dir.iterdir()):
            if entry.is_dir() and (entry / "rigid").is_dir():
                recs.append(entry.name)
            elif entry.suffix.lower() == ".pdbqt":
                recs.append(entry.stem)
        for r in recs:
            self.cmb_receptors.addItem(r)

        if current_sel and current_sel in recs:
            self.cmb_receptors.setCurrentText(current_sel)
        self.cmb_receptors.blockSignals(False)

    def _on_receptor_changed(self, idx: int) -> None:
        rec_name = self.cmb_receptors.currentText().strip()
        if not rec_name:
            return
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        candidates = [
            rec_dir / rec_name / "config.txt",
            rec_dir / f"{rec_name}.txt",
        ]
        from validators import parse_vina_config
        for cand in candidates:
            if cand.is_file():
                try:
                    cfg = parse_vina_config(cand)
                    if "center_x" in cfg:
                        self.spn_cx.setValue(float(cfg["center_x"]))
                        self.spn_cy.setValue(float(cfg["center_y"]))
                        self.spn_cz.setValue(float(cfg["center_z"]))
                    if "size_x" in cfg:
                        self.spn_sx.setValue(float(cfg["size_x"]))
                        self.spn_sy.setValue(float(cfg["size_y"]))
                        self.spn_sz.setValue(float(cfg["size_z"]))
                    self._update_grid_metrics()
                    break
                except Exception:
                    pass

    def _browse_coligand(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self, "Select Co-crystallized Ligand",
            str(getattr(self.config, "ligand_directory", Path.cwd())),
            "Molecular Files (*.pdb *.pdbqt *.sdf *.mol2);;All Files (*.*)"
        )
        if f:
            self.ent_coligand.setText(f)

    def _autofill_from_coligand(self) -> None:
        path_str = self.ent_coligand.text().strip()
        if not path_str or not Path(path_str).is_file():
            QMessageBox.warning(self, "File Required", "Please select a valid ligand file.")
            return

        from prepare import get_coligand_bbox
        try:
            bbox = get_coligand_bbox(Path(path_str), padding=5.0)
            if bbox:
                cx, cy, cz = bbox["center_x"], bbox["center_y"], bbox["center_z"]
                sx, sy, sz = bbox["size_x"], bbox["size_y"], bbox["size_z"]
                self.spn_cx.setValue(cx)
                self.spn_cy.setValue(cy)
                self.spn_cz.setValue(cz)
                self.spn_sx.setValue(sx)
                self.spn_sy.setValue(sy)
                self.spn_sz.setValue(sz)
                self._update_grid_metrics()
                QMessageBox.information(
                    self, "Grid Box Calculated",
                    f"Computed bounding box from {Path(path_str).name}:\n\n"
                    f"Center: ({cx:.2f}, {cy:.2f}, {cz:.2f})\n"
                    f"Size:   ({sx:.1f}, {sy:.1f}, {sz:.1f}) Å (with 5.0 Å padding)"
                )
            else:
                QMessageBox.warning(self, "Calculation Failed", "Could not extract atom coordinates.")
        except Exception as e:
            QMessageBox.critical(self, "Calculation Error", str(e))

    def _generate_vina_config(self) -> None:
        rec_name = self.cmb_receptors.currentText().strip()
        if not rec_name:
            QMessageBox.warning(self, "Selection Required", "Select a target receptor first.")
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        target_dir = rec_dir / rec_name
        target_dir.mkdir(parents=True, exist_ok=True)
        cfg_path = target_dir / "config.txt"

        content = (
            f"# AutoDock Vina Grid Configuration for {rec_name}\n"
            f"center_x = {self.spn_cx.value():.3f}\n"
            f"center_y = {self.spn_cy.value():.3f}\n"
            f"center_z = {self.spn_cz.value():.3f}\n\n"
            f"size_x = {self.spn_sx.value():.1f}\n"
            f"size_y = {self.spn_sy.value():.1f}\n"
            f"size_z = {self.spn_sz.value():.1f}\n\n"
            f"exhaustiveness = {getattr(self.config, 'exhaustiveness', 8)}\n"
            f"num_modes = {getattr(self.config, 'num_modes', 9)}\n"
            f"energy_range = {getattr(self.config, 'energy_range', 6.0):.1f}\n"
        )
        cfg_path.write_text(content, encoding="utf-8")
        QMessageBox.information(self, "Config Saved", f"Vina config.txt generated at:\n{cfg_path}")
        self.on_config_change()

    def _generate_ad4_gpf(self) -> None:
        rec_name = self.cmb_receptors.currentText().strip()
        if not rec_name:
            QMessageBox.warning(self, "Selection Required", "Select a target receptor first.")
            return

        from autodock4_workflow import generate_gpf_file
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        rec_pdbqt = rec_dir / rec_name / "rigid" / f"{rec_name}.pdbqt"
        if not rec_pdbqt.is_file():
            rec_pdbqt = rec_dir / rec_name / "rigid" / "receptor.pdbqt"
        if not rec_pdbqt.is_file():
            rec_pdbqt = rec_dir / f"{rec_name}.pdbqt"
        if not rec_pdbqt.is_file():
            cand = list((rec_dir / rec_name / "rigid").glob("*.pdbqt")) if (rec_dir / rec_name / "rigid").is_dir() else []
            if cand:
                rec_pdbqt = cand[0]

        if not rec_pdbqt.is_file():
            QMessageBox.warning(
                self, "Receptor Not Found",
                f"Could not locate a prepared PDBQT file for receptor '{rec_name}'.\n\n"
                f"Expected path:\n{rec_dir / rec_name / 'rigid' / f'{rec_name}.pdbqt'}\n\n"
                f"Please prepare the receptor in the Preparation Studio tab first."
            )
            return

        center = (self.spn_cx.value(), self.spn_cy.value(), self.spn_cz.value())
        size = (self.spn_sx.value(), self.spn_sy.value(), self.spn_sz.value())
        spacing = self.spn_spacing.value()

        gpf_out = rec_dir / rec_name / f"{rec_name}.gpf"
        try:
            generate_gpf_file(
                receptor_pdbqt=rec_pdbqt,
                output_gpf=gpf_out,
                center=center,
                size=size,
                spacing=spacing,
            )
            QMessageBox.information(self, "GPF Generated", f"AutoDock4 GPF generated at:\n{gpf_out}")
            self.on_config_change()
        except Exception as e:
            QMessageBox.critical(self, "GPF Generation Error", str(e))

    def _run_autogrid(self) -> None:
        rec_name = self.cmb_receptors.currentText().strip()
        if not rec_name:
            QMessageBox.warning(self, "Selection Required", "Select a target receptor first.")
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        gpf_cand = rec_dir / rec_name / f"{rec_name}.gpf"
        if not gpf_cand.is_file():
            QMessageBox.warning(self, "GPF Required", f"GPF file not found: {gpf_cand}\nPlease click 'Generate AutoDock4 GPF' first.")
            return

        ag_exe = self.config.autogrid4_executable
        if not ag_exe or not Path(ag_exe).is_file():
            QMessageBox.critical(self, "AutoGrid4 Missing", f"AutoGrid4 executable not found: {ag_exe}")
            return

        self.txt_autogrid_log.clear()
        self.txt_autogrid_log.append(f"▶ Starting AutoGrid4: {ag_exe.name} -p {gpf_cand.name}")
        self.btn_run_autogrid.setEnabled(False)

        import threading
        def _worker():
            try:
                proc = subprocess.Popen(
                    [str(ag_exe), "-p", str(gpf_cand), "-l", f"{rec_name}.glg"],
                    cwd=str(gpf_cand.parent),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                for line in iter(proc.stdout.readline, ""):
                    if line:
                        self._signals.log_signal.emit(line.rstrip())
                proc.stdout.close()
                proc.wait()
                self._signals.done_signal.emit(proc.returncode)
            except Exception as e:
                self._signals.log_signal.emit(f"Execution error: {e}")
                self._signals.done_signal.emit(-1)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_autogrid_log(self, text: str) -> None:
        self.txt_autogrid_log.append(text)

    def _on_autogrid_done(self, code: int) -> None:
        self.btn_run_autogrid.setEnabled(True)
        if code == 0:
            self.txt_autogrid_log.append("\n✔ AutoGrid4 grid map calculation completed successfully!")
            QMessageBox.information(self, "AutoGrid4 Complete", "Grid map calculation finished successfully!")
        else:
            self.txt_autogrid_log.append(f"\n✖ AutoGrid4 exited with code {code}")
            QMessageBox.warning(self, "AutoGrid4 Warning", f"AutoGrid4 exited with return code {code}.")
