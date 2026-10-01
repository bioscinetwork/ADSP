"""
AutoDock Suite Pro — PySide6 Scientific Platform Settings (gui_qt/settings_tab.py)
=================================================================================
Centralized settings management grouped into:
  - Appearance: Theme selection with palette swatch preview, UI scaling
  - Docking Engines: Executable resolution (Vina, AutoDock4, AutoGrid4, OpenBabel) with validation
  - Biophysical Analysis: DLG parsing options, ADMET descriptors, interaction cutoffs
  - Reporting & Export: CSV, XLSX, JSON, HTML, PyMOL session directories
  - Reproducibility & Integrity: Provenance records, SHA-256 input hashing, verified resume
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Callable, Optional

from PySide6.QtCore import Qt
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
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from config import resolve_executable
from executables import (
    get_autodock4_version,
    get_autogrid4_version,
    get_vina_split_version,
    get_vina_version,
    validate_executable,
)
from gui_qt.styles import THEME_NAMES, ThemePreviewWidget, get_theme
from models import __version__
from ad4_compatibility import PARAMETER_PROFILES, validate_profile_id

if TYPE_CHECKING:
    from config import ProjectConfig


class SettingsTab(QWidget):
    """Clean, grouped scientific settings and platform preferences."""

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, on_theme_change: Optional[Callable[[str], None]] = None, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self.on_theme_change = on_theme_change
        self._build_ui()
        self.load_from_config()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(12)

        # Header card
        header_card = QFrame()
        header_card.setObjectName("cardFrame")
        h_layout = QHBoxLayout(header_card)
        h_layout.setContentsMargins(14, 10, 14, 10)

        lbl_title = QLabel("⚙ Platform Settings & Scientific Configuration")
        lbl_title.setObjectName("sectionHeader")
        h_layout.addWidget(lbl_title)
        h_layout.addStretch()

        self.btn_save = QPushButton("💾 Save to project_config.toml")
        self.btn_save.setObjectName("accentButton")
        self.btn_save.clicked.connect(self._save_settings)
        h_layout.addWidget(self.btn_save)

        self.btn_reset = QPushButton("↺ Reload")
        self.btn_reset.clicked.connect(self.load_from_config)
        h_layout.addWidget(self.btn_reset)

        main_layout.addWidget(header_card)

        # Scroll area for grouped sections
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)

        container = QWidget()
        c_layout = QVBoxLayout(container)
        c_layout.setContentsMargins(0, 0, 0, 0)
        c_layout.setSpacing(16)

        # ── 1. Appearance & Theme ───────────────────────────────
        app_box = QGroupBox("1. Appearance & Visual Design System")
        app_l = QVBoxLayout(app_box)
        app_l.setSpacing(10)

        row_theme = QHBoxLayout()
        row_theme.addWidget(QLabel("Interface Theme:"))
        self.cmb_theme = QComboBox()
        self.cmb_theme.addItems(THEME_NAMES)
        self.cmb_theme.currentTextChanged.connect(self._on_theme_changed)
        row_theme.addWidget(self.cmb_theme)

        self.theme_swatch = ThemePreviewWidget(THEME_NAMES[0])
        row_theme.addWidget(self.theme_swatch)
        row_theme.addStretch()
        app_l.addLayout(row_theme)

        self.chk_high_dpi = QCheckBox("Enable High-DPI Automatic Geometry Scaling")
        self.chk_high_dpi.setChecked(True)
        app_l.addWidget(self.chk_high_dpi)

        c_layout.addWidget(app_box)

        # ── 2. Docking Engines & Binaries ───────────────────────
        eng_box = QGroupBox("2. Docking Engines & Executable Resolution")
        eng_l = QVBoxLayout(eng_box)
        eng_l.setSpacing(10)

        def make_path_row(label: str, attr_name: str, placeholder: str):
            box = QVBoxLayout()
            box.setSpacing(4)
            lbl = QLabel(label)
            lbl.setObjectName("subHeader")
            box.addWidget(lbl)

            h = QHBoxLayout()
            txt = QLineEdit()
            txt.setPlaceholderText(placeholder)
            btn = QPushButton("Browse...")
            lbl_ver = QLabel("—")
            lbl_ver.setStyleSheet("color: #94a3b8; font-size: 11px;")

            def on_browse():
                f, _ = QFileDialog.getOpenFileName(self, f"Select {label}", str(Path.cwd()), "Executables (*.exe);;All Files (*.*)")
                if f:
                    txt.setText(f)

            btn.clicked.connect(on_browse)
            h.addWidget(txt, 1)
            h.addWidget(btn)
            box.addLayout(h)
            box.addWidget(lbl_ver)
            eng_l.addLayout(box)
            return txt, lbl_ver

        self.txt_vina, self.lbl_vina_ver = make_path_row("AutoDock Vina Binary:", "vina_executable", "bin/vina.exe")
        self.txt_vina_split, self.lbl_vina_split_ver = make_path_row("Vina Splitter Binary:", "vina_split_executable", "bin/vina_split.exe")
        self.txt_ad4, self.lbl_ad4_ver = make_path_row("AutoDock 4.2.6 Binary:", "autodock4_executable", "bin/autodock4.exe")
        self.txt_autogrid, self.lbl_autogrid_ver = make_path_row("AutoGrid 4.2.6 Binary:", "autogrid4_executable", "bin/autogrid4.exe")

        profile_row = QHBoxLayout()
        profile_row.addWidget(QLabel("AutoDock4 Parameter Profile:"))
        self.cmb_ad4_profile = QComboBox()
        self.cmb_ad4_profile.addItems(sorted(PARAMETER_PROFILES))
        self.cmb_ad4_profile.currentTextChanged.connect(self._on_ad4_profile_changed)
        profile_row.addWidget(self.cmb_ad4_profile)
        self.lbl_ad4_profile_status = QLabel("—")
        self.lbl_ad4_profile_status.setStyleSheet("color: #94a3b8; font-size: 11px;")
        profile_row.addWidget(self.lbl_ad4_profile_status)
        profile_row.addStretch()
        eng_l.addLayout(profile_row)

        c_layout.addWidget(eng_box)

        # ── 3. Biophysical Analysis & Validation ────────────────
        ana_box = QGroupBox("3. Biophysical Analysis, Thermodynamics & ADMET")
        ana_l = QVBoxLayout(ana_box)
        ana_l.setSpacing(10)

        row_ana1 = QHBoxLayout()
        self.chk_enable_admet = QCheckBox("Calculate ADMET Physicochemical Descriptors (Lipinski / Veber)")
        self.chk_enable_admet.setChecked(True)
        row_ana1.addWidget(self.chk_enable_admet)

        self.chk_enable_thermo = QCheckBox("Extract Statistical Mechanics & Thermodynamic Quantities (Q, S, A)")
        self.chk_enable_thermo.setChecked(True)
        row_ana1.addWidget(self.chk_enable_thermo)
        ana_l.addLayout(row_ana1)

        row_ana2 = QHBoxLayout()
        row_ana2.addWidget(QLabel("Hydrogen Bond Distance Cutoff (Å):"))
        self.spn_hbond_dist = QDoubleSpinBox()
        self.spn_hbond_dist.setRange(2.0, 4.5)
        self.spn_hbond_dist.setValue(3.5)
        self.spn_hbond_dist.setSingleStep(0.1)
        row_ana2.addWidget(self.spn_hbond_dist)

        row_ana2.addWidget(QLabel("Cluster RMSD Tolerance (Å):"))
        self.spn_cluster_tol = QDoubleSpinBox()
        self.spn_cluster_tol.setRange(0.5, 5.0)
        self.spn_cluster_tol.setValue(2.0)
        self.spn_cluster_tol.setSingleStep(0.25)
        row_ana2.addWidget(self.spn_cluster_tol)
        row_ana2.addStretch()
        ana_l.addLayout(row_ana2)

        c_layout.addWidget(ana_box)

        # ── 4. Reporting & Publishing ───────────────────────────
        rep_box = QGroupBox("4. Publication Reports & File Exports")
        rep_l = QVBoxLayout(rep_box)
        rep_l.setSpacing(10)

        row_rep = QHBoxLayout()
        self.chk_rep_html = QCheckBox("HTML Research Report (with figures & provenance)")
        self.chk_rep_html.setChecked(True)
        row_rep.addWidget(self.chk_rep_html)

        self.chk_rep_excel = QCheckBox("Excel Workbook (.xlsx)")
        self.chk_rep_excel.setChecked(getattr(self.config, "export_excel", True))
        row_rep.addWidget(self.chk_rep_excel)

        self.chk_rep_csv = QCheckBox("Per-table CSV Files")
        self.chk_rep_csv.setChecked(getattr(self.config, "export_csv", True))
        row_rep.addWidget(self.chk_rep_csv)
        rep_l.addLayout(row_rep)

        c_layout.addWidget(rep_box)

        # ── 5. Reproducibility & Integrity ──────────────────────
        rep_box2 = QGroupBox("5. Provenance, Checksums & Checkpoint Verification")
        rep_l2 = QVBoxLayout(rep_box2)
        rep_l2.setSpacing(10)

        self.chk_strict_resume = QCheckBox("Strict Verified Resume (Verify SHA-256 hashes of input files before skipping)")
        self.chk_strict_resume.setChecked(True)
        rep_l2.addWidget(self.chk_strict_resume)

        self.chk_auto_checkpoint = QCheckBox("Save JSON Job Manifest Checkpoints after each job completes")
        self.chk_auto_checkpoint.setChecked(True)
        rep_l2.addWidget(self.chk_auto_checkpoint)

        c_layout.addWidget(rep_box2)
        c_layout.addStretch()

        scroll.setWidget(container)
        main_layout.addWidget(scroll, 1)

    def _on_theme_changed(self, theme_name: str) -> None:
        self.theme_swatch.set_theme(theme_name)
        if self.on_theme_change:
            self.on_theme_change(theme_name)

    def _on_ad4_profile_changed(self, profile_id: str) -> None:
        result = validate_profile_id(profile_id)
        self.lbl_ad4_profile_status.setText(
            f"{result.status}: {result.errors[0]}" if result.errors else result.status
        )
        self.lbl_ad4_profile_status.setStyleSheet(
            "color: #ef4444; font-size: 11px;" if result.errors else "color: #22c55e; font-size: 11px;"
        )

    def load_from_config(self) -> None:
        """Populate settings widgets from ProjectConfig."""
        saved_theme = getattr(self.config, "ui_theme", "Noir")
        if saved_theme in THEME_NAMES:
            self.cmb_theme.setCurrentText(saved_theme)
            self.theme_swatch.set_theme(saved_theme)

        self.txt_vina.setText(str(self.config.vina_executable or ""))
        self.txt_vina_split.setText(str(self.config.vina_split_executable or ""))
        self.txt_ad4.setText(str(self.config.autodock4_executable or ""))
        self.txt_autogrid.setText(str(self.config.autogrid4_executable or ""))
        profile_id = getattr(self.config, "ad4_parameter_profile", "ad4_standard_4.2")
        if profile_id in PARAMETER_PROFILES:
            self.cmb_ad4_profile.setCurrentText(profile_id)
        self._on_ad4_profile_changed(self.cmb_ad4_profile.currentText())

        # Check binary versions
        try:
            v_ver = get_vina_version(self.config.vina_executable)
            self.lbl_vina_ver.setText(f"Detected: {v_ver}" if v_ver else "Executable not found")
        except Exception:
            self.lbl_vina_ver.setText("Executable not found")

        try:
            a_ver = get_autodock4_version(self.config.autodock4_executable)
            self.lbl_ad4_ver.setText(f"Detected: {a_ver}" if a_ver else "Executable not found")
        except Exception:
            self.lbl_ad4_ver.setText("Executable not found")

        try:
            ag_ver = get_autogrid4_version(self.config.autogrid4_executable)
            self.lbl_autogrid_ver.setText(f"Detected: {ag_ver}" if ag_ver else "Executable not found")
        except Exception:
            self.lbl_autogrid_ver.setText("Executable not found")

    def _save_settings(self) -> None:
        """Write current settings back to configuration and project_config.toml."""
        from gui.settings import save_config_to_toml

        setattr(self.config, "ui_theme", self.cmb_theme.currentText())

        if self.txt_vina.text().strip():
            self.config.vina_executable = Path(self.txt_vina.text().strip())
        if self.txt_vina_split.text().strip():
            self.config.vina_split_executable = Path(self.txt_vina_split.text().strip())
        if self.txt_ad4.text().strip():
            self.config.autodock4_executable = Path(self.txt_ad4.text().strip())
        if self.txt_autogrid.text().strip():
            self.config.autogrid4_executable = Path(self.txt_autogrid.text().strip())
        self.config.ad4_parameter_profile = self.cmb_ad4_profile.currentText()

        self.config.export_excel = self.chk_rep_excel.isChecked()
        self.config.export_csv = self.chk_rep_csv.isChecked()

        ws = getattr(self.config, "project_root", Path.cwd())
        toml_path = Path(ws) / "project_config.toml"
        try:
            save_config_to_toml(self.config, toml_path)
            self.on_config_change()
            QMessageBox.information(self, "Settings Saved", f"Settings successfully updated and saved to:\n{toml_path}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save settings: {e}")
