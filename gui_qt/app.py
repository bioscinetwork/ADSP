"""
AutoDock Suite Pro — PySide6 Main Application Window (gui_qt/app.py)
====================================================================
The centralized desktop application window assembling:
  - Modern sidebar navigation architecture with QStackedWidget
  - Professional scientific header with project, engine, and status indicators
  - Live ThemePreviewWidget swatch & centralized design token themes
  - Non-blocking notification banner system
  - 8 core scientific workspaces:
      1. Workspace Studio (project, receptors, ligands)
      2. Preparation Studio (cleanup, protonation, flexible residues)
      3. Grid Box Studio (AutoGrid4 & Vina search space)
      4. Screening Library (ADMET & chemical properties)
      5. Docking Console (execution & live streaming log)
      6. Results & Profiler (poses, clusters, validation, thermodynamics, ADMET)
      7. Cross-Job Comparison (multi-engine experiment comparison)
      8. Platform Settings (paths, cutoffs, reporting, themes)
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional

from PySide6.QtCore import QTimer, Qt
from PySide6.QtGui import (
    QColor,
    QIcon,
    QKeySequence,
    QLinearGradient,
    QPainter,
    QPixmap,
    QShortcut,
)
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QComboBox,
    QDialog,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QStatusBar,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from gui_qt.comparison_tab import ExperimentComparisonTab
from gui_qt.grid_tab import GridBoxTab
from gui_qt.library_tab import LibraryTab
from gui_qt.preparation_tab import PreparationTab
from gui_qt.results_tab import ResultsTab
from gui_qt.runner_tab import RunnerTab
from gui_qt.settings_tab import SettingsTab
from gui_qt.styles import THEME_NAMES, ThemePreviewWidget, get_stylesheet, resolve_theme_name
from gui_qt.workspace_tab import WorkspaceTab
from models import SUITE_BANNER, SUITE_NAME, __version__

if TYPE_CHECKING:
    from config import ProjectConfig


class _GradientAccentBar(QWidget):
    """Subtle animated gradient accent strip under the titlebar."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(3)
        self._offset = 0.0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._step)
        self._timer.start(40)

    def _step(self):
        self._offset = (self._offset + 0.003) % 1.0
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        grad = QLinearGradient(0, 0, self.width(), 0)

        # 5 harmonic spectral colors
        colors = [
            QColor("#38bdf8"),  # Sky cyan
            QColor("#818cf8"),  # Indigo
            QColor("#c084fc"),  # Violet
            QColor("#f472b6"),  # Rose
            QColor("#38bdf8"),  # Loop back
        ]
        num = len(colors)
        for i, col in enumerate(colors):
            pos = (i / (num - 1) + self._offset) % 1.0
            grad.setColorAt(pos, col)

        painter.fillRect(self.rect(), grad)


class _NotificationBanner(QFrame):
    """Non-blocking transient notification banner for scientific events."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("notificationBanner")
        self.hide()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(16, 8, 16, 8)
        layout.setSpacing(10)

        self.lbl_icon = QLabel("ℹ")
        self.lbl_icon.setStyleSheet("font-size: 14px; font-weight: bold;")
        layout.addWidget(self.lbl_icon)

        self.lbl_msg = QLabel("")
        self.lbl_msg.setStyleSheet("font-size: 12px; font-weight: 600;")
        layout.addWidget(self.lbl_msg, 1)

        self.btn_dismiss = QPushButton("✕")
        self.btn_dismiss.setFixedSize(22, 22)
        self.btn_dismiss.setStyleSheet(
            "background: transparent; border: none; font-size: 12px; font-weight: bold; padding: 0;"
        )
        self.btn_dismiss.clicked.connect(self.hide)
        layout.addWidget(self.btn_dismiss)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self.hide)

    def show_notification(self, text: str, level: str = "info", duration_ms: int = 4500):
        self.lbl_msg.setText(text)
        level_lower = level.lower()
        if "success" in level_lower or "done" in level_lower or "complete" in level_lower:
            self.lbl_icon.setText("✔")
            self.lbl_icon.setStyleSheet("color: #10b981; font-size: 14px; font-weight: bold;")
            self.setStyleSheet(
                "background-color: rgba(16, 185, 129, 0.12); border-bottom: 1px solid rgba(16, 185, 129, 0.35);"
            )
        elif "warn" in level_lower:
            self.lbl_icon.setText("⚠")
            self.lbl_icon.setStyleSheet("color: #f59e0b; font-size: 14px; font-weight: bold;")
            self.setStyleSheet(
                "background-color: rgba(245, 158, 11, 0.12); border-bottom: 1px solid rgba(245, 158, 11, 0.35);"
            )
        elif "error" in level_lower or "fail" in level_lower:
            self.lbl_icon.setText("✖")
            self.lbl_icon.setStyleSheet("color: #ef4444; font-size: 14px; font-weight: bold;")
            self.setStyleSheet(
                "background-color: rgba(239, 68, 68, 0.12); border-bottom: 1px solid rgba(239, 68, 68, 0.35);"
            )
        else:
            self.lbl_icon.setText("ℹ")
            self.lbl_icon.setStyleSheet("color: #0284c7; font-size: 14px; font-weight: bold;")
            self.setStyleSheet(
                "background-color: rgba(2, 132, 199, 0.12); border-bottom: 1px solid rgba(2, 132, 199, 0.35);"
            )

        self.show()
        if duration_ms > 0:
            self._timer.start(duration_ms)


class HelpDialog(QDialog):
    """Interactive Scientific Reference & User Guide modal."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"About & Scientific Guide — {SUITE_NAME}")
        self.resize(780, 580)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        tb = QTextBrowser()
        tb.setOpenExternalLinks(True)
        tb.setHtml(f"""
        <div style="font-family: -apple-system, Segoe UI, sans-serif; line-height: 1.55;">
            <h2 style="color: #0284c7; margin-top: 0;">AutoDock Suite Pro <span style="font-size: 14px; color: #64748b;">v{__version__}</span></h2>
            <p><strong>AutoDock Suite Pro (ADSP)</strong> is an offline, publication-grade biophysical molecular docking platform unifying <strong>AutoDock Vina 1.2</strong> and <strong>AutoDock 4.2.6</strong> into an automated high-throughput virtual screening pipeline.</p>
            
            <hr style="border: none; border-top: 1px solid #cbd5e1; margin: 16px 0;" />

            <h3 style="color: #0284c7;">Workflow Architecture</h3>
            <ol>
                <li><strong>Workspace:</strong> Organize target macromolecules, screening ligands, and output directories. Direct drag-and-drop file import supported.</li>
                <li><strong>Preparation Studio:</strong> Clean up crystal waters, add polar hydrogens (pH 7.4), remove non-standard chains or cofactors, perform Alanine scanning / residue mutations, and split sidechains for flexible docking.</li>
                <li><strong>Grid Box:</strong> Set search space center and dimensions (Å). Auto-calculate bounding box from co-crystallized ligands or precompute AutoGrid 4.2 grid maps.</li>
                <li><strong>Screening Library:</strong> Small-molecule management with real-time ADMET physicochemical profiling (Lipinski Rule of 5 and Veber rules).</li>
                <li><strong>Docking Console:</strong> Execute AutoDock Vina, AutoDock 4, or dual-engine sequentially with live progress and colorized streaming console.</li>
                <li><strong>Results & Profiler:</strong> Analyze binding energies (ΔG), estimated inhibition constants (Ki / Kd), RMSD conformational clusters, 2D interaction maps, non-covalent residue contacts, and 3D WebGL / PyMOL inspection.</li>
                <li><strong>Cross-Job Comparison:</strong> Side-by-side multi-engine experiment comparison, delta-delta-G scoring analysis, and data export.</li>
                <li><strong>Settings:</strong> Centralized configuration for docking engine binaries, biophysical analysis thresholds, reporting formats, and visual theme selection.</li>
            </ol>

            <hr style="border: none; border-top: 1px solid #cbd5e1; margin: 16px 0;" />

            <h3 style="color: #0284c7;">Key Citations & Methodologies</h3>
            <ul>
                <li><em>AutoDock Vina:</em> Trott, O., & Olson, A. J. (2010). J. Comput. Chem., 31(2), 455–461.</li>
                <li><em>Vinardo:</em> Quiroga, R., & Villarreal, M. A. (2016). PLoS ONE, 11(5), e0155183.</li>
                <li><em>AutoDock 4.2:</em> Morris, G. M. et al. (2009). J. Comput. Chem., 30(16), 2785–2791.</li>
                <li><em>RDKit:</em> Landrum, G. et al. Open-source cheminformatics toolkit.</li>
                <li><em>OpenBabel:</em> O'Boyle, N. M. et al. (2011). J. Cheminform., 3, 33.</li>
            </ul>
        </div>
        """)
        layout.addWidget(tb)

        btn_close = QPushButton("Close")
        btn_close.setObjectName("accentButton")
        btn_close.clicked.connect(self.accept)
        layout.addWidget(btn_close, 0, Qt.AlignRight)


class MainWindow(QMainWindow):
    """Main application window for AutoDock Suite Pro in PySide6 with Sidebar Navigation."""

    NAV_ITEMS = [
        ("🗂", "Workspace", "Manage project files, macromolecules, and ligands"),
        ("🧪", "Preparation", "Receptor & ligand cleanup, protonation, and flexible splitting"),
        ("🎯", "Grid Box", "Search space definition & AutoGrid4 map precalculation"),
        ("📚", "Screening Library", "Compound library filtering and ADMET physicochemical profiling"),
        ("▶", "Docking Console", "Docking execution, stage monitor, and real-time streaming log"),
        ("📊", "Results & Profiler", "Poses, LGA clusters, interactions, and thermodynamic analysis"),
        ("⚖", "Comparison", "Cross-job side-by-side analysis and ΔΔG scoring comparison"),
        ("⚙", "Settings", "Docking engine paths, analysis thresholds, and appearance"),
    ]

    def __init__(self, config: Optional["ProjectConfig"] = None):
        super().__init__()
        from config import ProjectConfig, load_config

        self.config = config or self._load_default_config()

        self.setWindowTitle(f"{SUITE_NAME} v{__version__} — Molecular Docking & Virtual Screening")
        self.resize(1440, 960)
        self.setMinimumSize(1200, 820)

        # Assets
        self._assets_dir = Path(__file__).resolve().parent.parent / "gui" / "assets"
        ico = self._assets_dir / "app_icon.ico"
        if ico.is_file():
            self.setWindowIcon(QIcon(str(ico)))

        self._build_ui()
        self._setup_shortcuts()
        self._apply_theme(getattr(self.config, "ui_theme", "Noir"))

        # Clock timer for statusbar
        self._clock_timer = QTimer(self)
        self._clock_timer.timeout.connect(self._update_clock)
        self._clock_timer.start(1000)
        self._update_clock()
        self._update_engine_status()

    def _load_default_config(self):
        from config import ProjectConfig, load_config

        toml_path = Path("project_config.toml")
        if toml_path.is_file():
            try:
                return load_config(toml_path)
            except Exception:
                pass
        return ProjectConfig()

    def _build_ui(self) -> None:
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # ── Header Title Bar ──────────────────────────────────
        titlebar = QFrame()
        titlebar.setObjectName("titleBar")
        tb_layout = QHBoxLayout(titlebar)
        tb_layout.setContentsMargins(16, 8, 16, 8)
        tb_layout.setSpacing(12)

        logo_path = self._assets_dir / "logo_36.png"
        if not logo_path.is_file():
            logo_path = self._assets_dir / "logo_48.png"
        if logo_path.is_file():
            lbl_logo = QLabel()
            pix = QPixmap(str(logo_path)).scaled(32, 32, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            lbl_logo.setPixmap(pix)
            tb_layout.addWidget(lbl_logo)

        title_block = QVBoxLayout()
        title_block.setSpacing(1)
        lbl_app_title = QLabel(SUITE_NAME)
        lbl_app_title.setStyleSheet("font-size: 16px; font-weight: 800; letter-spacing: 0.4px;")
        title_block.addWidget(lbl_app_title)

        lbl_app_sub = QLabel(f"v{__version__}  ·  Molecular Docking & Virtual Screening")
        lbl_app_sub.setStyleSheet("font-size: 11px; color: #94a3b8;")
        title_block.addWidget(lbl_app_sub)
        tb_layout.addLayout(title_block)

        # Workspace Badge Pill
        self.lbl_ws_badge = QLabel("")
        self.lbl_ws_badge.setStyleSheet(
            "background-color: rgba(56, 189, 248, 0.12); color: #38bdf8; "
            "border: 1px solid rgba(56, 189, 248, 0.35); border-radius: 12px; "
            "padding: 3px 12px; font-size: 11px; font-weight: 600;"
        )
        tb_layout.addWidget(self.lbl_ws_badge)

        # Global Status Pill
        self.lbl_global_status = QLabel("● READY")
        self.lbl_global_status.setStyleSheet(
            "background-color: rgba(16, 185, 129, 0.12); color: #10b981; "
            "border: 1px solid rgba(16, 185, 129, 0.35); border-radius: 12px; "
            "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
        )
        tb_layout.addWidget(self.lbl_global_status)

        tb_layout.addStretch()

        # Engine Quick Indicators
        self.lbl_engine_tags = QLabel("")
        self.lbl_engine_tags.setStyleSheet("font-size: 11px; font-weight: 600; color: #94a3b8;")
        tb_layout.addWidget(self.lbl_engine_tags)

        # Live Theme Preview Swatch
        self.theme_preview = ThemePreviewWidget()
        tb_layout.addWidget(self.theme_preview)

        # Theme Selector Dropdown
        self.cmb_theme = QComboBox()
        self.cmb_theme.addItems(THEME_NAMES)
        cur_t = resolve_theme_name(getattr(self.config, "ui_theme", "Noir"))
        if cur_t in THEME_NAMES:
            self.cmb_theme.setCurrentText(cur_t)
        self.theme_preview.set_theme(self.cmb_theme.currentText())
        self.cmb_theme.currentTextChanged.connect(self._on_theme_dropdown_changed)
        tb_layout.addWidget(self.cmb_theme)

        # Quick Save & Open Config
        btn_save = QPushButton("💾 Save")
        btn_save.setToolTip("Save project configuration (Ctrl+S)")
        btn_save.clicked.connect(self._save_config)
        tb_layout.addWidget(btn_save)

        btn_open = QPushButton("📂 Open")
        btn_open.setToolTip("Open project configuration (Ctrl+O)")
        btn_open.clicked.connect(self._open_config)
        tb_layout.addWidget(btn_open)

        btn_help = QPushButton("❓ Guide")
        btn_help.setToolTip("Scientific guide & citations (F1)")
        btn_help.clicked.connect(self._show_help)
        tb_layout.addWidget(btn_help)

        main_layout.addWidget(titlebar)

        # ── Animated Gradient Accent Line ──────────────────────
        self.accent_bar = _GradientAccentBar()
        main_layout.addWidget(self.accent_bar)

        # ── Non-Blocking Notification Banner ───────────────────
        self.notification_banner = _NotificationBanner()
        main_layout.addWidget(self.notification_banner)

        # ── Main Content Area: Sidebar + Stacked Workspaces ───
        body_widget = QWidget()
        body_layout = QHBoxLayout(body_widget)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(0)

        # ── Left Navigation Sidebar ────────────────────────────
        self.sidebar_frame = QFrame()
        self.sidebar_frame.setObjectName("sidebarFrame")
        self.sidebar_frame.setFixedWidth(210)
        sidebar_layout = QVBoxLayout(self.sidebar_frame)
        sidebar_layout.setContentsMargins(8, 12, 8, 12)
        sidebar_layout.setSpacing(6)

        lbl_nav_title = QLabel("WORKSTATIONS")
        lbl_nav_title.setStyleSheet(
            "font-size: 10px; font-weight: 800; color: #64748b; letter-spacing: 1px; padding-left: 10px; margin-bottom: 4px;"
        )
        sidebar_layout.addWidget(lbl_nav_title)

        self.nav_button_group = QButtonGroup(self)
        self.nav_button_group.setExclusive(True)
        self.nav_buttons: List[QPushButton] = []

        for idx, (icon, text, tooltip) in enumerate(self.NAV_ITEMS):
            btn = QPushButton(f" {icon}  {text}")
            btn.setObjectName("navButton")
            btn.setCheckable(True)
            btn.setToolTip(tooltip)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda checked=False, i=idx: self.set_active_tab(i))
            self.nav_button_group.addButton(btn, idx)
            self.nav_buttons.append(btn)
            sidebar_layout.addWidget(btn)

        sidebar_layout.addStretch()

        # Sidebar Footer: Quick status info
        self.lbl_sidebar_info = QLabel("AutoDockSuite Pro\nLevel 3 Research Platform")
        self.lbl_sidebar_info.setStyleSheet(
            "font-size: 10px; color: #64748b; text-align: center; padding: 10px; line-height: 1.4;"
        )
        self.lbl_sidebar_info.setAlignment(Qt.AlignCenter)
        sidebar_layout.addWidget(self.lbl_sidebar_info)

        body_layout.addWidget(self.sidebar_frame)

        # ── Central Workspace Stack ────────────────────────────
        self.stack = QStackedWidget()

        # 1. Workspace Tab
        self.tab_workspace = WorkspaceTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_workspace)

        # 2. Preparation Studio
        self.tab_prep = PreparationTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_prep)

        # 3. Grid Box
        self.tab_grid = GridBoxTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_grid)

        # 4. Screening Library
        self.tab_library = LibraryTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_library)

        # 5. Run Docking
        self.tab_runner = RunnerTab(self.config, self._on_config_changed, self._on_run_complete)
        self.stack.addWidget(self.tab_runner)

        # 6. Results
        self.tab_results = ResultsTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_results)

        # 7. Comparison Studio
        self.tab_comparison = ExperimentComparisonTab(self.config, self._on_config_changed)
        self.stack.addWidget(self.tab_comparison)

        # 8. Settings
        self.tab_settings = SettingsTab(self.config, self._on_config_changed, self._apply_theme)
        self.stack.addWidget(self.tab_settings)

        # Backward compatibility alias
        self.tab_widget = self.stack

        body_layout.addWidget(self.stack, 1)
        main_layout.addWidget(body_widget, 1)

        # Select initial tab (Workspace)
        self.set_active_tab(0)

        # ── Status Bar ────────────────────────────────────────
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

        self.lbl_status_msg = QLabel("● Ready")
        self.lbl_status_msg.setStyleSheet("color: #10b981; font-weight: 600;")
        self.status_bar.addWidget(self.lbl_status_msg, 1)

        self.lbl_engine_info = QLabel("")
        self.status_bar.addPermanentWidget(self.lbl_engine_info)

        self.lbl_clock = QLabel("")
        self.lbl_clock.setStyleSheet("font-family: Consolas, monospace; color: #94a3b8; font-weight: 600;")
        self.status_bar.addPermanentWidget(self.lbl_clock)

        self._update_workspace_badge()

    def set_active_tab(self, index: int) -> None:
        """Switch active workspace by index and sync sidebar checked state."""
        if 0 <= index < self.stack.count():
            self.stack.setCurrentIndex(index)
            if index < len(self.nav_buttons):
                self.nav_buttons[index].setChecked(True)
            # Auto-refresh tabs that support it
            widget = self.stack.widget(index)
            if index == 6 and hasattr(self, "tab_results") and self.tab_results._jobs:
                self.tab_comparison.set_jobs(self.tab_results._jobs)
            elif hasattr(widget, "refresh"):
                try:
                    widget.refresh()
                except Exception:
                    pass

    def set_global_status(self, text: str, level: str = "ready") -> None:
        """Update the centralized top-level status pill."""
        self.lbl_global_status.setText(f"● {text.upper()}")
        lvl = level.lower()
        if "run" in lvl or "active" in lvl:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(2, 132, 199, 0.15); color: #0284c7; "
                "border: 1px solid rgba(2, 132, 199, 0.4); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )
        elif "analys" in lvl or "process" in lvl:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(124, 58, 237, 0.15); color: #8b5cf6; "
                "border: 1px solid rgba(124, 58, 237, 0.4); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )
        elif "success" in lvl or "done" in lvl or "complete" in lvl:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(16, 185, 129, 0.15); color: #10b981; "
                "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )
        elif "warn" in lvl or "cancel" in lvl:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(245, 158, 11, 0.15); color: #f59e0b; "
                "border: 1px solid rgba(245, 158, 11, 0.4); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )
        elif "error" in lvl or "fail" in lvl:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(239, 68, 68, 0.15); color: #ef4444; "
                "border: 1px solid rgba(239, 68, 68, 0.4); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )
        else:
            self.lbl_global_status.setStyleSheet(
                "background-color: rgba(16, 185, 129, 0.12); color: #10b981; "
                "border: 1px solid rgba(16, 185, 129, 0.35); border-radius: 12px; "
                "padding: 3px 12px; font-size: 11px; font-weight: 700; letter-spacing: 0.5px;"
            )

    def notify(self, message: str, level: str = "info", duration_ms: int = 4500) -> None:
        """Trigger a non-blocking in-app notification banner."""
        self.notification_banner.show_notification(message, level=level, duration_ms=duration_ms)

    def _update_workspace_badge(self) -> None:
        ws = getattr(self.config, "project_root", None) or Path.cwd()
        ws_name = Path(ws).name or str(ws)
        self.lbl_ws_badge.setText(f"📁 {ws_name}")
        self.lbl_ws_badge.setToolTip(f"Active project workspace directory:\n{ws}")

    def _setup_shortcuts(self) -> None:
        QShortcut(QKeySequence("F1"), self, self._show_help)
        QShortcut(QKeySequence.Save, self, self._save_config)
        QShortcut(QKeySequence.Open, self, self._open_config)
        QShortcut(QKeySequence("F5"), self, self._refresh_all)
        QShortcut(QKeySequence("Ctrl+1"), self, lambda: self.set_active_tab(0))
        QShortcut(QKeySequence("Ctrl+2"), self, lambda: self.set_active_tab(1))
        QShortcut(QKeySequence("Ctrl+3"), self, lambda: self.set_active_tab(2))
        QShortcut(QKeySequence("Ctrl+4"), self, lambda: self.set_active_tab(3))
        QShortcut(QKeySequence("Ctrl+5"), self, lambda: self.set_active_tab(4))
        QShortcut(QKeySequence("Ctrl+6"), self, lambda: self.set_active_tab(5))
        QShortcut(QKeySequence("Ctrl+7"), self, lambda: self.set_active_tab(6))
        QShortcut(QKeySequence("Ctrl+8"), self, lambda: self.set_active_tab(7))

    def _refresh_all(self) -> None:
        self._on_config_changed()
        self.lbl_status_msg.setText("Refreshed all workspaces and tabs.")
        self.notify("Refreshed all workspaces and experimental datasets.", level="info")

    def _show_help(self) -> None:
        dlg = HelpDialog(self)
        dlg.exec()

    def _on_theme_dropdown_changed(self, theme_name: str) -> None:
        self.theme_preview.set_theme(theme_name)
        self._apply_theme(theme_name)

    def _apply_theme(self, theme_name: str) -> None:
        resolved = resolve_theme_name(theme_name)
        qss = get_stylesheet(resolved)
        app = QApplication.instance()
        if app:
            app.setStyleSheet(qss)
        setattr(self.config, "ui_theme", resolved)
        if self.cmb_theme.currentText() != resolved:
            self.cmb_theme.blockSignals(True)
            self.cmb_theme.setCurrentText(resolved)
            self.cmb_theme.blockSignals(False)
        self.theme_preview.set_theme(resolved)

    def _update_clock(self) -> None:
        self.lbl_clock.setText(datetime.now().strftime("%H:%M:%S  "))

    def _update_engine_status(self) -> None:
        try:
            from executables import get_autodock4_version, get_vina_version

            v_ver = get_vina_version(self.config.vina_executable)
            a_ver = get_autodock4_version(self.config.autodock4_executable)
            v_ok = bool(v_ver and "not" not in v_ver.lower())
            a_ok = bool(a_ver and "not" not in a_ver.lower())
            v_dot = "🟢" if v_ok else "🔴"
            a_dot = "🟢" if a_ok else "🔴"
            self.lbl_engine_info.setText(f"{v_dot} Vina {v_ver}   |   {a_dot} AutoDock {a_ver}    ")
            self.lbl_engine_tags.setText(
                f"[Vina: {'OK' if v_ok else 'MISSING'}]  [AD4: {'OK' if a_ok else 'MISSING'}]"
            )
        except Exception:
            pass

    def _on_config_changed(self) -> None:
        self._update_workspace_badge()
        self._update_engine_status()
        # Propagate config updates across all tabs
        for tab in (
            self.tab_workspace,
            self.tab_prep,
            self.tab_grid,
            self.tab_library,
            self.tab_runner,
            self.tab_results,
            self.tab_comparison,
            self.tab_settings,
        ):
            if hasattr(tab, "config"):
                tab.config = self.config
            if hasattr(tab, "refresh"):
                try:
                    tab.refresh()
                except Exception:
                    pass

    def _on_run_complete(self, jobs: list) -> None:
        self.lbl_status_msg.setText(f"Docking complete — {len(jobs)} job(s) processed.")
        self.set_global_status("COMPLETED", level="success")
        self.notify(f"Docking workflow complete: {len(jobs)} job(s) processed.", level="success")
        self.tab_results.refresh()
        self.tab_comparison.set_jobs(self.tab_results._jobs if self.tab_results._jobs else jobs)
        # Automatically switch to Results tab (index 5)
        self.set_active_tab(5)

    def _save_config(self) -> None:
        from gui.settings import save_config_to_toml

        cfg_path = getattr(self.config, "project_root", Path.cwd()) / "project_config.toml"
        try:
            save_config_to_toml(self.config, cfg_path)
            self.lbl_status_msg.setText(f"Saved configuration to {cfg_path.name}")
            self.notify(f"Project configuration saved to {cfg_path.name}", level="success")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", str(e))

    def _open_config(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self,
            "Open project_config.toml",
            str(getattr(self.config, "project_root", Path.cwd())),
            "TOML Files (*.toml);;All Files (*.*)",
        )
        if f:
            from config import load_config

            try:
                self.config = load_config(Path(f))
                self._on_config_changed()
                self.lbl_status_msg.setText(f"Loaded configuration from {Path(f).name}")
                self.notify(f"Loaded configuration from {Path(f).name}", level="info")
            except Exception as e:
                QMessageBox.critical(self, "Load Error", str(e))


def run_gui(config: Optional["ProjectConfig"] = None) -> None:
    """Launch the PySide6 AutoDock Suite Pro application."""
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    window = MainWindow(config=config)
    window.show()
    app.exec()


launch_gui = run_gui
