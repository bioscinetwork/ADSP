"""
AutoDock Suite Pro — PySide6 Scientific Design System & Themes (gui_qt/styles.py)
================================================================================
Centralized scientific design tokens, typography, spacing, palettes, and QSS for:
  - Dark Studio      (Deep scientific workstation: obsidian/slate with cyan accent)
  - Scientific Light (Publication/journal academic: crisp white/slate with navy accent)
  - Molecular Plasma (Ultraviolet molecular visualization: deep indigo with magenta accent)
  - Bio-Neutral      (Restrained biological laboratory: dark forest with emerald/mint accent)
  - Noir             (Minimal high-contrast monochrome workstation)

Provides ThemePreviewWidget for visual palette preview, chart color palettes,
and unified component stylesheets.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from PySide6.QtCore import QRect, QSize, Qt
from PySide6.QtGui import QColor, QFont, QPainter, QPaintEvent
from PySide6.QtWidgets import QWidget


@dataclass
class ThemeDefinition:
    """Scientific design tokens for a complete application theme."""
    name: str
    # Surfaces
    bg_main: str
    bg_card: str
    bg_card_elevated: str
    bg_input: str
    # Borders
    border: str
    border_subtle: str
    # Typography
    text_primary: str
    text_secondary: str
    text_muted: str
    text_disabled: str
    # Accents & Interactions
    accent: str
    accent_hover: str
    accent_pressed: str
    accent_subtle: str
    selection_bg: str
    selection_text: str
    # Semantic Status
    status_ok: str
    status_warn: str
    status_fail: str
    status_info: str
    # Engine-Specific Identifiers
    engine_vina: str
    engine_ad4: str
    # Tables & Lists
    row_alt: str
    # Mode & Chart Palette
    is_dark: bool
    chart_palette: List[str] = field(default_factory=list)


THEMES: Dict[str, ThemeDefinition] = {
    "Dark Studio": ThemeDefinition(
        name="Dark Studio",
        bg_main="#080d1a",
        bg_card="#0f172a",
        bg_card_elevated="#1e293b",
        bg_input="#131d31",
        border="#24334a",
        border_subtle="#192437",
        text_primary="#f8fafc",
        text_secondary="#94a3b8",
        text_muted="#64748b",
        text_disabled="#475569",
        accent="#38bdf8",
        accent_hover="#0ea5e9",
        accent_pressed="#0284c7",
        accent_subtle="rgba(56, 189, 248, 0.14)",
        selection_bg="#0284c7",
        selection_text="#ffffff",
        status_ok="#10b981",
        status_warn="#f59e0b",
        status_fail="#ef4444",
        status_info="#38bdf8",
        engine_vina="#38bdf8",
        engine_ad4="#f59e0b",
        row_alt="#0c1424",
        is_dark=True,
        chart_palette=["#38bdf8", "#818cf8", "#c084fc", "#34d399", "#f59e0b", "#f472b6", "#2dd4bf"],
    ),
    "Scientific Light": ThemeDefinition(
        name="Scientific Light",
        bg_main="#f8fafc",
        bg_card="#ffffff",
        bg_card_elevated="#f1f5f9",
        bg_input="#f8fafc",
        border="#cbd5e1",
        border_subtle="#e2e8f0",
        text_primary="#0f172a",
        text_secondary="#475569",
        text_muted="#64748b",
        text_disabled="#94a3b8",
        accent="#0284c7",
        accent_hover="#0369a1",
        accent_pressed="#075985",
        accent_subtle="rgba(2, 132, 199, 0.10)",
        selection_bg="#0284c7",
        selection_text="#ffffff",
        status_ok="#16a34a",
        status_warn="#d97706",
        status_fail="#dc2626",
        status_info="#0284c7",
        engine_vina="#0284c7",
        engine_ad4="#d97706",
        row_alt="#f1f5f9",
        is_dark=False,
        chart_palette=["#0284c7", "#4f46e5", "#0d9488", "#16a34a", "#d97706", "#dc2626", "#7c3aed"],
    ),
    "Molecular Plasma": ThemeDefinition(
        name="Molecular Plasma",
        bg_main="#0a0117",
        bg_card="#150529",
        bg_card_elevated="#230a42",
        bg_input="#1e0638",
        border="#4c1d95",
        border_subtle="#321063",
        text_primary="#f5e8ff",
        text_secondary="#c084fc",
        text_muted="#9333ea",
        text_disabled="#7e22ce",
        accent="#e040fb",
        accent_hover="#ce93d8",
        accent_pressed="#ba68c8",
        accent_subtle="rgba(224, 64, 251, 0.14)",
        selection_bg="#7e22ce",
        selection_text="#ffffff",
        status_ok="#4ade80",
        status_warn="#facc15",
        status_fail="#f43f5e",
        status_info="#38bdf8",
        engine_vina="#38bdf8",
        engine_ad4="#e040fb",
        row_alt="#100320",
        is_dark=True,
        chart_palette=["#e040fb", "#818cf8", "#38bdf8", "#4ade80", "#fbbf24", "#f43f5e", "#a855f7"],
    ),
    "Bio-Neutral": ThemeDefinition(
        name="Bio-Neutral",
        bg_main="#021510",
        bg_card="#06241c",
        bg_card_elevated="#0d3b2e",
        bg_input="#093025",
        border="#065f46",
        border_subtle="#044331",
        text_primary="#ecfdf5",
        text_secondary="#6ee7b7",
        text_muted="#34d399",
        text_disabled="#047857",
        accent="#10b981",
        accent_hover="#059669",
        accent_pressed="#047857",
        accent_subtle="rgba(16, 185, 129, 0.14)",
        selection_bg="#059669",
        selection_text="#ffffff",
        status_ok="#34d399",
        status_warn="#fbbf24",
        status_fail="#f87171",
        status_info="#38bdf8",
        engine_vina="#38bdf8",
        engine_ad4="#10b981",
        row_alt="#041e17",
        is_dark=True,
        chart_palette=["#10b981", "#34d399", "#38bdf8", "#a78bfa", "#f59e0b", "#fb7185", "#06b6d4"],
    ),
    "Noir": ThemeDefinition(
        name="Noir",
        bg_main="#000000",
        bg_card="#0c0c0e",
        bg_card_elevated="#18181b",
        bg_input="#141416",
        border="#27272a",
        border_subtle="#18181b",
        text_primary="#ffffff",
        text_secondary="#a1a1aa",
        text_muted="#71717a",
        text_disabled="#52525b",
        accent="#f4f4f5",
        accent_hover="#e4e4e7",
        accent_pressed="#d4d4d8",
        accent_subtle="rgba(244, 244, 245, 0.14)",
        selection_bg="#27272a",
        selection_text="#ffffff",
        status_ok="#ffffff",
        status_warn="#e4e4e7",
        status_fail="#a1a1aa",
        status_info="#ffffff",
        engine_vina="#ffffff",
        engine_ad4="#a1a1aa",
        row_alt="#070708",
        is_dark=True,
        chart_palette=["#f4f4f5", "#a1a1aa", "#71717a", "#52525b", "#3f3f46", "#27272a", "#d4d4d8"],
    ),
}

THEME_DEFINITIONS: Dict[str, ThemeDefinition] = THEMES

# Backward compatibility alias mapping
THEME_ALIASES: Dict[str, str] = {
    "Dark": "Dark Studio",
    "Scientific": "Scientific Light",
    "Plasma": "Molecular Plasma",
    "Bio-Neon": "Bio-Neutral",
    "Noir": "Noir",
}

THEME_NAMES: List[str] = list(THEMES.keys())


def resolve_theme_name(name: str) -> str:
    """Normalize legacy or aliased theme names."""
    if name in THEMES:
        return name
    return THEME_ALIASES.get(name, "Dark Studio")


def get_theme(theme_name: str = "Dark Studio") -> ThemeDefinition:
    """Retrieve theme tokens by name, falling back gracefully to Dark Studio."""
    canonical = resolve_theme_name(theme_name)
    return THEMES.get(canonical, THEMES["Dark Studio"])


def get_chart_palette(theme_name: str = "Dark Studio") -> List[str]:
    """Retrieve publication-grade chart colors tailored for the theme."""
    t = get_theme(theme_name)
    return t.chart_palette or ["#38bdf8", "#818cf8", "#34d399", "#f59e0b", "#f472b6"]


# ═══════════════════════════════════════════════════════════════════════════════
# Theme Preview Swatch Widget
# ═══════════════════════════════════════════════════════════════════════════════

class ThemePreviewWidget(QWidget):
    """Compact swatch preview rendering background, card, and accent tokens."""

    def __init__(self, theme_name: str = "Dark Studio", parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._theme = get_theme(theme_name)
        self.setFixedSize(QSize(44, 22))
        self.setToolTip(f"Theme palette: {self._theme.name}")

    def set_theme(self, theme_name: str) -> None:
        self._theme = get_theme(theme_name)
        self.setToolTip(f"Theme palette: {self._theme.name}")
        self.update()

    def paintEvent(self, event: QPaintEvent) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)

        w = self.width()
        h = self.height()

        # Outer rounded border
        painter.setPen(QColor(self._theme.border))
        painter.setBrush(QColor(self._theme.bg_main))
        painter.drawRoundedRect(0, 0, w - 1, h - 1, 4, 4)

        # Card surface pill
        painter.setPen(Qt.NoPen)
        painter.setBrush(QColor(self._theme.bg_card))
        painter.drawRoundedRect(4, 3, w - 8, h - 6, 3, 3)

        # Accent dot
        painter.setBrush(QColor(self._theme.accent))
        painter.drawEllipse(w - 14, (h - 8) // 2, 8, 8)


# ═══════════════════════════════════════════════════════════════════════════════
# Centralized QSS Stylesheet Generator
# ═══════════════════════════════════════════════════════════════════════════════

def get_stylesheet(theme_name: str = "Dark Studio") -> str:
    """Generate a comprehensive, publication-grade Qt Style Sheet for the specified theme."""
    t = get_theme(theme_name)
    btn_text = "#000000" if t.name == "Noir" else ("#ffffff" if t.is_dark else "#ffffff")
    border_subtle = t.border_subtle or t.border
    bg_elevated = t.bg_card_elevated or t.bg_input

    return f"""
/* ── Publication-Grade Global Base ────────────────────────── */
QWidget {{
    background-color: {t.bg_main};
    color: {t.text_primary};
    font-family: "Segoe UI Variable Text", "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    font-size: 12px;
    letter-spacing: 0.15px;
    selection-background-color: {t.selection_bg};
    selection-color: {t.selection_text};
}}

QMainWindow, QDialog {{
    background-color: {t.bg_main};
}}

/* ── Typography & Section Headers ─────────────────────────── */
QLabel#appTitle {{
    font-size: 16px;
    font-weight: 800;
    color: {t.text_primary};
    letter-spacing: 0.5px;
}}

QLabel#appSubtitle {{
    font-size: 11px;
    color: {t.text_secondary};
    letter-spacing: 0.2px;
}}

QLabel#pageTitle {{
    font-size: 15px;
    font-weight: 800;
    color: {t.accent};
    letter-spacing: 0.4px;
}}

QLabel#sectionHeader {{
    font-size: 12px;
    font-weight: 800;
    color: {t.accent};
    letter-spacing: 0.6px;
    text-transform: uppercase;
    border-bottom: 2px solid {t.accent};
    padding-bottom: 4px;
    margin-bottom: 4px;
}}

QLabel#subHeader {{
    font-size: 11px;
    font-weight: 600;
    color: {t.text_secondary};
    letter-spacing: 0.2px;
}}

QLabel#kpiTitle {{
    font-size: 10px;
    font-weight: 700;
    color: {t.text_secondary};
    text-transform: uppercase;
    letter-spacing: 0.6px;
}}

QLabel#kpiValue {{
    font-size: 16px;
    font-weight: 800;
    color: {t.text_primary};
    font-family: "Consolas", monospace;
}}

/* ── Scientific KPI Summary Cards ─────────────────────────── */
QFrame#kpiCard {{
    background-color: {t.bg_card};
    border: 1px solid {t.border};
    border-radius: 8px;
    padding: 10px 14px;
}}

QFrame#kpiCard:hover {{
    border-color: {t.accent};
}}

QFrame#cardFrame {{
    background-color: {t.bg_card};
    border: 1px solid {t.border};
    border-radius: 8px;
}}

QFrame#cardFrameElevated {{
    background-color: {bg_elevated};
    border: 1px solid {t.border};
    border-radius: 8px;
}}

QFrame#titleBar {{
    background-color: {t.bg_card};
    border-bottom: 1px solid {t.border};
}}

/* ── Modern Navigation Sidebar ────────────────────────────── */
QFrame#sidebarFrame {{
    background-color: {t.bg_card};
    border-right: 1px solid {t.border};
}}

QPushButton#navButton {{
    background-color: transparent;
    color: {t.text_secondary};
    border: none;
    border-radius: 6px;
    padding: 9px 14px;
    font-size: 12px;
    font-weight: 600;
    text-align: left;
    margin: 2px 6px;
}}

QPushButton#navButton:hover {{
    background-color: {t.bg_card_elevated};
    color: {t.text_primary};
}}

QPushButton#navButton:checked {{
    background-color: {t.accent_subtle};
    color: {t.accent};
    font-weight: 700;
    border-left: 3px solid {t.accent};
}}

/* ── Group Boxes & Fieldsets ──────────────────────────────── */
QGroupBox {{
    background-color: {t.bg_card};
    border: 1px solid {t.border};
    border-radius: 8px;
    margin-top: 22px;
    padding: 16px 14px 14px 14px;
    font-weight: 700;
    font-size: 11px;
    color: {t.accent};
    letter-spacing: 0.5px;
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 14px;
    padding: 3px 10px;
    background-color: {t.bg_input};
    border: 1px solid {t.border};
    border-radius: 4px;
    color: {t.accent};
    font-weight: 700;
}}

/* ── Scientific Stepper / Workflow Progression ────────────── */
QFrame#stepperContainer {{
    background-color: {t.bg_card};
    border: 1px solid {t.border};
    border-radius: 8px;
    padding: 6px 12px;
}}

QLabel#stepperDone {{
    color: {t.status_ok};
    font-weight: 700;
    font-size: 11px;
}}

QLabel#stepperActive {{
    color: {t.accent};
    font-weight: 800;
    font-size: 11px;
    background-color: {t.accent_subtle};
    border: 1px solid {t.accent};
    border-radius: 4px;
    padding: 2px 8px;
}}

QLabel#stepperPending {{
    color: {t.text_muted};
    font-size: 11px;
}}

/* ── Status & Engine Badges ───────────────────────────────── */
QLabel#badgeVina {{
    background-color: rgba(56, 189, 248, 0.16);
    color: {t.engine_vina};
    border: 1px solid {t.engine_vina};
    border-radius: 10px;
    padding: 2px 10px;
    font-size: 11px;
    font-weight: 700;
}}

QLabel#badgeAD4 {{
    background-color: rgba(245, 158, 11, 0.16);
    color: {t.engine_ad4};
    border: 1px solid {t.engine_ad4};
    border-radius: 10px;
    padding: 2px 10px;
    font-size: 11px;
    font-weight: 700;
}}

QLabel#badgeRigid {{
    background-color: {t.bg_card_elevated};
    color: {t.text_secondary};
    border: 1px solid {t.border};
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 600;
}}

QLabel#badgeFlex {{
    background-color: rgba(168, 85, 247, 0.16);
    color: #c084fc;
    border: 1px solid #c084fc;
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 700;
}}

QLabel#statusBadgeOk {{
    background-color: rgba(16, 185, 129, 0.16);
    color: {t.status_ok};
    border: 1px solid {t.status_ok};
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 700;
}}

QLabel#statusBadgeWarn {{
    background-color: rgba(245, 158, 11, 0.16);
    color: {t.status_warn};
    border: 1px solid {t.status_warn};
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 700;
}}

QLabel#statusBadgeFail {{
    background-color: rgba(239, 68, 68, 0.16);
    color: {t.status_fail};
    border: 1px solid {t.status_fail};
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 700;
}}

QLabel#statusBadgeInfo {{
    background-color: {t.accent_subtle};
    color: {t.accent};
    border: 1px solid {t.accent};
    border-radius: 10px;
    padding: 2px 8px;
    font-size: 10px;
    font-weight: 700;
}}

/* ── Interactive Buttons & Primary Actions ────────────────── */
QPushButton {{
    background-color: {t.bg_input};
    color: {t.text_primary};
    border: 1px solid {t.border};
    border-radius: 6px;
    padding: 6px 14px;
    font-weight: 600;
    font-size: 12px;
    min-height: 24px;
}}

QPushButton:hover {{
    background-color: {t.bg_card_elevated};
    border-color: {t.accent};
    color: {t.text_primary};
}}

QPushButton:pressed {{
    background-color: {t.accent_pressed};
    color: #ffffff;
}}

QPushButton:disabled {{
    background-color: {t.bg_main};
    color: {t.text_disabled};
    border-color: {border_subtle};
}}

QPushButton#accentButton {{
    background-color: {t.accent};
    color: {btn_text};
    border: 1px solid {t.accent};
    font-weight: 700;
    border-radius: 6px;
}}

QPushButton#accentButton:hover {{
    background-color: {t.accent_hover};
    border-color: {t.accent_hover};
}}

QPushButton#accentButton:pressed {{
    background-color: {t.accent_pressed};
}}

QPushButton#dangerButton {{
    background-color: {t.status_fail};
    color: #ffffff;
    border: 1px solid {t.status_fail};
    font-weight: 700;
    border-radius: 6px;
}}

QPushButton#dangerButton:hover {{
    background-color: #dc2626;
    border-color: #b91c1c;
}}

QPushButton#secondaryButton {{
    background-color: {t.bg_card};
    color: {t.text_secondary};
    border: 1px solid {t.border};
    border-radius: 6px;
    font-weight: 600;
}}

QPushButton#secondaryButton:hover {{
    color: {t.text_primary};
    border-color: {t.accent};
}}

/* ── Text Inputs, SpinBoxes & Combos ──────────────────────── */
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background-color: {t.bg_input};
    color: {t.text_primary};
    border: 1px solid {t.border};
    border-radius: 6px;
    padding: 5px 10px;
    min-height: 24px;
    font-size: 12px;
    selection-background-color: {t.selection_bg};
    selection-color: {t.selection_text};
}}

QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover {{
    border-color: {t.accent_hover};
}}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border: 1px solid {t.accent};
    background-color: {t.bg_input};
}}

QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 24px;
    border-left: 1px solid {t.border};
    border-top-right-radius: 6px;
    border-bottom-right-radius: 6px;
}}

QComboBox QAbstractItemView {{
    background-color: {t.bg_card};
    color: {t.text_primary};
    border: 1px solid {t.border};
    selection-background-color: {t.selection_bg};
    selection-color: {t.selection_text};
    border-radius: 6px;
    padding: 4px;
    outline: none;
}}

/* ── Publication-Grade Tables & Tree Views ────────────────── */
QTableWidget, QTreeWidget, QTableView, QTreeView, QListWidget {{
    background-color: {t.bg_card};
    alternate-background-color: {t.row_alt};
    color: {t.text_primary};
    border: 1px solid {t.border};
    border-radius: 6px;
    gridline-color: {border_subtle};
    selection-background-color: {t.selection_bg};
    selection-color: {t.selection_text};
    font-size: 12px;
    outline: none;
}}

QHeaderView::section {{
    background-color: {t.bg_input};
    color: {t.accent};
    font-weight: 700;
    font-size: 11px;
    letter-spacing: 0.4px;
    padding: 7px 10px;
    border: 1px solid {t.border};
    border-top: none;
    border-left: none;
}}

QTableWidget::item {{
    padding: 5px 8px;
}}

QTableWidget::item:selected, QListWidget::item:selected {{
    background-color: {t.selection_bg};
    color: {t.selection_text};
}}

/* ── Scientific Terminal & Logs ───────────────────────────── */
QTextEdit, QPlainTextEdit {{
    background-color: {t.bg_input};
    color: {t.text_primary};
    border: 1px solid {t.border};
    border-radius: 6px;
    font-family: "Consolas", "Cascadia Code", "Fira Code", monospace;
    font-size: 10pt;
    padding: 8px;
    selection-background-color: {t.selection_bg};
    selection-color: {t.selection_text};
}}

QTextEdit#consoleBox {{
    background-color: {'#030712' if t.is_dark else '#0f172a'};
    color: #f1f5f9;
    border: 1px solid {t.border};
    font-family: "Consolas", "Cascadia Code", monospace;
    font-size: 10pt;
    border-radius: 6px;
}}

/* ── Minimal Scrollbars ───────────────────────────────────── */
QScrollBar:vertical {{
    background: {t.bg_main};
    width: 8px;
    margin: 0;
}}

QScrollBar::handle:vertical {{
    background: {t.border};
    border-radius: 4px;
    min-height: 24px;
}}

QScrollBar::handle:vertical:hover {{
    background: {t.accent};
}}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
    background: none;
}}

QScrollBar:horizontal {{
    background: {t.bg_main};
    height: 8px;
    margin: 0;
}}

QScrollBar::handle:horizontal {{
    background: {t.border};
    border-radius: 4px;
    min-width: 24px;
}}

QScrollBar::handle:horizontal:hover {{
    background: {t.accent};
}}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
    background: none;
}}

/* ── Progress Indicators ──────────────────────────────────── */
QProgressBar {{
    background-color: {t.bg_input};
    border: 1px solid {t.border};
    border-radius: 6px;
    text-align: center;
    color: {t.text_primary};
    font-weight: 700;
    font-size: 11px;
    height: 18px;
}}

QProgressBar::chunk {{
    background-color: {t.accent};
    border-radius: 5px;
}}

/* ── Tab Widgets (Used in Sub-Views) ──────────────────────── */
QTabWidget::pane {{
    border: 1px solid {t.border};
    border-radius: 6px;
    background-color: {t.bg_card};
    top: -1px;
}}

QTabBar::tab {{
    background-color: {t.bg_input};
    color: {t.text_secondary};
    border: 1px solid {t.border};
    border-bottom: none;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    min-width: 100px;
    padding: 7px 14px;
    margin-right: 3px;
    font-weight: 700;
    font-size: 11px;
    letter-spacing: 0.2px;
}}

QTabBar::tab:selected {{
    background-color: {t.accent};
    color: {btn_text};
    border-color: {t.accent};
}}

QTabBar::tab:hover:!selected {{
    background-color: {t.bg_card};
    color: {t.text_primary};
    border-color: {t.accent_hover};
}}

/* ── Status Bar ───────────────────────────────────────────── */
QStatusBar {{
    background-color: {t.bg_card};
    border-top: 1px solid {t.border};
    color: {t.text_secondary};
    font-size: 11px;
}}

QStatusBar QLabel {{
    background: transparent;
    color: {t.text_secondary};
    padding: 0 8px;
}}

/* ── Checkboxes & Radio Buttons ───────────────────────────── */
QCheckBox, QRadioButton {{
    spacing: 7px;
    font-size: 12px;
    color: {t.text_primary};
}}

QCheckBox::indicator, QRadioButton::indicator {{
    width: 16px;
    height: 16px;
    border: 1px solid {t.border};
    border-radius: 4px;
    background-color: {t.bg_input};
}}

QCheckBox::indicator:hover, QRadioButton::indicator:hover {{
    border-color: {t.accent};
}}

QCheckBox::indicator:checked, QRadioButton::indicator:checked {{
    background-color: {t.accent};
    border-color: {t.accent};
}}

/* ── Tooltips ─────────────────────────────────────────────── */
QToolTip {{
    background-color: {t.bg_card};
    color: {t.text_primary};
    border: 1px solid {t.accent};
    border-radius: 6px;
    padding: 6px 12px;
    font-size: 11px;
}}
"""
