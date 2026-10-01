"""
AutoDock Suite Pro — Theme & Typography System (gui/themes.py)
==============================================================
Defines four professional, high-clarity color palettes and a centralized
typography scale. Engineered for scientific visual excellence and optimal
readability across high-DPI displays.

Themes:
  Dark        — Deep obsidian & slate surfaces, vibrant cyan/blue accent (default)
  Light       — Clean light slate, navy accent, high contrast
  Noir        — High-contrast OLED/minimalist monochrome
  Scientific  — Pure white, crisp slate text, deep navy accent (journal-ready)
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple
import customtkinter as ctk


# ---------------------------------------------------------------------------
# Typography System (Centralized Font Scale)
# ---------------------------------------------------------------------------

class Typography:
    """Standardized typography scale.
    
    Guarantees readable, balanced font sizes across all monitors.
    No microscopic text (minimum 11pt, standard body 12pt, headings 13-18pt).
    """
    FONT_FAMILY = "Segoe UI"
    CODE_FAMILY = "Consolas"

    # Font tuples for Tkinter / ttk widgets
    TITLE_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 18, "bold")
    HEADER_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 15, "bold")
    SECTION_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 13, "bold")
    SUBSECTION_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 12, "bold")
    BODY_TUPLE: Tuple[str, int] = (FONT_FAMILY, 12)
    BODY_BOLD_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 12, "bold")
    BUTTON_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 12, "bold")
    TAB_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 13, "bold")
    CAPTION_TUPLE: Tuple[str, int] = (FONT_FAMILY, 11)
    CAPTION_BOLD_TUPLE: Tuple[str, int, str] = (FONT_FAMILY, 11, "bold")
    CODE_TUPLE: Tuple[str, int] = (CODE_FAMILY, 11)
    CODE_BOLD_TUPLE: Tuple[str, int, str] = (CODE_FAMILY, 11, "bold")


def font_title() -> ctk.CTkFont:
    """Window & main title font (18pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=18, weight="bold")


def font_header() -> ctk.CTkFont:
    """Major panel or dialog header font (15pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=15, weight="bold")


def font_section() -> ctk.CTkFont:
    """Card or group section header (13pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=13, weight="bold")


def font_body() -> ctk.CTkFont:
    """Primary readable body text for inputs and descriptions (12pt)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=12)


def font_body_bold() -> ctk.CTkFont:
    """Emphasized body text (12pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=12, weight="bold")


def font_button() -> ctk.CTkFont:
    """Interactive button font (12pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=12, weight="bold")


def font_tab() -> ctk.CTkFont:
    """Top-level tab navigation font (13pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=13, weight="bold")


def font_caption() -> ctk.CTkFont:
    """Secondary metadata, hints, small badges, and tooltips (11pt)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=11)


def font_caption_bold() -> ctk.CTkFont:
    """Badges and highlighted small indicators (11pt bold)."""
    return ctk.CTkFont(family=Typography.FONT_FAMILY, size=11, weight="bold")


def font_code() -> ctk.CTkFont:
    """Monospace code, coordinates, formulas, and log entries (11pt)."""
    return ctk.CTkFont(family=Typography.CODE_FAMILY, size=11)


def font_code_bold() -> ctk.CTkFont:
    """Highlighted monospace values (11pt bold)."""
    return ctk.CTkFont(family=Typography.CODE_FAMILY, size=11, weight="bold")


# ---------------------------------------------------------------------------
# Theme Dataclass
# ---------------------------------------------------------------------------

@dataclass
class ThemeColors:
    """Complete color token set for a theme."""
    name: str

    # Core surfaces
    bg_primary: str       # Main window background
    bg_secondary: str     # Cards, panels, frames
    bg_tertiary: str      # Hover states, subtle dividers, input fields

    # Text
    text_primary: str     # Main readable text
    text_secondary: str   # Labels, captions, hints
    text_disabled: str    # Greyed-out content

    # Accent (modern biophysical cyan/blues)
    accent: str           # Interactive elements, selected state
    accent_hover: str     # Hover on accent elements

    # Borders & separators
    border: str

    # Status indicators
    status_ok_text: str       # Text color for OK states
    status_warn_text: str     # Text color for warnings
    status_fail_text: str     # Text color for errors

    # Treeview / table row alternation
    row_even: str
    row_odd: str
    row_selected: str

    # Console / log area
    console_bg: str
    console_text: str

    # CTk appearance mode string ("dark" / "light")
    ctk_mode: str

    # Button states
    btn_fg: str
    btn_hover: str
    btn_text: str


# ---------------------------------------------------------------------------
# Theme Definitions
# ---------------------------------------------------------------------------

DARK = ThemeColors(
    name="Dark",
    bg_primary="#0f172a",       # Deep slate 900
    bg_secondary="#1e293b",     # Slate 800 card surface
    bg_tertiary="#334155",      # Slate 700 input / hover
    text_primary="#f8fafc",     # Slate 50 crisp high-contrast text
    text_secondary="#94a3b8",   # Slate 400 clear subtitle
    text_disabled="#64748b",    # Slate 500 disabled text
    accent="#38bdf8",           # Sky 400 — brighter, more vivid biophysical blue
    accent_hover="#0ea5e9",     # Sky 500 hover state
    border="#334155",           # Crisp boundary
    status_ok_text="#34d399",   # Emerald 400
    status_warn_text="#fbbf24", # Amber 400
    status_fail_text="#f87171", # Coral 400
    row_even="#1e293b",
    row_odd="#182234",
    row_selected="#0369a1",
    console_bg="#090d16",
    console_text="#e2e8f0",
    ctk_mode="dark",
    btn_fg="#0284c7",
    btn_hover="#0369a1",
    btn_text="#ffffff",
)

LIGHT = ThemeColors(
    name="Light",
    bg_primary="#f8fafc",       # Slate 50 soft background
    bg_secondary="#ffffff",     # Pure white elevated card
    bg_tertiary="#f1f5f9",      # Slate 100 input / hover
    text_primary="#0f172a",     # Slate 900 high contrast
    text_secondary="#475569",   # Slate 600 clear text
    text_disabled="#94a3b8",    # Slate 400
    accent="#0284c7",           # Vibrant biophysical blue
    accent_hover="#0369a1",     # Darker blue hover
    border="#cbd5e1",           # Slate 300 clean border
    status_ok_text="#059669",   # Emerald 600
    status_warn_text="#d97706", # Amber 600
    status_fail_text="#dc2626", # Red 600
    row_even="#ffffff",
    row_odd="#f8fafc",
    row_selected="#bae6fd",
    console_bg="#0f172a",
    console_text="#f1f5f9",
    ctk_mode="light",
    btn_fg="#0284c7",
    btn_hover="#0369a1",
    btn_text="#ffffff",
)

NOIR = ThemeColors(
    name="Noir",
    bg_primary="#000000",
    bg_secondary="#111111",
    bg_tertiary="#222222",
    text_primary="#ffffff",     # High-contrast pure white
    text_secondary="#a1a1aa",   # Zinc 400
    text_disabled="#52525b",    # Zinc 600
    accent="#e4e4e7",           # Zinc 200
    accent_hover="#d4d4d8",
    border="#27272a",
    status_ok_text="#ffffff",
    status_warn_text="#d4d4d8",
    status_fail_text="#a1a1aa",
    row_even="#080808",
    row_odd="#141414",
    row_selected="#2a2a2a",
    console_bg="#000000",
    console_text="#e4e4e7",
    ctk_mode="dark",
    btn_fg="#27272a",
    btn_hover="#3f3f46",
    btn_text="#ffffff",
)

SCIENTIFIC = ThemeColors(
    name="Scientific",
    bg_primary="#ffffff",
    bg_secondary="#f8fafc",
    bg_tertiary="#f1f5f9",
    text_primary="#1e293b",
    text_secondary="#475569",
    text_disabled="#94a3b8",
    accent="#0369a1",           # Academic journal navy blue
    accent_hover="#075985",
    border="#e2e8f0",
    status_ok_text="#047857",
    status_warn_text="#b45309",
    status_fail_text="#b91c1c",
    row_even="#ffffff",
    row_odd="#f8fafc",
    row_selected="#e0f2fe",
    console_bg="#0b1120",
    console_text="#f8fafc",
    ctk_mode="light",
    btn_fg="#0369a1",
    btn_hover="#075985",
    btn_text="#ffffff",
)

PLASMA = ThemeColors(
    name="Plasma",
    # Deep violet near-black inspired by VMD & Chimera molecular visualization
    bg_primary="#0d0221",       # Near-black violet
    bg_secondary="#1a0533",     # Deep purple card surface
    bg_tertiary="#2d1154",      # Medium purple inputs / hover
    text_primary="#f0e6ff",     # Lavender-white main text
    text_secondary="#b388ff",   # Violet 300 subtitle
    text_disabled="#6a3fa0",    # Muted deep violet
    accent="#e040fb",           # Electric magenta accent (Material AM1)
    accent_hover="#ce93d8",     # Softer purple hover
    border="#4a1a7a",           # Violet boundary
    status_ok_text="#69f0ae",   # Green A200
    status_warn_text="#ffeb3b", # Yellow A200
    status_fail_text="#ff5252", # Red A200
    row_even="#1a0533",
    row_odd="#150427",
    row_selected="#4a148c",
    console_bg="#060011",
    console_text="#e8d5ff",
    ctk_mode="dark",
    btn_fg="#7b1fa2",
    btn_hover="#9c27b0",
    btn_text="#f0e6ff",
)

BIO_NEON = ThemeColors(
    name="Bio-Neon",
    # Dark forest-green substrate with neon emerald accents —
    # inspired by bioinformatics terminal / sequencer display aesthetics
    bg_primary="#001a12",       # Near-black forest green
    bg_secondary="#00291c",     # Dark green card surface
    bg_tertiary="#003d28",      # Input / hover background
    text_primary="#ccffe6",     # Mint-white main text
    text_secondary="#66cc99",   # Emerald 400 subtitle
    text_disabled="#2d8a5e",    # Muted mid-green
    accent="#00e676",           # Neon green A400
    accent_hover="#00c853",     # Green A700 hover
    border="#005c3a",           # Forest boundary
    status_ok_text="#69f0ae",   # Green A200
    status_warn_text="#ffeb3b", # Bright amber warning
    status_fail_text="#ff5252", # Red fail
    row_even="#00291c",
    row_odd="#00221a",
    row_selected="#00701a",
    console_bg="#000d08",
    console_text="#b3ffda",
    ctk_mode="dark",
    btn_fg="#00701a",
    btn_hover="#00c853",
    btn_text="#ccffe6",
)

ALL_THEMES: Dict[str, ThemeColors] = {
    "Dark": DARK,
    "Light": LIGHT,
    "Noir": NOIR,
    "Scientific": SCIENTIFIC,
    "Plasma": PLASMA,
    "Bio-Neon": BIO_NEON,
}
THEME_NAMES = list(ALL_THEMES.keys())
DEFAULT_THEME = "Dark"


# ---------------------------------------------------------------------------
# Active Theme Registry
# ---------------------------------------------------------------------------

_active: ThemeColors = DARK


def get() -> ThemeColors:
    """Returns the currently active theme."""
    return _active


def set_theme(name: str) -> ThemeColors:
    """Switch to a named theme. Returns the new ThemeColors."""
    global _active
    theme = ALL_THEMES.get(name, DARK)
    _active = theme
    ctk.set_appearance_mode(theme.ctk_mode)
    return theme
