"""
AutoDock Suite Pro — Main Application Window (gui/app.py)
==========================================================
Clean, modular replacement for the legacy gui.py God Object.
The window assembles all tab modules and handles:
  - Scientific branding & application icon
  - Standardized typography scale (no microscopic text)
  - Theme switching (live, no restart)
  - Config persistence
  - Cross-tab communication
  - Status & engine monitoring

Tab order (matches scientific workflow):
  1. Workspace   — Receptor & Ligand management
  2. Preparation — Receptor / Ligand PDBQT preparation
  3. Grid Box    — Vina config / AutoDock4 GPF
  4. Library     — Screening library + ADMET filter
  5. Run Docking — Execution console
  6. Results     — Poses, interactions, 2D structure & PyMOL export

Author: AutoDock Suite Pro Development Team
Version: 0.3.0
"""

from __future__ import annotations

import logging
import math
import os
import sys
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import List, Optional, Tuple

import customtkinter as ctk
from PIL import Image

from config import ProjectConfig, load_config
from models import SUITE_NAME, __version__
from gui import themes
from gui.themes import ALL_THEMES, THEME_NAMES, set_theme

logger = logging.getLogger("docking_automation.gui")


# ---------------------------------------------------------------------------
# Ki helper — used by results_tab via direct import
# ---------------------------------------------------------------------------

def calculate_inhibition_constant(
    delta_g_kcal_mol: float, temp_kelvin: float = 298.15
) -> Tuple[Optional[float], str]:
    """Thermodynamic inference from Vina affinity is strictly prohibited in ADSP.
    AutoDock4 Ki values must be DLG-driven only.
    """
    raise NotImplementedError(
        "Thermodynamic inference is not supported: Vina does not report Ki, "
        "and AutoDock4 Ki must be parsed directly from the DLG."
    )


# ---------------------------------------------------------------------------
# Main Window
# ---------------------------------------------------------------------------

class AutoDockSuiteProApp(ctk.CTk):
    """AutoDock Suite Pro — top-level application window."""

    TAB_NAMES = [
        "  🗂  Workspace  ",
        "  🧪  Preparation  ",
        "  🎯  Grid Box  ",
        "  📚  Library  ",
        "  ▶  Run Docking  ",
        "  📊  Results  ",
    ]

    def __init__(self, config: Optional[ProjectConfig] = None) -> None:
        super().__init__()
        self.title(f"{SUITE_NAME} v{__version__} — Molecular Docking & Virtual Screening")
        self.geometry("1440x960")
        self.minsize(1200, 820)

        # Application Icon
        self._assets_dir = Path(__file__).parent / "assets"
        ico_path = self._assets_dir / "app_icon.ico"
        if ico_path.is_file():
            try:
                self.iconbitmap(str(ico_path))
            except Exception as e:
                logger.debug(f"iconbitmap failed: {e}")

        # Preload Logo Image
        self._logo_image = None
        logo_path = self._assets_dir / "logo_36.png"
        if not logo_path.is_file():
            logo_path = self._assets_dir / "logo_48.png"
        if logo_path.is_file():
            try:
                pil_logo = Image.open(logo_path)
                self._logo_image = ctk.CTkImage(
                    light_image=pil_logo,
                    dark_image=pil_logo,
                    size=(36, 36),
                )
            except Exception as e:
                logger.debug(f"Could not load logo: {e}")

        self._config = config if config else self._get_default_config()

        # Apply saved theme first
        saved_theme = getattr(self._config, "ui_theme", "Noir")
        set_theme(saved_theme)

        self._build_ui()
        self._apply_window_colors()
        self.after(100, self._tabs["  🗂  Workspace  "].refresh)

    # ─────────────────────────────────────────────────
    # UI Construction
    # ─────────────────────────────────────────────────

    def _build_ui(self) -> None:
        t = themes.get()
        self.configure(fg_color=t.bg_primary)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # ── Title bar ────────────────────────────────
        self._build_titlebar()

        # ── Tab view ─────────────────────────────────
        self.tabview = ctk.CTkTabview(
            self,
            fg_color=t.bg_secondary,
            segmented_button_fg_color=t.bg_tertiary,
            segmented_button_selected_color=t.accent,
            segmented_button_selected_hover_color=t.accent_hover,
            segmented_button_unselected_color=t.bg_tertiary,
            segmented_button_unselected_hover_color=t.bg_tertiary,
            text_color=t.text_secondary,
            text_color_disabled=t.text_disabled,
            border_width=0,
        )
        self.tabview.grid(row=1, column=0, sticky="nsew", padx=10, pady=(0, 6))

        for name in self.TAB_NAMES:
            self.tabview.add(name)

        # Apply standardized tab typography (13pt bold)
        try:
            self.tabview._segmented_button.configure(font=themes.font_tab())
        except Exception:
            pass

        self._tabs = {}
        self._build_tabs()

        # ── Status bar ───────────────────────────────
        self._build_statusbar()

    def _build_titlebar(self) -> None:
        t = themes.get()
        self.titlebar = ctk.CTkFrame(self, fg_color=t.bg_secondary, height=60, corner_radius=0)
        self.titlebar.grid(row=0, column=0, sticky="ew")
        self.titlebar.grid_columnconfigure(1, weight=1)

        # Left branding block: Logo + App Name + Subtitle
        brand_frame = ctk.CTkFrame(self.titlebar, fg_color="transparent")
        brand_frame.grid(row=0, column=0, padx=(14, 8), pady=8, sticky="w")

        if self._logo_image:
            logo_lbl = ctk.CTkLabel(brand_frame, image=self._logo_image, text="")
            logo_lbl.pack(side="left", padx=(0, 10))

        text_block = ctk.CTkFrame(brand_frame, fg_color="transparent")
        text_block.pack(side="left")

        self.lbl_title = ctk.CTkLabel(
            text_block,
            text=SUITE_NAME,
            font=themes.font_title(),
            text_color=t.text_primary,
            anchor="w",
        )
        self.lbl_title.pack(anchor="w")

        self.lbl_subtitle = ctk.CTkLabel(
            text_block,
            text=f"v{__version__}  ·  AutoDock Vina 1.2 & AutoDock 4.2  ·  Biophysical Virtual Screening",
            font=themes.font_caption(),
            text_color=t.text_secondary,
            anchor="w",
        )
        self.lbl_subtitle.pack(anchor="w")

        # Right actions block: Theme + Save/Open Project
        actions_frame = ctk.CTkFrame(self.titlebar, fg_color="transparent")
        actions_frame.grid(row=0, column=2, padx=14, pady=8, sticky="e")

        self.lbl_theme_prompt = ctk.CTkLabel(
            actions_frame,
            text="Theme:",
            text_color=t.text_secondary,
            font=themes.font_body(),
        )
        self.lbl_theme_prompt.pack(side="left", padx=(0, 6))

        self.cmb_theme = ctk.CTkOptionMenu(
            actions_frame,
            values=THEME_NAMES,
            width=120,
            height=32,
            fg_color=t.bg_tertiary,
            button_color=t.accent,
            text_color=t.text_primary,
            font=themes.font_body(),
            dropdown_font=themes.font_body(),
            command=self._change_theme,
        )
        self.cmb_theme.set(themes.get().name)
        self.cmb_theme.pack(side="left", padx=(0, 10))

        self.btn_save = ctk.CTkButton(
            actions_frame,
            text="💾 Save",
            command=self._save_config,
            width=76,
            height=32,
            fg_color="transparent",
            border_width=1,
            border_color=t.border,
            text_color=t.text_secondary,
            font=themes.font_button(),
        )
        self.btn_save.pack(side="left", padx=4)

        self.btn_open = ctk.CTkButton(
            actions_frame,
            text="📂 Open",
            command=self._load_config,
            width=76,
            height=32,
            fg_color="transparent",
            border_width=1,
            border_color=t.border,
            text_color=t.text_secondary,
            font=themes.font_button(),
        )
        self.btn_open.pack(side="left", padx=4)

        # ── Animated gradient accent strip (3px at bottom of titlebar) ──────
        import tkinter as tk
        self._accent_canvas = tk.Canvas(
            self.titlebar, height=3, highlightthickness=0, bd=0,
        )
        self._accent_canvas.grid(row=1, column=0, columnspan=3, sticky="ew")
        self._gradient_phase = 0.0
        self._animate_gradient()

    def _animate_gradient(self) -> None:
        """Pulse a sweeping gradient on the accent strip below the titlebar."""
        import tkinter as tk
        try:
            canvas = self._accent_canvas
            t = themes.get()
            w = self.winfo_width() or 1440
            canvas.configure(bg=t.bg_secondary, width=w)
            canvas.delete("all")

            # Sweep position: 0..1 sinusoid over time
            import math
            self._gradient_phase = (self._gradient_phase + 0.015) % (2 * math.pi)
            pos = (math.sin(self._gradient_phase) + 1) / 2  # 0..1

            # Parse accent hex to RGB
            acc = t.accent.lstrip("#")
            r1, g1, b1 = int(acc[0:2], 16), int(acc[2:4], 16), int(acc[4:6], 16)
            # Complementary: blend with bg_primary
            bg = t.bg_primary.lstrip("#")
            r2, g2, b2 = int(bg[0:2], 16), int(bg[2:4], 16), int(bg[4:6], 16)

            # Draw horizontal gradient segments
            segments = max(20, w // 6)
            seg_w = max(1, w // segments)
            for i in range(segments):
                frac = i / segments
                mix = max(0.0, min(1.0, 1.0 - abs(frac - pos) * 3))
                r = int(r1 * mix + r2 * (1 - mix))
                g = int(g1 * mix + g2 * (1 - mix))
                b = int(b1 * mix + b2 * (1 - mix))
                color = f"#{r:02x}{g:02x}{b:02x}"
                x0, x1 = i * seg_w, min((i + 1) * seg_w + 1, w)
                canvas.create_rectangle(x0, 0, x1, 3, fill=color, outline="")

            self.after(40, self._animate_gradient)   # ~25 fps, lightweight
        except Exception:
            pass  # Window destroyed or canvas gone

    def _build_tabs(self) -> None:
        """Instantiate each tab module and mount it inside its CTkTabview frame."""
        from gui.workspace_tab import WorkspaceTab
        from gui.library_tab import LibraryTab
        from gui.runner_tab import RunnerTab
        from gui.results_tab import ResultsTab
        from gui.preparation_tab import PreparationTab
        from gui.grid_tab import GridBoxTab

        tab_frames = {name: self.tabview.tab(name) for name in self.TAB_NAMES}

        # Workspace
        ws_tab = WorkspaceTab(
            tab_frames["  🗂  Workspace  "],
            config=self._config,
            on_config_change=self._on_config_change,
        )
        ws_tab.pack(fill="both", expand=True)
        self._tabs["  🗂  Workspace  "] = ws_tab

        # Preparation
        prep_tab = PreparationTab(
            tab_frames["  🧪  Preparation  "],
            config=self._config,
            on_config_change=self._on_config_change,
        )
        prep_tab.pack(fill="both", expand=True)
        self._tabs["  🧪  Preparation  "] = prep_tab

        # Grid Box
        grid_tab = GridBoxTab(
            tab_frames["  🎯  Grid Box  "],
            config=self._config,
            on_config_change=self._on_config_change,
        )
        grid_tab.pack(fill="both", expand=True)
        self._tabs["  🎯  Grid Box  "] = grid_tab

        # Library
        lib_tab = LibraryTab(
            tab_frames["  📚  Library  "],
            config=self._config,
            on_config_change=self._on_config_change,
        )
        lib_tab.pack(fill="both", expand=True)
        self._tabs["  📚  Library  "] = lib_tab

        # Runner
        run_tab = RunnerTab(
            tab_frames["  ▶  Run Docking  "],
            config=self._config,
            on_config_change=self._on_config_change,
            on_run_complete=self._on_run_complete,
        )
        run_tab.pack(fill="both", expand=True)
        self._tabs["  ▶  Run Docking  "] = run_tab

        # Results
        res_tab = ResultsTab(
            tab_frames["  📊  Results  "],
            config=self._config,
            on_config_change=self._on_config_change,
        )
        res_tab.pack(fill="both", expand=True)
        self._tabs["  📊  Results  "] = res_tab

    def _build_statusbar(self) -> None:
        t = themes.get()
        self.statusbar = ctk.CTkFrame(self, fg_color=t.bg_secondary, height=28, corner_radius=0)
        self.statusbar.grid(row=2, column=0, sticky="ew")
        self.statusbar.grid_columnconfigure(1, weight=1)

        self.lbl_status = ctk.CTkLabel(
            self.statusbar,
            text="  Ready",
            text_color=t.text_secondary,
            font=themes.font_caption(),
            anchor="w",
        )
        self.lbl_status.grid(row=0, column=0, sticky="w", padx=12, pady=4)

        self.lbl_engine_status = ctk.CTkLabel(
            self.statusbar,
            text="",
            text_color=t.text_secondary,
            font=themes.font_caption(),
            anchor="e",
        )
        self.lbl_engine_status.grid(row=0, column=1, sticky="e", padx=8, pady=4)

        self.lbl_clock = ctk.CTkLabel(
            self.statusbar,
            text="",
            text_color=t.text_disabled,
            font=themes.font_caption(),
            anchor="e",
        )
        self.lbl_clock.grid(row=0, column=2, sticky="e", padx=12, pady=4)

        self.after(500, self._update_engine_status)
        self._tick_clock()

    # ─────────────────────────────────────────────────
    # Theme switching (live, no restart)
    # ─────────────────────────────────────────────────

    def _change_theme(self, name: str) -> None:
        set_theme(name)
        setattr(self._config, "ui_theme", name)
        self._apply_window_colors()
        self._on_config_change()
        # Save silently
        try:
            cfg_path = self._config.project_root / "project_config.toml"
            if cfg_path.is_file():
                from gui.settings import save_config_to_toml
                save_config_to_toml(self._config, cfg_path)
        except Exception:
            pass

    def _apply_window_colors(self) -> None:
        t = themes.get()
        self.configure(fg_color=t.bg_primary)
        if hasattr(self, "titlebar"):
            self.titlebar.configure(fg_color=t.bg_secondary)
            self.lbl_title.configure(text_color=t.text_primary)
            self.lbl_subtitle.configure(text_color=t.text_secondary)
            self.lbl_theme_prompt.configure(text_color=t.text_secondary)
            self.cmb_theme.configure(
                fg_color=t.bg_tertiary,
                button_color=t.accent,
                button_hover_color=t.accent_hover,
                text_color=t.text_primary,
            )
            self.btn_save.configure(border_color=t.border, text_color=t.text_secondary)
            self.btn_open.configure(border_color=t.border, text_color=t.text_secondary)

        if hasattr(self, "statusbar"):
            self.statusbar.configure(fg_color=t.bg_secondary)
            self.lbl_status.configure(text_color=t.text_secondary)
            self.lbl_engine_status.configure(text_color=t.text_secondary)
            if hasattr(self, "lbl_clock"):
                self.lbl_clock.configure(text_color=t.text_disabled)

        if hasattr(self, "tabview"):
            self.tabview.configure(
                fg_color=t.bg_secondary,
                segmented_button_fg_color=t.bg_tertiary,
                segmented_button_selected_color=t.accent,
                segmented_button_selected_hover_color=t.accent_hover,
                segmented_button_unselected_color=t.bg_tertiary,
                segmented_button_unselected_hover_color=t.bg_tertiary,
                text_color=t.text_secondary,
            )
            try:
                self.tabview._segmented_button.configure(
                    selected_color=t.accent,
                    selected_hover_color=t.accent_hover,
                    unselected_color=t.bg_tertiary,
                    unselected_hover_color=t.bg_tertiary,
                    text_color=t.text_primary,
                )
            except Exception:
                pass

    # ─────────────────────────────────────────────────
    # Config persistence
    # ─────────────────────────────────────────────────

    def _get_default_config(self) -> ProjectConfig:
        toml_path = Path("project_config.toml")
        if toml_path.is_file():
            try:
                return load_config(toml_path)
            except Exception:
                pass
        return ProjectConfig()

    def _save_config(self) -> None:
        cfg_path = self._config.project_root / "project_config.toml"
        try:
            from gui.settings import save_config_to_toml
            save_config_to_toml(self._config, cfg_path)
            self._set_status(f"Saved: {cfg_path.name}")
        except Exception as e:
            messagebox.showerror("Save Failed", str(e))

    def _load_config(self) -> None:
        path = filedialog.askopenfilename(
            title="Open project_config.toml",
            filetypes=[("TOML Config", "*.toml"), ("All Files", "*.*")],
        )
        if not path:
            return
        try:
            new_config = load_config(Path(path))
            self._config = new_config
            self._on_config_change()
            self._set_status(f"Loaded: {Path(path).name}")
        except Exception as e:
            messagebox.showerror("Load Failed", str(e))

    def _on_config_change(self) -> None:
        """Propagate config changes to all tab modules that hold a reference."""
        for tab in self._tabs.values():
            if hasattr(tab, "config"):
                tab.config = self._config
            if hasattr(tab, "refresh"):
                try:
                    tab.refresh()
                except Exception:
                    pass

    # ─────────────────────────────────────────────────
    # Run pipeline callbacks
    # ─────────────────────────────────────────────────

    def _on_run_complete(self, jobs: list) -> None:
        """Called when the docking pipeline finishes — auto-refresh results."""
        self._set_status(f"Docking complete — {len(jobs)} job(s) processed.")
        res_tab = self._tabs.get("  📊  Results  ")
        if res_tab and hasattr(res_tab, "refresh"):
            res_tab.refresh()
        # Switch to Results tab automatically
        try:
            self.tabview.set("  📊  Results  ")
        except Exception:
            pass

    # ─────────────────────────────────────────────────
    # Status bar helpers
    # ─────────────────────────────────────────────────

    def _set_status(self, msg: str) -> None:
        self.lbl_status.configure(text=f"  {msg}")

    def _update_engine_status(self) -> None:
        try:
            from executables import get_vina_version, get_autodock4_version
            vina_ver = get_vina_version(self._config.vina_executable)
            ad4_ver = get_autodock4_version(self._config.autodock4_executable)
            self.lbl_engine_status.configure(
                text=f"Vina {vina_ver}  |  AutoDock4 {ad4_ver}  "
            )
        except Exception:
            pass

    def _tick_clock(self) -> None:
        try:
            from datetime import datetime
            now_str = datetime.now().strftime("%H:%M:%S")
            if hasattr(self, "lbl_clock"):
                self.lbl_clock.configure(text=f"{now_str}  ")
            self.after(1000, self._tick_clock)
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_gui(config: Optional[ProjectConfig] = None) -> None:
    """Launch the AutoDock Suite Pro GUI application."""
    app = AutoDockSuiteProApp(config=config)
    app.mainloop()
