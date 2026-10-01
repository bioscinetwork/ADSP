"""
AutoDock Suite Pro — Runner Tab (gui/runner_tab.py)
====================================================
Docking console with per-job progress reporting, engine selection,
full parameter control panel (Vina & AutoDock4/AutoDockTools options),
and live log streaming.
"""

from __future__ import annotations

import os
import queue
import subprocess
import sys
import threading
from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Optional

import customtkinter as ctk

from gui import themes, widgets
from gui.widgets import TextRedirector

if TYPE_CHECKING:
    from config import ProjectConfig
    from models import DockingJob


class RunnerTab(ctk.CTkFrame):
    """Docking execution console tab with comprehensive parameter controls."""

    def __init__(self, parent: ctk.CTkTabview, config: "ProjectConfig",
                 on_config_change: Callable,
                 on_run_complete: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self.on_run_complete = on_run_complete
        self._running = False
        self._jobs_completed = 0
        self._jobs_total = 0
        self._params_visible = True
        self._build()

    # ─────────────────────────────────────────────────
    # Build
    # ─────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(4, weight=1)

        widgets.section_header(self, "Docking Console", row=0)

        # Controls row
        ctrl = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        ctrl.grid(row=1, column=0, sticky="ew", padx=12, pady=6)
        ctrl.grid_columnconfigure(5, weight=1)

        ctk.CTkLabel(ctrl, text="Engine:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11, weight="bold")).grid(row=0, column=0, padx=(10, 4), pady=8)
        engine_val = self.config.engine.value if hasattr(self.config.engine, "value") else str(self.config.engine)
        self.cmb_engine = ctk.CTkOptionMenu(
            ctrl, values=["VINA", "AUTODOCK4", "BOTH"],
            width=120, fg_color=t.bg_tertiary, button_color=t.accent,
            text_color=t.text_primary, command=self._on_engine_change,
        )
        if engine_val in ["VINA", "AUTODOCK4", "BOTH"]:
            self.cmb_engine.set(engine_val)
        self.cmb_engine.grid(row=0, column=1, padx=4, pady=8)

        ctk.CTkLabel(ctrl, text="Mode:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11, weight="bold")).grid(row=0, column=2, padx=(12, 4))
        mode_val = self.config.docking_mode.value if hasattr(self.config.docking_mode, "value") else str(self.config.docking_mode)
        self.cmb_mode = ctk.CTkOptionMenu(
            ctrl, values=["AUTO", "RIGID", "FLEXIBLE", "BOTH"],
            width=110, fg_color=t.bg_tertiary, button_color=t.accent,
            text_color=t.text_primary,
        )
        if mode_val in ["AUTO", "RIGID", "FLEXIBLE", "BOTH"]:
            self.cmb_mode.set(mode_val)
        self.cmb_mode.grid(row=0, column=3, padx=4)

        self.btn_toggle_params = ctk.CTkButton(
            ctrl, text="⚙ Parameters ▲", command=self._toggle_params,
            width=125, fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.accent_hover, font=ctk.CTkFont(size=11),
        )
        self.btn_toggle_params.grid(row=0, column=4, padx=6)

        self.btn_start = widgets.accent_button(ctrl, "▶  Start Docking",
                                               self._start_docking, width=140)
        self.btn_start.grid(row=0, column=6, padx=6, pady=8)

        self.btn_stop = ctk.CTkButton(
            ctrl, text="⏹  Stop", command=self._stop_docking,
            width=80, state="disabled",
            fg_color=t.bg_tertiary, text_color=t.text_disabled,
        )
        self.btn_stop.grid(row=0, column=7, padx=(4, 10), pady=8)

        # ── Parameter Panel (collapsible) ───────────────
        self.params_panel = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        self.params_panel.grid(row=2, column=0, sticky="ew", padx=12, pady=(0, 6))
        self.params_panel.grid_columnconfigure(0, weight=1)

        # Panel header
        panel_hdr = ctk.CTkFrame(self.params_panel, fg_color="transparent")
        panel_hdr.grid(row=0, column=0, sticky="ew", padx=10, pady=(6, 2))
        panel_hdr.grid_columnconfigure(0, weight=1)

        self.lbl_panel_title = ctk.CTkLabel(
            panel_hdr, text="⚙ DOCKING PARAMETERS",
            font=themes.font_section(),
            text_color=t.accent, anchor="w",
        )
        self.lbl_panel_title.grid(row=0, column=0, sticky="w")

        ctk.CTkButton(
            panel_hdr, text="↺ Reset Defaults", command=self._reset_params_to_defaults,
            width=120, height=28, fg_color=t.bg_tertiary, text_color=t.text_secondary,
            hover_color=t.border, font=themes.font_button(),
        ).grid(row=0, column=1, sticky="e")

        # Container for engine-specific frames
        self._build_vina_panel()
        self._build_ad4_panel()
        self._update_param_view()

        # ── Progress bar + label ────────────────────────
        prog_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        prog_frame.grid(row=3, column=0, sticky="ew", padx=12, pady=(0, 6))
        prog_frame.grid_columnconfigure(0, weight=1)

        self.lbl_progress = ctk.CTkLabel(
            prog_frame, text="Ready — configure workspace & parameters, then click Start Docking",
            text_color=t.text_secondary, font=themes.font_body(), anchor="w",
        )
        self.lbl_progress.grid(row=0, column=0, sticky="w", padx=10, pady=(6, 2))

        self.progress_bar = ctk.CTkProgressBar(
            prog_frame, height=8, fg_color=t.bg_tertiary, progress_color=t.accent,
        )
        self.progress_bar.set(0)
        self.progress_bar.grid(row=1, column=0, sticky="ew", padx=10, pady=(2, 8))

        # ── Console output ──────────────────────────────
        console_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        console_frame.grid(row=4, column=0, sticky="nsew", padx=12, pady=(0, 8))
        console_frame.grid_rowconfigure(0, weight=1)
        console_frame.grid_columnconfigure(0, weight=1)

        self._console_font_size = 13
        self.console = ctk.CTkTextbox(
            console_frame, state="disabled",
            fg_color=t.console_bg, text_color=t.console_text,
            font=ctk.CTkFont(family="Consolas", size=self._console_font_size),
            wrap="word",
        )
        self.console.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)

        console_tb = ctk.CTkFrame(console_frame, fg_color="transparent")
        console_tb.grid(row=1, column=0, sticky="ew", padx=6, pady=4)
        console_tb.grid_columnconfigure(3, weight=1)

        # Left side: Font size zoom controls
        ctk.CTkLabel(console_tb, text="Font Size:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=(2, 4), sticky="w")

        ctk.CTkButton(
            console_tb, text="A -", width=36, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=ctk.CTkFont(size=11, weight="bold"),
            command=self._zoom_out,
        ).grid(row=0, column=1, padx=2, sticky="w")

        self.lbl_font_size = ctk.CTkLabel(
            console_tb, text=f"{self._console_font_size} pt", text_color=t.text_primary,
            font=ctk.CTkFont(size=11, weight="bold"), width=44,
        )
        self.lbl_font_size.grid(row=0, column=2, padx=2, sticky="w")

        ctk.CTkButton(
            console_tb, text="A +", width=36, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=ctk.CTkFont(size=11, weight="bold"),
            command=self._zoom_in,
        ).grid(row=0, column=3, padx=2, sticky="w")

        # Right side: Folder quick access and clear
        actions_frame = ctk.CTkFrame(console_tb, fg_color="transparent")
        actions_frame.grid(row=0, column=4, sticky="e")

        ctk.CTkButton(
            actions_frame, text="📁 Open DLG Folder", width=130, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=ctk.CTkFont(size=11),
            command=self._open_dlg_folder,
        ).pack(side="left", padx=4)

        ctk.CTkButton(
            actions_frame, text="📁 Open Results", width=115, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=ctk.CTkFont(size=11),
            command=self._open_results_folder,
        ).pack(side="left", padx=4)

        ctk.CTkButton(
            actions_frame, text="🗑 Clear Log", width=90, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_secondary,
            hover_color=t.border, font=ctk.CTkFont(size=11),
            command=self._clear_console,
        ).pack(side="left", padx=4)

    # ─────────────────────────────────────────────────
    # Vina Parameter Panel
    # ─────────────────────────────────────────────────

    def _build_vina_panel(self) -> None:
        t = themes.get()
        self.vina_frame = ctk.CTkFrame(self.params_panel, fg_color="transparent")
        for c in range(8):
            self.vina_frame.grid_columnconfigure(c, weight=1 if c in (1, 3, 5, 7) else 0)

        # Row 0: Scoring function, Exhaustiveness, Num modes, Energy range
        ctk.CTkLabel(self.vina_frame, text="Scoring Function:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=(10, 4), pady=4, sticky="w")
        self.cmb_vina_scoring = ctk.CTkOptionMenu(
            self.vina_frame, values=["vina", "vinardo", "ad4"],
            width=90, fg_color=t.bg_tertiary, button_color=t.accent,
            text_color=t.text_primary,
        )
        self.cmb_vina_scoring.set(getattr(self.config, "vina_scoring", "vina"))
        self.cmb_vina_scoring.grid(row=0, column=1, padx=4, pady=4, sticky="ew")

        ctk.CTkLabel(self.vina_frame, text="Exhaustiveness:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=2, padx=(10, 4), pady=4, sticky="w")
        self.ent_vina_exhaust = ctk.CTkEntry(self.vina_frame, width=60, fg_color=t.bg_tertiary,
                                             border_color=t.border, text_color=t.text_primary)
        self.ent_vina_exhaust.insert(0, str(getattr(self.config, "exhaustiveness", 8)))
        self.ent_vina_exhaust.grid(row=0, column=3, padx=4, pady=4, sticky="ew")

        ctk.CTkLabel(self.vina_frame, text="Num Modes:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=4, padx=(10, 4), pady=4, sticky="w")
        self.ent_vina_modes = ctk.CTkEntry(self.vina_frame, width=60, fg_color=t.bg_tertiary,
                                           border_color=t.border, text_color=t.text_primary)
        self.ent_vina_modes.insert(0, str(getattr(self.config, "num_modes", 9)))
        self.ent_vina_modes.grid(row=0, column=5, padx=4, pady=4, sticky="ew")

        ctk.CTkLabel(self.vina_frame, text="Energy Range (kcal):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=6, padx=(10, 4), pady=4, sticky="w")
        self.ent_vina_energy = ctk.CTkEntry(self.vina_frame, width=60, fg_color=t.bg_tertiary,
                                            border_color=t.border, text_color=t.text_primary)
        self.ent_vina_energy.insert(0, str(getattr(self.config, "energy_range", 3.0)))
        self.ent_vina_energy.grid(row=0, column=7, padx=(4, 10), pady=4, sticky="ew")

        # Row 1: Min RMSD, CPUs, Seed
        ctk.CTkLabel(self.vina_frame, text="Min RMSD (Å):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=0, padx=(10, 4), pady=(0, 6), sticky="w")
        self.ent_vina_min_rmsd = ctk.CTkEntry(self.vina_frame, width=90, fg_color=t.bg_tertiary,
                                              border_color=t.border, text_color=t.text_primary)
        self.ent_vina_min_rmsd.insert(0, str(getattr(self.config, "vina_min_rmsd", 1.0)))
        self.ent_vina_min_rmsd.grid(row=1, column=1, padx=4, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(self.vina_frame, text="CPU Workers:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=2, padx=(10, 4), pady=(0, 6), sticky="w")
        self.ent_vina_cpu = ctk.CTkEntry(self.vina_frame, width=60, fg_color=t.bg_tertiary,
                                         border_color=t.border, text_color=t.text_primary)
        self.ent_vina_cpu.insert(0, str(getattr(self.config, "cpu", "AUTO")))
        self.ent_vina_cpu.grid(row=1, column=3, padx=4, pady=(0, 6), sticky="ew")

        ctk.CTkLabel(self.vina_frame, text="Random Seed:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=4, padx=(10, 4), pady=(0, 6), sticky="w")
        self.ent_vina_seed = ctk.CTkEntry(self.vina_frame, width=60, fg_color=t.bg_tertiary,
                                          border_color=t.border, text_color=t.text_primary)
        self.ent_vina_seed.insert(0, str(getattr(self.config, "seed", "AUTO")))
        self.ent_vina_seed.grid(row=1, column=5, padx=4, pady=(0, 6), sticky="ew")

    # ─────────────────────────────────────────────────
    # AutoDock4 / ADT Parameter Panel
    # ─────────────────────────────────────────────────

    def _build_ad4_panel(self) -> None:
        t = themes.get()
        self.ad4_frame = ctk.CTkFrame(self.params_panel, fg_color="transparent")
        self.ad4_frame.grid_columnconfigure(0, weight=1)

        self.ad4_tabs = ctk.CTkTabview(
            self.ad4_frame, height=130,
            fg_color=t.bg_tertiary,
            segmented_button_fg_color=t.bg_primary,
            segmented_button_selected_color=t.accent,
            segmented_button_unselected_color=t.bg_tertiary,
            text_color=t.text_primary,
        )
        self.ad4_tabs.grid(row=0, column=0, sticky="nsew", padx=6, pady=(0, 4))

        tab_ga = self.ad4_tabs.add("🧬 Search & GA (DPF)")
        tab_sw = self.ad4_tabs.add("🔍 Solis-Wets Local Search")
        tab_out = self.ad4_tabs.add("📊 Analysis & Maps")

        # ── Tab 1: GA ────────────────────────────────
        for c in range(8):
            tab_ga.grid_columnconfigure(c, weight=1 if c in (1, 3, 5, 7) else 0)

        # Row 0
        ctk.CTkLabel(tab_ga, text="Algorithm:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=4, pady=3, sticky="w")
        self.cmb_ad4_alg = ctk.CTkOptionMenu(
            tab_ga, values=["LGA", "GA", "LS"], width=80,
            fg_color=t.bg_primary, button_color=t.accent, text_color=t.text_primary,
            command=self._on_algorithm_change,
        )
        self.cmb_ad4_alg.set(getattr(self.config, "ad4_algorithm", "LGA"))
        self.cmb_ad4_alg.grid(row=0, column=1, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_run (Runs):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=2, padx=4, pady=3, sticky="w")
        self.ent_ad4_run = ctk.CTkEntry(tab_ga, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_run.insert(0, str(getattr(self.config, "ga_run", 100)))
        self.ent_ad4_run.grid(row=0, column=3, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_pop_size:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=4, padx=4, pady=3, sticky="w")
        self.ent_ad4_pop = ctk.CTkEntry(tab_ga, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_pop.insert(0, str(getattr(self.config, "ga_pop_size", 150)))
        self.ent_ad4_pop.grid(row=0, column=5, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_num_evals:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=6, padx=4, pady=3, sticky="w")
        self.ent_ad4_evals = ctk.CTkEntry(tab_ga, width=70, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_evals.insert(0, str(getattr(self.config, "ga_num_evals", 2500000)))
        self.ent_ad4_evals.grid(row=0, column=7, padx=4, pady=3, sticky="ew")

        # Row 1
        ctk.CTkLabel(tab_ga, text="ga_num_generations:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=0, padx=4, pady=3, sticky="w")
        self.ent_ad4_gens = ctk.CTkEntry(tab_ga, width=70, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_gens.insert(0, str(getattr(self.config, "ga_num_generations", 27000)))
        self.ent_ad4_gens.grid(row=1, column=1, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_elitism:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=2, padx=4, pady=3, sticky="w")
        self.ent_ad4_elitism = ctk.CTkEntry(tab_ga, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_elitism.insert(0, str(getattr(self.config, "ga_elitism", 1)))
        self.ent_ad4_elitism.grid(row=1, column=3, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_mutation_rate:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=4, padx=4, pady=3, sticky="w")
        self.ent_ad4_mut = ctk.CTkEntry(tab_ga, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_mut.insert(0, str(getattr(self.config, "ga_mutation_rate", 0.02)))
        self.ent_ad4_mut.grid(row=1, column=5, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_ga, text="ga_crossover_rate:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=6, padx=4, pady=3, sticky="w")
        self.ent_ad4_cross = ctk.CTkEntry(tab_ga, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_cross.insert(0, str(getattr(self.config, "ga_crossover_rate", 0.80)))
        self.ent_ad4_cross.grid(row=1, column=7, padx=4, pady=3, sticky="ew")

        # ── Tab 2: Solis-Wets Local Search ────────────
        for c in range(6):
            tab_sw.grid_columnconfigure(c, weight=1 if c in (1, 3, 5) else 0)

        ctk.CTkLabel(tab_sw, text="sw_max_its (Max Iters):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=4, pady=3, sticky="w")
        self.ent_sw_its = ctk.CTkEntry(tab_sw, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_sw_its.insert(0, str(getattr(self.config, "sw_max_its", 300)))
        self.ent_sw_its.grid(row=0, column=1, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_sw, text="sw_max_succ:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=2, padx=4, pady=3, sticky="w")
        self.ent_sw_succ = ctk.CTkEntry(tab_sw, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_sw_succ.insert(0, str(getattr(self.config, "sw_max_succ", 4)))
        self.ent_sw_succ.grid(row=0, column=3, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_sw, text="sw_max_fail:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=4, padx=4, pady=3, sticky="w")
        self.ent_sw_fail = ctk.CTkEntry(tab_sw, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_sw_fail.insert(0, str(getattr(self.config, "sw_max_fail", 4)))
        self.ent_sw_fail.grid(row=0, column=5, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_sw, text="sw_rho (Step Size):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=0, padx=4, pady=3, sticky="w")
        self.ent_sw_rho = ctk.CTkEntry(tab_sw, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_sw_rho.insert(0, str(getattr(self.config, "sw_rho", 1.0)))
        self.ent_sw_rho.grid(row=1, column=1, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_sw, text="sw_lb_rho (Lower Bound):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=2, padx=4, pady=3, sticky="w")
        self.ent_sw_lb_rho = ctk.CTkEntry(tab_sw, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_sw_lb_rho.insert(0, str(getattr(self.config, "sw_lb_rho", 0.01)))
        self.ent_sw_lb_rho.grid(row=1, column=3, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_sw, text="ls_search_freq:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=4, padx=4, pady=3, sticky="w")
        self.ent_ls_freq = ctk.CTkEntry(tab_sw, width=50, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ls_freq.insert(0, str(getattr(self.config, "ls_search_freq", 0.06)))
        self.ent_ls_freq.grid(row=1, column=5, padx=4, pady=3, sticky="ew")

        # ── Tab 3: Output & Maps ──────────────────────
        for c in range(6):
            tab_out.grid_columnconfigure(c, weight=1 if c in (1, 3, 5) else 0)

        ctk.CTkLabel(tab_out, text="rmstol (RMSD Tol Å):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=4, pady=3, sticky="w")
        self.ent_ad4_rmstol = ctk.CTkEntry(tab_out, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_rmstol.insert(0, str(getattr(self.config, "rmstol", 2.0)))
        self.ent_ad4_rmstol.grid(row=0, column=1, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_out, text="outlev (Log Level):", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=2, padx=4, pady=3, sticky="w")
        self.cmb_ad4_outlev = ctk.CTkOptionMenu(
            tab_out, values=["0", "1", "2"], width=70,
            fg_color=t.bg_primary, button_color=t.accent, text_color=t.text_primary,
        )
        self.cmb_ad4_outlev.set(str(getattr(self.config, "outlev", 1)))
        self.cmb_ad4_outlev.grid(row=0, column=3, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_out, text="unbound_model:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=4, padx=4, pady=3, sticky="w")
        self.cmb_ad4_unbound = ctk.CTkOptionMenu(
            tab_out, values=["bound", "extended", "compact"], width=90,
            fg_color=t.bg_primary, button_color=t.accent, text_color=t.text_primary,
        )
        self.cmb_ad4_unbound.set(str(getattr(self.config, "unbound_model", "bound")))
        self.cmb_ad4_unbound.grid(row=0, column=5, padx=4, pady=3, sticky="ew")

        ctk.CTkLabel(tab_out, text="Random Seed:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=1, column=0, padx=4, pady=3, sticky="w")
        self.ent_ad4_seed = ctk.CTkEntry(tab_out, width=60, fg_color=t.bg_primary, text_color=t.text_primary)
        self.ent_ad4_seed.insert(0, str(getattr(self.config, "ad4_seed", "AUTO")))
        self.ent_ad4_seed.grid(row=1, column=1, padx=4, pady=3, sticky="ew")

        self.chk_ad4_reuse_maps = ctk.CTkCheckBox(
            tab_out, text="Reuse Existing Maps (.map)",
            text_color=t.text_secondary, font=ctk.CTkFont(size=11),
        )
        if getattr(self.config, "reuse_existing_maps", False):
            self.chk_ad4_reuse_maps.select()
        self.chk_ad4_reuse_maps.grid(row=1, column=2, columnspan=2, padx=4, pady=3, sticky="w")

    # ─────────────────────────────────────────────────
    # Panel View Switching & Controls
    # ─────────────────────────────────────────────────

    def _toggle_params(self) -> None:
        self._params_visible = not self._params_visible
        if self._params_visible:
            self.params_panel.grid()
            self.btn_toggle_params.configure(text="⚙ Parameters ▲")
        else:
            self.params_panel.grid_remove()
            self.btn_toggle_params.configure(text="⚙ Parameters ▼")

    def _update_param_view(self) -> None:
        engine = self.cmb_engine.get().upper()
        if engine == "VINA":
            self.ad4_frame.grid_remove()
            self.vina_frame.grid(row=1, column=0, sticky="ew", padx=6, pady=4)
            self.lbl_panel_title.configure(text="⚙ AUTODOCK VINA PARAMETERS")
        elif engine == "BOTH":
            # Show BOTH panels stacked: Vina on top, AD4 below
            self.vina_frame.grid(row=1, column=0, sticky="ew", padx=6, pady=(4, 0))
            self.ad4_frame.grid(row=2, column=0, sticky="ew", padx=6, pady=(0, 4))
            self.lbl_panel_title.configure(text="⚙ DUAL-ENGINE: VINA (above)  +  AUTODOCK 4 (below)")
        else:
            self.vina_frame.grid_remove()
            self.ad4_frame.grid(row=1, column=0, sticky="ew", padx=6, pady=4)
            self.lbl_panel_title.configure(text="⚙ AUTODOCK 4.2.6 (ADT) PARAMETERS")

    def _on_engine_change(self, value: str) -> None:
        from models import Engine
        try:
            self.config.engine = Engine(value.upper())
        except ValueError:
            pass
        self._update_param_view()
        self.on_config_change()

    def _on_algorithm_change(self, value: str) -> None:
        self.config.ad4_algorithm = value.upper()

    def _reset_params_to_defaults(self) -> None:
        """Reset parameter input fields to default recommended values."""
        engine = self.cmb_engine.get().upper()
        if engine == "VINA":
            self.cmb_vina_scoring.set("vina")
            self.ent_vina_exhaust.delete(0, "end"); self.ent_vina_exhaust.insert(0, "8")
            self.ent_vina_modes.delete(0, "end"); self.ent_vina_modes.insert(0, "9")
            self.ent_vina_energy.delete(0, "end"); self.ent_vina_energy.insert(0, "3.0")
            self.ent_vina_min_rmsd.delete(0, "end"); self.ent_vina_min_rmsd.insert(0, "1.0")
            self.ent_vina_cpu.delete(0, "end"); self.ent_vina_cpu.insert(0, "AUTO")
            self.ent_vina_seed.delete(0, "end"); self.ent_vina_seed.insert(0, "AUTO")
        else:
            self.cmb_ad4_alg.set("LGA")
            self.ent_ad4_run.delete(0, "end"); self.ent_ad4_run.insert(0, "100")
            self.ent_ad4_pop.delete(0, "end"); self.ent_ad4_pop.insert(0, "150")
            self.ent_ad4_evals.delete(0, "end"); self.ent_ad4_evals.insert(0, "2500000")
            self.ent_ad4_gens.delete(0, "end"); self.ent_ad4_gens.insert(0, "27000")
            self.ent_ad4_elitism.delete(0, "end"); self.ent_ad4_elitism.insert(0, "1")
            self.ent_ad4_mut.delete(0, "end"); self.ent_ad4_mut.insert(0, "0.02")
            self.ent_ad4_cross.delete(0, "end"); self.ent_ad4_cross.insert(0, "0.80")
            self.ent_sw_its.delete(0, "end"); self.ent_sw_its.insert(0, "300")
            self.ent_sw_succ.delete(0, "end"); self.ent_sw_succ.insert(0, "4")
            self.ent_sw_fail.delete(0, "end"); self.ent_sw_fail.insert(0, "4")
            self.ent_sw_rho.delete(0, "end"); self.ent_sw_rho.insert(0, "1.0")
            self.ent_sw_lb_rho.delete(0, "end"); self.ent_sw_lb_rho.insert(0, "0.01")
            self.ent_ls_freq.delete(0, "end"); self.ent_ls_freq.insert(0, "0.06")
            self.ent_ad4_rmstol.delete(0, "end"); self.ent_ad4_rmstol.insert(0, "2.0")
            self.cmb_ad4_outlev.set("1")
            self.cmb_ad4_unbound.set("bound")
            self.ent_ad4_seed.delete(0, "end"); self.ent_ad4_seed.insert(0, "AUTO")
            self.chk_ad4_reuse_maps.deselect()

    def _sync_params_to_config(self) -> None:
        """Pulls all parameters from UI widgets into ProjectConfig."""
        from models import Engine, DockingMode
        try:
            self.config.engine = Engine(self.cmb_engine.get().upper())
        except Exception:
            pass

        try:
            self.config.docking_mode = DockingMode(self.cmb_mode.get().upper())
        except Exception:
            pass

        # Vina parameters
        try:
            self.config.vina_scoring = self.cmb_vina_scoring.get().strip().lower()
            self.config.exhaustiveness = int(self.ent_vina_exhaust.get().strip() or "8")
            self.config.num_modes = int(self.ent_vina_modes.get().strip() or "9")
            self.config.energy_range = float(self.ent_vina_energy.get().strip() or "3.0")
            self.config.vina_min_rmsd = float(self.ent_vina_min_rmsd.get().strip() or "1.0")
            self.config.cpu = self.ent_vina_cpu.get().strip() or "AUTO"
            self.config.seed = self.ent_vina_seed.get().strip() or "AUTO"
        except Exception as e:
            print(f"[Warning] Error parsing Vina parameters: {e}")

        # AutoDock4 parameters
        try:
            self.config.ad4_algorithm = self.cmb_ad4_alg.get().strip().upper()
            self.config.ga_run = int(self.ent_ad4_run.get().strip() or "100")
            self.config.ga_pop_size = int(self.ent_ad4_pop.get().strip() or "150")
            self.config.ga_num_evals = int(self.ent_ad4_evals.get().strip() or "2500000")
            self.config.ga_num_generations = int(self.ent_ad4_gens.get().strip() or "27000")
            self.config.ga_elitism = int(self.ent_ad4_elitism.get().strip() or "1")
            self.config.ga_mutation_rate = float(self.ent_ad4_mut.get().strip() or "0.02")
            self.config.ga_crossover_rate = float(self.ent_ad4_cross.get().strip() or "0.80")

            self.config.sw_max_its = int(self.ent_sw_its.get().strip() or "300")
            self.config.sw_max_succ = int(self.ent_sw_succ.get().strip() or "4")
            self.config.sw_max_fail = int(self.ent_sw_fail.get().strip() or "4")
            self.config.sw_rho = float(self.ent_sw_rho.get().strip() or "1.0")
            self.config.sw_lb_rho = float(self.ent_sw_lb_rho.get().strip() or "0.01")
            self.config.ls_search_freq = float(self.ent_ls_freq.get().strip() or "0.06")

            self.config.rmstol = float(self.ent_ad4_rmstol.get().strip() or "2.0")
            self.config.outlev = int(self.cmb_ad4_outlev.get().strip() or "1")
            self.config.unbound_model = self.cmb_ad4_unbound.get().strip().lower()
            self.config.ad4_seed = self.ent_ad4_seed.get().strip() or "AUTO"
            self.config.reuse_existing_maps = bool(self.chk_ad4_reuse_maps.get())
        except Exception as e:
            print(f"[Warning] Error parsing AutoDock4 parameters: {e}")

        self.on_config_change()

    # ─────────────────────────────────────────────────
    # Docking Pipeline Execution
    # ─────────────────────────────────────────────────

    def _start_docking(self) -> None:
        if self._running:
            return
        self._running = True
        self.btn_start.configure(state="disabled")
        self.btn_stop.configure(state="normal", fg_color=themes.get().bg_tertiary,
                                text_color=themes.get().text_primary)
        self.progress_bar.set(0)
        self.lbl_progress.configure(text="Initialising docking pipeline …")

        # Sync all parameters directly from controls
        self._sync_params_to_config()

        # Redirect stdout to console
        self._orig_stdout = sys.stdout
        self._redirector = TextRedirector(self.console)
        sys.stdout = self._redirector

        thread = threading.Thread(target=self._run_pipeline, daemon=True)
        thread.start()

    def _run_pipeline(self) -> None:
        """Runs the full docking workflow in a background thread."""
        import copy
        from models import ResumeMode, Engine
        try:
            if self.config.engine == Engine.VINA:
                from vina_workflow import run_vina_workflow
                jobs = run_vina_workflow(
                    self.config,
                    resume_mode=ResumeMode.RESUME,
                    progress_callback=self._progress_callback,
                )
            elif self.config.engine == Engine.AUTODOCK4:
                from autodock4_workflow import run_autodock4_workflow
                jobs = run_autodock4_workflow(
                    self.config,
                    progress_callback=self._progress_callback,
                )
            elif self.config.engine == Engine.BOTH:
                # ── Sequential: Vina first, then AutoDock4 ──────────────────
                print("\n" + "=" * 60)
                print("  DUAL-ENGINE — Phase 1 of 2: AutoDock Vina")
                print("=" * 60)
                vina_cfg = copy.copy(self.config)
                vina_cfg.engine = Engine.VINA
                from vina_workflow import run_vina_workflow
                vina_jobs = run_vina_workflow(
                    vina_cfg,
                    resume_mode=ResumeMode.RESUME,
                    progress_callback=self._progress_callback,
                ) or []

                print("\n" + "=" * 60)
                print("  DUAL-ENGINE — Phase 2 of 2: AutoDock 4.2.6")
                print("=" * 60)
                ad4_cfg = copy.copy(self.config)
                ad4_cfg.engine = Engine.AUTODOCK4
                from autodock4_workflow import run_autodock4_workflow
                ad4_jobs = run_autodock4_workflow(
                    ad4_cfg,
                    progress_callback=self._progress_callback,
                ) or []

                jobs = vina_jobs + ad4_jobs
                if jobs:
                    from reporting import generate_reports
                    generate_reports(jobs, self.config)
                print("\n✔  Dual-engine run complete.")
            else:
                jobs = []
            self.after(0, self._on_done, jobs, None)
        except Exception as e:
            self.after(0, self._on_done, [], e)

    def _progress_callback(self, completed: int, total: int, current_label: str = "") -> None:
        """Called from the docking thread after each job finishes."""
        self.after(0, self._update_progress, completed, total, current_label)

    def _update_progress(self, completed: int, total: int, label: str) -> None:
        self._jobs_completed = completed
        self._jobs_total = total
        if total > 0:
            self.progress_bar.set(completed / total)
        pct = int(100 * completed / total) if total > 0 else 0
        self.lbl_progress.configure(
            text=f"Job {completed} / {total} — {label}  [{pct}%]"
        )

    def _on_done(self, jobs: List, error: Optional[Exception]) -> None:
        self._running = False
        sys.stdout = self._orig_stdout
        self.btn_start.configure(state="normal")
        self.btn_stop.configure(state="disabled",
                                fg_color=themes.get().bg_tertiary,
                                text_color=themes.get().text_disabled)
        if error:
            self.progress_bar.set(0)
            self.lbl_progress.configure(
                text=f"Pipeline error: {error}", text_color=themes.get().status_fail_text
            )
        else:
            self.progress_bar.set(1.0)
            self.lbl_progress.configure(
                text=f"✔ Completed — {len(jobs)} job(s) processed",
                text_color=themes.get().status_ok_text,
            )
        self.on_run_complete(jobs)

    def _stop_docking(self) -> None:
        """Signal a graceful stop (marks flag — actual stop depends on workflow)."""
        self._running = False
        self.lbl_progress.configure(text="Stop requested — waiting for current job …")

    def _zoom_in(self) -> None:
        if self._console_font_size < 26:
            self._console_font_size += 1
            self.console.configure(font=ctk.CTkFont(family="Consolas", size=self._console_font_size))
            self.lbl_font_size.configure(text=f"{self._console_font_size} pt")

    def _zoom_out(self) -> None:
        if self._console_font_size > 9:
            self._console_font_size -= 1
            self.console.configure(font=ctk.CTkFont(family="Consolas", size=self._console_font_size))
            self.lbl_font_size.configure(text=f"{self._console_font_size} pt")

    def _open_dlg_folder(self) -> None:
        dlg_dir = getattr(self.config, "result_directory", None)
        target = None
        if dlg_dir and (Path(dlg_dir) / "DLG").is_dir():
            target = Path(dlg_dir) / "DLG"
        elif getattr(self.config, "project_root", None) and (Path(self.config.project_root) / "DLG").is_dir():
            target = Path(self.config.project_root) / "DLG"
        elif dlg_dir:
            target = Path(dlg_dir) / "DLG"
        else:
            target = Path(getattr(self.config, "project_root", Path.cwd())) / "results" / "DLG"
        target.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(target))
        else:
            subprocess.run(["xdg-open", str(target)])

    def _open_results_folder(self) -> None:
        res_dir = getattr(self.config, "result_directory", None)
        if not res_dir:
            res_dir = Path(getattr(self.config, "project_root", Path.cwd())) / "results"
        res_dir = Path(res_dir)
        res_dir.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(res_dir))
        else:
            subprocess.run(["xdg-open", str(res_dir)])

    def _clear_console(self) -> None:
        self.console.configure(state="normal")
        self.console.delete("1.0", "end")
        self.console.configure(state="disabled")

    def log(self, text: str) -> None:
        """Write text directly to the console (called from other tabs)."""
        self._redirector.write(text + "\n") if hasattr(self, "_redirector") else None
