"""
AutoDock Suite Pro — Results Tab (gui/results_tab.py)
=====================================================
Results browser with:
  - Job list with receptor/ligand/status
  - Docked Poses table (pose-responsive — clicking a pose updates interactions)
  - Protein–Ligand Interactions table (responds to selected pose)
  - 2D Ligand Structure depiction panel
  - Compound ADMET physicochemical profile card
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import threading
from io import BytesIO
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List, Optional, Tuple

import customtkinter as ctk
from tkinter import ttk

from gui import themes, widgets

if TYPE_CHECKING:
    from config import ProjectConfig
    from models import DockingJob


class ResultsTab(ctk.CTkFrame):
    """Full results and interaction analysis tab."""

    POSE_COLS = ("Pose", "Affinity", "Ki", "RMSD_lb", "RMSD_ub")
    INTER_COLS = ("Pose", "Type", "Residue", "Rec.Atom", "Lig.Atom", "Dist.(Å)", "Details")
    JOB_COLS = ("Job ID", "Receptor", "Ligand", "Mode", "Status", "ΔG (best)")

    def __init__(self, parent: ctk.CTkTabview, config: "ProjectConfig",
                 on_config_change: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self._jobs: List["DockingJob"] = []
        self._current_job: Optional["DockingJob"] = None
        self._structure_image = None  # PhotoImage reference kept alive
        self._build()

    # ─────────────────────────────────────────────────
    # Build
    # ─────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(2, weight=1)

        # Header controls
        hdr = ctk.CTkFrame(self, fg_color="transparent")
        hdr.grid(row=0, column=0, sticky="ew", padx=12, pady=(10, 0))
        hdr.grid_columnconfigure(6, weight=1)

        widgets.section_header(hdr, "Results & Screening", row=0, col=0)
        self.btn_refresh = widgets.accent_button(hdr, "⟳ Refresh Results",
                              self.refresh, width=130, small=True
                              )
        self.btn_refresh.grid(row=0, column=1, padx=(16, 4))
        self.btn_dlg_extract = widgets.accent_button(hdr, "🔬 Run DLG Extract",
                              self._run_dlg_extract, width=140, small=True
                              )
        self.btn_dlg_extract.grid(row=0, column=2, padx=4)
        self.btn_export = widgets.accent_button(hdr, "📤 Export Report",
                              self._export_report, width=120, small=True
                              )
        self.btn_export.grid(row=0, column=3, padx=4)
        self.btn_dlg_folder = ctk.CTkButton(hdr, text="📁 DLG Folder",
                      command=self._open_dlg_folder, width=110, height=26,
                      fg_color=t.bg_secondary, text_color=t.text_primary,
                      hover_color=t.border, font=ctk.CTkFont(size=11)
                      )
        self.btn_dlg_folder.grid(row=0, column=4, padx=4)
        self.btn_results_folder = ctk.CTkButton(hdr, text="📁 Results Folder",
                      command=self._open_results_folder, width=120, height=26,
                      fg_color=t.bg_secondary, text_color=t.text_primary,
                      hover_color=t.border, font=ctk.CTkFont(size=11)
                      )
        self.btn_results_folder.grid(row=0, column=5, padx=4)

        self.lbl_action_status = ctk.CTkLabel(
            hdr, text="", font=ctk.CTkFont(size=11),
            text_color=t.accent_primary, anchor="w"
        )
        self.lbl_action_status.grid(row=0, column=6, padx=(12, 4), sticky="w")

        # ── KPI bar ─────────────────────────────────
        kpi = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        kpi.grid(row=1, column=0, sticky="ew", padx=12, pady=6)
        for col in range(4):
            kpi.grid_columnconfigure(col, weight=1)

        self.lbl_kpi_total = self._kpi_cell(kpi, "Total Jobs", "—", 0)
        self.lbl_kpi_success = self._kpi_cell(kpi, "Completed", "—", 1)
        self.lbl_kpi_lead = self._kpi_cell(kpi, "Top Lead", "Run docking to see results", 2, wide=True)
        self.lbl_kpi_failed = self._kpi_cell(kpi, "Failed", "—", 3)

        # ── Main split: job list left, detail right ──
        body = ctk.CTkFrame(self, fg_color="transparent")
        body.grid(row=2, column=0, sticky="nsew", padx=12, pady=(0, 8))
        body.grid_columnconfigure(0, weight=3)
        body.grid_columnconfigure(1, weight=2)
        body.grid_rowconfigure(0, weight=1)

        # Left: Job list
        left = ctk.CTkFrame(body, fg_color=t.bg_secondary, corner_radius=8)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 6))
        left.grid_rowconfigure(1, weight=1)
        left.grid_columnconfigure(0, weight=1)

        ctk.CTkLabel(left, text="Docking Jobs", text_color=t.accent,
                     font=themes.font_section(), anchor="w"
                     ).grid(row=0, column=0, sticky="w", padx=8, pady=(6, 2))

        self.tree_jobs = widgets.styled_treeview(left, self.JOB_COLS, height=12)
        for col, w in zip(self.JOB_COLS, (160, 90, 130, 70, 80, 100)):
            self.tree_jobs.heading(col, text=col)
            self.tree_jobs.column(col, width=w, anchor="center")
        self.tree_jobs.grid(row=1, column=0, sticky="nsew", padx=4, pady=(2, 4))
        self.tree_jobs.bind("<<TreeviewSelect>>", self._on_job_select)

        vsb_jobs = ttk.Scrollbar(left, orient="vertical", command=self.tree_jobs.yview)
        vsb_jobs.grid(row=1, column=1, sticky="ns")
        self.tree_jobs.configure(yscrollcommand=vsb_jobs.set)

        # Poses sub-table (below job list)
        ctk.CTkLabel(left, text="Docked Poses", text_color=t.accent,
                     font=themes.font_section(), anchor="w"
                     ).grid(row=2, column=0, sticky="w", padx=8, pady=(8, 2))

        self.tree_poses = widgets.styled_treeview(left, self.POSE_COLS, height=5)
        for col, w in zip(self.POSE_COLS, (80, 130, 150, 120, 120)):
            self.tree_poses.heading(col, text={
                "Affinity": "ΔG (kcal/mol)", "RMSD_lb": "RMSD l.b. (Å)",
                "RMSD_ub": "RMSD u.b. (Å)", "Ki": "Est. Ki / Kd"
            }.get(col, col))
            self.tree_poses.column(col, width=w, anchor="center")
        self.tree_poses.grid(row=3, column=0, columnspan=2, sticky="ew", padx=4, pady=(2, 6))
        self.tree_poses.bind("<<TreeviewSelect>>", self._on_pose_select)

        # Right panel: 2D structure + ADMET card + interactions
        right = ctk.CTkFrame(body, fg_color="transparent")
        right.grid(row=0, column=1, sticky="nsew")
        right.grid_rowconfigure(2, weight=1)
        right.grid_columnconfigure(0, weight=1)

        # 2D Interaction Diagram card (replaces old 2D ligand PNG panel)
        struct_card = ctk.CTkFrame(right, fg_color=t.bg_secondary, corner_radius=8)
        struct_card.grid(row=0, column=0, sticky="nsew", pady=(0, 6))
        struct_card.grid_columnconfigure(0, weight=1)
        struct_card.grid_rowconfigure(1, weight=1)
        right.grid_rowconfigure(0, weight=2)   # give the diagram more vertical space

        ctk.CTkLabel(struct_card, text="2D Interaction Diagram", text_color=t.accent,
                     font=themes.font_section(), anchor="w"
                     ).grid(row=0, column=0, sticky="w", padx=8, pady=(6, 0))

        # Try to embed tkinterweb for rich SVG interaction diagram
        self._has_html_frame = False
        self._render_generation = 0  # increments each selection; stale callbacks are dropped
        try:
            import tkinterweb  # type: ignore
            self._html_frame = tkinterweb.HtmlFrame(
                struct_card, messages_enabled=False,
                vertical_scrollbar=False, horizontal_scrollbar=False,
            )
            self._html_frame.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)
            self._has_html_frame = True
            # Seed with a "select a job" placeholder
            self._html_frame.load_html(self._placeholder_html("Select a job to view"))
        except Exception:
            # Graceful fallback to plain text label + old PIL rendering
            self.lbl_structure = ctk.CTkLabel(
                struct_card, text="Select a job to view",
                text_color=t.text_secondary, font=themes.font_caption(),
            )
            self.lbl_structure.grid(row=1, column=0, padx=8, pady=6, sticky="nsew")
            self._structure_image = None

        # ADMET card
        admet_card = ctk.CTkFrame(right, fg_color=t.bg_secondary, corner_radius=8)
        admet_card.grid(row=1, column=0, sticky="ew", pady=(0, 6))
        admet_card.grid_columnconfigure(0, weight=1)
        admet_card.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(admet_card, text="Compound Profile",
                     text_color=t.accent,
                     font=themes.font_section(), anchor="w"
                     ).grid(row=0, column=0, columnspan=4, sticky="w", padx=8, pady=(6, 2))

        self._admet_labels: Dict[str, ctk.CTkLabel] = {}
        admet_fields = [
            ("MW", "mw", "Da"), ("LogP", "logp", ""),
            ("PSA", "psa", "Å²"), ("HBD", "hbd", ""),
            ("HBA", "hba", ""), ("RotBonds", "rot_bonds", ""),
        ]
        for i, (label, key, unit) in enumerate(admet_fields):
            row, col_base = divmod(i, 3)
            ctk.CTkLabel(admet_card, text=f"{label}:",
                         text_color=t.text_secondary,
                         font=themes.font_body()).grid(
                row=row + 1, column=col_base * 2, sticky="w", padx=(8, 2), pady=2)
            lbl = ctk.CTkLabel(admet_card, text="—", text_color=t.text_primary,
                               font=themes.font_body_bold())
            lbl.grid(row=row + 1, column=col_base * 2 + 1, sticky="w", padx=(0, 10))
            self._admet_labels[key] = lbl

        # RO5 / Veber badges
        badge_frame = ctk.CTkFrame(admet_card, fg_color="transparent")
        badge_frame.grid(row=3, column=0, columnspan=4, sticky="w", padx=8, pady=(4, 6))
        ctk.CTkLabel(badge_frame, text="RO5:", text_color=t.text_secondary,
                     font=themes.font_body()).pack(side="left", padx=(0, 4))
        self.lbl_ro5_badge = ctk.CTkLabel(badge_frame, text="—",
                                          font=themes.font_body_bold(),
                                          text_color=t.text_secondary)
        self.lbl_ro5_badge.pack(side="left", padx=(0, 12))
        ctk.CTkLabel(badge_frame, text="Veber:", text_color=t.text_secondary,
                     font=themes.font_body()).pack(side="left", padx=(0, 4))
        self.lbl_veber_badge = ctk.CTkLabel(badge_frame, text="—",
                                            font=themes.font_body_bold(),
                                            text_color=t.text_secondary)
        self.lbl_veber_badge.pack(side="left")

        # Interactions table
        inter_card = ctk.CTkFrame(right, fg_color=t.bg_secondary, corner_radius=8)
        inter_card.grid(row=2, column=0, sticky="nsew")
        inter_card.grid_rowconfigure(1, weight=1)
        inter_card.grid_columnconfigure(0, weight=1)

        inter_hdr = ctk.CTkFrame(inter_card, fg_color="transparent")
        inter_hdr.grid(row=0, column=0, sticky="ew", padx=8, pady=(6, 2))

        ctk.CTkLabel(inter_hdr, text="Protein–Ligand Interactions",
                     text_color=t.accent,
                     font=themes.font_section(), anchor="w"
                     ).pack(side="left")

        btn_bar = ctk.CTkFrame(inter_hdr, fg_color="transparent")
        btn_bar.pack(side="right")

        ctk.CTkButton(btn_bar, text="🌐 WebGL", width=80, height=26,
                      command=self._launch_3d_inspector,
                      fg_color=t.accent, hover_color=t.accent_hover,
                      font=themes.font_button()
                      ).pack(side="left", padx=2)

        ctk.CTkButton(btn_bar, text="🔬 Open in PyMOL", width=125, height=26,
                      command=self._launch_pymol_inspector,
                      fg_color="#107C41", hover_color="#0b5c30",
                      font=themes.font_button()
                      ).pack(side="left", padx=2)

        ctk.CTkButton(btn_bar, text="💾 Export PyMOL", width=115, height=26,
                      command=self._export_pymol_session_ui,
                      fg_color=t.bg_tertiary, hover_color=t.accent_hover,
                      font=themes.font_button()
                      ).pack(side="left", padx=2)

        self.tree_inter = widgets.styled_treeview(inter_card, self.INTER_COLS, height=10)
        for col, w in zip(self.INTER_COLS, (50, 160, 100, 80, 70, 75, 200)):
            self.tree_inter.heading(col, text=col)
            self.tree_inter.column(col, width=w, anchor="center" if col not in ("Details",) else "w")
        self.tree_inter.grid(row=1, column=0, sticky="nsew", padx=4, pady=(2, 4))

        vsb_inter = ttk.Scrollbar(inter_card, orient="vertical", command=self.tree_inter.yview)
        vsb_inter.grid(row=1, column=1, sticky="ns")
        self.tree_inter.configure(yscrollcommand=vsb_inter.set)

    # ─────────────────────────────────────────────────
    # KPI helper
    # ─────────────────────────────────────────────────

    def _kpi_cell(self, parent: ctk.CTkFrame, title: str, value: str,
                  col: int, wide: bool = False) -> ctk.CTkLabel:
        t = themes.get()
        f = ctk.CTkFrame(parent, fg_color=t.bg_tertiary, corner_radius=6)
        f.grid(row=0, column=col, sticky="ew", padx=4, pady=6)
        ctk.CTkLabel(f, text=title, text_color=t.text_secondary,
                     font=themes.font_caption()).pack(pady=(4, 0))
        lbl = ctk.CTkLabel(f, text=value, text_color=t.text_primary,
                           font=ctk.CTkFont(family="Segoe UI", size=11 if wide else 14, weight="bold"),
                           wraplength=200 if wide else 0)
        lbl.pack(pady=(0, 4))
        return lbl

    # ─────────────────────────────────────────────────
    # Data loading
    # ─────────────────────────────────────────────────

    def refresh(self) -> None:
        """Re-scan results directory and populate job list."""
        from job_manager import load_docking_jobs
        try:
            # FIX: pass Path, not ProjectConfig
            self._jobs = load_docking_jobs(self.config.result_directory)
            # Resolve any relative paths stored in job_status.json
            project_root = self.config.project_root
            for job in self._jobs:
                for attr in ("receptor_path", "ligand_path", "output_path",
                             "log_path", "dlg_path"):
                    val = getattr(job, attr, None)
                    if val and not Path(val).is_absolute():
                        setattr(job, attr, project_root / val)
        except Exception as e:
            self._jobs = []
            import logging
            logging.getLogger("docking_automation.gui").warning(
                f"ResultsTab.refresh error: {e}")
        self._populate_job_table()

    def _populate_job_table(self) -> None:
        for item in self.tree_jobs.get_children():
            self.tree_jobs.delete(item)

        from models import JobStatus
        success = 0
        failed = 0
        best_affinity = None
        best_lead = ""

        for i, job in enumerate(self._jobs):
            status_str = job.status.value if hasattr(job.status, "value") else str(job.status)
            best_dg = "—"
            if job.vina_results:
                dg = job.vina_results[0].binding_affinity
                best_dg = f"{dg:.2f}"
                if best_affinity is None or dg < best_affinity:
                    best_affinity = dg
                    best_lead = f"{job.ligand_name} ({dg:.2f} kcal/mol)"
            else:
                ad4_dg = self._get_ad4_affinity(job)
                if ad4_dg is not None:
                    best_dg = f"{ad4_dg:.2f}"
                    if best_affinity is None or ad4_dg < best_affinity:
                        best_affinity = ad4_dg
                        best_lead = f"{job.ligand_name} ({ad4_dg:.2f} kcal/mol)"

            if "SUCCESS" in status_str.upper():
                success += 1
            elif "FAIL" in status_str.upper() or "ERROR" in status_str.upper():
                failed += 1

            mode_str = job.docking_mode.value if hasattr(job.docking_mode, "value") else "?"
            row_tag = "even" if i % 2 == 0 else "odd"
            self.tree_jobs.insert("", "end", iid=str(i), tags=(row_tag,), values=(
                job.job_id, job.receptor_name, job.ligand_name,
                mode_str, status_str, best_dg,
            ))

        total = len(self._jobs)
        self.lbl_kpi_total.configure(text=str(total))
        self.lbl_kpi_success.configure(text=str(success))
        self.lbl_kpi_failed.configure(text=str(failed))
        self.lbl_kpi_lead.configure(
            text=best_lead if best_lead else "No completed jobs yet"
        )

    # ─────────────────────────────────────────────────
    # Event handlers
    # ─────────────────────────────────────────────────

    def _on_job_select(self, _event=None) -> None:
        sel = self.tree_jobs.selection()
        if not sel:
            return
        idx = int(sel[0])
        if idx >= len(self._jobs):
            return
        job = self._jobs[idx]
        self._current_job = job
        self._populate_poses(job)
        # Start diagram with empty interactions immediately; will refresh once
        # _populate_interactions completes and passes the real interaction list.
        self._update_interaction_diagram(job, interactions=[], pose_index=1)
        self._update_admet(job)

    def _on_pose_select(self, _event=None) -> None:
        """Called when the user clicks on a pose row — updates interactions for THAT pose."""
        sel = self.tree_poses.selection()
        if not sel or self._current_job is None:
            return
        row_values = self.tree_poses.item(sel[0], "values")
        try:
            pose_idx = int(row_values[0])
        except (ValueError, IndexError):
            pose_idx = 1
        self._populate_interactions(self._current_job, pose_index=pose_idx)

    def _populate_poses(self, job: "DockingJob") -> None:
        for item in self.tree_poses.get_children():
            self.tree_poses.delete(item)

        if not job.vina_results:
            ad4_poses = self._get_ad4_poses(job)
            for i, p in enumerate(ad4_poses):
                tag = "even" if i % 2 == 0 else "odd"
                self.tree_poses.insert("", "end", tags=(tag,), values=(
                    p.get("pose", i + 1),
                    p.get("affinity", "—"),
                    p.get("ki", "—"),
                    p.get("rmsd_lb", "—"),
                    p.get("rmsd_ub", "—"),
                ))
            first = self.tree_poses.get_children()
            if first:
                self.tree_poses.selection_set(first[0])
            return

        from gui.app import calculate_inhibition_constant  # imported from app util
        for i, res in enumerate(job.vina_results):
            try:
                _, ki_str = calculate_inhibition_constant(res.binding_affinity)
            except Exception:
                ki_str = "—"
            tag = "even" if i % 2 == 0 else "odd"
            self.tree_poses.insert("", "end", tags=(tag,), values=(
                res.pose,
                f"{res.binding_affinity:.2f}",
                ki_str,
                f"{res.rmsd_lower_bound:.3f}",
                f"{res.rmsd_upper_bound:.3f}",
            ))

        # Auto-select first pose and show its interactions
        first = self.tree_poses.get_children()
        if first:
            self.tree_poses.selection_set(first[0])
            self._populate_interactions(job, pose_index=1)

    def _populate_interactions(self, job: "DockingJob", pose_index: int = 1) -> None:
        """Populate interactions for the specified pose number (1-based)."""
        for item in self.tree_inter.get_children():
            self.tree_inter.delete(item)

        receptor_path = job.receptor_path
        output_pdbqt = job.output_pdbqt

        # Find receptor if missing or relative
        if not receptor_path or not Path(receptor_path).is_file():
            rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
            for cand in [rec_dir / f"{job.receptor_name}.pdbqt",
                         rec_dir / job.receptor_name / f"{job.receptor_name}.pdbqt",
                         rec_dir / job.receptor_name / "rigid" / f"{job.receptor_name}.pdbqt"]:
                if cand.is_file():
                    receptor_path = cand
                    break

        # Fallback for AutoDock4 or when output_pdbqt is not directly set
        if not output_pdbqt or not Path(output_pdbqt).is_file():
            res_dir = getattr(self.config, "result_directory", None)
            root = getattr(self.config, "project_root", None)
            candidates: List[Path] = []
            if getattr(job, "dlg_path", None) and Path(job.dlg_path).is_file():
                candidates.append(Path(job.dlg_path))
            if res_dir:
                candidates.append(Path(res_dir) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                candidates.append(Path(res_dir) / f"{job.receptor_name}_{job.ligand_name}.dlg")
                candidates.extend(list(Path(res_dir).rglob(f"*{job.ligand_name}*pose*.pdbqt")))
                candidates.extend(list(Path(res_dir).rglob(f"*{job.ligand_name}*run*.pdbqt")))
            if root:
                candidates.append(Path(root) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                candidates.append(Path(root) / "results" / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                candidates.extend(list(Path(root).rglob(f"*{job.ligand_name}*pose*.pdbqt")))
                candidates.extend(list(Path(root).rglob(f"*{job.ligand_name}*run*.pdbqt")))
            for cand in candidates:
                if cand and cand.is_file():
                    output_pdbqt = cand
                    break

        self._current_selected_job = job
        self._current_selected_pose = pose_index
        self._current_resolved_receptor = receptor_path
        self._current_resolved_ligand = output_pdbqt

        if not receptor_path or not Path(receptor_path).is_file():
            return
        if not output_pdbqt or not Path(output_pdbqt).is_file():
            return

        def _worker() -> None:
            try:
                from interactions import profile_docking_job
                interactions = profile_docking_job(
                    receptor_path, output_pdbqt,
                    pose_index=pose_index,
                )
                self.after(0, self._fill_interactions, interactions)
            except Exception:
                pass

        threading.Thread(target=_worker, daemon=True).start()

    def _get_current_grid_box(self) -> Tuple[Optional[Tuple[float, float, float]], Optional[Tuple[float, float, float]]]:
        """Resolves (center, size) coordinates for currently selected job."""
        job = getattr(self, "_current_selected_job", None)
        if not job:
            return None, None

        from validators import parse_vina_config
        candidates: List[Path] = []
        if getattr(job, "config_path", None) and Path(job.config_path).is_file():
            candidates.append(Path(job.config_path))

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        candidates.append(rec_dir / job.receptor_name / "config.txt")
        candidates.append(rec_dir / f"{job.receptor_name}_clean" / "config.txt")
        if job.receptor_name.endswith("_clean"):
            candidates.append(rec_dir / job.receptor_name.replace("_clean", "") / "config.txt")

        if getattr(job, "output_dir", None):
            candidates.append(Path(job.output_dir) / "input" / "config.txt")
            candidates.append(Path(job.output_dir) / "config.txt")

        for cand in candidates:
            if cand.is_file():
                try:
                    cfg = parse_vina_config(cand)
                    if "center_x" in cfg and "size_x" in cfg:
                        center = (float(cfg["center_x"]), float(cfg["center_y"]), float(cfg["center_z"]))
                        size = (float(cfg["size_x"]), float(cfg["size_y"]), float(cfg["size_z"]))
                        return center, size
                except Exception:
                    continue
        return None, None

    def _launch_3d_inspector(self) -> None:
        """Launch the standalone 3D WebGL viewer for the currently selected job and pose."""
        from tkinter import messagebox
        rec_path = getattr(self, "_current_resolved_receptor", None)
        lig_path = getattr(self, "_current_resolved_ligand", None)

        if not rec_path or not Path(rec_path).is_file():
            messagebox.showwarning("3D Inspector", "Select a job with an available receptor first.")
            return

        try:
            from viewer_3d import launch_3d_viewer
            # Extract target residues from the interaction table if populated
            target_res = []
            for item in self.tree_inter.get_children():
                vals = self.tree_inter.item(item, "values")
                if len(vals) > 2 and vals[2] and vals[2] not in target_res:
                    target_res.append(vals[2])

            grid_center, grid_size = self._get_current_grid_box()

            html_path = launch_3d_viewer(
                receptor_path=rec_path,
                ligand_path=lig_path,
                grid_center=grid_center,
                grid_size=grid_size,
                target_residues=target_res[:10],
            )
            messagebox.showinfo(
                "3D WebGL Inspector Launched",
                f"Opened interactive 3D WebGL viewer in your web browser!\n\n"
                f"Receptor: {Path(rec_path).name}\n"
                f"Ligand: {Path(lig_path).name if lig_path else 'None'}\n"
                f"Interactions highlighted: {len(target_res)} residues",
            )
        except Exception as e:
            messagebox.showerror("3D Inspector Error", f"Failed to launch 3D WebGL viewer:\n{e}")

    def _launch_pymol_inspector(self) -> None:
        """Launch native PyMOL desktop application for current job and pose.

        Automatically exports the session to a deterministic directory under
        results/pymol_sessions/ without prompting the user for a location.
        Works for every receptor/ligand/pose combination independently.
        """
        from tkinter import filedialog, messagebox
        from pymol_exporter import export_pymol_session, find_pymol_executable, launch_pymol

        rec_path = getattr(self, "_current_resolved_receptor", None)
        lig_path = getattr(self, "_current_resolved_ligand", None)
        pose_index = getattr(self, "_current_selected_pose", 1)

        if not rec_path or not Path(rec_path).is_file():
            messagebox.showwarning("PyMOL Inspector", "Select a job with an available receptor first.")
            return

        if not lig_path or not Path(lig_path).is_file():
            messagebox.showwarning("PyMOL Inspector", "No docked ligand file found for this job.")
            return

        # Check for PyMOL executable (prompt once, remember choice)
        custom_pymol = getattr(self, "_custom_pymol_path", None)
        pymol_exe = find_pymol_executable(custom_pymol)
        if not pymol_exe:
            resp = messagebox.askyesno(
                "PyMOL Not Detected",
                "PyMOL executable was not detected automatically on your system.\n\n"
                "Would you like to browse and locate your 'PyMOL.exe' / 'PyMOLWin.exe' now?\n\n"
                "(Choose 'No' to export the PyMOL session files (.pml / .pse) instead.)"
            )
            if resp:
                chosen = filedialog.askopenfilename(
                    title="Select PyMOL Executable",
                    filetypes=[("Executable Files", "*.exe;*.bat;*.cmd"), ("All Files", "*.*")]
                )
                if chosen and Path(chosen).is_file():
                    self._custom_pymol_path = chosen  # remember for future calls
                    pymol_exe = Path(chosen)
                else:
                    return
            else:
                # User chose to export only — no-dialog export
                self._auto_export_pymol_session()
                return

        # Auto-export to deterministic path (no dialog) then launch
        bundle = self._auto_export_pymol_session(silent=True)
        if bundle is None:
            return

        target = bundle.get("pse") or bundle.get("pml")
        if not target or not Path(target).is_file():
            messagebox.showerror("PyMOL Error", "Session file was not created. Check logs for details.")
            return

        try:
            success = launch_pymol(target, custom_pymol_path=pymol_exe)
            if success:
                inter_count = bundle.get("interactions_count", 0)
                messagebox.showinfo(
                    "PyMOL Launched",
                    f"PyMOL 3D Session launched successfully!\n\n"
                    f"Receptor: {Path(rec_path).name}\n"
                    f"Ligand:   {Path(lig_path).name} (Pose {pose_index})\n"
                    f"Session:  {Path(target).name}\n"
                    f"Interactions: {inter_count} non-covalent contacts modeled\n"
                    f"{'Grid Box: Active' if self._get_current_grid_box()[0] else ''}"
                )
            else:
                messagebox.showerror("Launch Failed", "Could not start PyMOL process.")
        except Exception as e:
            messagebox.showerror("PyMOL Error", f"Failed to launch PyMOL:\n{e}")

    def _auto_export_pymol_session(self, silent: bool = False) -> Optional[Dict]:
        """Export PyMOL session to a deterministic path automatically (no dialog).

        The output directory is always:
            <result_dir>/pymol_sessions/<rec_stem>_<lig_stem>_pose<N>/

        If the directory already exists it is overwritten cleanly so that
        every ligand/pose combination gets its own independent folder.

        Args:
            silent: If True, suppresses the success dialog (used by _launch_pymol_inspector).

        Returns:
            The export bundle dict, or None on failure/missing inputs.
        """
        from tkinter import messagebox
        from pymol_exporter import export_pymol_session

        rec_path = getattr(self, "_current_resolved_receptor", None)
        lig_path = getattr(self, "_current_resolved_ligand", None)
        pose_index = getattr(self, "_current_selected_pose", 1)

        if not rec_path or not Path(rec_path).is_file():
            if not silent:
                messagebox.showwarning("Export PyMOL", "Select a job with an available receptor first.")
            return None

        if not lig_path or not Path(lig_path).is_file():
            if not silent:
                messagebox.showwarning("Export PyMOL", "No docked ligand file found for this job.")
            return None

        rec_stem = Path(rec_path).stem
        lig_stem = Path(lig_path).stem

        # Deterministic, unique directory per receptor/ligand/pose combination
        base_sessions_dir = (
            Path(getattr(self.config, "result_directory", "results")).resolve() / "pymol_sessions"
        )
        dest_dir = base_sessions_dir / f"{rec_stem}_{lig_stem}_pose{pose_index}"
        # exist_ok=True so re-exporting the same pose always works
        dest_dir.mkdir(parents=True, exist_ok=True)

        grid_center, grid_size = self._get_current_grid_box()

        try:
            bundle = export_pymol_session(
                receptor_path=rec_path,
                ligand_path=lig_path,
                output_dir=dest_dir,
                pose_index=pose_index,
                grid_center=grid_center,
                grid_size=grid_size,
                compile_pse=True,
                custom_pymol_path=getattr(self, "_custom_pymol_path", None),
            )

            if not silent:
                pml_name = bundle["pml"].name
                pse_name = (
                    bundle["pse"].name
                    if bundle.get("pse")
                    else "(no local PyMOL found — open .pml in PyMOL)"
                )
                inter_count = bundle.get("interactions_count", 0)
                msg = (
                    f"Publication-Grade PyMOL Session Exported!\n\n"
                    f"Directory:\n  {dest_dir}\n\n"
                    f"PyMOL Script (.pml):\n  {pml_name}\n"
                    f"PyMOL Session (.pse):\n  {pse_name}\n"
                    f"Interactions Modeled:\n  {inter_count} contacts (H-bonds, Salt bridges, etc.)\n"
                    f"Publication Figures:\n  Auto 300 DPI ray-traced renders (front, side, back, top)\n\n"
                    f"Would you like to open the destination folder in File Explorer?"
                )
                if messagebox.askyesno("Export Complete", msg):
                    if sys.platform == "win32":
                        os.startfile(str(dest_dir))
                    else:
                        subprocess.run(["xdg-open", str(dest_dir)])

            return bundle
        except Exception as e:
            if not silent:
                messagebox.showerror("Export Error", f"Failed to export PyMOL session:\n{e}")
            else:
                import logging
                logging.getLogger("docking_automation.gui").error(
                    f"Auto-export PyMOL session failed: {e}"
                )
            return None

    def _export_pymol_session_ui(self) -> None:
        """Export PyMOL session — automatic no-dialog export to deterministic path.

        Previously this opened a directory picker dialog every time. It now
        automatically exports to results/pymol_sessions/<rec>_<lig>_pose<N>/
        so the workflow is seamless and repeatable for every ligand and pose.
        """
        self._auto_export_pymol_session(silent=False)

    def _fill_interactions(self, interactions) -> None:
        t = themes.get()
        # Color-code by interaction type
        type_colors = {
            "Hydrogen Bond": t.status_ok_text,
            "Salt Bridge": t.accent,
            "Pi-Stacking (Face-to-face)": t.text_secondary,
            "Pi-Stacking (T-shaped)": t.text_secondary,
            "Pi-Cation": t.text_secondary,
            "Halogen Bond": t.status_warn_text,
        }
        for i, inter in enumerate(interactions):
            tag = "even" if i % 2 == 0 else "odd"
            color = type_colors.get(inter.interaction_type, t.text_primary)
            iid = self.tree_inter.insert("", "end", tags=(tag,), values=(
                inter.pose_index,
                inter.interaction_type,
                inter.receptor_residue,
                inter.receptor_atom,
                inter.ligand_atom,
                f"{inter.distance_angstrom:.2f}",
                inter.details,
            ))
            self.tree_inter.tag_configure(
                f"type_{i}", foreground=color
            )

        # Now re-render the interaction diagram with the real contacts
        job = getattr(self, "_current_selected_job", None)
        if job is not None:
            pose_idx = getattr(self, "_current_selected_pose", 1)
            self._update_interaction_diagram(job, interactions=interactions, pose_index=pose_idx)

    def _placeholder_html(self, message: str) -> str:
        """Generate a minimal placeholder HTML string for the diagram panel."""
        t = themes.get()
        bg = t.bg_secondary
        fg = t.text_secondary
        return (
            f'<!DOCTYPE html><html><body style="margin:0;padding:0;background:{bg};">'
            f'<div style="display:flex;align-items:center;justify-content:center;'
            f'height:100vh;color:{fg};font-family:Segoe UI,Arial,sans-serif;font-size:13px;">'
            f'{message}</div></body></html>'
        )

    def _update_interaction_diagram(
        self, job: "DockingJob",
        interactions=None,  # already-computed list from _fill_interactions
        pose_index: int = 1,
    ) -> None:
        """Render and display the 2D protein-ligand interaction diagram.

        Uses a generation counter so that rapid ligand switching never shows
        stale results — any background callback whose generation no longer
        matches the current one is silently discarded.
        """
        self._render_generation += 1
        gen = self._render_generation

        ligand_path = self._resolve_ligand_path(job)

        # Show "Rendering…" immediately
        if self._has_html_frame:
            try:
                self._html_frame.load_html(self._placeholder_html("Rendering interaction diagram…"))
            except Exception:
                pass
        else:
            try:
                self.lbl_structure.configure(image=None, text="Rendering …")
            except Exception:
                pass

        def _worker() -> None:
            try:
                from gui.interaction_diagram import render_interaction_diagram_html
                t = themes.get()
                diagram_theme = "light" if t.ctk_mode == "light" else "dark"
                html = render_interaction_diagram_html(
                    ligand_path,
                    interactions or [],
                    obabel_exe=self.config.obabel_executable,
                    width=390, height=340,
                    theme=diagram_theme,
                )
                # Only apply if generation still matches (no newer selection arrived)
                if gen == self._render_generation:
                    if html and self._has_html_frame:
                        self.after(0, self._show_interaction_diagram, html)
                    elif html:
                        # Fallback: render the plain 2D ligand PNG instead
                        self.after(0, self._fallback_render_png, ligand_path)
                    else:
                        if self._has_html_frame:
                            self.after(0, self._html_frame.load_html,
                                       self._placeholder_html("2D diagram unavailable"))
                        else:
                            self.after(0, self.lbl_structure.configure,
                                       {"text": "2D diagram unavailable"})
            except Exception as exc:
                if gen == self._render_generation:
                    msg = "2D diagram unavailable"
                    if self._has_html_frame:
                        self.after(0, self._html_frame.load_html,
                                   self._placeholder_html(msg))
                    else:
                        try:
                            self.after(0, self.lbl_structure.configure, {"text": msg})
                        except Exception:
                            pass

        threading.Thread(target=_worker, daemon=True).start()

    def _show_interaction_diagram(self, html: str) -> None:
        """Display rendered HTML in the embedded HtmlFrame."""
        try:
            self._html_frame.load_html(html)
        except Exception:
            pass

    def _fallback_render_png(self, ligand_path) -> None:
        """Fallback PIL-based 2D ligand render when tkinterweb is unavailable."""
        if not ligand_path or not Path(ligand_path).is_file():
            return
        try:
            from depiction import render_ligand_2d
            png_bytes = render_ligand_2d(
                ligand_path, size=(280, 200),
                obabel_exe=self.config.obabel_executable,
            )
            if png_bytes:
                from io import BytesIO
                from PIL import Image, ImageTk
                img = Image.open(BytesIO(png_bytes))
                photo = ImageTk.PhotoImage(img)
                self._structure_image = photo
                self.lbl_structure.configure(image=photo, text="")
            else:
                self.lbl_structure.configure(text="2D structure unavailable")
        except Exception:
            try:
                self.lbl_structure.configure(text="2D structure unavailable")
            except Exception:
                pass

    # Keep old method name as a shim so any internal code that still calls it works
    def _update_2d_structure(self, job: "DockingJob") -> None:
        """Shim → delegates to _update_interaction_diagram (no pre-computed interactions)."""
        self._update_interaction_diagram(job, interactions=[])

    def _show_structure(self, png_bytes) -> None:
        """Legacy shim — kept for backwards compatibility."""
        if not png_bytes or not hasattr(self, "lbl_structure"):
            return
        try:
            from io import BytesIO
            from PIL import Image, ImageTk
            img = Image.open(BytesIO(png_bytes))
            photo = ImageTk.PhotoImage(img)
            self._structure_image = photo
            self.lbl_structure.configure(image=photo, text="")
        except Exception:
            try:
                self.lbl_structure.configure(text="PIL not available — install Pillow")
            except Exception:
                pass

    def _update_admet(self, job: "DockingJob") -> None:
        """Compute and display ADMET properties for the selected compound.\n        (unchanged)"""
        t = themes.get()
        ligand_path = self._resolve_ligand_path(job)
        if not ligand_path or not Path(ligand_path).is_file():
            return

        def _worker() -> None:
            try:
                from prepare import calculate_molecule_descriptors
                info = calculate_molecule_descriptors(ligand_path)
                self.after(0, self._fill_admet, info)
            except Exception:
                pass

        threading.Thread(target=_worker, daemon=True).start()

    def _resolve_ligand_path(self, job: "DockingJob") -> Optional[Path]:
        """Resolve the ligand source file (.sdf, .mol2, .pdbqt, or .pdb)."""
        if getattr(job, "ligand_path", None) and Path(job.ligand_path).is_file():
            return Path(job.ligand_path)

        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands")).resolve()
        root = getattr(self.config, "project_root", None)
        candidates = [
            lig_dir / f"{job.ligand_name}.sdf",
            lig_dir / f"{job.ligand_name}.mol2",
            lig_dir / f"{job.ligand_name}.pdbqt",
            lig_dir / f"{job.ligand_name}.pdb",
        ]
        if root:
            candidates.extend([
                Path(root) / "ligands" / f"{job.ligand_name}.sdf",
                Path(root) / "ligands" / f"{job.ligand_name}.mol2",
                Path(root) / "ligands" / f"{job.ligand_name}.pdbqt",
            ])
        for c in candidates:
            if c.is_file():
                return c
        return None

    def _fill_admet(self, info: Dict) -> None:
        t = themes.get()
        field_map = {
            "mw": f"{info.get('mw', 0):.2f} Da",
            "logp": f"{info.get('logp', 0):.2f}",
            "psa": f"{info.get('psa', 0):.1f} Å²",
            "hbd": str(info.get("hbd", 0)),
            "hba": str(info.get("hba", 0)),
            "rot_bonds": str(info.get("rot_bonds", 0)),
        }
        for key, val in field_map.items():
            if key in self._admet_labels:
                self._admet_labels[key].configure(text=val)

        ro5 = info.get("lipinski", "Pass")
        veber = info.get("veber_pass", True)

        self.lbl_ro5_badge.configure(
            text=f"✔ {ro5}" if ro5 == "Pass" else f"✖ {ro5}",
            text_color=t.status_ok_text if ro5 == "Pass" else t.status_fail_text,
        )
        self.lbl_veber_badge.configure(
            text="✔ Pass" if veber else "✖ Fail",
            text_color=t.status_ok_text if veber else t.status_fail_text,
        )

    def _set_action_status(self, text: str) -> None:
        """Update asynchronous action status text safely."""
        if hasattr(self, "lbl_action_status"):
            self.lbl_action_status.configure(text=text)

    def _export_report(self) -> None:
        """Trigger report generation for all completed jobs asynchronously in a worker thread."""
        from tkinter import messagebox
        try:
            if hasattr(self, "btn_export") and str(self.btn_export.cget("state")) == "disabled":
                return

            if hasattr(self, "btn_export"):
                self.btn_export.configure(state="disabled", text="⏳ Exporting...")
            if hasattr(self, "lbl_action_status"):
                self.lbl_action_status.configure(
                    text="⏳ Preparing reports...",
                    text_color=themes.get().accent_primary,
                )

            def _progress(msg: str) -> None:
                self.after(0, lambda m=msg: self._set_action_status(m))

            def _worker() -> None:
                try:
                    if sys.platform == "win32":
                        for _s in (sys.stdout, sys.stderr):
                            if hasattr(_s, "reconfigure"):
                                try:
                                    _s.reconfigure(encoding="utf-8", errors="replace")
                                except Exception:
                                    pass

                    from models import JobStatus
                    jobs = list(self._jobs) if self._jobs else []
                    if not jobs and getattr(self.config, "result_directory", None):
                        _progress("Scanning workspace for completed jobs...")
                        from job_manager import load_job_status
                        saved_status = load_job_status(self.config.result_directory)
                        if saved_status:
                            from job_manager import build_job_queue
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
                        self.after(0, self._on_export_no_jobs)
                        return

                    from reporting import generate_reports
                    res = generate_reports(jobs, self.config, progress_callback=_progress)
                    self.after(0, self._on_export_done, res, None)
                except Exception as ex:
                    self.after(0, self._on_export_done, None, ex)

            threading.Thread(target=_worker, daemon=True).start()

        except Exception as e:
            if hasattr(self, "btn_export"):
                self.btn_export.configure(state="normal", text="📤 Export Report")
            if hasattr(self, "lbl_action_status"):
                self.lbl_action_status.configure(text="")
            messagebox.showerror("Export Error", f"Failed to start report export:\n{e}")

    def _on_export_no_jobs(self) -> None:
        """Handle case where no completed jobs are available for report export."""
        from tkinter import messagebox
        if hasattr(self, "btn_export"):
            self.btn_export.configure(state="normal", text="📤 Export Report")
        if hasattr(self, "lbl_action_status"):
            self.lbl_action_status.configure(text="")
        messagebox.showwarning(
            "No Docking Results",
            "There are no docking results or completed jobs available to export.\n\n"
            "Please run docking first before exporting a report.",
        )

    def _on_export_done(self, res: Optional[Dict], error: Optional[Exception]) -> None:
        """Handle asynchronous report generation completion."""
        from tkinter import messagebox
        if hasattr(self, "btn_export"):
            self.btn_export.configure(state="normal", text="📤 Export Report")
        if hasattr(self, "lbl_action_status"):
            self.lbl_action_status.configure(text="")

        if error:
            messagebox.showerror("Export Error", f"Failed to export report:\n{error}")
            return

        xlsx_path = res.get("xlsx_path") if isinstance(res, dict) else None
        report_dir = res.get("report_dir") if isinstance(res, dict) else self.config.report_directory

        msg = "Docking report generated successfully!\n\n"
        if xlsx_path and Path(xlsx_path).is_file():
            msg += f"Excel Report:\n{xlsx_path}\n\n"
        msg += f"Reports Directory:\n{report_dir}\n\nWould you like to open the reports folder?"

        if messagebox.askyesno("Export Complete", msg):
            if report_dir and Path(report_dir).is_dir():
                if sys.platform == "win32":
                    os.startfile(str(report_dir))
                else:
                    subprocess.run(["xdg-open", str(report_dir)])

    def _run_dlg_extract(self) -> None:
        """Trigger BSNDVP DLG post-processing extraction on the current workspace."""
        if hasattr(self, "btn_dlg_extract") and str(self.btn_dlg_extract.cget("state")) == "disabled":
            return

        root = getattr(self.config, "project_root", None)
        res_dir = getattr(self.config, "result_directory", None)

        # Locate DLG folder (prioritizing nested results/AD4 without needing flat copies)
        dlg_dir = None
        candidates = []
        if res_dir:
            candidates.extend([Path(res_dir) / "AD4", Path(res_dir) / "DLG"])
        if root:
            candidates.extend([Path(root) / "results" / "AD4", Path(root) / "DLG", Path(root) / "results" / "DLG"])

        for cand in candidates:
            if cand.is_dir() and any(cand.rglob("*.dlg")):
                dlg_dir = cand
                break

        if not dlg_dir and root:
            try:
                import importlib.util
                dlg_extract_path = Path(__file__).resolve().parent.parent / "dlg_extract.py"
                spec = importlib.util.spec_from_file_location("dlg_extract", str(dlg_extract_path))
                mod = importlib.util.module_from_spec(spec)
                sys.modules["dlg_extract"] = mod
                spec.loader.exec_module(mod)
                dlg_dir = mod._find_dlg_dir(Path(root))
            except Exception:
                pass

        if not dlg_dir:
            from tkinter import messagebox
            messagebox.showwarning(
                "No DLG Files Found",
                "No AutoDock 4 .dlg files were found in the current workspace.\n\n"
                f"Searched in:\n"
                + "\n".join(f"• {c}" for c in candidates)
                + "\n\nPlease run AutoDock 4 docking first.",
            )
            return

        if hasattr(self, "btn_dlg_extract"):
            self.btn_dlg_extract.configure(state="disabled", text="⏳ Extracting...")
        if hasattr(self, "lbl_action_status"):
            self.lbl_action_status.configure(
                text="⏳ Extracting DLG results...",
                text_color=themes.get().accent_primary,
            )

        # Always save BSNDVP_RESULTS to workspace root, keeping results/ clean (VINA and AD4 only)
        out_dir = (Path(root) / "BSNDVP_RESULTS") if root else (Path(res_dir).parent / "BSNDVP_RESULTS" if res_dir else Path("BSNDVP_RESULTS"))
        out_dir.mkdir(parents=True, exist_ok=True)

        def _worker():
            try:
                if sys.platform == "win32":
                    for _s in (sys.stdout, sys.stderr):
                        if hasattr(_s, "reconfigure"):
                            try:
                                _s.reconfigure(encoding="utf-8", errors="replace")
                            except Exception:
                                pass

                import importlib.util
                dlg_extract_path = Path(__file__).resolve().parent.parent / "dlg_extract.py"
                spec = importlib.util.spec_from_file_location("dlg_extract", str(dlg_extract_path))
                mod = importlib.util.module_from_spec(spec)
                sys.modules["dlg_extract"] = mod
                spec.loader.exec_module(mod)

                pipeline = mod.DLGPipeline(dlg_dir, out_dir, root or Path.cwd())
                pipeline.run()

                self.after(0, self._on_dlg_extract_done, None, out_dir)
            except Exception as e:
                self.after(0, self._on_dlg_extract_done, e, out_dir)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_dlg_extract_done(self, error: Optional[Exception], out_dir: Path) -> None:
        from tkinter import messagebox
        if hasattr(self, "btn_dlg_extract"):
            self.btn_dlg_extract.configure(state="normal", text="🔬 Run DLG Extract")
        if hasattr(self, "lbl_action_status"):
            self.lbl_action_status.configure(text="")

        if error:
            messagebox.showerror("DLG Extract Error", f"Failed to extract DLG results:\n{error}")
        else:
            self.refresh()
            msg = (
                f"BSNDVP™ DLG analysis finished successfully!\n\n"
                f"Results saved to:\n{out_dir}\n\n"
                f"Would you like to open the results folder?"
            )
            if messagebox.askyesno("DLG Extract Complete", msg):
                if sys.platform == "win32":
                    os.startfile(str(out_dir))
                else:
                    subprocess.run(["xdg-open", str(out_dir)])

    def _open_dlg_folder(self) -> None:
        root = getattr(self.config, "project_root", Path.cwd())
        res_dir = getattr(self.config, "result_directory", None)
        target = None
        if root and (Path(root) / "BSNDVP_RESULTS").is_dir():
            target = Path(root) / "BSNDVP_RESULTS"
        elif res_dir and (Path(res_dir) / "AD4").is_dir():
            target = Path(res_dir) / "AD4"
        elif root and (Path(root) / "results" / "AD4").is_dir():
            target = Path(root) / "results" / "AD4"
        elif root:
            target = Path(root) / "BSNDVP_RESULTS"
        else:
            target = Path("BSNDVP_RESULTS")
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

    def _get_ad4_affinity(self, job: "DockingJob") -> Optional[float]:
        """Extract best binding affinity for an AutoDock4 job."""
        res_dir = getattr(self.config, "result_directory", None)
        root = getattr(self.config, "project_root", None)

        csv_candidates = []
        if res_dir:
            csv_candidates.append(Path(res_dir) / "BSNDVP_RESULTS" / "03_Ligand_Summary.csv")
        if root:
            csv_candidates.append(Path(root) / "BSNDVP_RESULTS" / "03_Ligand_Summary.csv")

        for sum_csv in csv_candidates:
            if sum_csv.is_file():
                try:
                    lines = sum_csv.read_text(encoding="utf-8", errors="replace").splitlines()
                    for line in lines:
                        if not line.strip() or line.startswith("#") or line.startswith("---"):
                            continue
                        parts = [p.strip() for p in line.split(",")]
                        if len(parts) > 4 and parts[1].lower() == job.ligand_name.lower():
                            try:
                                return float(parts[4])
                            except ValueError:
                                pass
                except Exception:
                    pass

        dlg_p = getattr(job, "dlg_path", None) or (job.ad4_results.get("dlg_path") if getattr(job, "ad4_results", None) else None)
        if not dlg_p:
            dlg_candidates = []
            if res_dir:
                dlg_candidates.append(Path(res_dir) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                dlg_candidates.append(Path(res_dir) / f"{job.receptor_name}_{job.ligand_name}.dlg")
            if root:
                dlg_candidates.append(Path(root) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                dlg_candidates.append(Path(root) / "results" / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
            for candidate in dlg_candidates:
                if candidate.is_file():
                    dlg_p = candidate
                    break

        if dlg_p and Path(dlg_p).is_file():
            try:
                lines = Path(dlg_p).read_text(encoding="utf-8", errors="replace").splitlines()
                for line in lines:
                    if "Estimated Free Energy of Binding" in line:
                        m = re.search(r"=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)\s*kcal/mol", line)
                        if m:
                            return float(m.group(1))
            except Exception:
                pass
        return None

    def _get_ad4_poses(self, job: "DockingJob") -> List[Dict]:
        """Extract pose/cluster rows for an AutoDock4 job."""
        poses = []
        res_dir = getattr(self.config, "result_directory", None)
        root = getattr(self.config, "project_root", None)

        csv_candidates = []
        if res_dir:
            csv_candidates.append(Path(res_dir) / "BSNDVP_RESULTS" / "01_Docking_Master.csv")
        if root:
            csv_candidates.append(Path(root) / "BSNDVP_RESULTS" / "01_Docking_Master.csv")

        for master_csv in csv_candidates:
            if master_csv.is_file():
                try:
                    lines = master_csv.read_text(encoding="utf-8", errors="replace").splitlines()
                    for line in lines[5:]:
                        if not line.strip() or line.startswith("---"):
                            continue
                        parts = line.split(",")
                        if len(parts) > 15 and parts[1].strip().lower() == job.ligand_name.lower():
                            poses.append({
                                "pose": f"Run {parts[3].strip()}",
                                "affinity": f"{float(parts[4].strip()):.2f}" if parts[4].strip() not in ("N/A", "") else "—",
                                "ki": parts[5].strip(),
                                "rmsd_lb": f"{float(parts[13].strip()):.3f}" if parts[13].strip() not in ("N/A", "") else "—",
                                "rmsd_ub": f"{float(parts[14].strip()):.3f}" if parts[14].strip() not in ("N/A", "") else "—",
                            })
                    if poses:
                        return poses
                except Exception:
                    pass

        if not poses:
            # Fallback to direct DLG parse
            dlg_p = getattr(job, "dlg_path", None) or (job.ad4_results.get("dlg_path") if getattr(job, "ad4_results", None) else None)
            if not dlg_p:
                dlg_candidates = []
                if res_dir:
                    dlg_candidates.append(Path(res_dir) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                if root:
                    dlg_candidates.append(Path(root) / "DLG" / f"{job.receptor_name}_{job.ligand_name}.dlg")
                for candidate in dlg_candidates:
                    if candidate.is_file():
                        dlg_p = candidate
                        break

            if dlg_p and Path(dlg_p).is_file():
                try:
                    lines = Path(dlg_p).read_text(encoding="utf-8", errors="replace").splitlines()
                    run_id = 0
                    cur_energy = None
                    cur_ki = "—"
                    for line in lines:
                        m_run = re.search(r"\[\s*Run\s+(\d+)\s+of\s+(\d+)", line)
                        if m_run:
                            run_id = int(m_run.group(1))
                        if "Estimated Free Energy of Binding" in line:
                            m_eng = re.search(r"=\s*([+-]?\d+\.?\d*)\s*kcal/mol", line)
                            if m_eng:
                                cur_energy = m_eng.group(1)
                        if "Estimated Inhibition Constant" in line:
                            m_ki = re.search(r"=\s*([^\s\[]+(?:\s*[µunm]?M)?)", line)
                            if m_ki:
                                cur_ki = m_ki.group(1)
                            if cur_energy is not None:
                                poses.append({
                                    "pose": f"Run {run_id}",
                                    "affinity": cur_energy,
                                    "ki": cur_ki,
                                    "rmsd_lb": "—",
                                    "rmsd_ub": "—",
                                })
                                cur_energy = None
                                cur_ki = "—"
                except Exception:
                    pass
        return poses
