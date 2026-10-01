"""
AutoDock Suite Pro — Workspace Tab (gui/workspace_tab.py)
=========================================================
First tab the user sees. Manages receptors, ligands, and
the project workspace folder. Provides Add/Remove buttons
that handle folder conventions automatically.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import TYPE_CHECKING, Callable, List, Optional

import customtkinter as ctk

from gui import themes, widgets

if TYPE_CHECKING:
    from config import ProjectConfig


class WorkspaceTab(ctk.CTkFrame):
    """Workspace management tab: Add/Remove receptors and ligands."""

    def __init__(self, parent: ctk.CTkTabview, config: "ProjectConfig",
                 on_config_change: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self._build()
        self.refresh()

    # ─────────────────────────────────────────────────
    # Build
    # ─────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        # ── Workspace path ──────────────────────────
        widgets.section_header(self, "Project Workspace", row=0, colspan=2)

        ws_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        ws_frame.grid(row=1, column=0, columnspan=2, sticky="ew", padx=12, pady=(0, 8))
        ws_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(ws_frame, text="Workspace folder:",
                     text_color=t.text_secondary, font=themes.font_body()).grid(row=0, column=0, padx=10, pady=8, sticky="w")
        self.lbl_workspace = ctk.CTkLabel(
            ws_frame, text="", text_color=t.text_primary,
            font=themes.font_body(), anchor="w",
        )
        self.lbl_workspace.grid(row=0, column=1, sticky="ew", padx=4)
        widgets.accent_button(ws_frame, "Scaffold Workspace",
                              self._scaffold_workspace, width=170, small=True
                              ).grid(row=0, column=2, padx=(0, 4), pady=6)
        widgets.accent_button(ws_frame, "Change Workspace",
                              self._change_workspace, width=160, small=True
                              ).grid(row=0, column=3, padx=(0, 10), pady=6)

        widgets.divider(self, row=2, colspan=2)

        # ── Receptors column ────────────────────────
        widgets.section_header(self, "Receptors", row=3, col=0, colspan=1)

        rec_btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        rec_btn_frame.grid(row=4, column=0, sticky="w", padx=12, pady=(0, 4))

        widgets.accent_button(rec_btn_frame, "➕ Add PDBQT",
                              self._add_receptor, width=120, small=True
                              ).grid(row=0, column=0, padx=(0, 4))
        widgets.accent_button(rec_btn_frame, "🗑 Remove",
                              self._remove_receptor, width=100, small=True
                              ).grid(row=0, column=1)

        rec_list_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        rec_list_frame.grid(row=5, column=0, sticky="nsew", padx=12, pady=(0, 12))
        rec_list_frame.grid_rowconfigure(0, weight=1)
        rec_list_frame.grid_columnconfigure(0, weight=1)

        self.lb_receptors = ctk.CTkTextbox(
            rec_list_frame, state="disabled", height=200,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            font=themes.font_code(),
        )
        self.lb_receptors.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)

        # ── Ligands column ──────────────────────────
        widgets.section_header(self, "Ligands", row=3, col=1, colspan=1)

        lig_btn_frame = ctk.CTkFrame(self, fg_color="transparent")
        lig_btn_frame.grid(row=4, column=1, sticky="w", padx=12, pady=(0, 4))

        widgets.accent_button(lig_btn_frame, "➕ Add PDBQTs",
                              self._add_ligands, width=120, small=True
                              ).grid(row=0, column=0, padx=(0, 4))
        widgets.accent_button(lig_btn_frame, "🗑 Clear All",
                              self._clear_ligands, width=100, small=True
                              ).grid(row=0, column=1)
        ctk.CTkButton(
            lig_btn_frame, text="📁 Open Folder",
            command=self._open_ligand_folder,
            width=110, fg_color="transparent",
            border_width=1, border_color=themes.get().border,
            text_color=themes.get().text_secondary,
            font=themes.font_button(),
        ).grid(row=0, column=2, padx=(4, 0))

        lig_list_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        lig_list_frame.grid(row=5, column=1, sticky="nsew", padx=12, pady=(0, 12))
        lig_list_frame.grid_rowconfigure(0, weight=1)
        lig_list_frame.grid_columnconfigure(0, weight=1)

        self.lb_ligands = ctk.CTkTextbox(
            lig_list_frame, state="disabled", height=200,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            font=themes.font_code(),
        )
        self.lb_ligands.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)

        self.grid_rowconfigure(5, weight=1)

        # ── Status bar ──────────────────────────────
        widgets.divider(self, row=6, colspan=2)
        self.lbl_status = ctk.CTkLabel(
            self, text="", text_color=t.text_secondary,
            font=themes.font_caption(), anchor="w",
        )
        self.lbl_status.grid(row=7, column=0, columnspan=2, sticky="w", padx=12, pady=4)

    # ─────────────────────────────────────────────────
    # Refresh
    # ─────────────────────────────────────────────────

    def refresh(self) -> None:
        """Re-scan workspace and update receptor/ligand lists."""
        self.lbl_workspace.configure(text=str(self.config.project_root.resolve()))

        rec_dir = self.config.receptor_directory
        lig_dir = self.config.ligand_directory

        # Receptors
        rec_lines: List[str] = []
        if rec_dir.exists():
            for item in sorted(rec_dir.iterdir()):
                if not item.is_dir():
                    continue
                has_rigid = any((item / "rigid").glob("*.pdbqt")) if (item / "rigid").is_dir() else False
                has_flex = any((item / "flex").glob("*.pdbqt")) if (item / "flex").is_dir() else False
                has_config = (item / "config.txt").is_file()

                if has_rigid and has_flex:
                    mode = "[RIGID+FLEX]"
                elif has_flex:
                    mode = "[FLEX]     "
                elif has_rigid:
                    mode = "[RIGID]    "
                else:
                    mode = "[NO PDBQT] "

                cfg_flag = "✔ grid" if has_config else "⚠ no grid"
                ok_sym = "✔" if (has_rigid or has_flex) else "✖"
                rec_lines.append(f" {ok_sym}  {item.name:<24s}  {mode}  {cfg_flag}")

        rec_text = "\n".join(rec_lines) if rec_lines else " (no receptors found)"
        self._set_textbox(self.lb_receptors, rec_text)

        # Ligands
        lig_count = 0
        lig_lines: List[str] = []
        if lig_dir.exists():
            pdbqts = sorted(lig_dir.glob("*.pdbqt"))
            lig_count = len(pdbqts)
            for p in pdbqts[:50]:
                lig_lines.append(f"  {p.name}")
            if lig_count > 50:
                lig_lines.append(f"  ... and {lig_count - 50} more")

        lig_text = "\n".join(lig_lines) if lig_lines else " (no ligands found)"
        self._set_textbox(self.lb_ligands, lig_text)

        # Status bar
        n_rec = len(rec_lines)
        self.lbl_status.configure(
            text=f"  {n_rec} receptor(s)  |  {lig_count} ligand(s) ready  |  "
                 f"Project: {self.config.project_name}"
        )

    def _set_textbox(self, tb: ctk.CTkTextbox, text: str) -> None:
        tb.configure(state="normal")
        tb.delete("1.0", "end")
        tb.insert("end", text)
        tb.configure(state="disabled")

    # ─────────────────────────────────────────────────
    # Actions
    # ─────────────────────────────────────────────────

    def _change_workspace(self) -> None:
        folder = filedialog.askdirectory(title="Select Project Workspace Folder")
        if not folder:
            return
        from config import set_workspace_directory
        set_workspace_directory(self.config, Path(folder))
        self.on_config_change()
        self.refresh()

    def _scaffold_workspace(self) -> None:
        """Create the standard folder scaffold in a selected or current workspace."""
        root = self.config.project_root
        suite_root = Path(__file__).resolve().parent.parent
        is_suite_dir = root and (root.resolve() == suite_root.resolve())

        if not root or not root.exists() or is_suite_dir:
            folder = filedialog.askdirectory(title="Select Folder to Create New Workspace Scaffold")
            if not folder:
                return
            root = Path(folder)
        else:
            ans = messagebox.askyesnocancel(
                "Scaffold Workspace",
                f"Create scaffold in the current workspace folder:\n{root}?\n\n"
                f"• Click 'Yes' to scaffold in {root.name}\n"
                f"• Click 'No' to choose a DIFFERENT folder\n"
                f"• Click 'Cancel' to abort",
            )
            if ans is None:
                return
            if not ans:
                folder = filedialog.askdirectory(title="Select Folder to Create New Workspace Scaffold")
                if not folder:
                    return
                root = Path(folder)

        from config import set_workspace_directory
        from prepare import create_scaffold

        create_scaffold(root)
        set_workspace_directory(self.config, root)

        self.on_config_change()
        self.refresh()
        self.lbl_status.configure(text=f"  Scaffold created and active: {root}")

        has_bin = (root / "bin").is_dir()
        bin_msg = "Binaries: bin/ (AutoGrid4, AutoDock4, Vina, OpenBabel copied)\n" if has_bin else ""
        messagebox.showinfo(
            "Scaffold Created & Activated",
            f"Workspace structure created and ACTIVATED in:\n{root}\n\n"
            f"Folders:  receptors/  ligands/  results/ (with DLG/, DPF/, VINA/)  DLG/  BSNDVP_RESULTS/  logs/  reports/\n"
            f"{bin_msg}"
            f"Config:   project_config.toml\n"
            f"Post-processing: dlg_extract.py (ready in workspace)\n\n"
            "All docking results, DLGs, and DLG extraction reports will be saved directly in this folder.",
        )

        # Open Explorer at the new workspace
        try:
            import subprocess as sp
            sp.Popen(["explorer", str(root)])
        except Exception:
            pass

    def _add_receptor(self) -> None:
        """Browse to a receptor PDBQT and place it in receptors/<name>/rigid/."""
        files = filedialog.askopenfilenames(
            title="Select Rigid Receptor PDBQT file(s)",
            filetypes=[("PDBQT Files", "*.pdbqt"), ("All Files", "*.*")],
        )
        if not files:
            return
        for src in files:
            src_path = Path(src)
            name = src_path.stem
            dest_dir = self.config.receptor_directory / name / "rigid"
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, dest_dir / src_path.name)

            # Create stub config.txt if missing
            config_path = self.config.receptor_directory / name / "config.txt"
            if not config_path.is_file():
                config_path.write_text(
                    "# EDIT THIS FILE: set the correct grid center and size.\n"
                    "center_x = 0.000\ncenter_y = 0.000\ncenter_z = 0.000\n"
                    "size_x = 25.000\nsize_y = 25.000\nsize_z = 25.000\n"
                    "exhaustiveness = 8\nnum_modes = 9\nenergy_range = 3.0\n",
                    encoding="utf-8",
                )
        self.refresh()
        messagebox.showinfo(
            "Receptors Added",
            f"{len(files)} receptor(s) added.\n"
            "Remember to set the Grid Box in the Grid Box tab before docking.",
        )

    def _remove_receptor(self) -> None:
        """Remove a receptor folder after confirmation."""
        rec_dir = self.config.receptor_directory
        if not rec_dir.exists():
            return
        receptors = [d.name for d in rec_dir.iterdir() if d.is_dir()]
        if not receptors:
            messagebox.showinfo("No Receptors", "No receptor folders found.")
            return

        # Simple dialog: ask user to type the name
        from tkinter.simpledialog import askstring
        name = askstring(
            "Remove Receptor",
            f"Available receptors:\n{chr(10).join(receptors)}\n\nEnter name to remove:",
        )
        if not name or name not in receptors:
            return
        if not messagebox.askyesno("Confirm", f"Delete receptor folder '{name}' and all its files?"):
            return
        shutil.rmtree(rec_dir / name)
        self.refresh()

    def _add_ligands(self) -> None:
        """Browse for multiple PDBQT files and copy them to ligands/."""
        files = filedialog.askopenfilenames(
            title="Select Ligand PDBQT file(s)",
            filetypes=[("PDBQT Files", "*.pdbqt"), ("All Files", "*.*")],
        )
        if not files:
            return
        lig_dir = self.config.ligand_directory
        lig_dir.mkdir(parents=True, exist_ok=True)
        for src in files:
            shutil.copy2(src, lig_dir / Path(src).name)
        self.refresh()

    def _clear_ligands(self) -> None:
        """Delete all ligands from the ligands/ folder after confirmation."""
        lig_dir = self.config.ligand_directory
        count = len(list(lig_dir.glob("*.pdbqt"))) if lig_dir.exists() else 0
        if count == 0:
            messagebox.showinfo("No Ligands", "Ligands folder is already empty.")
            return
        if not messagebox.askyesno("Confirm", f"Delete all {count} ligand PDBQT file(s)?"):
            return
        for f in lig_dir.glob("*.pdbqt"):
            f.unlink()
        self.refresh()

    def _open_ligand_folder(self) -> None:
        """Open the ligands/ folder in Windows Explorer."""
        lig_dir = self.config.ligand_directory
        lig_dir.mkdir(parents=True, exist_ok=True)
        try:
            os.startfile(str(lig_dir))
        except AttributeError:
            subprocess.Popen(["explorer", str(lig_dir)])
