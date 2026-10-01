"""
AutoDock Suite Pro -- Grid Box Tab (gui/grid_tab.py)
====================================================
Provides:
  1. Grid center / size input (6 spinboxes with +/- step controls)
  2. Auto-fill from co-crystallized ligand bounding box
  3. Generate Vina config.txt
  4. Generate AutoDock4 GPF
  5. Run AutoGrid4 with live scrolling log
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import TYPE_CHECKING, Callable, Dict, List, Optional

import customtkinter as ctk

from gui import themes, widgets

if TYPE_CHECKING:
    from config import ProjectConfig


def _make_spinbox(parent, label: str, default: float,
                  step: float = 1.0, lo: float = -999, hi: float = 999,
                  width: int = 90) -> ctk.CTkEntry:
    """Return a labelled entry that behaves like a spinbox."""
    t = themes.get()
    row_f = ctk.CTkFrame(parent, fg_color="transparent")
    row_f.pack(side="left", padx=6)
    ctk.CTkLabel(row_f, text=label, text_color=t.text_secondary,
                 font=themes.font_caption()).pack()
    ent = ctk.CTkEntry(row_f, width=width, justify="center",
                       font=themes.font_code_bold())
    ent.insert(0, str(default))
    ent.pack()
    btn_row = ctk.CTkFrame(row_f, fg_color="transparent")
    btn_row.pack()
    for sym, delta in [("-", -step), ("+", step)]:
        def _cb(d=delta, e=ent):
            try:
                v = float(e.get()) + d
                v = max(lo, min(hi, round(v, 3)))
                e.delete(0, "end")
                e.insert(0, str(v))
            except ValueError:
                pass
        ctk.CTkButton(btn_row, text=sym, width=32, height=22,
                      fg_color=themes.get().bg_tertiary,
                      text_color=themes.get().text_primary,
                      font=themes.font_button(),
                      command=_cb).pack(side="left", padx=1)
    return ent


class GridBoxTab(ctk.CTkFrame):
    """Grid Box configuration, file generation, and AutoGrid4 runner."""

    def __init__(self, parent, config: "ProjectConfig",
                 on_config_change: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self._receptor_name: Optional[str] = None
        self._build()

    # ── Build ────────────────────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        widgets.section_header(self, "Grid Box Configuration", row=0, colspan=2)

        # Left: grid inputs + file generators
        left = ctk.CTkFrame(self, fg_color="transparent")
        left.grid(row=1, column=0, sticky="nsew", padx=(12, 6), pady=4)
        left.grid_columnconfigure(0, weight=1)
        self._build_grid_inputs(left)
        self._build_coligand_section(left)
        self._build_generate_section(left)

        # Right: AutoGrid runner
        right = ctk.CTkFrame(self, fg_color="transparent")
        right.grid(row=1, column=1, sticky="nsew", padx=(6, 12), pady=4)
        right.grid_columnconfigure(0, weight=1)
        right.grid_rowconfigure(1, weight=1)
        self._build_autogrid_section(right)

        widgets.divider(self, row=2, colspan=2)
        self.lbl_status = ctk.CTkLabel(
            self, text="  Select a receptor and set grid parameters.",
            text_color=t.text_secondary, font=themes.font_caption(), anchor="w")
        self.lbl_status.grid(row=3, column=0, columnspan=2, sticky="w", padx=12, pady=4)

    # ── Grid inputs section ───────────────────────────────────────────────────

    def _build_grid_inputs(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.pack(fill="x", pady=(0, 8))

        ctk.CTkLabel(card, text="Grid Box Parameters",
                     text_color=t.accent, font=themes.font_section(),
                     anchor="w").pack(fill="x", padx=10, pady=(8, 4))

        # Receptor picker
        rec_row = ctk.CTkFrame(card, fg_color="transparent")
        rec_row.pack(fill="x", padx=8, pady=4)
        ctk.CTkLabel(rec_row, text="Receptor:", text_color=t.text_secondary,
                     font=themes.font_body()).pack(side="left", padx=(0, 6))
        self.var_receptor = ctk.StringVar(value="(select receptor)")
        self.opt_receptor = ctk.CTkOptionMenu(
            rec_row, variable=self.var_receptor,
            values=self._list_receptors(),
            command=self._on_receptor_change,
            width=220, height=30,
            font=themes.font_body(), dropdown_font=themes.font_body())
        self.opt_receptor.pack(side="left")
        ctk.CTkButton(rec_row, text="Refresh", width=75, height=28,
                      fg_color="transparent", border_width=1,
                      border_color=t.border, text_color=t.text_secondary,
                      font=themes.font_button(),
                      command=self._refresh_receptor_list
                      ).pack(side="left", padx=(6, 0))
        ctk.CTkButton(rec_row, text="📁 Folder", width=75, height=28,
                      fg_color="transparent", border_width=1,
                      border_color=t.border, text_color=t.text_secondary,
                      font=themes.font_button(),
                      command=self._open_receptor_folder
                      ).pack(side="left", padx=(4, 0))

        # Center row
        ctr = ctk.CTkFrame(card, fg_color="transparent")
        ctr.pack(pady=(4, 0))
        ctk.CTkLabel(ctr, text="Grid Center (Angstroms):",
                     text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).pack(anchor="w")
        center_row = ctk.CTkFrame(ctr, fg_color="transparent")
        center_row.pack()
        self.ent_cx = _make_spinbox(center_row, "center_x", 0.0, step=0.5)
        self.ent_cy = _make_spinbox(center_row, "center_y", 0.0, step=0.5)
        self.ent_cz = _make_spinbox(center_row, "center_z", 0.0, step=0.5)

        # Size row
        szf = ctk.CTkFrame(card, fg_color="transparent")
        szf.pack(pady=(8, 0))
        ctk.CTkLabel(szf, text="Grid Size (Angstroms):",
                     text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).pack(anchor="w")
        size_row = ctk.CTkFrame(szf, fg_color="transparent")
        size_row.pack()
        self.ent_sx = _make_spinbox(size_row, "size_x", 25.0, step=1.0, lo=1)
        self.ent_sy = _make_spinbox(size_row, "size_y", 25.0, step=1.0, lo=1)
        self.ent_sz = _make_spinbox(size_row, "size_z", 25.0, step=1.0, lo=1)

        # GPF extra params
        gpf_row = ctk.CTkFrame(card, fg_color="transparent")
        gpf_row.pack(fill="x", padx=8, pady=(8, 8))
        for label, attr, default, w in [
            ("spacing (A)", "ent_spacing", "0.375", 60),
            ("smooth",      "ent_smooth",  "0.5",   60),
            ("dielectric",  "ent_dielec",  "-0.1465", 80),
        ]:
            f = ctk.CTkFrame(gpf_row, fg_color="transparent")
            f.pack(side="left", padx=6)
            ctk.CTkLabel(f, text=label, text_color=t.text_secondary,
                         font=themes.font_caption()).pack()
            ent = ctk.CTkEntry(f, width=w, justify="center",
                               font=themes.font_code())
            ent.insert(0, default)
            ent.pack()
            setattr(self, attr, ent)

    # ── Co-ligand section ─────────────────────────────────────────────────────

    def _build_coligand_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.pack(fill="x", pady=(0, 8))

        ctk.CTkLabel(card, text="Auto-fill from Co-crystallized Ligand",
                     text_color=t.accent, font=themes.font_section(),
                     anchor="w").pack(fill="x", padx=10, pady=(8, 4))

        inner = ctk.CTkFrame(card, fg_color="transparent")
        inner.pack(fill="x", padx=8, pady=(0, 8))
        inner.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(inner, text="Ligand file:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=0, padx=(0, 6))
        self.lbl_colig = ctk.CTkLabel(inner, text="(none selected)",
                                      text_color=t.text_secondary,
                                      font=themes.font_code(),
                                      anchor="w")
        self.lbl_colig.grid(row=0, column=1, sticky="ew")
        widgets.accent_button(inner, "Browse", self._browse_coligand,
                              width=80, small=True).grid(row=0, column=2, padx=(8, 0))

        pad_row = ctk.CTkFrame(card, fg_color="transparent")
        pad_row.pack(fill="x", padx=8, pady=(0, 4))
        ctk.CTkLabel(pad_row, text="Padding (A):", text_color=t.text_secondary,
                     font=themes.font_body()).pack(side="left", padx=(0, 6))
        self.ent_padding = ctk.CTkEntry(pad_row, width=60, justify="center",
                                         font=themes.font_code())
        self.ent_padding.insert(0, "5.0")
        self.ent_padding.pack(side="left")
        ctk.CTkLabel(pad_row, text="added on each side",
                     text_color=t.text_secondary,
                     font=themes.font_caption()).pack(side="left", padx=(6, 0))

        widgets.accent_button(card, "Get Dimensions from Co-ligand",
                              self._autofill_from_coligand, width=250
                              ).pack(pady=(0, 10))

    # ── Generate config files section ─────────────────────────────────────────

    def _build_generate_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.pack(fill="x")

        ctk.CTkLabel(card, text="Generate Configuration Files",
                     text_color=t.accent, font=themes.font_section(),
                     anchor="w").pack(fill="x", padx=10, pady=(8, 4))

        btn_row = ctk.CTkFrame(card, fg_color="transparent")
        btn_row.pack(pady=(0, 10))
        widgets.accent_button(btn_row, "Generate config.txt  (Vina)",
                              self._gen_vina_config, width=175
                              ).pack(side="left", padx=4)
        widgets.accent_button(btn_row, "Generate GPF  (AutoDock4)",
                              self._gen_gpf, width=175
                              ).pack(side="left", padx=4)
        widgets.accent_button(btn_row, "Generate DPF  (AutoDock4)",
                              self._gen_dpf, width=175
                              ).pack(side="left", padx=4)

    # ── AutoGrid runner section ────────────────────────────────────────────────

    def _build_autogrid_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.grid(row=0, column=0, sticky="nsew")
        card.grid_columnconfigure(0, weight=1)
        card.grid_rowconfigure(2, weight=1)
        parent.grid_rowconfigure(0, weight=1)

        ctk.CTkLabel(card, text="Run AutoGrid4",
                     text_color=t.accent, font=ctk.CTkFont(size=12, weight="bold"),
                     anchor="w").grid(row=0, column=0, sticky="w", padx=10, pady=(8, 4))

        btn_row = ctk.CTkFrame(card, fg_color="transparent")
        btn_row.grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 4))
        btn_row.grid_columnconfigure(1, weight=1)
        widgets.accent_button(btn_row, "Run AutoGrid4",
                              self._run_autogrid, width=140
                              ).grid(row=0, column=0)
        self.lbl_ag_status = ctk.CTkLabel(
            btn_row, text="",
            text_color=t.text_secondary, font=themes.font_caption(), anchor="w")
        self.lbl_ag_status.grid(row=0, column=1, sticky="w", padx=(10, 0))

        # Font zoom controls on right
        ag_zoom_frame = ctk.CTkFrame(btn_row, fg_color="transparent")
        ag_zoom_frame.grid(row=0, column=2, sticky="e")
        ctk.CTkLabel(ag_zoom_frame, text="Font:", text_color=t.text_secondary,
                     font=themes.font_caption()).pack(side="left", padx=(2, 2))
        ctk.CTkButton(
            ag_zoom_frame, text="A -", width=34, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=themes.font_button(),
            command=self._zoom_out_ag,
        ).pack(side="left", padx=1)
        self.lbl_ag_font = ctk.CTkLabel(
            ag_zoom_frame, text="13 pt", text_color=t.text_primary,
            font=themes.font_caption_bold(), width=40)
        self.lbl_ag_font.pack(side="left", padx=1)
        ctk.CTkButton(
            ag_zoom_frame, text="A +", width=34, height=24,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            hover_color=t.border, font=themes.font_button(),
            command=self._zoom_in_ag,
        ).pack(side="left", padx=1)

        self._ag_font_size = 13
        self.tb_autogrid = ctk.CTkTextbox(
            card, state="disabled", height=400,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            font=ctk.CTkFont(family="Consolas", size=self._ag_font_size))
        self.tb_autogrid.grid(row=2, column=0, sticky="nsew", padx=8, pady=(0, 8))

    # ── Receptor helpers ──────────────────────────────────────────────────────

    def _list_receptors(self) -> List[str]:
        try:
            candidates = set()
            rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
            if rec_dir.is_dir():
                for item in rec_dir.iterdir():
                    if item.is_dir() and not item.name.startswith("."):
                        candidates.add(item.name)
                    elif item.is_file() and item.suffix.lower() == ".pdbqt":
                        candidates.add(item.stem)

            proj_root = getattr(self.config, "project_root", None)
            if proj_root:
                macro_dir = (Path(proj_root) / "Macromolecules").resolve()
                if macro_dir.is_dir() and macro_dir != rec_dir:
                    for item in macro_dir.iterdir():
                        if item.is_dir() and not item.name.startswith("."):
                            candidates.add(item.name)
                        elif item.is_file() and item.suffix.lower() == ".pdbqt":
                            candidates.add(item.stem)

            real_candidates = [c for c in sorted(candidates) if c != "example_receptor"]
            if real_candidates:
                return real_candidates
            if candidates:
                return sorted(candidates)
        except Exception:
            pass
        return ["(no receptors found)"]

    def _open_receptor_folder(self) -> None:
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        name = self.var_receptor.get().strip()
        target = (rec_dir / name) if (name and not name.startswith("(")) else rec_dir
        target.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(target))
        else:
            subprocess.run(["xdg-open", str(target)])

    def _find_receptor_pdbqt(self, name: str) -> Optional[Path]:
        """Find receptor PDBQT across standard, rigid, alias, and Macromolecules locations."""
        if not name or name.startswith("("):
            return None

        aliases = [name]
        if "_clean" in name:
            aliases.append(name.replace("_clean", ""))
        else:
            aliases.append(f"{name}_clean")

        # Lowercase / uppercase variants
        for a in list(aliases):
            if a.lower() not in aliases:
                aliases.append(a.lower())
            if a.upper() not in aliases:
                aliases.append(a.upper())

        search_dirs: List[Path] = []
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        search_dirs.append(rec_dir)

        proj_root = getattr(self.config, "project_root", None)
        if proj_root:
            m_dir = (Path(proj_root) / "Macromolecules").resolve()
            if m_dir.is_dir() and m_dir not in search_dirs:
                search_dirs.append(m_dir)

        found_path: Optional[Path] = None

        # 1. Search subdirectories
        for s_dir in search_dirs:
            if not s_dir.is_dir():
                continue
            for alias in aliases:
                sub = s_dir / alias
                if sub.is_dir():
                    rigid_sub = sub / "rigid"
                    if rigid_sub.is_dir():
                        for a in aliases:
                            p = rigid_sub / f"{a}.pdbqt"
                            if p.is_file():
                                found_path = p
                                break
                        if not found_path:
                            all_p = list(rigid_sub.glob("*.pdbqt"))
                            if all_p:
                                found_path = all_p[0]

                    if not found_path:
                        for a in aliases:
                            p = sub / f"{a}.pdbqt"
                            if p.is_file():
                                found_path = p
                                break
                    if not found_path:
                        all_p = list(sub.glob("*.pdbqt"))
                        if all_p:
                            found_path = all_p[0]

                if found_path:
                    break
            if found_path:
                break

        # 2. Search loose .pdbqt in search_dirs
        if not found_path:
            for s_dir in search_dirs:
                if not s_dir.is_dir():
                    continue
                for alias in aliases:
                    p = s_dir / f"{alias}.pdbqt"
                    if p.is_file():
                        found_path = p
                        break
                if found_path:
                    break

        if not found_path or not found_path.is_file():
            return None

        # Ensure mirrored to standard location: rec_dir / name / [rigid/]
        target_sub = rec_dir / name
        target_sub.mkdir(parents=True, exist_ok=True)
        (target_sub / "rigid").mkdir(parents=True, exist_ok=True)

        std_pdbqt = target_sub / f"{name}.pdbqt"
        std_rigid_pdbqt = target_sub / "rigid" / f"{name}.pdbqt"

        if not std_pdbqt.is_file():
            try:
                shutil.copy2(found_path, std_pdbqt)
            except Exception:
                pass
        if not std_rigid_pdbqt.is_file():
            try:
                shutil.copy2(found_path, std_rigid_pdbqt)
            except Exception:
                pass

        return std_pdbqt if std_pdbqt.is_file() else found_path

    def _refresh_receptor_list(self) -> None:
        names = self._list_receptors()
        self.opt_receptor.configure(values=names)
        cur = self.var_receptor.get()
        if (not cur or cur.startswith("(") or cur not in names) and names:
            valid = [n for n in names if not n.startswith("(")]
            if valid:
                self.var_receptor.set(valid[0])
                self._on_receptor_change(valid[0])
            else:
                self.var_receptor.set(names[0])
        self._set_status("Receptor list refreshed.")

    def _on_receptor_change(self, name: str) -> None:
        self._receptor_name = name
        self._load_existing_config(name)

    def _load_existing_config(self, name: str) -> None:
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        config_path = rec_dir / name / "config.txt"
        if not config_path.is_file():
            aliases = [name]
            if "_clean" in name:
                aliases.append(name.replace("_clean", ""))
            else:
                aliases.append(f"{name}_clean")
            for a in aliases:
                alt = rec_dir / a / "config.txt"
                if alt.is_file():
                    config_path = alt
                    break
            proj_root = getattr(self.config, "project_root", None)
            if not config_path.is_file() and proj_root:
                m_dir = (Path(proj_root) / "Macromolecules").resolve()
                for a in aliases:
                    alt = m_dir / a / "config.txt"
                    if alt.is_file():
                        config_path = alt
                        break

        if not config_path.is_file():
            return
        params = {}
        for line in config_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, _, v = line.partition("=")
                params[k.strip()] = v.strip()
        mapping = [
            ("center_x", self.ent_cx), ("center_y", self.ent_cy),
            ("center_z", self.ent_cz), ("size_x", self.ent_sx),
            ("size_y", self.ent_sy), ("size_z", self.ent_sz),
        ]
        for key, ent in mapping:
            if key in params:
                ent.delete(0, "end")
                ent.insert(0, params[key])
        self._set_status(f"Loaded existing config.txt for {name}.")

    # ── Co-ligand ─────────────────────────────────────────────────────────────

    def _browse_coligand(self) -> None:
        path = filedialog.askopenfilename(
            title="Select co-crystallized ligand file",
            filetypes=[("Structure Files", "*.sdf *.mol2 *.pdb *.mol"),
                       ("All Files", "*.*")])
        if path:
            self._coligand_path = Path(path)
            self.lbl_colig.configure(text=self._coligand_path.name,
                                     text_color=themes.get().text_primary)

    def _autofill_from_coligand(self) -> None:
        colig = getattr(self, "_coligand_path", None)
        if not colig or not colig.is_file():
            messagebox.showwarning("No File", "Please browse a ligand file first.")
            return
        try:
            padding = float(self.ent_padding.get())
        except ValueError:
            padding = 5.0
        self._set_status("Computing bounding box ...")
        threading.Thread(target=self._calc_bbox,
                         args=(colig, padding), daemon=True).start()

    def _calc_bbox(self, colig: Path, padding: float) -> None:
        try:
            from prepare import get_coligand_bbox
            bbox = get_coligand_bbox(colig, padding)
        except Exception as e:
            self.after(0, self._set_status, f"Bbox error: {e}")
            return
        if not bbox:
            self.after(0, self._set_status, "Could not compute bounding box.")
            return
        self.after(0, self._apply_bbox, bbox)

    def _apply_bbox(self, bbox: dict) -> None:
        for key, ent in [
            ("center_x", self.ent_cx), ("center_y", self.ent_cy),
            ("center_z", self.ent_cz), ("size_x", self.ent_sx),
            ("size_y", self.ent_sy), ("size_z", self.ent_sz),
        ]:
            ent.delete(0, "end")
            ent.insert(0, str(bbox[key]))
        self._set_status(
            f"Auto-filled: center ({bbox['center_x']}, {bbox['center_y']}, "
            f"{bbox['center_z']})  size ({bbox['size_x']}, {bbox['size_y']}, {bbox['size_z']})")

    # ── Grid coordinate helpers ────────────────────────────────────────────────

    def _get_grid_params(self) -> Optional[dict]:
        try:
            return {
                "center_x": float(self.ent_cx.get()),
                "center_y": float(self.ent_cy.get()),
                "center_z": float(self.ent_cz.get()),
                "size_x": float(self.ent_sx.get()),
                "size_y": float(self.ent_sy.get()),
                "size_z": float(self.ent_sz.get()),
            }
        except ValueError as e:
            messagebox.showerror("Invalid Value", f"Grid parameter error: {e}")
            return None

    # ── Generate Vina config.txt ──────────────────────────────────────────────

    def _gen_vina_config(self) -> None:
        name = self.var_receptor.get()
        if not name or name.startswith("("):
            messagebox.showwarning("No Receptor", "Select a receptor first.")
            return
        params = self._get_grid_params()
        if not params:
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve() / name
        rec_dir.mkdir(parents=True, exist_ok=True)
        config_path = rec_dir / "config.txt"
        lines = [
            f"# Vina config — {name}",
            f"center_x = {params['center_x']:.3f}",
            f"center_y = {params['center_y']:.3f}",
            f"center_z = {params['center_z']:.3f}",
            f"size_x = {params['size_x']:.3f}",
            f"size_y = {params['size_y']:.3f}",
            f"size_z = {params['size_z']:.3f}",
            f"exhaustiveness = {self.config.exhaustiveness}",
            f"num_modes = {self.config.num_modes}",
            f"energy_range = {self.config.energy_range}",
        ]
        config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        # Mirror config.txt to alias folder and Macromolecules if present
        aliases = []
        if "_clean" in name:
            aliases.append(name.replace("_clean", ""))
        else:
            aliases.append(f"{name}_clean")
        for a in aliases:
            alt_dir = rec_dir.parent / a
            if alt_dir.is_dir():
                shutil.copy2(config_path, alt_dir / "config.txt")

        proj_root = getattr(self.config, "project_root", None)
        if proj_root:
            m_dir = (Path(proj_root) / "Macromolecules").resolve()
            if m_dir.is_dir() and m_dir != rec_dir.parent:
                for a in [name] + aliases:
                    m_sub = m_dir / a
                    if m_sub.is_dir():
                        shutil.copy2(config_path, m_sub / "config.txt")

        # Ensure receptor PDBQT exists in this folder as well
        self._find_receptor_pdbqt(name)

        self._set_status(f"Written: {config_path}")
        self.on_config_change()
        messagebox.showinfo("Config Written",
                            f"Vina config.txt saved:\n{config_path}")

    # ── Generate AutoDock4 GPF ────────────────────────────────────────────────

    def _gen_gpf(self) -> None:
        name = self.var_receptor.get()
        if not name or name.startswith("("):
            messagebox.showwarning("No Receptor", "Select a receptor first.")
            return
        params = self._get_grid_params()
        if not params:
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve() / name
        rec_dir.mkdir(parents=True, exist_ok=True)

        rigid_pdbqt = self._find_receptor_pdbqt(name)
        if not rigid_pdbqt or not rigid_pdbqt.is_file():
            messagebox.showerror(
                "Missing PDBQT",
                f"No PDBQT found for '{name}' in receptors/ or Macromolecules/.\n\n"
                "Please prepare the receptor in the Preparation tab first."
            )
            return

        # Ensure receptor PDBQT is directly in rec_dir so AutoGrid can find it
        dest_pdbqt = rec_dir / f"{name}.pdbqt"
        if not dest_pdbqt.is_file():
            shutil.copy2(rigid_pdbqt, dest_pdbqt)

        # Find all available ligand PDBQTs for atom types
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands")).resolve()
        lig_files = sorted(list(lig_dir.glob("*.pdbqt"))) if lig_dir.is_dir() else []
        if not lig_files:
            messagebox.showerror("No Ligands",
                                  "No ligand PDBQTs found in ligands/.\n"
                                  "Prepare at least one ligand first.")
            return

        try:
            spacing = float(self.ent_spacing.get())
        except ValueError:
            spacing = 0.375
        try:
            smooth = float(self.ent_smooth.get())
        except ValueError:
            smooth = 0.5
        try:
            dielectric = float(self.ent_dielec.get())
        except ValueError:
            dielectric = -0.1465

        gpf_path = rec_dir / f"{name}.gpf"
        self._set_status(f"Generating GPF for {len(lig_files)} ligand(s) ...")

        threading.Thread(
            target=self._run_gen_gpf,
            args=(dest_pdbqt, lig_files, gpf_path, params,
                  spacing, smooth, dielectric),
            daemon=True).start()

    def _run_gen_gpf(self, rigid_pdbqt, lig_files, gpf_path,
                      params, spacing, smooth, dielectric) -> None:
        try:
            from autodock4_workflow import generate_gpf
            generate_gpf(
                receptor_path=rigid_pdbqt,
                ligand_path=lig_files,
                output_path=gpf_path,
                config=self.config,
                grid_center={k.replace("center_", ""): params[k]
                              for k in ["center_x", "center_y", "center_z"]},
                grid_size={k.replace("size_", ""): params[k]
                           for k in ["size_x", "size_y", "size_z"]},
                spacing=spacing,
                smooth=smooth,
                dielectric=dielectric,
            )
            self.after(0, self._on_gpf_done, gpf_path, None)
        except Exception as e:
            self.after(0, self._on_gpf_done, gpf_path, e)

    def _on_gpf_done(self, gpf_path: Path, error: Optional[Exception]) -> None:
        if error:
            self._set_status(f"GPF generation failed: {error}")
            messagebox.showerror("GPF Failed", str(error))
            return
        self._set_status(f"GPF written: {gpf_path}")
        self.on_config_change()
        messagebox.showinfo("GPF Ready",
                            f"AutoGrid4 GPF written:\n{gpf_path}\n\n"
                            "This GPF covers all atom types in your ligand library.\n"
                            "Click 'Run AutoGrid4' to compute affinity maps.")

    # ── Generate AutoDock4 DPF ────────────────────────────────────────────────

    def _gen_dpf(self) -> None:
        name = self.var_receptor.get()
        if not name or name.startswith("("):
            messagebox.showwarning("No Receptor", "Select a receptor first.")
            return
        params = self._get_grid_params()
        if not params:
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve() / name
        rec_dir.mkdir(parents=True, exist_ok=True)

        rigid_pdbqt = self._find_receptor_pdbqt(name)
        if not rigid_pdbqt or not rigid_pdbqt.is_file():
            messagebox.showerror(
                "Missing PDBQT",
                f"No PDBQT found for '{name}' in receptors/ or Macromolecules/.\n\n"
                "Please prepare the receptor in the Preparation tab first."
            )
            return

        dest_pdbqt = rec_dir / f"{name}.pdbqt"
        if not dest_pdbqt.is_file():
            shutil.copy2(rigid_pdbqt, dest_pdbqt)

        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands")).resolve()
        lig_files = sorted(list(lig_dir.glob("*.pdbqt"))) if lig_dir.is_dir() else []
        if not lig_files:
            messagebox.showerror("No Ligands",
                                  "No ligand PDBQTs found in ligands/.\n"
                                  "Prepare at least one ligand first.")
            return

        self._set_status(f"Generating DPF for {len(lig_files)} ligand(s) ...")

        threading.Thread(
            target=self._run_gen_dpf,
            args=(dest_pdbqt, lig_files, rec_dir, params),
            daemon=True).start()

    def _run_gen_dpf(self, rigid_pdbqt, lig_files: List[Path], rec_dir: Path, params) -> None:
        try:
            from autodock4_workflow import generate_dpf
            center = {k.replace("center_", ""): params[k]
                      for k in ["center_x", "center_y", "center_z"]}

            generated_dpfs = []
            for lig in lig_files:
                dpf_path = rec_dir / f"{rigid_pdbqt.stem}_{lig.stem}.dpf"
                generate_dpf(
                    receptor_path=rigid_pdbqt,
                    ligand_path=lig,
                    output_path=dpf_path,
                    config=self.config,
                    grid_center=center,
                    algorithm=getattr(self.config, "ad4_algorithm", "LGA"),
                )
                generated_dpfs.append(dpf_path)

            self.after(0, self._on_dpf_done, generated_dpfs, None)
        except Exception as e:
            self.after(0, self._on_dpf_done, [], e)

    def _on_dpf_done(self, dpf_paths: List[Path], error: Optional[Exception]) -> None:
        if error:
            self._set_status(f"DPF generation failed: {error}")
            messagebox.showerror("DPF Failed", str(error))
            return
        n = len(dpf_paths)
        preview = "\n".join(f"  • {p.name}" for p in dpf_paths[:6])
        if n > 6:
            preview += f"\n  ... and {n - 6} more"
        self._set_status(f"Generated {n} DPF file(s)")
        self.on_config_change()
        messagebox.showinfo(
            "DPFs Ready",
            f"Successfully generated {n} AutoDock DPF file(s) for all ligands in ligands/:\n\n"
            f"{preview}\n\n"
            f"Saved to receptors/{self.var_receptor.get()}/ and results/AD4/DPF/.\n"
            "You can inspect them or go to Docking Console to run AutoDock4."
        )

    def _zoom_in_ag(self) -> None:
        if self._ag_font_size < 26:
            self._ag_font_size += 1
            self.tb_autogrid.configure(font=ctk.CTkFont(family="Consolas", size=self._ag_font_size))
            self.lbl_ag_font.configure(text=f"{self._ag_font_size} pt")

    def _zoom_out_ag(self) -> None:
        if self._ag_font_size > 9:
            self._ag_font_size -= 1
            self.tb_autogrid.configure(font=ctk.CTkFont(family="Consolas", size=self._ag_font_size))
            self.lbl_ag_font.configure(text=f"{self._ag_font_size} pt")

    # ── Run AutoGrid4 ─────────────────────────────────────────────────────────

    def _run_autogrid(self) -> None:
        name = self.var_receptor.get()
        if not name or name.startswith("("):
            messagebox.showwarning("No Receptor", "Select a receptor first.")
            return
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve() / name
        rec_dir.mkdir(parents=True, exist_ok=True)
        gpf_path = rec_dir / f"{name}.gpf"

        if not gpf_path.is_file():
            # Check if an alias GPF exists
            aliases = [name]
            if "_clean" in name:
                aliases.append(name.replace("_clean", ""))
            else:
                aliases.append(f"{name}_clean")
            for a in aliases:
                alt_gpf = (rec_dir.parent / a / f"{a}.gpf")
                if alt_gpf.is_file():
                    shutil.copy2(alt_gpf, gpf_path)
                    break

        if not gpf_path.is_file():
            messagebox.showerror("No GPF",
                                  f"No GPF found:\n{gpf_path}\n\n"
                                  "Generate the GPF first.")
            return

        ag_exe = self.config.autogrid4_executable
        if not ag_exe or not Path(ag_exe).is_file():
            from executables import resolve_executable
            ag_exe = resolve_executable(None, "autogrid4.exe")

        if not ag_exe or not Path(ag_exe).is_file():
            messagebox.showerror(
                "AutoGrid4 Not Found",
                f"AutoGrid4 executable was not found:\n{ag_exe}\n\n"
                "Please ensure bin/autogrid4.exe is present in the workspace or configure its path in Settings."
            )
            return

        # Ensure receptor file is present
        dest_pdbqt = self._find_receptor_pdbqt(name)

        # Check what receptor filename the GPF file explicitly references
        expected_rec_name = f"{name}.pdbqt"
        try:
            for line in gpf_path.read_text(encoding="utf-8", errors="replace").splitlines():
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0].lower() == "receptor":
                    expected_rec_name = parts[1]
                    break
        except Exception:
            pass

        target_expected = rec_dir / expected_rec_name
        if not target_expected.is_file() and dest_pdbqt and dest_pdbqt.is_file():
            shutil.copy2(dest_pdbqt, target_expected)

        if not target_expected.is_file():
            # Search for any candidate pdbqt in rec_dir or rec_dir / rigid
            candidates = []
            if (rec_dir / "rigid").is_dir():
                candidates += list((rec_dir / "rigid").glob("*.pdbqt"))
            candidates += [p for p in rec_dir.glob("*.pdbqt") if p.name != expected_rec_name]
            if candidates:
                shutil.copy2(candidates[0], target_expected)

        if not target_expected.is_file():
            messagebox.showerror(
                "Receptor PDBQT Not Found",
                f"AutoGrid requires the receptor file:\n{target_expected}\n\n"
                f"Please prepare the receptor PDBQT in the Preparation tab first."
            )
            return

        try:
            from prepare import ensure_receptor_has_charges
            ensure_receptor_has_charges(dest_pdbqt)
        except Exception:
            pass

        glg_path = rec_dir / f"{name}.glg"

        self._set_textbox(self.tb_autogrid, f"Starting AutoGrid4 ({ag_exe})...\n")
        self.lbl_ag_status.configure(text="Running...", text_color=themes.get().accent)
        threading.Thread(
            target=self._stream_autogrid,
            args=(ag_exe, gpf_path, glg_path, rec_dir),
            daemon=True).start()

    def _stream_autogrid(self, ag_exe, gpf_path, glg_path, cwd) -> None:
        try:
            proc = subprocess.Popen(
                [str(ag_exe), "-p", gpf_path.name, "-l", glg_path.name],
                cwd=str(cwd),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace")
            for line in proc.stdout:
                self.after(0, self._append_ag_log, line)
            proc.wait()
            rc = proc.returncode
            self.after(0, self._on_autogrid_done, cwd, rc)
        except Exception as e:
            self.after(0, self._append_ag_log, f"\nError: {e}\n")
            self.after(0, self.lbl_ag_status.configure,
                       {"text": "Error", "text_color": "red"})

    def _append_ag_log(self, text: str) -> None:
        self.tb_autogrid.configure(state="normal")
        self.tb_autogrid.insert("end", text)
        self.tb_autogrid.see("end")
        self.tb_autogrid.configure(state="disabled")

    def _on_autogrid_done(self, cwd: Path, rc: int) -> None:
        t = themes.get()
        if rc == 0:
            maps = list(cwd.glob("*.map"))
            self.lbl_ag_status.configure(text=f"Done — {len(maps)} map(s)",
                                          text_color=t.success if hasattr(t, "success") else "green")
            self._set_status(f"AutoGrid4 done. {len(maps)} affinity map(s) generated in {cwd}.")
            self.on_config_change()
        else:
            self.lbl_ag_status.configure(text=f"Failed (rc={rc})", text_color="red")
            self._set_status(f"AutoGrid4 failed (exit code {rc}).")

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _set_status(self, msg: str) -> None:
        self.lbl_status.configure(text=f"  {msg}")

    def _set_textbox(self, tb: ctk.CTkTextbox, text: str) -> None:
        tb.configure(state="normal")
        tb.delete("1.0", "end")
        tb.insert("end", text)
        tb.configure(state="disabled")

    def refresh(self) -> None:
        self._refresh_receptor_list()
