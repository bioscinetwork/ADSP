"""
AutoDock Suite Pro -- Preparation Tab (gui/preparation_tab.py)
Full molecular preparation studio:
  1. Receptor Preparation (PDB/CIF -> PDBQT) with HETATM multi-removal checklist
  2. Flexible Receptor Splitting with searchable residue browser
  3. Ligand Batch Preparation (SDF/MOL2/PDB -> PDBQT)
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import TYPE_CHECKING, Callable, Dict, List, Optional, Set, Tuple

import customtkinter as ctk

from gui import themes, widgets

if TYPE_CHECKING:
    from config import ProjectConfig


# ── Scrollable Checklist ─────────────────────────────────────────────────────

class _ChecklistFrame(ctk.CTkScrollableFrame):
    """Scrollable frame of labeled checkboxes with live search filter."""

    def __init__(self, parent, **kwargs):
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_tertiary, **kwargs)
        self._vars: Dict[str, ctk.BooleanVar] = {}
        self._rows: Dict[str, ctk.CTkCheckBox] = {}
        self._data: list = []
        self._on_change_cb = None

    def set_change_callback(self, cb) -> None:
        self._on_change_cb = cb

    def _on_change(self) -> None:
        if self._on_change_cb:
            self._on_change_cb()

    def populate(self, records: list, label_key: str = "label",
                 pre_checked_pred=None) -> None:
        t = themes.get()
        for w in self.winfo_children():
            w.destroy()
        self._vars.clear()
        self._rows.clear()
        self._data = records
        for rec in records:
            label = str(rec.get(label_key, ""))
            var = ctk.BooleanVar(
                value=bool(pre_checked_pred and pre_checked_pred(rec))
            )
            cb = ctk.CTkCheckBox(
                self, text=label, variable=var,
                text_color=t.text_primary,
                font=ctk.CTkFont(family="Consolas", size=11),
                checkbox_width=16, checkbox_height=16,
                command=self._on_change,
            )
            cb.pack(anchor="w", padx=6, pady=1)
            self._vars[label] = var
            self._rows[label] = cb

    def filter(self, query: str) -> None:
        q = query.lower()
        for label, cb in self._rows.items():
            if q in label.lower():
                cb.pack(anchor="w", padx=6, pady=1)
            else:
                cb.pack_forget()

    def get_checked_records(self) -> list:
        return [
            rec for rec in self._data
            if self._vars.get(str(rec.get("label", "")),
                              ctk.BooleanVar()).get()
        ]

    def check_all(self) -> None:
        for v in self._vars.values():
            v.set(True)

    def uncheck_all(self) -> None:
        for v in self._vars.values():
            v.set(False)


# ── Main Tab ─────────────────────────────────────────────────────────────────

class PreparationTab(ctk.CTkFrame):
    """Full molecular preparation studio tab."""

    def __init__(self, parent, config: "ProjectConfig",
                 on_config_change: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self._pdb_path: Optional[Path] = None
        self._receptor_pdbqt: Optional[Path] = None
        self._hetatm_records: list = []
        self._chain_records: list = []
        self._pdb_residues: list = []
        self._mutations: Dict[Tuple[str, int], str] = {}
        self._deleted_residues: Set[Tuple[str, int]] = set()
        self._atom_residues: list = []
        self._lig_files: List[Path] = []
        self._build()

    # ── Build ────────────────────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        widgets.section_header(self, "Molecular Preparation Studio", row=0, colspan=2)

        # Left: receptor + flex
        left = ctk.CTkFrame(self, fg_color="transparent")
        left.grid(row=1, column=0, sticky="nsew", padx=(12, 6), pady=4)
        left.grid_columnconfigure(0, weight=1)
        left.grid_rowconfigure(0, weight=2)
        left.grid_rowconfigure(1, weight=1)
        self._build_receptor_section(left)
        self._build_flex_section(left)

        # Right: ligands
        right = ctk.CTkFrame(self, fg_color="transparent")
        right.grid(row=1, column=1, sticky="nsew", padx=(6, 12), pady=4)
        right.grid_columnconfigure(0, weight=1)
        right.grid_rowconfigure(0, weight=1)
        self._build_ligand_section(right)

        widgets.divider(self, row=2, colspan=2)
        self.lbl_status = ctk.CTkLabel(
            self, text="  Ready -- browse a PDB/CIF file to begin.",
            text_color=t.text_secondary, font=themes.font_caption(), anchor="w",
        )
        self.lbl_status.grid(row=3, column=0, columnspan=2, sticky="w", padx=12, pady=4)

    # ── Receptor section ─────────────────────────────────────────────────────

    def _build_receptor_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.grid(row=0, column=0, sticky="nsew", pady=(0, 8))
        card.grid_columnconfigure(0, weight=1)
        card.grid_rowconfigure(3, weight=1)

        ctk.CTkLabel(
            card, text="Receptor Preparation  (PDB / CIF -> PDBQT)",
            text_color=t.accent, font=themes.font_section(),
            anchor="w",
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=10, pady=(8, 4))

        # File picker
        fr = ctk.CTkFrame(card, fg_color="transparent")
        fr.grid(row=1, column=0, sticky="ew", padx=8, pady=4)
        fr.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(fr, text="Source:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=0, padx=(0, 6))
        self.lbl_pdb = ctk.CTkLabel(fr, text="(none selected)",
                                    text_color=t.text_secondary,
                                    font=themes.font_code(),
                                    anchor="w")
        self.lbl_pdb.grid(row=0, column=1, sticky="ew")
        widgets.accent_button(fr, "Browse PDB/CIF",
                              self._browse_pdb, width=130, small=True
                              ).grid(row=0, column=2, padx=(8, 0))

        # Basic cleanup options
        opt = ctk.CTkFrame(card, fg_color="transparent")
        opt.grid(row=2, column=0, sticky="w", padx=8, pady=(0, 4))
        self.var_water = ctk.BooleanVar(value=True)
        self.var_hydrogens = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(opt, text="Remove waters", variable=self.var_water,
                        text_color=t.text_primary, font=ctk.CTkFont(size=11)
                        ).pack(side="left", padx=(0, 12))
        ctk.CTkCheckBox(opt, text="Add polar hydrogens", variable=self.var_hydrogens,
                        text_color=t.text_primary, font=ctk.CTkFont(size=11)
                        ).pack(side="left")

        # Studio tabview: Tab 1 = Chains & Cofactors (Remove), Tab 2 = Residue Editor & Mutations
        self.rec_tabview = ctk.CTkTabview(
            card, height=195,
            fg_color=t.bg_tertiary,
            segmented_button_fg_color=t.bg_primary,
            segmented_button_selected_color=t.accent,
            segmented_button_unselected_color=t.bg_tertiary,
            text_color=t.text_primary,
        )
        self.rec_tabview.grid(row=3, column=0, sticky="nsew", padx=8, pady=(0, 4))

        tab_clean = self.rec_tabview.add("🧬 Chains & Cofactors")
        tab_mut = self.rec_tabview.add("✏ Residue Editor & Mutations")

        # ── Tab 1: Chains & Cofactors ──────────────────────
        tab_clean.grid_columnconfigure(0, weight=1)
        tab_clean.grid_columnconfigure(1, weight=1)
        tab_clean.grid_rowconfigure(1, weight=1)

        # Left: Chains checklist
        chain_hdr = ctk.CTkFrame(tab_clean, fg_color="transparent")
        chain_hdr.grid(row=0, column=0, sticky="ew", padx=(4, 6), pady=(2, 2))
        chain_hdr.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(
            chain_hdr, text="Chains to REMOVE:",
            text_color=t.accent, font=themes.font_body_bold(),
        ).grid(row=0, column=0, sticky="w")
        ch_bf = ctk.CTkFrame(chain_hdr, fg_color="transparent")
        ch_bf.grid(row=0, column=1, sticky="e")
        for txt, cmd in [("All", lambda: self.chains_list.check_all()),
                         ("None", lambda: self.chains_list.uncheck_all())]:
            ctk.CTkButton(ch_bf, text=txt, width=44, height=24,
                          fg_color="transparent", border_width=1,
                          border_color=t.border, text_color=t.text_secondary,
                          font=themes.font_caption_bold(), command=cmd
                          ).pack(side="left", padx=1)

        self.chains_list = _ChecklistFrame(tab_clean, height=110)
        self.chains_list.grid(row=1, column=0, sticky="nsew", padx=(4, 6), pady=(0, 4))
        self.chains_list.populate([{"label": "(load a PDB to see chains)"}])

        # Right: HETATM checklist
        het_hdr = ctk.CTkFrame(tab_clean, fg_color="transparent")
        het_hdr.grid(row=0, column=1, sticky="ew", padx=(6, 4), pady=(2, 2))
        het_hdr.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(
            het_hdr, text="Cofactors/Ligands to REMOVE:",
            text_color=t.accent, font=themes.font_body_bold(),
        ).grid(row=0, column=0, sticky="w")
        het_bf = ctk.CTkFrame(het_hdr, fg_color="transparent")
        het_bf.grid(row=0, column=1, sticky="e")
        for txt, cmd in [("All", lambda: self.hetatm_list.check_all()),
                         ("None", lambda: self.hetatm_list.uncheck_all())]:
            ctk.CTkButton(het_bf, text=txt, width=44, height=24,
                          fg_color="transparent", border_width=1,
                          border_color=t.border, text_color=t.text_secondary,
                          font=themes.font_caption_bold(), command=cmd
                          ).pack(side="left", padx=1)

        self.hetatm_list = _ChecklistFrame(tab_clean, height=110)
        self.hetatm_list.grid(row=1, column=1, sticky="nsew", padx=(6, 4), pady=(0, 4))
        self.hetatm_list.populate([{"label": "(load a PDB to see HETATM entries)"}])

        # ── Tab 2: Residue Editor & Mutations ─────────────
        tab_mut.grid_columnconfigure(0, weight=3)
        tab_mut.grid_columnconfigure(1, weight=2)
        tab_mut.grid_rowconfigure(1, weight=1)

        # Top controls
        mut_ctrl = ctk.CTkFrame(tab_mut, fg_color="transparent")
        mut_ctrl.grid(row=0, column=0, columnspan=2, sticky="ew", padx=4, pady=(2, 4))
        mut_ctrl.grid_columnconfigure(4, weight=1)

        ctk.CTkLabel(mut_ctrl, text="Chain:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=0, padx=(0, 4))
        self.cmb_res_chain = ctk.CTkOptionMenu(
            mut_ctrl, values=["All Chains"], width=95, height=28,
            fg_color=t.bg_secondary, button_color=t.accent,
            font=themes.font_body(), dropdown_font=themes.font_body(),
            command=lambda _: self._filter_pdb_residues(),
        )
        self.cmb_res_chain.grid(row=0, column=1, padx=(0, 8))

        ctk.CTkLabel(mut_ctrl, text="Search:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=2, padx=(0, 4))
        self.ent_res_search = ctk.CTkEntry(
            mut_ctrl, placeholder_text="seq / name...", width=110, height=28,
            font=themes.font_body(),
        )
        self.ent_res_search.grid(row=0, column=3, padx=(0, 8))
        self.ent_res_search.bind("<KeyRelease>", lambda _: self._filter_pdb_residues())

        ctk.CTkLabel(mut_ctrl, text="Target:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=4, sticky="e", padx=(0, 4))
        self.cmb_mutate_to = ctk.CTkOptionMenu(
            mut_ctrl,
            values=[
                "ALA (Ala-scan)", "GLY (Gly-scan)", "VAL", "LEU", "ILE",
                "SER", "THR", "CYS", "MET", "ASP", "GLU", "ASN", "GLN",
                "LYS", "ARG", "HIS", "PHE", "TYR", "TRP", "PRO",
                "HID (neutral δ)", "HIE (neutral ε)", "HIP (protonated)",
                "CYX (disulfide)", "ASH (neutral Asp)", "GLH (neutral Glu)",
            ],
            width=125, height=28,
            fg_color=t.bg_secondary, button_color=t.accent,
            font=themes.font_body(), dropdown_font=themes.font_body(),
        )
        self.cmb_mutate_to.grid(row=0, column=5, padx=(0, 4))

        ctk.CTkButton(
            mut_ctrl, text="Mutate", width=65, height=28,
            fg_color=t.accent, hover_color=t.accent_hover,
            font=themes.font_button(),
            command=self._apply_mutation,
        ).grid(row=0, column=6, padx=2)

        ctk.CTkButton(
            mut_ctrl, text="Delete", width=60, height=28,
            fg_color=t.bg_secondary, text_color=t.status_fail_text,
            hover_color=t.border, font=themes.font_button(),
            command=self._apply_deletion,
        ).grid(row=0, column=7, padx=2)

        ctk.CTkButton(
            mut_ctrl, text="Revert", width=60, height=28,
            fg_color=t.bg_secondary, text_color=t.text_secondary,
            hover_color=t.border, font=themes.font_button(),
            command=self._revert_selected_residue,
        ).grid(row=0, column=8, padx=2)

        # Left: Residues treeview
        res_tree_frame = ctk.CTkFrame(tab_mut, fg_color="transparent")
        res_tree_frame.grid(row=1, column=0, sticky="nsew", padx=(4, 4), pady=(0, 4))
        res_tree_frame.grid_columnconfigure(0, weight=1)
        res_tree_frame.grid_rowconfigure(0, weight=1)

        self.tree_pdb_res = widgets.styled_treeview(
            res_tree_frame, ("Chain", "Seq", "Name", "Atoms", "Status"), height=5,
        )
        for col, w in zip(("Chain", "Seq", "Name", "Atoms", "Status"), (45, 45, 55, 45, 100)):
            self.tree_pdb_res.heading(col, text=col)
            self.tree_pdb_res.column(col, width=w, anchor="center")
        self.tree_pdb_res.grid(row=0, column=0, sticky="nsew")

        vsb_res = ttk.Scrollbar(res_tree_frame, orient="vertical", command=self.tree_pdb_res.yview)
        vsb_res.grid(row=0, column=1, sticky="ns")
        self.tree_pdb_res.configure(yscrollcommand=vsb_res.set)

        # Right: Pending modifications list
        mod_frame = ctk.CTkFrame(tab_mut, fg_color=t.bg_secondary, corner_radius=6)
        mod_frame.grid(row=1, column=1, sticky="nsew", padx=(4, 4), pady=(0, 4))
        mod_frame.grid_columnconfigure(0, weight=1)
        mod_frame.grid_rowconfigure(1, weight=1)

        mod_hdr = ctk.CTkFrame(mod_frame, fg_color="transparent")
        mod_hdr.grid(row=0, column=0, sticky="ew", padx=6, pady=(4, 2))
        mod_hdr.grid_columnconfigure(0, weight=1)
        self.lbl_mod_hdr = ctk.CTkLabel(
            mod_hdr, text="Pending Edits (0):",
            text_color=t.text_secondary, font=themes.font_body_bold(),
        )
        self.lbl_mod_hdr.grid(row=0, column=0, sticky="w")
        ctk.CTkButton(
            mod_hdr, text="Clear All", width=65, height=22,
            fg_color="transparent", text_color=t.text_secondary,
            font=themes.font_caption_bold(), command=self._clear_all_edits,
        ).grid(row=0, column=1, sticky="e")

        self.tb_pending_edits = ctk.CTkTextbox(
            mod_frame, state="disabled", fg_color=t.bg_primary,
            text_color=t.text_primary, font=themes.font_code(),
        )
        self.tb_pending_edits.grid(row=1, column=0, sticky="nsew", padx=6, pady=(0, 6))

        # Bottom row of card: Name + Edits Badge + Prepare Receptor button
        nr = ctk.CTkFrame(card, fg_color="transparent")
        nr.grid(row=4, column=0, sticky="ew", padx=8, pady=(4, 8))
        nr.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(nr, text="Name:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=0, padx=(0, 6))
        self.ent_rec_name = ctk.CTkEntry(
            nr, placeholder_text="receptor_name", font=themes.font_body(), height=32)
        self.ent_rec_name.grid(row=0, column=1, sticky="ew", padx=(0, 8))

        self.lbl_edits_badge = ctk.CTkLabel(
            nr, text="0 chains removed | 0 mutations",
            text_color=t.text_secondary, font=themes.font_caption(),
        )
        self.lbl_edits_badge.grid(row=0, column=2, padx=(0, 6))

        widgets.accent_button(nr, "Prepare Receptor",
                              self._prepare_receptor, width=140
                              ).grid(row=0, column=3, padx=(0, 4))

        ctk.CTkButton(nr, text="📁 Open Folder",
                      command=self._open_receptor_folder, width=95, height=28,
                      fg_color=t.bg_secondary, text_color=t.text_primary,
                      hover_color=t.border, font=ctk.CTkFont(size=11)
                      ).grid(row=0, column=4)

    # ── Flex section ─────────────────────────────────────────────────────────

    def _build_flex_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.grid(row=1, column=0, sticky="nsew")
        card.grid_columnconfigure(0, weight=1)
        card.grid_rowconfigure(2, weight=1)

        ctk.CTkLabel(
            card, text="Flexible Receptor Splitting",
            text_color=t.accent, font=ctk.CTkFont(size=12, weight="bold"),
            anchor="w",
        ).grid(row=0, column=0, columnspan=3, sticky="w", padx=10, pady=(8, 4))

        # Search filter
        sr = ctk.CTkFrame(card, fg_color="transparent")
        sr.grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 4))
        sr.grid_columnconfigure(1, weight=1)
        ctk.CTkLabel(sr, text="Filter:", text_color=t.text_secondary,
                     font=ctk.CTkFont(size=11)).grid(row=0, column=0, padx=(0, 6))
        self.ent_flex_filter = ctk.CTkEntry(
            sr, placeholder_text="residue name or number...",
            font=ctk.CTkFont(size=11))
        self.ent_flex_filter.grid(row=0, column=1, sticky="ew")
        self.ent_flex_filter.bind(
            "<KeyRelease>",
            lambda e: self.flex_list.filter(self.ent_flex_filter.get()))

        # Residue checklist
        self.flex_list = _ChecklistFrame(card, height=120)
        self.flex_list.grid(row=2, column=0, sticky="nsew", padx=8, pady=(0, 4))
        self.flex_list.populate([{"label": "(prepare a receptor first)"}])
        self.flex_list.set_change_callback(self._update_flex_summary)

        # Summary + split button
        splitr = ctk.CTkFrame(card, fg_color="transparent")
        splitr.grid(row=3, column=0, sticky="ew", padx=8, pady=(2, 8))
        splitr.grid_columnconfigure(0, weight=1)
        self.lbl_flex_sel = ctk.CTkLabel(
            splitr, text="Selected: (none)",
            text_color=t.text_secondary,
            font=themes.font_code(), anchor="w")
        self.lbl_flex_sel.grid(row=0, column=0, sticky="w")
        widgets.accent_button(splitr, "Split Flex/Rigid",
                              self._split_flex, width=150
                              ).grid(row=0, column=1, padx=(8, 0))

    # ── Ligand section ────────────────────────────────────────────────────────

    def _build_ligand_section(self, parent) -> None:
        t = themes.get()
        card = ctk.CTkFrame(parent, fg_color=t.bg_secondary, corner_radius=8)
        card.grid(row=0, column=0, sticky="nsew")
        card.grid_columnconfigure(0, weight=1)
        card.grid_rowconfigure(3, weight=1)

        ctk.CTkLabel(
            card, text="Ligand Preparation  (SDF / MOL2 / PDB -> PDBQT)",
            text_color=t.accent, font=themes.font_section(),
            anchor="w",
        ).grid(row=0, column=0, columnspan=2, sticky="w", padx=10, pady=(8, 4))

        br = ctk.CTkFrame(card, fg_color="transparent")
        br.grid(row=1, column=0, sticky="ew", padx=8, pady=(0, 4))
        widgets.accent_button(br, "Add Files",
                              self._add_lig_files, width=110, small=True
                              ).grid(row=0, column=0, padx=(0, 4))
        ctk.CTkButton(
            br, text="Clear List", command=self._clear_lig_list,
            width=100, height=28, fg_color="transparent",
            border_width=1, border_color=t.border,
            text_color=t.text_secondary, font=themes.font_button(),
        ).grid(row=0, column=1, padx=4)
        self.var_gen3d = ctk.BooleanVar(value=True)
        ctk.CTkCheckBox(br, text="Generate 3D conformer", variable=self.var_gen3d,
                        text_color=t.text_primary, font=themes.font_body()
                        ).grid(row=0, column=2, padx=(12, 0))

        self.lbl_lig_count = ctk.CTkLabel(
            card, text="0 files selected",
            text_color=t.text_secondary, font=themes.font_caption(), anchor="w")
        self.lbl_lig_count.grid(row=2, column=0, sticky="w", padx=10)

        self.tb_lig_status = ctk.CTkTextbox(
            card, state="disabled", height=260,
            fg_color=t.bg_tertiary, text_color=t.text_primary,
            font=themes.font_code())
        self.tb_lig_status.grid(row=3, column=0, sticky="nsew", padx=8, pady=4)

        self.prog_lig = ctk.CTkProgressBar(
            card, height=8, fg_color=t.bg_tertiary, progress_color=t.accent)
        self.prog_lig.set(0)
        self.prog_lig.grid(row=4, column=0, sticky="ew", padx=8, pady=(0, 4))

        widgets.accent_button(card, "Prepare All Ligands",
                              self._prepare_ligands, width=180
                              ).grid(row=5, column=0, pady=(2, 10))

    # ── Receptor callbacks ────────────────────────────────────────────────────

    def _browse_pdb(self) -> None:
        path = filedialog.askopenfilename(
            title="Select receptor PDB/CIF file",
            filetypes=[("Structure Files", "*.pdb *.cif *.ent"),
                       ("PDB", "*.pdb"), ("CIF", "*.cif"),
                       ("All Files", "*.*")])
        if not path:
            return
        self._pdb_path = Path(path)
        self.lbl_pdb.configure(text=self._pdb_path.name,
                               text_color=themes.get().text_primary)
        
        # Always update receptor name with newly selected file stem
        self.ent_rec_name.delete(0, "end")
        self.ent_rec_name.insert(0, self._pdb_path.stem)

        # Reset previous receptor preparation state
        self._receptor_pdbqt = None
        self._mutations.clear()
        self._deleted_residues.clear()
        self._pdb_residues = []
        self._chain_records = []
        self._hetatm_records = []
        self._filter_pdb_residues()
        self._update_edits_display()
        self.flex_list.populate([{"label": "(prepare a receptor first)"}])

        self._set_status(f"Parsing {self._pdb_path.name} ...")
        threading.Thread(target=self._load_pdb_metadata, daemon=True).start()

    def _load_pdb_metadata(self) -> None:
        try:
            from prepare import parse_hetatm_records, parse_pdb_chains, parse_pdb_residues
            hetatms = parse_hetatm_records(self._pdb_path)
            chains = parse_pdb_chains(self._pdb_path)
            residues = parse_pdb_residues(self._pdb_path)
        except Exception as e:
            self.after(0, self._set_status, f"Error parsing PDB: {e}")
            return
        self.after(0, self._populate_pdb_metadata, chains, hetatms, residues)

    def _populate_pdb_metadata(self, chains: list, hetatms: list, residues: list) -> None:
        self._chain_records = chains
        self._hetatm_records = hetatms
        self._pdb_residues = residues
        self._mutations.clear()
        self._deleted_residues.clear()

        # Populate chains
        if not chains:
            self.chains_list.populate([{"label": "(no chain records found)"}])
        else:
            self.chains_list.populate(chains, label_key="label")

        # Populate hetatms
        if not hetatms:
            self.hetatm_list.populate([{"label": "(no HETATM records found)"}])
        else:
            self.hetatm_list.populate(
                hetatms, label_key="label",
                pre_checked_pred=lambda r: r.get("is_water", False))

        # Populate chain dropdown in residue editor
        chain_ids = sorted(list(set(r["chain"] for r in residues)))
        self.cmb_res_chain.configure(values=["All Chains"] + chain_ids)
        self.cmb_res_chain.set("All Chains")

        self._filter_pdb_residues()
        self._update_edits_display()

        n_lig = sum(1 for r in hetatms if not r.get("is_water", False))
        n_wat = sum(1 for r in hetatms if r.get("is_water", False))
        self._set_status(
            f"Loaded {self._pdb_path.name}: {len(chains)} chain(s), {len(residues)} residues, "
            f"{n_lig} cofactors, {n_wat} waters.")

    # ── Residue Editor & Mutation Methods ─────────────────────────────────────

    def _filter_pdb_residues(self) -> None:
        for item in self.tree_pdb_res.get_children():
            self.tree_pdb_res.delete(item)

        sel_chain = self.cmb_res_chain.get().strip()
        q = self.ent_res_search.get().strip().upper()

        idx = 0
        for r in self._pdb_residues:
            c = r["chain"]
            s = r["seq"]
            name = r["name"]
            atoms = r.get("atoms", 0)

            if sel_chain != "All Chains" and c != sel_chain:
                continue
            if q and (q not in str(s) and q not in name):
                continue

            key = (c, s)
            if key in self._mutations:
                status = f"MUTATE -> {self._mutations[key]}"
            elif key in self._deleted_residues:
                status = "DELETE"
            else:
                status = "—"

            tag = "even" if idx % 2 == 0 else "odd"
            self.tree_pdb_res.insert(
                "", "end", iid=f"{c}_{s}", tags=(tag,),
                values=(c, s, name, atoms, status),
            )
            idx += 1

    def _apply_mutation(self) -> None:
        sel = self.tree_pdb_res.selection()
        if not sel:
            messagebox.showwarning("Select Residue", "Please select one or more residues from the table.")
            return

        target = self.cmb_mutate_to.get().split()[0].upper()
        for item in sel:
            vals = self.tree_pdb_res.item(item, "values")
            chain = str(vals[0])
            seq = int(vals[1])
            key = (chain, seq)
            self._mutations[key] = target
            self._deleted_residues.discard(key)

        self._filter_pdb_residues()
        self._update_edits_display()

    def _apply_deletion(self) -> None:
        sel = self.tree_pdb_res.selection()
        if not sel:
            messagebox.showwarning("Select Residue", "Please select one or more residues from the table.")
            return

        for item in sel:
            vals = self.tree_pdb_res.item(item, "values")
            chain = str(vals[0])
            seq = int(vals[1])
            key = (chain, seq)
            self._deleted_residues.add(key)
            self._mutations.pop(key, None)

        self._filter_pdb_residues()
        self._update_edits_display()

    def _revert_selected_residue(self) -> None:
        sel = self.tree_pdb_res.selection()
        if not sel:
            return

        for item in sel:
            vals = self.tree_pdb_res.item(item, "values")
            chain = str(vals[0])
            seq = int(vals[1])
            key = (chain, seq)
            self._mutations.pop(key, None)
            self._deleted_residues.discard(key)

        self._filter_pdb_residues()
        self._update_edits_display()

    def _clear_all_edits(self) -> None:
        self._mutations.clear()
        self._deleted_residues.clear()
        self._filter_pdb_residues()
        self._update_edits_display()

    def _update_edits_display(self) -> None:
        self.tb_pending_edits.configure(state="normal")
        self.tb_pending_edits.delete("1.0", "end")

        lines = []
        for (c, s), target in sorted(self._mutations.items()):
            lines.append(f"• [MUTATE] Chain {c}:{s} -> {target}")
        for (c, s) in sorted(self._deleted_residues):
            lines.append(f"• [DELETE] Chain {c}:{s}")

        if not lines:
            self.tb_pending_edits.insert("end", "(no modifications queued)\n")
        else:
            self.tb_pending_edits.insert("end", "\n".join(lines) + "\n")
        self.tb_pending_edits.configure(state="disabled")

        n_mut = len(self._mutations)
        n_del = len(self._deleted_residues)
        self.lbl_mod_hdr.configure(text=f"Pending Edits ({n_mut + n_del}):")

        checked_chains = [r["chain"] for r in self.chains_list.get_checked_records() if "chain" in r]
        self.lbl_edits_badge.configure(
            text=f"{len(checked_chains)} chains removed | {n_mut} mutations | {n_del} deleted"
        )

    # ── Preparation Execution ─────────────────────────────────────────────────

    def _prepare_receptor(self) -> None:
        if not self._pdb_path or not self._pdb_path.is_file():
            messagebox.showwarning("No File", "Please browse a PDB/CIF file first.")
            return
        name = self.ent_rec_name.get().strip()
        if not name:
            messagebox.showwarning("No Name", "Please enter a receptor name.")
            return

        checked_chains = self.chains_list.get_checked_records()
        remove_chains = {r["chain"] for r in checked_chains if "chain" in r}

        checked_het = self.hetatm_list.get_checked_records()
        remove_keys = {r["key"] for r in checked_het if "key" in r}

        mutations = self._mutations.copy()
        deleted_residues = self._deleted_residues.copy()

        out_dir = self.config.receptor_directory / name / "rigid"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / f"{name}.pdbqt"

        n_mut = len(mutations)
        n_del = len(deleted_residues)
        self._set_status(f"Preparing {name} (Chains removed: {len(remove_chains)}, Mutations: {n_mut}, Deleted: {n_del}) ...")

        threading.Thread(
            target=self._run_receptor_prep,
            args=(name, out_pdbqt, remove_keys, remove_chains, mutations, deleted_residues),
            daemon=True).start()

    def _run_receptor_prep(self, name: str, out_pdbqt: Path,
                           remove_keys: set, remove_chains: set,
                           mutations: dict, deleted_residues: set) -> None:
        try:
            from prepare import prepare_receptor_pdbqt
            prepare_receptor_pdbqt(
                self._pdb_path, out_pdbqt,
                cleanup_water=self.var_water.get(),
                add_hydrogens=self.var_hydrogens.get(),
                remove_residue_keys=remove_keys,
                remove_chains=remove_chains,
                mutations=mutations,
                deleted_residues=deleted_residues,
            )
            self._receptor_pdbqt = out_pdbqt
            self.after(0, self._on_receptor_done, name, out_pdbqt, remove_chains, mutations, deleted_residues, None)
        except Exception as e:
            self.after(0, self._on_receptor_done, name, out_pdbqt, remove_chains, mutations, deleted_residues, e)

    def _open_receptor_folder(self) -> None:
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        name = self.ent_rec_name.get().strip() if hasattr(self, "ent_rec_name") else ""
        target = (rec_dir / name) if name else rec_dir
        target.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(target))
        else:
            subprocess.run(["xdg-open", str(target)])

    def refresh(self) -> None:
        """Called when project config or workspace changes."""
        pass

    def _on_receptor_done(self, name: str, out_pdbqt: Path,
                           remove_chains: set, mutations: dict,
                           deleted_residues: set, error: Optional[Exception]) -> None:
        if error:
            self._set_status(f"Failed: {error}")
            messagebox.showerror("Preparation Failed", str(error))
            return

        rec_dir = Path(self.config.receptor_directory).resolve()
        rec_dir.mkdir(parents=True, exist_ok=True)
        rec_folder = rec_dir / name
        rec_folder.mkdir(parents=True, exist_ok=True)

        config_path = rec_folder / "config.txt"
        if not config_path.is_file():
            config_path.write_text(
                "# EDIT THIS FILE: set the correct grid center and size.\n"
                "center_x = 0.000\ncenter_y = 0.000\ncenter_z = 0.000\n"
                "size_x = 25.000\nsize_y = 25.000\nsize_z = 25.000\n",
                encoding="utf-8",
            )

        # Ensure receptor PDBQT exists in both rigid/ subfolder and receptor root folder
        rec_root_pdbqt = rec_folder / f"{name}.pdbqt"
        if out_pdbqt.is_file():
            import shutil
            shutil.copy2(out_pdbqt, rec_root_pdbqt)

        # Ensure receptor files are cleanly populated in rec_folder
        (rec_folder / "rigid").mkdir(parents=True, exist_ok=True)
        if out_pdbqt.is_file():
            dest_root = rec_folder / f"{name}.pdbqt"
            if dest_root.resolve() != out_pdbqt.resolve():
                shutil.copy2(out_pdbqt, dest_root)
            dest_rigid = rec_folder / "rigid" / f"{name}.pdbqt"
            if dest_rigid.resolve() != out_pdbqt.resolve():
                shutil.copy2(out_pdbqt, dest_rigid)
        rec_cfg = rec_folder / "config.txt"
        if not rec_cfg.is_file() and config_path.is_file() and rec_cfg.resolve() != config_path.resolve():
            shutil.copy2(config_path, rec_cfg)

        # Mirror to Macromolecules folder only if it exists in the workspace
        if getattr(self.config, "project_root", None):
            macro_dir = (Path(self.config.project_root) / "Macromolecules").resolve()
            if macro_dir.is_dir() and macro_dir != rec_dir:
                m_dest = macro_dir / name
                m_dest.mkdir(parents=True, exist_ok=True)
                (m_dest / "rigid").mkdir(parents=True, exist_ok=True)
                if out_pdbqt.is_file():
                    mr = m_dest / f"{name}.pdbqt"
                    if mr.resolve() != out_pdbqt.resolve():
                        shutil.copy2(out_pdbqt, mr)
                    mrg = m_dest / "rigid" / f"{name}.pdbqt"
                    if mrg.resolve() != out_pdbqt.resolve():
                        shutil.copy2(out_pdbqt, mrg)
                m_cfg = m_dest / "config.txt"
                if config_path.is_file() and not m_cfg.is_file() and m_cfg.resolve() != config_path.resolve():
                    shutil.copy2(config_path, m_cfg)

        edited_pdb = rec_folder / f"{name}_edited.pdb"
        edited_msg = f"\n• Edited PDB: {edited_pdb.name}" if edited_pdb.is_file() else ""
        chain_msg = f"\n• Chains removed: {', '.join(sorted(remove_chains))}" if remove_chains else ""
        mut_msg = f"\n• Mutations applied: {len(mutations)}" if mutations else ""
        del_msg = f"\n• Residues deleted: {len(deleted_residues)}" if deleted_residues else ""

        self._set_status(f"Prepared: {out_pdbqt}")
        threading.Thread(target=self._load_atom_residues, daemon=True).start()
        self.on_config_change()
        messagebox.showinfo(
            "Receptor Ready",
            f"Saved receptor PDBQT:\n{out_pdbqt}\n\n"
            f"Folder:\n{rec_folder}\n"
            f"{edited_msg}{chain_msg}{mut_msg}{del_msg}\n\n"
            "Grid config.txt generated.\n"
            "Proceed to the Grid Box tab to configure coordinates.")

    # ── Flex callbacks ────────────────────────────────────────────────────────

    def _load_atom_residues(self) -> None:
        if not self._receptor_pdbqt or not self._receptor_pdbqt.is_file():
            return
        try:
            from prepare import inspect_receptor_residues
            residues = inspect_receptor_residues(self._receptor_pdbqt)
        except Exception:
            residues = []
        self.after(0, self._populate_flex_list, residues)

    def _populate_flex_list(self, residues: list) -> None:
        if not residues:
            self.flex_list.populate([{"label": "(no ATOM records)"}])
            return
        for r in residues:
            r["label"] = f"{r['chain']}  {r['name']:<4s}  {r['seq']}"
        self.flex_list.populate(residues, label_key="label")
        self._atom_residues = residues
        self._set_status(f"Flex browser: {len(residues)} residues loaded.")

    def _update_flex_summary(self) -> None:
        checked = self.flex_list.get_checked_records()
        if not checked:
            self.lbl_flex_sel.configure(text="Selected: (none)")
        else:
            labels = ", ".join(
                f"{r.get('chain', '')}:{r.get('name', '')}{r.get('seq', '')}"
                for r in checked[:6])
            if len(checked) > 6:
                labels += f" ... +{len(checked)-6} more"
            self.lbl_flex_sel.configure(text=f"Selected: {labels}")

    def _split_flex(self) -> None:
        if not self._receptor_pdbqt or not self._receptor_pdbqt.is_file():
            messagebox.showwarning("No Receptor", "Prepare a receptor first.")
            return
        checked = self.flex_list.get_checked_records()
        if not checked:
            messagebox.showwarning("No Selection",
                                   "Tick at least one residue for flexible docking.")
            return
        flex_ids = [f"{r.get('chain', 'A')}:{r.get('seq', 0)}" for r in checked]
        name = self._receptor_pdbqt.parent.parent.name
        flex_dir = self.config.receptor_directory / name / "flex"
        flex_dir.mkdir(parents=True, exist_ok=True)
        out_rigid = flex_dir / f"{name}_rigid.pdbqt"
        out_flex = flex_dir / f"{name}_flex.pdbqt"
        self._set_status(f"Splitting {len(checked)} residue(s) ...")
        threading.Thread(
            target=self._run_flex_split,
            args=(flex_ids, out_rigid, out_flex), daemon=True).start()

    def _run_flex_split(self, flex_ids: List[str],
                         out_rigid: Path, out_flex: Path) -> None:
        try:
            from prepare import split_flexible_receptor
            split_flexible_receptor(self._receptor_pdbqt, flex_ids, out_rigid, out_flex)
            self.after(0, self._on_flex_done, out_rigid, out_flex, None)
        except Exception as e:
            self.after(0, self._on_flex_done, out_rigid, out_flex, e)

    def _on_flex_done(self, out_rigid: Path, out_flex: Path,
                       error: Optional[Exception]) -> None:
        if error:
            self._set_status(f"Split failed: {error}")
            messagebox.showerror("Split Failed", str(error))
            return
        self._set_status(f"Split done: {out_rigid.name} + {out_flex.name}")
        self.on_config_change()
        messagebox.showinfo(
            "Flex Split Complete",
            f"Rigid backbone:  {out_rigid}\n"
            f"Flex sidechains: {out_flex}\n\n"
            "The Workspace tab will now show [FLEX] for this receptor.")

    # ── Ligand callbacks ──────────────────────────────────────────────────────

    def _add_lig_files(self) -> None:
        files = filedialog.askopenfilenames(
            title="Select ligand structure files",
            filetypes=[("Structure Files", "*.sdf *.mol2 *.pdb *.mol *.smi"),
                       ("SDF", "*.sdf"), ("MOL2", "*.mol2"),
                       ("PDB", "*.pdb"), ("SMILES", "*.smi"),
                       ("All Files", "*.*")])
        if not files:
            return
        for f in files:
            p = Path(f)
            if p not in self._lig_files:
                self._lig_files.append(p)
        self._refresh_lig_label()

    def _clear_lig_list(self) -> None:
        self._lig_files.clear()
        self._refresh_lig_label()
        self._set_textbox(self.tb_lig_status, "")
        self.prog_lig.set(0)

    def _refresh_lig_label(self) -> None:
        n = len(self._lig_files)
        self.lbl_lig_count.configure(text=f"{n} file(s) selected")
        preview = "\n".join(f"  {p.name}" for p in self._lig_files[:20])
        if n > 20:
            preview += f"\n  ... and {n-20} more"
        self._set_textbox(self.tb_lig_status, preview or "")

    def _prepare_ligands(self) -> None:
        if not self._lig_files:
            messagebox.showwarning("No Files", "Add ligand files first.")
            return
        lig_dir = self.config.ligand_directory
        lig_dir.mkdir(parents=True, exist_ok=True)
        self.prog_lig.set(0)
        self._set_textbox(self.tb_lig_status, "Starting ligand preparation...\n")
        threading.Thread(
            target=self._run_ligand_prep,
            args=(list(self._lig_files), lig_dir), daemon=True).start()

    def _run_ligand_prep(self, files: List[Path], lig_dir: Path) -> None:
        from prepare import prepare_ligand_pdbqt
        total = len(files)
        ok = fail = 0
        for i, src in enumerate(files):
            out = lig_dir / (src.stem + ".pdbqt")
            try:
                prepare_ligand_pdbqt(src, out, generate_3d=self.var_gen3d.get())
                ok += 1
                msg = f"  OK  {src.name} -> {out.name}\n"
            except Exception as e:
                fail += 1
                logger.error(f"Failed to prepare ligand {src.name}: {e}", exc_info=True)
                msg = f"  ERR {src.name}: {e}\n"
            self.after(0, self._append_lig_log, msg, (i + 1) / total)
        summary = f"\n  Done: {ok} prepared, {fail} failed.\n"
        self.after(0, self._append_lig_log, summary, 1.0)
        self.after(0, self.on_config_change)

    def _append_lig_log(self, msg: str, progress: float) -> None:
        self.tb_lig_status.configure(state="normal")
        self.tb_lig_status.insert("end", msg)
        self.tb_lig_status.see("end")
        self.tb_lig_status.configure(state="disabled")
        self.prog_lig.set(progress)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _set_status(self, msg: str) -> None:
        self.lbl_status.configure(text=f"  {msg}")

    def _set_textbox(self, tb: ctk.CTkTextbox, text: str) -> None:
        tb.configure(state="normal")
        tb.delete("1.0", "end")
        tb.insert("end", text)
        tb.configure(state="disabled")

    def refresh(self) -> None:
        pass
