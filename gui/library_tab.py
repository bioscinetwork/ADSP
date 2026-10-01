"""
AutoDock Suite Pro — Screening Library Tab (gui/library_tab.py)
===============================================================
Displays all ligands in the library with full physicochemical
properties computed by prepare.calculate_molecule_descriptors().

Columns (all defensible in publication):
  Name, Formula, MW, LogP (Crippen), PSA (Å²), HBD, HBA,
  RotBonds, Arom.Rings, RO5, Veber

Row color coding:
  ✔ accent row   — passes BOTH RO5 and Veber
  ✔ normal row   — passes RO5 only
  — dim row      — fails both
"""

from __future__ import annotations

import threading
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List, Optional

import customtkinter as ctk
from tkinter import ttk

from gui import themes, widgets

if TYPE_CHECKING:
    from config import ProjectConfig


class LibraryTab(ctk.CTkFrame):
    """Screening Library Explorer with ADMET physicochemical properties."""

    COLS = (
        "Name", "Formula", "MW", "LogP", "PSA (Å²)",
        "HBD", "HBA", "RotBonds", "Arom.Rings", "RO5", "Veber",
    )
    COL_WIDTHS = (170, 110, 70, 70, 80, 50, 50, 75, 80, 60, 60)

    def __init__(self, parent: ctk.CTkTabview, config: "ProjectConfig",
                 on_config_change: Callable) -> None:
        t = themes.get()
        super().__init__(parent, fg_color=t.bg_primary)
        self.config = config
        self.on_config_change = on_config_change
        self._data: List[Dict] = []
        self._build()

    # ─────────────────────────────────────────────────
    # Build UI
    # ─────────────────────────────────────────────────

    def _build(self) -> None:
        t = themes.get()
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(0, weight=1)

        # Header row
        hdr = ctk.CTkFrame(self, fg_color="transparent")
        hdr.grid(row=0, column=0, sticky="ew", padx=12, pady=(10, 0))
        hdr.grid_columnconfigure(1, weight=1)

        widgets.section_header(hdr, "Screening Library", row=0, col=0, colspan=1)

        btn_frame = ctk.CTkFrame(hdr, fg_color="transparent")
        btn_frame.grid(row=0, column=1, sticky="e")

        widgets.accent_button(btn_frame, "⟳ Scan Library",
                              self._scan_library, width=130, small=True
                              ).grid(row=0, column=0, padx=4)

        # Filter bar
        filter_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        filter_frame.grid(row=1, column=0, sticky="ew", padx=12, pady=6)
        filter_frame.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(filter_frame, text="Filter:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=0, padx=8, pady=6)
        self.entry_filter = ctk.CTkEntry(
            filter_frame, placeholder_text="Name / Formula ...",
            fg_color=t.bg_tertiary, border_color=t.border,
            text_color=t.text_primary, width=260,
            font=themes.font_body(),
        )
        self.entry_filter.grid(row=0, column=1, sticky="w", padx=4, pady=6)
        self.entry_filter.bind("<KeyRelease>", lambda _: self._apply_filter())

        ctk.CTkLabel(filter_frame, text="Show:", text_color=t.text_secondary,
                     font=themes.font_body()).grid(row=0, column=2, padx=(16, 4))
        self.cmb_filter_ro5 = ctk.CTkOptionMenu(
            filter_frame, values=["All", "RO5 Pass", "RO5 Fail", "Veber Pass", "Both Pass"],
            width=130, command=lambda _: self._apply_filter(),
            fg_color=t.bg_tertiary, button_color=t.accent,
            text_color=t.text_primary,
            font=themes.font_body(),
            dropdown_font=themes.font_body(),
        )
        self.cmb_filter_ro5.grid(row=0, column=3, padx=4, pady=6)
        self.lbl_count = ctk.CTkLabel(filter_frame, text="0 compounds",
                                      text_color=t.text_secondary,
                                      font=themes.font_caption_bold())
        self.lbl_count.grid(row=0, column=4, padx=12)

        # Treeview
        tree_frame = ctk.CTkFrame(self, fg_color=t.bg_secondary, corner_radius=8)
        tree_frame.grid(row=2, column=0, sticky="nsew", padx=12, pady=(0, 8))
        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        self.tree = widgets.styled_treeview(tree_frame, self.COLS, height=18)
        for col, w in zip(self.COLS, self.COL_WIDTHS):
            self.tree.heading(col, text=col, command=lambda c=col: self._sort_by(c))
            anchor = "center" if col not in ("Name", "Formula") else "w"
            self.tree.column(col, width=w, anchor=anchor, minwidth=40)
        self.tree.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        vsb.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=vsb.set)

        # Legend
        legend = ctk.CTkFrame(self, fg_color="transparent")
        legend.grid(row=3, column=0, sticky="w", padx=12, pady=(0, 8))
        for text, ok in [("✔ Passes RO5 + Veber", True),
                         ("  Passes RO5 only", None),
                         ("✖ Fails RO5", False)]:
            widgets.status_badge(legend, text, ok).pack(side="left", padx=8)

        # Property info footer
        ctk.CTkLabel(
            self,
            text=("LogP: Crippen–Wildman (RDKit)  |  PSA: Ertl et al. 2000  |  "
                  "RO5: Lipinski et al. 1997  |  Veber: Veber et al. 2002"),
            text_color=t.text_disabled, font=themes.font_caption(), anchor="w",
        ).grid(row=4, column=0, sticky="w", padx=12, pady=(0, 6))

    # ─────────────────────────────────────────────────
    # Data
    # ─────────────────────────────────────────────────

    def _scan_library(self) -> None:
        """Scan ligands/ directory and compute descriptors in background thread."""
        from prepare import calculate_molecule_descriptors
        lig_dir = self.config.ligand_directory
        if not lig_dir.exists():
            return

        def _worker() -> None:
            pdbqts = list(lig_dir.glob("*.pdbqt"))
            data = []
            for p in pdbqts:
                info = calculate_molecule_descriptors(p)
                data.append(info)
            self._data = data
            self.after(0, self._populate_tree, data)

        threading.Thread(target=_worker, daemon=True).start()
        self.lbl_count.configure(text="Scanning …")

    def _populate_tree(self, data: List[Dict]) -> None:
        for item in self.tree.get_children():
            self.tree.delete(item)
        t = themes.get()
        for i, info in enumerate(data):
            veber_str = "✔" if info.get("veber_pass") else "✖"
            ro5_str = info.get("lipinski", "Pass")
            passes_both = (ro5_str == "Pass") and info.get("veber_pass", False)
            passes_ro5 = ro5_str == "Pass"

            if passes_both:
                tag = "both_pass"
            elif passes_ro5:
                tag = "ro5_pass"
            else:
                tag = "fails"

            row_tag = "even" if i % 2 == 0 else "odd"
            self.tree.insert("", "end", tags=(row_tag, tag), values=(
                info.get("stem", info.get("file_name", "")),
                info.get("formula", ""),
                f"{info.get('mw', 0):.2f}",
                f"{info.get('logp', 0):.2f}",
                f"{info.get('psa', 0):.1f}",
                info.get("hbd", 0),
                info.get("hba", 0),
                info.get("rot_bonds", 0),
                info.get("arom_rings", 0),
                ro5_str,
                veber_str,
            ))

        self.tree.tag_configure("both_pass", foreground=t.status_ok_text)
        self.tree.tag_configure("ro5_pass", foreground=t.text_primary)
        self.tree.tag_configure("fails", foreground=t.text_disabled)

        self.lbl_count.configure(text=f"{len(data)} compound(s)")

    def _apply_filter(self) -> None:
        query = self.entry_filter.get().strip().lower()
        mode = self.cmb_filter_ro5.get()

        filtered = self._data
        if query:
            filtered = [
                d for d in filtered
                if query in d.get("stem", "").lower()
                or query in d.get("formula", "").lower()
            ]
        if mode == "RO5 Pass":
            filtered = [d for d in filtered if d.get("lipinski", "") == "Pass"]
        elif mode == "RO5 Fail":
            filtered = [d for d in filtered if d.get("lipinski", "") != "Pass"]
        elif mode == "Veber Pass":
            filtered = [d for d in filtered if d.get("veber_pass")]
        elif mode == "Both Pass":
            filtered = [d for d in filtered
                        if d.get("lipinski", "") == "Pass" and d.get("veber_pass")]

        self._populate_tree(filtered)

    def _sort_by(self, col: str) -> None:
        """Sort treeview by column header click."""
        items = [(self.tree.set(child, col), child)
                 for child in self.tree.get_children("")]
        try:
            items.sort(key=lambda x: float(x[0]))
        except (ValueError, TypeError):
            items.sort(key=lambda x: x[0].lower())
        for i, (_, child) in enumerate(items):
            self.tree.move(child, "", i)
