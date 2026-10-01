"""
AutoDock Suite Pro — PySide6 Molecular Preparation Studio (gui_qt/preparation_tab.py)
=====================================================================================
Comprehensive molecular preparation studio:
  1. Receptor Preparation (PDB/CIF -> PDBQT):
     - Chains to remove checklist
     - Non-standard residues / cofactors / ligands (HETATM) removal checklist
     - Full interactive Residue Editor (Alanine & Glycine scanning, mutations, deletions)
     - Polar hydrogens, Gasteiger charges, water cleanup
  2. Flexible Receptor Splitting (Rigid + Flex PDBQT)
  3. Ligand Batch Preparation (SDF/MOL2/PDB/SMI -> PDBQT) with 3D conformer generation
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional, Set, Tuple

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig


class _PrepWorkerSignals(QObject):
    log_signal = Signal(str)
    done_signal = Signal(str, object)  # message, error
    flex_log_signal = Signal(str)
    flex_done_signal = Signal(str, object)
    lig_log_signal = Signal(str)
    lig_done_signal = Signal(str, object)


class PreparationTab(QWidget):
    """Full molecular preparation studio in PySide6."""

    MUTATION_TARGETS = [
        "ALA (Ala-scan)", "GLY (Gly-scan)", "VAL", "LEU", "ILE",
        "SER", "THR", "CYS", "MET", "ASP", "GLU", "ASN", "GLN",
        "LYS", "ARG", "HIS", "PHE", "TYR", "TRP", "PRO",
        "HID (neutral δ)", "HIE (neutral ε)", "HIP (protonated)",
        "CYX (disulfide)", "ASH (neutral Asp)", "GLH (neutral Glu)",
    ]

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change

        self._pdb_path: Optional[Path] = None
        self._parsed_chains: List[Dict[str, Any]] = []
        self._parsed_hetatms: List[Dict[str, Any]] = []
        self._parsed_residues: List[Dict[str, Any]] = []

        self._mutations: Dict[Tuple[str, int], str] = {}
        self._deleted_residues: Set[Tuple[str, int]] = set()

        self._signals = _PrepWorkerSignals()
        self._signals.log_signal.connect(self._on_worker_log)
        self._signals.done_signal.connect(self._on_worker_done)
        self._signals.flex_log_signal.connect(self._on_flex_log)
        self._signals.flex_done_signal.connect(self._on_flex_done)
        self._signals.lig_log_signal.connect(self._on_lig_log)
        self._signals.lig_done_signal.connect(self._on_lig_done)

        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        # ── Scientific Workflow Progression Stepper ────────────
        self.prep_stepper = self._build_stepper_bar()
        main_layout.addWidget(self.prep_stepper)
        self._update_stepper_stage(0)

        # Main sub-tab container
        self.studio_tabs = QTabWidget()

        # ── 1. Receptor Studio Tab ──────────────────────────────────
        tab_rec = QWidget()
        rec_layout = QVBoxLayout(tab_rec)
        rec_layout.setContentsMargins(10, 10, 10, 10)
        rec_layout.setSpacing(8)

        # File selection bar
        file_card = QFrame()
        file_card.setObjectName("cardFrame")
        fc_l = QHBoxLayout(file_card)
        fc_l.setContentsMargins(10, 8, 10, 8)
        fc_l.setSpacing(8)

        fc_l.addWidget(QLabel("Source Receptor:"))
        self.ent_rec_source = QLineEdit()
        self.ent_rec_source.setPlaceholderText("Select raw receptor file (.pdb, .cif, .ent, .mol2)...")
        self.ent_rec_source.textChanged.connect(self._on_rec_path_changed)
        fc_l.addWidget(self.ent_rec_source, 1)

        self.btn_browse_rec = QPushButton("📂 Browse")
        self.btn_browse_rec.setObjectName("accentButton")
        self.btn_browse_rec.clicked.connect(self._browse_source_receptor)
        fc_l.addWidget(self.btn_browse_rec)

        fc_l.addWidget(QLabel("Receptor Name:"))
        self.ent_rec_name = QLineEdit()
        self.ent_rec_name.setPlaceholderText("receptor_name")
        self.ent_rec_name.setMinimumWidth(160)
        fc_l.addWidget(self.ent_rec_name)

        self.btn_open_rec_folder = QPushButton("📁 Folder")
        self.btn_open_rec_folder.clicked.connect(self._open_receptor_folder)
        fc_l.addWidget(self.btn_open_rec_folder)

        rec_layout.addWidget(file_card)

        # Basic prep options bar
        opt_card = QFrame()
        opt_card.setObjectName("cardFrame")
        opt_l = QHBoxLayout(opt_card)
        opt_l.setContentsMargins(10, 6, 10, 6)
        opt_l.setSpacing(16)

        self.chk_add_h = QCheckBox("Add Polar Hydrogens (pH 7.4)")
        self.chk_add_h.setChecked(True)
        opt_l.addWidget(self.chk_add_h)

        self.chk_remove_water = QCheckBox("Remove Waters (HOH / TIP3 / WAT)")
        self.chk_remove_water.setChecked(True)
        opt_l.addWidget(self.chk_remove_water)

        opt_l.addWidget(QLabel("Charges:"))
        self.cmb_charge_model = QComboBox()
        self.cmb_charge_model.addItems([
            "Kollman United-Atom (Standard AD4 / Vina)",
            "Gasteiger (PEOE Empirical)",
        ])
        self.cmb_charge_model.setToolTip(
            "Kollman United-Atom charges are standard AMBER-derived partial charges for protein receptors.\n"
            "Gasteiger charges compute empirical partial charges via Partial Equalization of Orbital Electronegativity."
        )
        opt_l.addWidget(self.cmb_charge_model)

        opt_l.addStretch()
        rec_layout.addWidget(opt_card)

        # Studio inner tabs: Tab A (Chains & HETATMs), Tab B (Residue Editor & Alanine Scanning)
        self.rec_inner_tabs = QTabWidget()

        # ── Tab A: Chains & Cofactors / Ligands ───────────────────
        tab_chains_het = QWidget()
        ch_layout = QHBoxLayout(tab_chains_het)
        ch_layout.setContentsMargins(8, 8, 8, 8)
        ch_layout.setSpacing(12)

        # Left Column: Chains Checklist
        box_chains = QGroupBox("Chains to REMOVE (Unchecked = Retain)")
        bc_l = QVBoxLayout(box_chains)
        bc_l.setContentsMargins(8, 8, 8, 8)

        ch_btns = QHBoxLayout()
        btn_chain_all = QPushButton("Select All")
        btn_chain_all.clicked.connect(self._check_all_chains)
        ch_btns.addWidget(btn_chain_all)
        btn_chain_none = QPushButton("Deselect All")
        btn_chain_none.clicked.connect(self._uncheck_all_chains)
        ch_btns.addWidget(btn_chain_none)
        ch_btns.addStretch()
        bc_l.addLayout(ch_btns)

        self.list_chains = QListWidget()
        self.list_chains.setAlternatingRowColors(True)
        bc_l.addWidget(self.list_chains)
        ch_layout.addWidget(box_chains, 1)

        # Right Column: HETATM Checklist
        box_het = QGroupBox("Non-Standard Residues / Cofactors / Ligands (HETATM to REMOVE)")
        bh_l = QVBoxLayout(box_het)
        bh_l.setContentsMargins(8, 8, 8, 8)

        het_btns = QHBoxLayout()
        btn_het_all = QPushButton("Select All")
        btn_het_all.clicked.connect(self._check_all_hetatms)
        het_btns.addWidget(btn_het_all)
        btn_het_none = QPushButton("Deselect All")
        btn_het_none.clicked.connect(self._uncheck_all_hetatms)
        het_btns.addWidget(btn_het_none)
        het_btns.addStretch()
        bh_l.addLayout(het_btns)

        self.list_hetatms = QListWidget()
        self.list_hetatms.setAlternatingRowColors(True)
        bh_l.addWidget(self.list_hetatms)
        ch_layout.addWidget(box_het, 1)

        self.rec_inner_tabs.addTab(tab_chains_het, "🧬 Chains & Cofactors (Remove)")

        # ── Tab B: Residue Editor & Alanine Scanning ──────────────
        tab_res_editor = QWidget()
        re_layout = QVBoxLayout(tab_res_editor)
        re_layout.setContentsMargins(8, 8, 8, 8)
        re_layout.setSpacing(6)

        # Top mutation control bar
        mut_bar = QHBoxLayout()
        mut_bar.setSpacing(8)

        mut_bar.addWidget(QLabel("Chain:"))
        self.cmb_res_chain = QComboBox()
        self.cmb_res_chain.addItem("All Chains")
        self.cmb_res_chain.currentIndexChanged.connect(self._filter_residues_table)
        mut_bar.addWidget(self.cmb_res_chain)

        mut_bar.addWidget(QLabel("Search:"))
        self.ent_res_search = QLineEdit()
        self.ent_res_search.setPlaceholderText("seq or residue name (e.g. 119, TRP)...")
        self.ent_res_search.textChanged.connect(self._filter_residues_table)
        mut_bar.addWidget(self.ent_res_search, 1)

        mut_bar.addWidget(QLabel("Target:"))
        self.cmb_mutate_to = QComboBox()
        self.cmb_mutate_to.addItems(self.MUTATION_TARGETS)
        mut_bar.addWidget(self.cmb_mutate_to)

        self.btn_apply_mutate = QPushButton("⚡ Mutate")
        self.btn_apply_mutate.setObjectName("accentButton")
        self.btn_apply_mutate.clicked.connect(self._apply_mutation)
        mut_bar.addWidget(self.btn_apply_mutate)

        self.btn_apply_delete = QPushButton("🗑 Delete")
        self.btn_apply_delete.clicked.connect(self._apply_deletion)
        mut_bar.addWidget(self.btn_apply_delete)

        self.btn_revert_res = QPushButton("↩ Revert")
        self.btn_revert_res.clicked.connect(self._revert_selected_residue)
        mut_bar.addWidget(self.btn_revert_res)

        re_layout.addLayout(mut_bar)

        # Residues Table + Pending Modifications Splitter
        res_splitter = QSplitter(Qt.Horizontal)

        # Table of residues
        self.table_residues = QTableWidget()
        self.table_residues.setColumnCount(5)
        self.table_residues.setHorizontalHeaderLabels(["Chain", "Seq", "Name", "Atoms", "Status"])
        self.table_residues.setAlternatingRowColors(True)
        self.table_residues.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_residues.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table_residues.horizontalHeader().setStretchLastSection(True)
        res_splitter.addWidget(self.table_residues)

        # Pending Edits Box
        pending_card = QFrame()
        pending_card.setObjectName("cardFrame")
        pc_l = QVBoxLayout(pending_card)
        pc_l.setContentsMargins(8, 8, 8, 8)
        pc_l.setSpacing(6)

        pc_hdr = QHBoxLayout()
        self.lbl_pending_hdr = QLabel("Pending Edits (0):")
        self.lbl_pending_hdr.setStyleSheet("font-weight: bold; color: #38bdf8;")
        pc_hdr.addWidget(self.lbl_pending_hdr)
        pc_hdr.addStretch()

        self.btn_clear_edits = QPushButton("Clear All")
        self.btn_clear_edits.clicked.connect(self._clear_all_edits)
        pc_hdr.addWidget(self.btn_clear_edits)
        pc_l.addLayout(pc_hdr)

        self.txt_pending_edits = QTextEdit()
        self.txt_pending_edits.setReadOnly(True)
        self.txt_pending_edits.setPlaceholderText("(No modifications queued)")
        pc_l.addWidget(self.txt_pending_edits)

        res_splitter.addWidget(pending_card)
        res_splitter.setStretchFactor(0, 3)
        res_splitter.setStretchFactor(1, 2)

        re_layout.addWidget(res_splitter)

        self.rec_inner_tabs.addTab(tab_res_editor, "✏ Residue Editor & Alanine Scanning")

        rec_layout.addWidget(self.rec_inner_tabs, 1)

        # Action execution footer for Receptor
        act_card = QFrame()
        act_card.setObjectName("cardFrame")
        ac_l = QHBoxLayout(act_card)
        ac_l.setContentsMargins(10, 8, 10, 8)
        ac_l.setSpacing(12)

        self.lbl_edits_badge = QLabel("0 chains removed | 0 mutations | 0 deleted")
        self.lbl_edits_badge.setStyleSheet("color: #94a3b8; font-weight: bold;")
        ac_l.addWidget(self.lbl_edits_badge)
        ac_l.addStretch()

        self.btn_run_rec_prep = QPushButton("⚡ Prepare Receptor (PDBQT)")
        self.btn_run_rec_prep.setObjectName("accentButton")
        self.btn_run_rec_prep.setMinimumHeight(34)
        self.btn_run_rec_prep.setMinimumWidth(200)
        self.btn_run_rec_prep.clicked.connect(self._run_receptor_prep)
        ac_l.addWidget(self.btn_run_rec_prep)

        rec_layout.addWidget(act_card)

        # Receptor Console Log
        self.rec_log = QTextEdit()
        self.rec_log.setObjectName("consoleBox")
        self.rec_log.setReadOnly(True)
        self.rec_log.setMaximumHeight(140)
        self.rec_log.setPlaceholderText("Receptor preparation output log will appear here...")
        rec_layout.addWidget(self.rec_log)

        self.studio_tabs.addTab(tab_rec, "Receptor Prep (Rigid)")

        # ── 2. Flexible Receptor Splitting Tab ──────────────────────
        tab_flex = QWidget()
        flex_layout = QVBoxLayout(tab_flex)
        flex_layout.setContentsMargins(10, 10, 10, 10)
        flex_layout.setSpacing(10)

        flex_bar = QHBoxLayout()
        flex_bar.addWidget(QLabel("Prepared Receptor PDBQT:"))
        self.ent_flex_source = QLineEdit()
        self.ent_flex_source.setPlaceholderText("Select prepared receptor .pdbqt...")
        flex_bar.addWidget(self.ent_flex_source, 1)
        self.btn_browse_flex = QPushButton("📂 Browse")
        self.btn_browse_flex.clicked.connect(self._browse_flex_receptor)
        flex_bar.addWidget(self.btn_browse_flex)
        flex_layout.addLayout(flex_bar)

        flex_body = QHBoxLayout()
        res_box = QGroupBox("Select Flexible Residues (Sidechains)")
        res_layout = QVBoxLayout(res_box)
        self.ent_search_res = QLineEdit()
        self.ent_search_res.setPlaceholderText("Filter residues (e.g. TRP, PHE, TYR, 119)...")
        self.ent_search_res.textChanged.connect(self._filter_flex_residues)
        res_layout.addWidget(self.ent_search_res)

        self.list_flex_residues = QListWidget()
        res_layout.addWidget(self.list_flex_residues)
        flex_body.addWidget(res_box, 1)

        action_box = QGroupBox("Flexible Docking Output")
        act_layout = QVBoxLayout(action_box)
        self.lbl_flex_summary = QLabel("0 residue(s) selected for sidechain flexibility.")
        self.lbl_flex_summary.setStyleSheet("color: #38bdf8; font-weight: bold;")
        act_layout.addWidget(self.lbl_flex_summary)

        self.btn_split_flex = QPushButton("✂ Split into Rigid + Flex PDBQT")
        self.btn_split_flex.setObjectName("accentButton")
        self.btn_split_flex.setMinimumHeight(34)
        self.btn_split_flex.clicked.connect(self._run_flex_split)
        act_layout.addWidget(self.btn_split_flex)

        self.flex_log = QTextEdit()
        self.flex_log.setObjectName("consoleBox")
        self.flex_log.setReadOnly(True)
        act_layout.addWidget(self.flex_log)

        flex_body.addWidget(action_box, 1)
        flex_layout.addLayout(flex_body)

        self.studio_tabs.addTab(tab_flex, "Flexible Receptor Splitting")

        # ── 3. Ligand Batch Preparation Tab ─────────────────────────
        tab_lig = QWidget()
        lig_layout = QVBoxLayout(tab_lig)
        lig_layout.setContentsMargins(10, 10, 10, 10)
        lig_layout.setSpacing(10)

        lig_file_bar = QHBoxLayout()
        lig_file_bar.addWidget(QLabel("Source Ligand File(s):"))
        self.ent_lig_source = QLineEdit()
        self.ent_lig_source.setPlaceholderText("Select SDF, MOL2, SMILES, or PDB file(s)...")
        lig_file_bar.addWidget(self.ent_lig_source, 1)
        self.btn_browse_lig = QPushButton("📂 Browse Files")
        self.btn_browse_lig.setObjectName("accentButton")
        self.btn_browse_lig.clicked.connect(self._browse_source_ligands)
        lig_file_bar.addWidget(self.btn_browse_lig)
        lig_layout.addLayout(lig_file_bar)

        lig_opts = QGroupBox("Ligand Preparation Parameters")
        opts_l = QHBoxLayout(lig_opts)

        self.chk_gen3d = QCheckBox("Generate 3D Conformers (if 2D or SMILES)")
        self.chk_gen3d.setChecked(True)
        opts_l.addWidget(self.chk_gen3d)

        self.chk_assign_lig_charges = QCheckBox("Assign Gasteiger Charges & Rotatable Bonds")
        self.chk_assign_lig_charges.setChecked(True)
        opts_l.addWidget(self.chk_assign_lig_charges)
        opts_l.addStretch()
        lig_layout.addWidget(lig_opts)

        self.btn_run_lig_prep = QPushButton("🚀 Run Ligand Batch Preparation")
        self.btn_run_lig_prep.setObjectName("accentButton")
        self.btn_run_lig_prep.setMinimumHeight(34)
        self.btn_run_lig_prep.clicked.connect(self._run_ligand_prep)
        lig_layout.addWidget(self.btn_run_lig_prep)

        self.lig_log = QTextEdit()
        self.lig_log.setObjectName("consoleBox")
        self.lig_log.setReadOnly(True)
        lig_layout.addWidget(self.lig_log)

        self.studio_tabs.addTab(tab_lig, "Ligand Batch Preparation")
        main_layout.addWidget(self.studio_tabs)

    def _build_stepper_bar(self) -> QFrame:
        frame = QFrame()
        frame.setObjectName("cardFrame")
        fl = QHBoxLayout(frame)
        fl.setContentsMargins(12, 8, 12, 8)
        fl.setSpacing(6)

        lbl_flow = QLabel("WORKFLOW PROGRESSION:")
        lbl_flow.setStyleSheet("font-size: 10px; font-weight: 800; color: #64748b; letter-spacing: 0.5px; margin-right: 6px;")
        fl.addWidget(lbl_flow)

        self.stepper_labels = []
        stages = [
            "1. INPUT",
            "2. INSPECTION",
            "3. CLEANING",
            "4. CHARGES / H",
            "5. FLEXIBLE RESIDUES",
            "6. VALIDATION",
            "7. READY",
        ]
        for i, s in enumerate(stages):
            lbl = QLabel(s)
            lbl.setStyleSheet(
                "background-color: #1e293b; color: #94a3b8; font-size: 10px; font-weight: 700; "
                "border: 1px solid #334155; border-radius: 4px; padding: 4px 8px;"
            )
            fl.addWidget(lbl)
            self.stepper_labels.append(lbl)
            if i < len(stages) - 1:
                arrow = QLabel("→")
                arrow.setStyleSheet("color: #64748b; font-weight: bold;")
                fl.addWidget(arrow)

        fl.addStretch()
        return frame

    def _update_stepper_stage(self, current_stage: int) -> None:
        """Update stepper badges: 0=input, 1=inspection, 2=cleaning, 3=charges, 4=flex, 5=validation, 6=ready."""
        if not hasattr(self, "stepper_labels"):
            return
        stages_base = [
            "1. INPUT",
            "2. INSPECTION",
            "3. CLEANING",
            "4. CHARGES / H",
            "5. FLEXIBLE RESIDUES",
            "6. VALIDATION",
            "7. READY",
        ]
        for i, (lbl, base) in enumerate(zip(self.stepper_labels, stages_base)):
            if i < current_stage:
                lbl.setText(f"✓ {base}")
                lbl.setStyleSheet(
                    "background-color: rgba(16, 185, 129, 0.2); color: #10b981; font-size: 10px; font-weight: 700; "
                    "border: 1px solid rgba(16, 185, 129, 0.4); border-radius: 4px; padding: 4px 8px;"
                )
            elif i == current_stage:
                lbl.setText(base)
                lbl.setStyleSheet(
                    "background-color: rgba(2, 132, 199, 0.2); color: #0284c7; font-size: 10px; font-weight: 700; "
                    "border: 1px solid rgba(2, 132, 199, 0.4); border-radius: 4px; padding: 4px 8px;"
                )
            else:
                lbl.setText(base)
                lbl.setStyleSheet(
                    "background-color: #1e293b; color: #94a3b8; font-size: 10px; font-weight: 700; "
                    "border: 1px solid #334155; border-radius: 4px; padding: 4px 8px;"
                )

    # ─────────────────────────────────────────────────────────
    # Receptor Studio Scanning & Actions
    # ─────────────────────────────────────────────────────────

    def _browse_source_receptor(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self, "Select Receptor File",
            str(getattr(self.config, "receptor_directory", Path.cwd())),
            "Structural Files (*.pdb *.cif *.ent *.mol2);;All Files (*.*)"
        )
        if f:
            self.ent_rec_name.setText(Path(f).stem)
            self.ent_rec_source.setText(f)

    def _on_rec_path_changed(self, text: str) -> None:
        p = Path(text.strip())
        if p.is_file():
            self._pdb_path = p
            self.ent_rec_name.setText(p.stem)
            self._scan_receptor_details(p)

    def _scan_receptor_details(self, path: Path) -> None:
        self.list_chains.clear()
        self.list_hetatms.clear()
        self._parsed_chains.clear()
        self._parsed_hetatms.clear()
        self._parsed_residues.clear()
        self._mutations.clear()
        self._deleted_residues.clear()

        from prepare import parse_hetatm_records, parse_pdb_chains, parse_pdb_residues

        # 1. Parse chains
        try:
            self._parsed_chains = parse_pdb_chains(path)
            for ch in self._parsed_chains:
                label = ch.get("label", f"Chain {ch.get('chain')}")
                item = QListWidgetItem(label)
                item.setData(Qt.UserRole, ch.get("chain"))
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Unchecked)  # Retain by default
                self.list_chains.addItem(item)
            self.rec_log.append(f"✔ Parsed {len(self._parsed_chains)} chain(s) from {path.name}")
        except Exception as e:
            self.rec_log.append(f"⚠ Chain parsing warning: {e}")

        # 2. Parse HETATMs
        try:
            self._parsed_hetatms = parse_hetatm_records(path)
            for r in self._parsed_hetatms:
                name = r.get("resname", "UNK")
                seq = r.get("resseq", "")
                chain = r.get("chain", "")
                label = r.get("label", f"[{chain}:{seq}] {name}")
                item = QListWidgetItem(label)
                item.setData(Qt.UserRole, r.get("key"))
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                is_water = r.get("is_water", False)
                # Preserve non-water components by default.  Removal of a
                # cofactor, ligand, or metal-containing group must be an
                # explicit user decision; the separate water option controls
                # routine crystallographic solvent cleanup.
                item.setCheckState(Qt.Unchecked)
                self.list_hetatms.addItem(item)
            self.rec_log.append(f"✔ Found {len(self._parsed_hetatms)} HETATM groups in {path.name}")
        except Exception as e:
            self.rec_log.append(f"⚠ HETATM parsing warning: {e}")

        # 3. Parse Residues
        try:
            self._parsed_residues = parse_pdb_residues(path)
            chain_ids = sorted(list(set(r["chain"] for r in self._parsed_residues)))
            self.cmb_res_chain.clear()
            self.cmb_res_chain.addItem("All Chains")
            for cid in chain_ids:
                self.cmb_res_chain.addItem(cid)
            self._filter_residues_table()
            self.rec_log.append(f"✔ Loaded {len(self._parsed_residues)} amino acid residues for mutation/scanning.")
            self._update_stepper_stage(2)
        except Exception as e:
            self.rec_log.append(f"⚠ Residue parsing warning: {e}")

        self._update_edits_display()

    def _check_all_chains(self) -> None:
        for i in range(self.list_chains.count()):
            self.list_chains.item(i).setCheckState(Qt.Checked)
        self._update_edits_display()

    def _uncheck_all_chains(self) -> None:
        for i in range(self.list_chains.count()):
            self.list_chains.item(i).setCheckState(Qt.Unchecked)
        self._update_edits_display()

    def _check_all_hetatms(self) -> None:
        for i in range(self.list_hetatms.count()):
            self.list_hetatms.item(i).setCheckState(Qt.Checked)

    def _uncheck_all_hetatms(self) -> None:
        for i in range(self.list_hetatms.count()):
            self.list_hetatms.item(i).setCheckState(Qt.Unchecked)

    # ─────────────────────────────────────────────────────────
    # Residue Table & Mutation Actions
    # ─────────────────────────────────────────────────────────

    def _filter_residues_table(self) -> None:
        sel_chain = self.cmb_res_chain.currentText().strip()
        q = self.ent_res_search.text().strip().upper()

        matching = []
        for r in self._parsed_residues:
            c = r["chain"]
            s = r["seq"]
            name = r["name"]
            if sel_chain != "All Chains" and c != sel_chain:
                continue
            if q and (q not in str(s) and q not in name):
                continue
            matching.append(r)

        self.table_residues.setRowCount(len(matching))
        for row, r in enumerate(matching):
            c = r["chain"]
            s = r["seq"]
            name = r["name"]
            atoms = r.get("atoms", 0)

            key = (c, s)
            if key in self._mutations:
                status = f"MUTATE -> {self._mutations[key]}"
                status_color = QColor("#38bdf8")
            elif key in self._deleted_residues:
                status = "DELETE"
                status_color = QColor("#f87171")
            else:
                status = "—"
                status_color = QColor("#94a3b8")

            vals = [c, str(s), name, str(atoms), status]
            for col, val in enumerate(vals):
                item = QTableWidgetItem(val)
                item.setTextAlignment(Qt.AlignCenter)
                if col == 4:
                    item.setForeground(status_color)
                self.table_residues.setItem(row, col, item)

    def _get_selected_residue_keys(self) -> List[Tuple[str, int]]:
        sel_rows = sorted(set(idx.row() for idx in self.table_residues.selectedIndexes()))
        keys = []
        for r in sel_rows:
            chain_item = self.table_residues.item(r, 0)
            seq_item = self.table_residues.item(r, 1)
            if chain_item and seq_item:
                keys.append((chain_item.text(), int(seq_item.text())))
        return keys

    def _apply_mutation(self) -> None:
        keys = self._get_selected_residue_keys()
        if not keys:
            QMessageBox.information(self, "Selection Required", "Select one or more residues from the table.")
            return
        target = self.cmb_mutate_to.currentText().split()[0].upper()
        for k in keys:
            self._mutations[k] = target
            self._deleted_residues.discard(k)
        self._filter_residues_table()
        self._update_edits_display()

    def _apply_deletion(self) -> None:
        keys = self._get_selected_residue_keys()
        if not keys:
            QMessageBox.information(self, "Selection Required", "Select one or more residues from the table.")
            return
        for k in keys:
            self._deleted_residues.add(k)
            self._mutations.pop(k, None)
        self._filter_residues_table()
        self._update_edits_display()

    def _revert_selected_residue(self) -> None:
        keys = self._get_selected_residue_keys()
        if not keys:
            return
        for k in keys:
            self._mutations.pop(k, None)
            self._deleted_residues.discard(k)
        self._filter_residues_table()
        self._update_edits_display()

    def _clear_all_edits(self) -> None:
        self._mutations.clear()
        self._deleted_residues.clear()
        self._filter_residues_table()
        self._update_edits_display()

    def _update_edits_display(self) -> None:
        lines = []
        for (c, s), target in sorted(self._mutations.items()):
            lines.append(f"• [MUTATE] Chain {c}:{s} -> {target}")
        for (c, s) in sorted(self._deleted_residues):
            lines.append(f"• [DELETE] Chain {c}:{s}")

        if not lines:
            self.txt_pending_edits.setPlainText("(No modifications queued)")
        else:
            self.txt_pending_edits.setPlainText("\n".join(lines))

        n_mut = len(self._mutations)
        n_del = len(self._deleted_residues)
        self.lbl_pending_hdr.setText(f"Pending Edits ({n_mut + n_del}):")

        checked_chains = sum(1 for i in range(self.list_chains.count()) if self.list_chains.item(i).checkState() == Qt.Checked)
        self.lbl_edits_badge.setText(
            f"{checked_chains} chains removed | {n_mut} mutations | {n_del} deleted"
        )

    def _open_receptor_folder(self) -> None:
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors")).resolve()
        name = self.ent_rec_name.text().strip()
        target = (rec_dir / name) if name else rec_dir
        target.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(target))
        else:
            subprocess.run(["xdg-open", str(target)])

    # ─────────────────────────────────────────────────────────
    # Receptor Preparation Execution
    # ─────────────────────────────────────────────────────────

    def _run_receptor_prep(self) -> None:
        if not self._pdb_path or not self._pdb_path.is_file():
            QMessageBox.warning(self, "Input Required", "Please select a valid receptor source file.")
            return

        name = self.ent_rec_name.text().strip()
        if not name:
            QMessageBox.warning(self, "Name Required", "Please enter a name for the prepared receptor.")
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        out_sub = rec_dir / name / "rigid"
        out_sub.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_sub / f"{name}.pdbqt"

        # Determine checked chains to remove
        remove_chains = set()
        for i in range(self.list_chains.count()):
            it = self.list_chains.item(i)
            if it.checkState() == Qt.Checked:
                ch = it.data(Qt.UserRole)
                if ch:
                    remove_chains.add(ch)

        # Determine checked HETATM keys to remove
        remove_keys = set()
        for i in range(self.list_hetatms.count()):
            it = self.list_hetatms.item(i)
            if it.checkState() == Qt.Checked:
                k = it.data(Qt.UserRole)
                if k:
                    remove_keys.add(k)

        mutations = self._mutations.copy()
        deleted_residues = self._deleted_residues.copy()
        charge_choice = "kollman" if "kollman" in self.cmb_charge_model.currentText().lower() else "gasteiger"

        self.rec_log.append(f"\n▶ Starting receptor preparation for: {self._pdb_path.name}")
        self.rec_log.append(f"  Target: {out_pdbqt}")
        self.rec_log.append(f"  Charge Model: {charge_choice.upper()} (United-Atom: {charge_choice == 'kollman'})")
        self.rec_log.append(f"  Chains to remove: {len(remove_chains)} | HETATMs to remove: {len(remove_keys)}")
        self.rec_log.append(f"  Mutations: {len(mutations)} | Deletions: {len(deleted_residues)}")

        self.btn_run_rec_prep.setEnabled(False)

        import threading
        def _worker():
            try:
                from prepare import prepare_receptor_pdbqt
                prepare_receptor_pdbqt(
                    input_file=self._pdb_path,
                    output_pdbqt=out_pdbqt,
                    cleanup_water=self.chk_remove_water.isChecked(),
                    add_hydrogens=self.chk_add_h.isChecked(),
                    charge_model=charge_choice,
                    remove_residue_keys=remove_keys,
                    remove_chains=remove_chains,
                    mutations=mutations,
                    deleted_residues=deleted_residues,
                )
                if out_pdbqt.is_file():
                    self._signals.log_signal.emit(f"✔ Receptor prepared successfully: {out_pdbqt.name}")
                    self._signals.done_signal.emit("Receptor preparation completed.", None)
                else:
                    self._signals.done_signal.emit("Output PDBQT was not generated.", RuntimeError("File missing"))
            except Exception as e:
                self._signals.done_signal.emit(f"Error: {e}", e)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_worker_log(self, msg: str) -> None:
        self.rec_log.append(msg)

    def _on_worker_done(self, msg: str, error: Optional[Exception]) -> None:
        self.btn_run_rec_prep.setEnabled(True)
        if error:
            self.rec_log.append(f"✖ Failed: {error}")
            QMessageBox.critical(self, "Preparation Failed", f"Receptor preparation failed:\n{error}")
        else:
            self.rec_log.append(f"✔ {msg}")
            self._update_stepper_stage(6)
            QMessageBox.information(self, "Preparation Complete", "Receptor prepared successfully!")
            self.on_config_change()

    # ─────────────────────────────────────────────────────────
    # Flex Split Actions
    # ─────────────────────────────────────────────────────────

    def _browse_flex_receptor(self) -> None:
        f, _ = QFileDialog.getOpenFileName(
            self, "Select Prepared Receptor PDBQT",
            str(getattr(self.config, "receptor_directory", Path.cwd())),
            "PDBQT Files (*.pdbqt);;All Files (*.*)"
        )
        if f:
            self.ent_flex_source.setText(f)
            self._load_flex_residues(Path(f))

    def _load_flex_residues(self, path: Path) -> None:
        self.list_flex_residues.clear()
        if not path.is_file():
            return
        from prepare import inspect_receptor_residues
        try:
            res_list = inspect_receptor_residues(path)
            for r in res_list:
                chain = r.get("chain", "A")
                seq = r.get("seq", 0)
                name = r.get("name", "RES")
                item = QListWidgetItem(f"[{chain}:{seq}] {name}")
                item.setData(Qt.UserRole, f"{chain}:{seq}")
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Unchecked)
                self.list_flex_residues.addItem(item)
            self.list_flex_residues.itemChanged.connect(self._update_flex_count)
            self.flex_log.append(f"✔ Loaded {len(res_list)} residues from {path.name}")
        except Exception as e:
            self.flex_log.append(f"⚠ Could not parse residues: {e}")

    def _filter_flex_residues(self, query: str) -> None:
        q = query.strip().lower()
        for i in range(self.list_flex_residues.count()):
            it = self.list_flex_residues.item(i)
            it.setHidden(q not in it.text().lower())

    def _update_flex_count(self, item=None) -> None:
        checked = sum(1 for i in range(self.list_flex_residues.count())
                      if self.list_flex_residues.item(i).checkState() == Qt.Checked)
        self.lbl_flex_summary.setText(f"{checked} residue(s) selected for sidechain flexibility.")

    def _run_flex_split(self) -> None:
        src = self.ent_flex_source.text().strip()
        if not src or not Path(src).is_file():
            QMessageBox.warning(self, "Input Required", "Select a valid prepared receptor PDBQT.")
            return

        selected_res = []
        for i in range(self.list_flex_residues.count()):
            it = self.list_flex_residues.item(i)
            if it.checkState() == Qt.Checked:
                selected_res.append(it.data(Qt.UserRole))

        if not selected_res:
            QMessageBox.warning(self, "Residues Required", "Please check at least one flexible residue.")
            return

        src_path = Path(src)
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        name = src_path.stem.replace("_rigid", "").replace("_flex", "")
        out_flex_dir = rec_dir / name / "flex"
        out_flex_dir.mkdir(parents=True, exist_ok=True)
        out_rigid = out_flex_dir / f"{name}_rigid.pdbqt"
        out_flex = out_flex_dir / f"{name}_flex.pdbqt"

        self.flex_log.append(f"\n▶ Splitting flexible residues for: {src_path.name}")
        self.flex_log.append(f"  Selected: {len(selected_res)} residues")

        import threading
        def _worker():
            try:
                from prepare import split_flexible_receptor
                split_flexible_receptor(src_path, selected_res, out_rigid, out_flex)
                self._signals.flex_log_signal.emit(f"✔ Flexible receptor split complete in {out_flex_dir}")
                self._signals.flex_done_signal.emit("Flexible split completed successfully!", None)
            except Exception as e:
                self._signals.flex_log_signal.emit(f"✖ Split error: {e}")
                self._signals.flex_done_signal.emit(f"Split error: {e}", e)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_flex_log(self, msg: str) -> None:
        self.flex_log.append(msg)

    def _on_flex_done(self, msg: str, error: Optional[Exception]) -> None:
        if error:
            QMessageBox.critical(self, "Flex Split Failed", f"Flexible splitting failed:\n{error}")
        else:
            self._update_stepper_stage(6)
            QMessageBox.information(self, "Flex Split Complete", msg)
            self.on_config_change()

    # ─────────────────────────────────────────────────────────
    # Ligand Batch Prep Actions
    # ─────────────────────────────────────────────────────────

    def _browse_source_ligands(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Ligand Files",
            str(getattr(self.config, "ligand_directory", Path.cwd())),
            "Molecular Structures (*.sdf *.mol2 *.smi *.pdb *.mol);;All Files (*.*)"
        )
        if files:
            self.ent_lig_source.setText("; ".join(files))

    def _run_ligand_prep(self) -> None:
        raw_text = self.ent_lig_source.text().strip()
        if not raw_text:
            QMessageBox.warning(self, "Input Required", "Please select source ligand file(s).")
            return

        paths = [Path(p.strip()) for p in raw_text.split(";") if p.strip()]
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        lig_dir.mkdir(parents=True, exist_ok=True)

        gen_3d = self.chk_gen3d.isChecked()

        self.lig_log.append(f"\n▶ Starting ligand batch preparation ({len(paths)} file(s))")

        import threading
        def _worker():
            from prepare import prepare_ligand_pdbqt
            count = 0
            for p in paths:
                if not p.is_file():
                    continue
                out_pdbqt = lig_dir / f"{p.stem}.pdbqt"
                try:
                    res = prepare_ligand_pdbqt(
                        input_file=p,
                        output_pdbqt=out_pdbqt,
                        generate_3d=gen_3d,
                    )
                    if out_pdbqt.is_file():
                        self._signals.lig_log_signal.emit(f"  ✔ {p.name} → {out_pdbqt.name}")
                        count += 1
                    else:
                        self._signals.lig_log_signal.emit(f"  ✖ {p.name} failed to generate PDBQT")
                except Exception as e:
                    self._signals.lig_log_signal.emit(f"  ✖ {p.name} error: {e}")
            self._signals.lig_log_signal.emit(f"\n✔ Finished preparing {count} of {len(paths)} ligands.")
            self._signals.lig_done_signal.emit(f"Prepared {count} of {len(paths)} ligands.", None)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_lig_log(self, msg: str) -> None:
        self.lig_log.append(msg)

    def _on_lig_done(self, msg: str, error: Optional[Exception]) -> None:
        if error:
            QMessageBox.critical(self, "Ligand Prep Failed", f"Ligand preparation failed:\n{error}")
        else:
            QMessageBox.information(self, "Ligand Prep Complete", msg)
            self.on_config_change()
