"""
AutoDock Suite Pro — PySide6 Workspace Tab (gui_qt/workspace_tab.py)
===================================================================
Project directory management, drag-and-drop receptor and ligand import/removal,
workspace scaffolding, and multi-file inspection.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Set

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig


class _DropListWidget(QListWidget):
    """QListWidget supporting drag-and-drop file imports from Windows Explorer."""
    files_dropped = Signal(list)

    def __init__(self, allowed_exts: Set[str], parent=None):
        super().__init__(parent)
        self.allowed_exts = {ext.lower() for ext in allowed_exts}
        self.setAcceptDrops(True)
        self.setAlternatingRowColors(True)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            for u in event.mimeData().urls():
                p = Path(u.toLocalFile())
                if p.is_file() and p.suffix.lower() in self.allowed_exts:
                    event.acceptProposedAction()
                    return
        event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            paths = []
            for u in event.mimeData().urls():
                p = Path(u.toLocalFile())
                if p.is_file() and p.suffix.lower() in self.allowed_exts:
                    paths.append(p)
            if paths:
                event.acceptProposedAction()
                self.files_dropped.emit(paths)
                return
        event.ignore()


class WorkspaceTab(QWidget):
    """PySide6 Workspace management tab with Drag-and-Drop file import."""

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self._build_ui()
        self.refresh()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(14)

        # ── Workspace Folder Card ──────────────────────────────
        ws_card = QFrame()
        ws_card.setObjectName("cardFrame")
        ws_layout = QHBoxLayout(ws_card)
        ws_layout.setContentsMargins(16, 12, 16, 12)
        ws_layout.setSpacing(12)

        lbl_ws_tag = QLabel("Active Workspace:")
        lbl_ws_tag.setObjectName("subHeader")
        ws_layout.addWidget(lbl_ws_tag)

        self.lbl_workspace = QLabel("")
        self.lbl_workspace.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.lbl_workspace.setStyleSheet("font-family: Consolas, monospace; font-size: 12px; font-weight: 600;")
        ws_layout.addWidget(self.lbl_workspace, 1)

        self.btn_scaffold = QPushButton("🏗 Scaffold Workspace")
        self.btn_scaffold.setToolTip("Create standard folder structure (receptors/, ligands/, results/, etc.)")
        self.btn_scaffold.clicked.connect(self._scaffold_workspace)
        ws_layout.addWidget(self.btn_scaffold)

        self.btn_change_ws = QPushButton("📂 Change Workspace")
        self.btn_change_ws.setObjectName("accentButton")
        self.btn_change_ws.clicked.connect(self._change_workspace)
        ws_layout.addWidget(self.btn_change_ws)

        main_layout.addWidget(ws_card)

        # ── Laboratory KPI Summary Row ─────────────────────────
        kpi_row = QHBoxLayout()
        kpi_row.setSpacing(12)

        self.lbl_kpi_receptors = self._create_kpi_card(kpi_row, "Receptors", "0", "#0284c7")
        self.lbl_kpi_ligands = self._create_kpi_card(kpi_row, "Ligands", "0", "#0ea5e9")
        self.lbl_kpi_completed = self._create_kpi_card(kpi_row, "Completed Jobs", "0", "#10b981")
        self.lbl_kpi_status = self._create_kpi_card(kpi_row, "Workspace Status", "Ready", "#8b5cf6")

        main_layout.addLayout(kpi_row)

        # ── Receptors & Ligands Columns ────────────────────────
        columns_layout = QHBoxLayout()
        columns_layout.setSpacing(14)

        # Left Column: Receptors
        rec_box = QGroupBox("RECEPTORS (Macromolecules — Drag & Drop PDB/PDBQT/CIF here)")
        rec_layout = QVBoxLayout(rec_box)
        rec_layout.setContentsMargins(12, 16, 12, 12)
        rec_layout.setSpacing(8)

        rec_actions = QHBoxLayout()
        self.btn_add_rec = QPushButton("➕ Add Receptor")
        self.btn_add_rec.setObjectName("accentButton")
        self.btn_add_rec.clicked.connect(self._add_receptor)
        rec_actions.addWidget(self.btn_add_rec)

        self.btn_rem_rec = QPushButton("🗑 Remove Selected")
        self.btn_rem_rec.clicked.connect(self._remove_receptor)
        rec_actions.addWidget(self.btn_rem_rec)
        rec_actions.addStretch()

        self.lbl_rec_count = QLabel("0 receptors")
        self.lbl_rec_count.setObjectName("subHeader")
        rec_actions.addWidget(self.lbl_rec_count)
        rec_layout.addLayout(rec_actions)

        self.list_receptors = _DropListWidget({".pdbqt", ".pdb", ".cif", ".ent"})
        self.list_receptors.files_dropped.connect(self._handle_receptors_dropped)
        self.list_receptors.itemDoubleClicked.connect(self._on_receptor_double_clicked)
        rec_layout.addWidget(self.list_receptors)

        columns_layout.addWidget(rec_box, 1)

        # Right Column: Ligands
        lig_box = QGroupBox("LIGANDS (Small Molecules — Drag & Drop PDBQT/SDF/MOL2 here)")
        lig_layout = QVBoxLayout(lig_box)
        lig_layout.setContentsMargins(12, 16, 12, 12)
        lig_layout.setSpacing(8)

        lig_actions = QHBoxLayout()
        self.btn_add_lig = QPushButton("➕ Add Ligand(s)")
        self.btn_add_lig.setObjectName("accentButton")
        self.btn_add_lig.clicked.connect(self._add_ligand)
        lig_actions.addWidget(self.btn_add_lig)

        self.btn_rem_lig = QPushButton("🗑 Remove Selected")
        self.btn_rem_lig.clicked.connect(self._remove_ligand)
        lig_actions.addWidget(self.btn_rem_lig)
        lig_actions.addStretch()

        self.lbl_lig_count = QLabel("0 ligands")
        self.lbl_lig_count.setObjectName("subHeader")
        lig_actions.addWidget(self.lbl_lig_count)
        lig_layout.addLayout(lig_actions)

        self.list_ligands = _DropListWidget({".pdbqt", ".sdf", ".mol2", ".mol", ".smi"})
        self.list_ligands.files_dropped.connect(self._handle_ligands_dropped)
        self.list_ligands.itemDoubleClicked.connect(self._on_ligand_double_clicked)
        lig_layout.addWidget(self.list_ligands)

        columns_layout.addWidget(lig_box, 1)
        main_layout.addLayout(columns_layout, 1)

        # ── Quick Open Folder Toolbar ──────────────────────────
        footer_layout = QHBoxLayout()
        self.btn_open_rec = QPushButton("📁 Open Receptors Folder")
        self.btn_open_rec.clicked.connect(lambda: self._open_dir(self.config.receptor_directory))
        footer_layout.addWidget(self.btn_open_rec)

        self.btn_open_lig = QPushButton("📁 Open Ligands Folder")
        self.btn_open_lig.clicked.connect(lambda: self._open_dir(self.config.ligand_directory))
        footer_layout.addWidget(self.btn_open_lig)

        self.btn_open_root = QPushButton("📁 Open Project Root")
        self.btn_open_root.clicked.connect(lambda: self._open_dir(self.config.project_root))
        footer_layout.addWidget(self.btn_open_root)

        footer_layout.addStretch()

        self.btn_refresh = QPushButton("⟳ Refresh")
        self.btn_refresh.clicked.connect(self.refresh)
        footer_layout.addWidget(self.btn_refresh)

        main_layout.addLayout(footer_layout)

    def _create_kpi_card(self, layout: QHBoxLayout, label: str, default: str, accent_color: str = "") -> QLabel:
        card = QFrame()
        card.setObjectName("cardFrame")
        cl = QVBoxLayout(card)
        cl.setContentsMargins(14, 10, 14, 10)
        cl.setSpacing(2)

        lbl_title = QLabel(label.upper())
        lbl_title.setStyleSheet("font-size: 10px; font-weight: 700; color: #64748b; letter-spacing: 0.5px;")
        cl.addWidget(lbl_title)

        lbl_val = QLabel(default)
        lbl_val.setStyleSheet(
            f"font-size: 18px; font-weight: 800; font-family: Consolas, monospace; "
            f"color: {accent_color if accent_color else '#0f172a'};"
        )
        cl.addWidget(lbl_val)
        layout.addWidget(card, 1)
        return lbl_val

    # ─────────────────────────────────────────────────────────
    # Logic & Refresh
    # ─────────────────────────────────────────────────────────

    def refresh(self) -> None:
        """Scan directories and refresh the lists and KPI summary cards."""
        p_root = getattr(self.config, "project_root", Path.cwd())
        self.lbl_workspace.setText(str(p_root))

        # Refresh receptors
        self.list_receptors.clear()
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        rec_count = 0
        if rec_dir.is_dir():
            for entry in sorted(rec_dir.iterdir()):
                if entry.is_dir():
                    rigid = entry / "rigid"
                    rigid_files = list(rigid.glob("*.pdbqt")) if rigid.is_dir() else []
                    flex = entry / "flex"
                    flex_files = list(flex.glob("*.pdbqt")) if flex.is_dir() else []
                    cfg = entry / "config.txt"
                    cfg_badge = " [Grid: ✔]" if cfg.is_file() else " [Grid: ✖]"
                    flex_badge = f" [Flex: {len(flex_files)}]" if flex_files else ""
                    if rigid_files:
                        self.list_receptors.addItem(f"📁 {entry.name} ({len(rigid_files)} rigid){flex_badge}{cfg_badge}")
                        rec_count += 1
                elif entry.suffix.lower() == ".pdbqt":
                    self.list_receptors.addItem(f"📄 {entry.name}")
                    rec_count += 1
        self.lbl_rec_count.setText(f"{rec_count} receptor(s)")
        self.lbl_kpi_receptors.setText(str(rec_count))

        # Refresh ligands
        self.list_ligands.clear()
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        lig_count = 0
        if lig_dir.is_dir():
            for p in sorted(lig_dir.glob("*.pdbqt")):
                size_kb = p.stat().st_size / 1024
                self.list_ligands.addItem(f"💊 {p.name}  ({size_kb:.1f} KB)")
                lig_count += 1
        self.lbl_lig_count.setText(f"{lig_count} ligand(s)")
        self.lbl_kpi_ligands.setText(str(lig_count))

        # Scan completed results
        res_dir = Path(getattr(self.config, "result_directory", "results"))
        completed_count = 0
        if res_dir.is_dir():
            try:
                from job_manager import load_docking_jobs
                jobs = load_docking_jobs(res_dir)
                completed_count = len([j for j in jobs if "SUCCESS" in str(getattr(j, "status", "")).upper()])
            except Exception:
                completed_count = len(list(res_dir.glob("**/*.dlg"))) + len(list(res_dir.glob("**/*_out.pdbqt")))
        self.lbl_kpi_completed.setText(str(completed_count))

        # Status
        if rec_count > 0 and lig_count > 0:
            self.lbl_kpi_status.setText("Ready to Dock")
            self.lbl_kpi_status.setStyleSheet("font-size: 14px; font-weight: 800; color: #10b981;")
        elif rec_count > 0:
            self.lbl_kpi_status.setText("Needs Ligands")
            self.lbl_kpi_status.setStyleSheet("font-size: 14px; font-weight: 800; color: #f59e0b;")
        else:
            self.lbl_kpi_status.setText("Needs Receptors")
            self.lbl_kpi_status.setStyleSheet("font-size: 14px; font-weight: 800; color: #ef4444;")

    def _scaffold_workspace(self) -> None:
        from prepare import create_scaffold
        root = getattr(self.config, "project_root", Path.cwd())
        try:
            create_scaffold(root)
            QMessageBox.information(
                self, "Scaffold Created",
                f"Standard project folder scaffold created at:\n{root}\n\n"
                "• receptors/\n• ligands/\n• results/\n• bin/\n• logs/\n• reports/"
            )
            self.refresh()
            self.on_config_change()
        except Exception as e:
            QMessageBox.critical(self, "Scaffold Error", f"Failed to scaffold:\n{e}")

    def _change_workspace(self) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self, "Select Project Workspace Directory",
            str(getattr(self.config, "project_root", Path.cwd()))
        )
        if chosen:
            from config import set_workspace_directory
            try:
                set_workspace_directory(self.config, chosen)
                self.refresh()
                self.on_config_change()
            except Exception as e:
                QMessageBox.critical(self, "Workspace Error", f"Failed to activate workspace:\n{e}")

    def _add_receptor(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Receptor File(s)",
            str(getattr(self.config, "receptor_directory", Path.cwd())),
            "Receptor Files (*.pdbqt *.pdb *.cif *.ent);;All Files (*.*)"
        )
        if files:
            self._import_receptors([Path(f) for f in files])

    def _handle_receptors_dropped(self, paths: List[Path]) -> None:
        self._import_receptors(paths)

    def _import_receptors(self, paths: List[Path]) -> None:
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        rec_dir.mkdir(parents=True, exist_ok=True)
        for p in paths:
            if not p.is_file():
                continue
            if p.suffix.lower() == ".pdbqt":
                target_sub = rec_dir / p.stem / "rigid"
                target_sub.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, target_sub / p.name)
            else:
                # Raw PDB/CIF, copy directly into receptor root for preparation
                shutil.copy2(p, rec_dir / p.name)
        self.refresh()
        self.on_config_change()

    def _remove_receptor(self) -> None:
        items = self.list_receptors.selectedItems()
        if not items:
            QMessageBox.information(self, "Selection", "Select one or more receptors to remove.")
            return

        res = QMessageBox.question(
            self, "Confirm Delete",
            f"Are you sure you want to remove {len(items)} selected receptor(s)?",
            QMessageBox.Yes | QMessageBox.No
        )
        if res != QMessageBox.Yes:
            return

        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        for it in items:
            name = it.text().split()[1]
            target = rec_dir / name
            if not target.exists():
                target = rec_dir / f"{name}.pdbqt"
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()

        self.refresh()
        self.on_config_change()

    def _add_ligand(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Ligand File(s)",
            str(getattr(self.config, "ligand_directory", Path.cwd())),
            "Ligand Files (*.pdbqt *.sdf *.mol2 *.mol *.smi);;All Files (*.*)"
        )
        if files:
            self._import_ligands([Path(f) for f in files])

    def _handle_ligands_dropped(self, paths: List[Path]) -> None:
        self._import_ligands(paths)

    def _import_ligands(self, paths: List[Path]) -> None:
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        lig_dir.mkdir(parents=True, exist_ok=True)
        for p in paths:
            if p.is_file():
                shutil.copy2(p, lig_dir / p.name)
        self.refresh()
        self.on_config_change()

    def _remove_ligand(self) -> None:
        items = self.list_ligands.selectedItems()
        if not items:
            QMessageBox.information(self, "Selection", "Select one or more ligands to remove.")
            return

        res = QMessageBox.question(
            self, "Confirm Delete",
            f"Are you sure you want to remove {len(items)} selected ligand(s)?",
            QMessageBox.Yes | QMessageBox.No
        )
        if res != QMessageBox.Yes:
            return

        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        for it in items:
            name = it.text().split()[1]
            target = lig_dir / name
            if target.is_file():
                target.unlink()

        self.refresh()
        self.on_config_change()

    def _on_receptor_double_clicked(self, item) -> None:
        name = item.text().split()[1]
        rec_dir = Path(getattr(self.config, "receptor_directory", "receptors"))
        target = rec_dir / name
        if not target.exists():
            target = rec_dir / f"{name}.pdbqt"
        if target.exists():
            self._open_dir(target if target.is_dir() else target.parent)

    def _on_ligand_double_clicked(self, item) -> None:
        name = item.text().split()[1]
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        target = lig_dir / name
        if target.is_file():
            self._open_dir(lig_dir)

    def _open_dir(self, directory: Path | str) -> None:
        p = Path(directory)
        p.mkdir(parents=True, exist_ok=True)
        if sys.platform == "win32":
            os.startfile(str(p))
        else:
            subprocess.run(["xdg-open", str(p)])
