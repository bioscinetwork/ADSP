"""
AutoDock Suite Pro — PySide6 Screening Library Tab (gui_qt/library_tab.py)
=========================================================================
Displays all ligands in the library with computed ADMET physicochemical
properties, Lipinski Rule of 5 and Veber rule compliance, summary KPI metrics,
interactive numeric sorting, and CSV export.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Dict, List

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig


class _NumericTableItem(QTableWidgetItem):
    """QTableWidgetItem supporting numeric sorting."""
    def __init__(self, text: str, sort_val: float):
        super().__init__(text)
        self.sort_val = sort_val

    def __lt__(self, other):
        if isinstance(other, _NumericTableItem):
            return self.sort_val < other.sort_val
        return super().__lt__(other)


class _ScanSignals(QObject):
    data_ready = Signal(list)
    progress_msg = Signal(str)


class LibraryTab(QWidget):
    """Screening Library Explorer in PySide6 with physicochemical ADMET profiling."""

    COLS = (
        "Name", "Formula", "MW (Da)", "LogP", "PSA (Å²)",
        "HBD", "HBA", "RotBonds", "Arom.Rings", "RO5", "Veber",
    )

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent=None):
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self._data: List[Dict] = []
        self._filtered_data: List[Dict] = []

        self._signals = _ScanSignals()
        self._signals.data_ready.connect(self._on_scan_ready)
        self._signals.progress_msg.connect(self._on_progress_msg)

        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(10)

        # Header Bar
        hdr_layout = QHBoxLayout()
        lbl_title = QLabel("Virtual Screening Library & Physicochemical Profiler (ADMET)")
        lbl_title.setObjectName("sectionHeader")
        hdr_layout.addWidget(lbl_title)
        hdr_layout.addStretch()

        self.btn_export_csv = QPushButton("📥 Export Filtered CSV")
        self.btn_export_csv.clicked.connect(self._export_csv)
        hdr_layout.addWidget(self.btn_export_csv)

        self.btn_scan = QPushButton("⟳ Scan Library")
        self.btn_scan.setObjectName("accentButton")
        self.btn_scan.clicked.connect(self._scan_library)
        hdr_layout.addWidget(self.btn_scan)
        main_layout.addLayout(hdr_layout)

        # ── KPI Summary Cards ──────────────────────────────────
        kpi_card = QFrame()
        kpi_card.setObjectName("cardFrame")
        kpi_l = QHBoxLayout(kpi_card)
        kpi_l.setContentsMargins(16, 10, 16, 10)

        self.lbl_kpi_total = self._create_kpi_cell(kpi_l, "Total Compounds", "0")
        self.lbl_kpi_ro5 = self._create_kpi_cell(kpi_l, "RO5 Compliant", "0%", "#10b981")
        self.lbl_kpi_veber = self._create_kpi_cell(kpi_l, "Veber Compliant", "0%", "#38bdf8")
        self.lbl_kpi_mw = self._create_kpi_cell(kpi_l, "Avg Molecular Weight", "0.0 Da")

        main_layout.addWidget(kpi_card)

        # Filter Bar
        filter_card = QFrame()
        filter_card.setObjectName("cardFrame")
        f_layout = QHBoxLayout(filter_card)
        f_layout.setContentsMargins(12, 8, 12, 8)
        f_layout.setSpacing(10)

        f_layout.addWidget(QLabel("Filter:"))
        self.ent_search = QLineEdit()
        self.ent_search.setPlaceholderText("Filter by name, formula, or substructure...")
        self.ent_search.setMinimumWidth(240)
        self.ent_search.textChanged.connect(self._apply_filter)
        f_layout.addWidget(self.ent_search)

        f_layout.addWidget(QLabel("Rule Compliance:"))
        self.cmb_filter_rule = QComboBox()
        self.cmb_filter_rule.addItems(["All", "RO5 Pass", "RO5 Fail", "Veber Pass", "Both Pass"])
        self.cmb_filter_rule.currentIndexChanged.connect(self._apply_filter)
        f_layout.addWidget(self.cmb_filter_rule)

        f_layout.addStretch()
        self.lbl_count = QLabel("0 compound(s)")
        self.lbl_count.setObjectName("subHeader")
        f_layout.addWidget(self.lbl_count)

        main_layout.addWidget(filter_card)

        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(len(self.COLS))
        self.table.setHorizontalHeaderLabels(self.COLS)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSortingEnabled(True)
        main_layout.addWidget(self.table, 1)

    def _create_kpi_cell(self, layout, title: str, default: str, val_color: str = "") -> QLabel:
        cell = QVBoxLayout()
        lbl_t = QLabel(title.upper())
        lbl_t.setObjectName("kpiTitle")
        lbl_v = QLabel(default)
        lbl_v.setObjectName("kpiValue")
        if val_color:
            lbl_v.setStyleSheet(f"color: {val_color};")
        cell.addWidget(lbl_t)
        cell.addWidget(lbl_v)
        layout.addLayout(cell)
        layout.addStretch()
        return lbl_v

    def refresh(self) -> None:
        self._scan_library()

    def _scan_library(self) -> None:
        lig_dir = Path(getattr(self.config, "ligand_directory", "ligands"))
        if not lig_dir.is_dir():
            self.lbl_count.setText("Ligand directory not found")
            return

        self.btn_scan.setEnabled(False)
        self.btn_scan.setText("Scanning...")

        import threading
        def _worker():
            from prepare import calculate_molecule_descriptors
            data = []
            files = sorted(lig_dir.glob("*.pdbqt"))
            for p in files:
                desc = calculate_molecule_descriptors(p, obabel_exe=self.config.obabel_executable)
                if desc:
                    data.append(desc)
                else:
                    data.append({
                        "name": p.stem, "formula": "—", "mw": 0.0, "logp": 0.0,
                        "psa": 0.0, "hbd": 0, "hba": 0, "rot_bonds": 0,
                        "aromatic_rings": 0, "ro5_pass": False, "veber_pass": False,
                    })
            try:
                self._signals.data_ready.emit(data)
            except RuntimeError:
                pass

        threading.Thread(target=_worker, daemon=True).start()

    def _on_progress_msg(self, msg: str) -> None:
        self.lbl_count.setText(msg)

    def _on_scan_ready(self, data: List[Dict]) -> None:
        self._data = data
        self.btn_scan.setEnabled(True)
        self.btn_scan.setText("⟳ Scan Library")
        self._apply_filter()

    def _apply_filter(self) -> None:
        query = self.ent_search.text().strip().lower()
        rule_sel = self.cmb_filter_rule.currentText()

        filtered = []
        for d in self._data:
            name = str(d.get("name", "")).lower()
            formula = str(d.get("formula", "")).lower()
            if query and (query not in name and query not in formula):
                continue
            ro5 = d.get("ro5_pass", False)
            veber = d.get("veber_pass", False)
            if rule_sel == "RO5 Pass" and not ro5:
                continue
            if rule_sel == "RO5 Fail" and ro5:
                continue
            if rule_sel == "Veber Pass" and not veber:
                continue
            if rule_sel == "Both Pass" and not (ro5 and veber):
                continue
            filtered.append(d)

        self._filtered_data = filtered

        # Update KPI Summary
        total = len(self._data)
        if total > 0:
            ro5_count = sum(1 for d in self._data if d.get("ro5_pass", False))
            veber_count = sum(1 for d in self._data if d.get("veber_pass", False))
            avg_mw = sum(d.get("mw", 0.0) for d in self._data) / total
            self.lbl_kpi_total.setText(str(total))
            self.lbl_kpi_ro5.setText(f"{ro5_count} ({ro5_count * 100 / total:.1f}%)")
            self.lbl_kpi_veber.setText(f"{veber_count} ({veber_count * 100 / total:.1f}%)")
            self.lbl_kpi_mw.setText(f"{avg_mw:.1f} Da")
        else:
            self.lbl_kpi_total.setText("0")
            self.lbl_kpi_ro5.setText("0%")
            self.lbl_kpi_veber.setText("0%")
            self.lbl_kpi_mw.setText("0.0 Da")

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(filtered))

        for row, d in enumerate(filtered):
            ro5_pass = d.get("ro5_pass", False)
            veber_pass = d.get("veber_pass", False)

            item_name = QTableWidgetItem(str(d.get("name", "")))
            item_formula = QTableWidgetItem(str(d.get("formula", "")))

            mw_val = float(d.get("mw", 0.0))
            item_mw = _NumericTableItem(f"{mw_val:.1f}", mw_val)
            item_mw.setTextAlignment(Qt.AlignCenter)

            logp_val = float(d.get("logp", 0.0))
            item_logp = _NumericTableItem(f"{logp_val:.2f}", logp_val)
            item_logp.setTextAlignment(Qt.AlignCenter)

            psa_val = float(d.get("psa", 0.0))
            item_psa = _NumericTableItem(f"{psa_val:.1f}", psa_val)
            item_psa.setTextAlignment(Qt.AlignCenter)

            hbd_val = int(d.get("hbd", 0))
            item_hbd = _NumericTableItem(str(hbd_val), hbd_val)
            item_hbd.setTextAlignment(Qt.AlignCenter)

            hba_val = int(d.get("hba", 0))
            item_hba = _NumericTableItem(str(hba_val), hba_val)
            item_hba.setTextAlignment(Qt.AlignCenter)

            rot_val = int(d.get("rot_bonds", 0))
            item_rot = _NumericTableItem(str(rot_val), rot_val)
            item_rot.setTextAlignment(Qt.AlignCenter)

            arom_val = int(d.get("aromatic_rings", 0))
            item_arom = _NumericTableItem(str(arom_val), arom_val)
            item_arom.setTextAlignment(Qt.AlignCenter)

            item_ro5 = QTableWidgetItem("✔ PASS" if ro5_pass else "✖ FAIL")
            item_ro5.setTextAlignment(Qt.AlignCenter)
            item_ro5.setForeground(QColor("#34d399") if ro5_pass else QColor("#f87171"))

            item_veber = QTableWidgetItem("✔ PASS" if veber_pass else "✖ FAIL")
            item_veber.setTextAlignment(Qt.AlignCenter)
            item_veber.setForeground(QColor("#34d399") if veber_pass else QColor("#f87171"))

            row_items = [
                item_name, item_formula, item_mw, item_logp, item_psa,
                item_hbd, item_hba, item_rot, item_arom, item_ro5, item_veber
            ]
            for col, item in enumerate(row_items):
                self.table.setItem(row, col, item)

        self.table.setSortingEnabled(True)
        self.lbl_count.setText(f"{len(filtered)} / {len(self._data)} compound(s)")

    def _export_csv(self) -> None:
        if not self._filtered_data:
            QMessageBox.information(self, "No Data", "There are no compounds to export.")
            return

        f, _ = QFileDialog.getSaveFileName(
            self, "Export Screening Library CSV",
            str(getattr(self.config, "project_root", Path.cwd()) / "screening_library.csv"),
            "CSV Files (*.csv);;All Files (*.*)"
        )
        if f:
            try:
                with open(f, "w", newline="", encoding="utf-8") as fp:
                    writer = csv.writer(fp)
                    writer.writerow(self.COLS)
                    for d in self._filtered_data:
                        writer.writerow([
                            d.get("name", ""),
                            d.get("formula", ""),
                            f"{d.get('mw', 0.0):.2f}",
                            f"{d.get('logp', 0.0):.2f}",
                            f"{d.get('psa', 0.0):.2f}",
                            d.get("hbd", 0),
                            d.get("hba", 0),
                            d.get("rot_bonds", 0),
                            d.get("aromatic_rings", 0),
                            "PASS" if d.get("ro5_pass", False) else "FAIL",
                            "PASS" if d.get("veber_pass", False) else "FAIL",
                        ])
                QMessageBox.information(self, "Export Complete", f"Saved {len(self._filtered_data)} compounds to:\n{f}")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to save CSV:\n{e}")
