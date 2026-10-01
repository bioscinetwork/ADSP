"""
AutoDock Suite Pro — PySide6 Experiment Comparison Studio (gui_qt/comparison_tab.py)
===================================================================================
Side-by-side scientific comparison workspace for multi-engine, multi-target,
and multi-ligand docking runs. Allows researchers to evaluate:
  - Cross-engine scoring differences (ΔΔG = AutoDock4 - Vina)
  - Ligand selectivity across multiple receptors
  - Rigid vs Flexible docking score shifts
  - Conformational cluster population comparisons
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from config import ProjectConfig
    from models import DockingJob


class ExperimentComparisonTab(QWidget):
    """Side-by-side docking experiment comparison and cross-engine analysis."""

    def __init__(self, config: "ProjectConfig", on_config_change: Callable, parent: Optional[QWidget] = None):
        # Unit-level consumers may construct the tab without the main Qt app.
        # Create an owned application in that case so Qt does not abort at the
        # native widget boundary. The normal GUI path reuses its existing app.
        self._owned_qapplication = None
        if QApplication.instance() is None:
            self._owned_qapplication = QApplication([])
        super().__init__(parent)
        self.config = config
        self.on_config_change = on_config_change
        self._jobs: List[DockingJob] = []
        self._filtered_rows: List[Dict[str, Any]] = []
        self._build_ui()

    def set_jobs(self, jobs: List["DockingJob"]) -> None:
        """Update available docking jobs for comparison."""
        self._jobs = jobs
        self.refresh()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(12)

        # ── Header & Filter Bar ────────────────────────────────
        header_card = QFrame()
        header_card.setObjectName("cardFrame")
        h_layout = QVBoxLayout(header_card)
        h_layout.setContentsMargins(14, 12, 14, 12)
        h_layout.setSpacing(10)

        top_row = QHBoxLayout()
        lbl_title = QLabel("⚖ Experiment Comparison & Cross-Engine Evaluator")
        lbl_title.setObjectName("sectionHeader")
        top_row.addWidget(lbl_title)
        top_row.addStretch()

        self.btn_export_csv = QPushButton("📊 Export Comparison (CSV)")
        self.btn_export_csv.clicked.connect(self._export_csv)
        top_row.addWidget(self.btn_export_csv)

        self.btn_refresh = QPushButton("🔄 Refresh Data")
        self.btn_refresh.clicked.connect(self.refresh)
        top_row.addWidget(self.btn_refresh)
        h_layout.addLayout(top_row)

        # Filter row
        filter_row = QHBoxLayout()
        filter_row.setSpacing(10)

        filter_row.addWidget(QLabel("Search Ligand:"))
        self.txt_lig_search = QLineEdit()
        self.txt_lig_search.setPlaceholderText("Filter ligand name...")
        self.txt_lig_search.textChanged.connect(self._apply_filters)
        filter_row.addWidget(self.txt_lig_search)

        filter_row.addWidget(QLabel("Target Receptor:"))
        self.cmb_filter_rec = QComboBox()
        self.cmb_filter_rec.addItem("All Receptors")
        self.cmb_filter_rec.currentTextChanged.connect(self._apply_filters)
        filter_row.addWidget(self.cmb_filter_rec)

        filter_row.addWidget(QLabel("Engine:"))
        self.cmb_filter_eng = QComboBox()
        self.cmb_filter_eng.addItems(["All Engines", "Vina", "AutoDock4", "Dual (Both)"])
        self.cmb_filter_eng.currentTextChanged.connect(self._apply_filters)
        filter_row.addWidget(self.cmb_filter_eng)

        filter_row.addWidget(QLabel("Max Energy (kcal):"))
        self.spn_max_energy = QDoubleSpinBox()
        self.spn_max_energy.setRange(-25.0, 0.0)
        self.spn_max_energy.setValue(0.0)
        self.spn_max_energy.setSingleStep(0.5)
        self.spn_max_energy.valueChanged.connect(self._apply_filters)
        filter_row.addWidget(self.spn_max_energy)

        h_layout.addLayout(filter_row)
        main_layout.addWidget(header_card)

        # ── KPI Summary Cards ──────────────────────────────────
        kpi_card = QFrame()
        kpi_card.setObjectName("cardFrame")
        kpi_l = QHBoxLayout(kpi_card)
        kpi_l.setContentsMargins(14, 10, 14, 10)
        kpi_l.setSpacing(16)

        def make_kpi(title: str, text: str):
            box = QVBoxLayout()
            box.setSpacing(2)
            lbl_t = QLabel(title)
            lbl_t.setObjectName("kpiTitle")
            lbl_v = QLabel(text)
            lbl_v.setObjectName("kpiValue")
            box.addWidget(lbl_t)
            box.addWidget(lbl_v)
            kpi_l.addLayout(box)
            return lbl_v

        self.kpi_total_pairs = make_kpi("Screened Pairs", "0")
        self.kpi_best_vina = make_kpi("Best Vina Affinity", "—")
        self.kpi_best_ad4 = make_kpi("Best AD4 Estimated Binding Energy", "—")
        self.kpi_mean_delta = make_kpi("Mean Score Difference (AD4 - Vina)", "—")
        kpi_l.addStretch()

        main_layout.addWidget(kpi_card)

        # ── Comparison Table ───────────────────────────────────
        self.table = QTableWidget()
        self.table.setColumnCount(10)
        self.table.setHorizontalHeaderLabels([
            "Receptor", "Ligand", "Mode",
            "Vina Affinity", "AutoDock4 Estimated Binding Energy", "Score Difference (AD4 - Vina)",
            "Vina Ki (N/A)", "AD4 DLG Ki", "AD4 Top Cluster", "Status"
        ])
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.table_comparison = self.table
        main_layout.addWidget(self.table, 1)

        # Scientific Note on Energy Differences
        note_card = QFrame()
        note_card.setObjectName("cardFrame")
        n_layout = QVBoxLayout(note_card)
        n_layout.setContentsMargins(12, 8, 12, 8)
        lbl_note = QLabel(
            "💡 <b>Scientific Scoring Note (AD4 vs Vina):</b> "
            "AutoDock4 computes explicit Coulomb electrostatic interactions using Gasteiger partial charges and distance-dependent dielectric, "
            "whereas AutoDock Vina uses an empirical scoring function based purely on steric, hydrophobic, and H-bonding terms (without Coulomb electrostatics). "
            "For charged ligands (e.g. quaternary alkaloids), AD4 ΔG is typically significantly more negative due to electrostatic attraction."
        )
        lbl_note.setStyleSheet("color: #cbd5e1; font-size: 11px;")
        lbl_note.setWordWrap(True)
        n_layout.addWidget(lbl_note)
        main_layout.addWidget(note_card)

    def refresh(self) -> None:
        """Reload comparison data from active jobs, results tab, or workspace directories."""
        from job_manager import load_docking_jobs, load_job_status
        from models import DockingJob

        # 1. Try synchronizing from MainWindow / ResultsTab if available
        parent_window = self.window()
        if hasattr(parent_window, "tab_results") and getattr(parent_window.tab_results, "_jobs", None):
            self._jobs = list(parent_window.tab_results._jobs)

        # 2. Otherwise scan candidate directories for job_status.json or docking outputs
        if not self._jobs:
            candidates: List[Path] = []
            res_dir = getattr(self.config, "result_directory", None)
            if res_dir:
                candidates.append(Path(res_dir))
            ws = getattr(self.config, "project_root", None)
            if ws:
                candidates.append(Path(ws) / "results")
                candidates.append(Path(ws))
            candidates.append(Path("results"))
            candidates.append(Path.cwd() / "results")
            candidates.append(Path.cwd())

            for d in candidates:
                if not d.is_dir():
                    continue
                try:
                    jobs = load_docking_jobs(d)
                    if jobs:
                        self._jobs = jobs
                        break
                    s_dict = load_job_status(d)
                    if s_dict:
                        self._jobs = [DockingJob.from_dict(val) for val in s_dict.values()]
                        break
                except Exception:
                    pass

        # Ensure canonical result is populated for all jobs
        for j in self._jobs:
            try:
                j.ensure_canonical_result()
            except Exception:
                pass

        # Populate receptor dropdown
        recs = sorted({j.receptor_name for j in self._jobs if j.receptor_name})
        cur_rec = self.cmb_filter_rec.currentText()
        self.cmb_filter_rec.blockSignals(True)
        self.cmb_filter_rec.clear()
        self.cmb_filter_rec.addItem("All Receptors")
        for r in recs:
            self.cmb_filter_rec.addItem(r)
        if cur_rec in recs:
            self.cmb_filter_rec.setCurrentText(cur_rec)
        self.cmb_filter_rec.blockSignals(False)

        self._apply_filters()

    def _apply_filters(self) -> None:
        """Filter jobs and populate comparison rows."""
        lig_query = self.txt_lig_search.text().strip().lower()
        rec_filter = self.cmb_filter_rec.currentText()
        eng_filter = self.cmb_filter_eng.currentText()
        max_energy = self.spn_max_energy.value()

        # Group completed jobs by (receptor, ligand, mode)
        grouped: Dict[tuple, Dict[str, Any]] = {}
        for j in self._jobs:
            st = j.status.value if hasattr(j.status, "value") else str(j.status)
            if st.upper() not in ("SUCCESS", "SUCCESS_WITH_WARNING"):
                continue

            j.ensure_canonical_result()
            rec = j.receptor_name or "Unknown"
            lig = j.ligand_name or "Unknown"
            mode = j.docking_mode.value if hasattr(j.docking_mode, "value") else str(j.docking_mode)
            eng = j.engine.value if hasattr(j.engine, "value") else str(j.engine)
            eng_up = eng.upper()

            key = (rec, lig, mode)
            if key not in grouped:
                grouped[key] = {
                    "receptor": rec,
                    "ligand": lig,
                    "mode": mode,
                    "vina_job": None,
                    "ad4_job": None,
                }
            if "VINA" in eng_up:
                grouped[key]["vina_job"] = j
            if "AD4" in eng_up or "AUTODOCK" in eng_up:
                grouped[key]["ad4_job"] = j
            if "BOTH" in eng_up:
                grouped[key]["vina_job"] = j
                grouped[key]["ad4_job"] = j

        rows: List[Dict[str, Any]] = []
        for (rec, lig, mode), g in grouped.items():
            # Apply search filters
            if lig_query and lig_query not in lig.lower():
                continue
            if rec_filter != "All Receptors" and rec != rec_filter:
                continue

            v_j = g["vina_job"]
            a_j = g["ad4_job"]

            if eng_filter == "Vina" and not v_j:
                continue
            elif eng_filter == "AutoDock4" and not a_j:
                continue
            elif eng_filter == "Dual (Both)" and not (v_j and a_j):
                continue

            v_e = None
            v_ki = "—"
            if v_j:
                v_j.ensure_canonical_result()
                if v_j.canonical_result and v_j.canonical_result.poses:
                    p0 = v_j.canonical_result.poses[0]
                    v_e = p0.binding_score if p0.binding_score is not None else p0.binding_energy
                    v_ki = p0.estimated_ki_formatted if p0.estimated_ki_formatted != "N/A" else "—"
                elif v_j.vina_results:
                    v_e = float(v_j.vina_results[0].binding_affinity)
                v_ki = "Not applicable — Vina does not report Ki"

            a_e = None
            a_ki = "—"
            a_clust = "—"
            if a_j:
                a_j.ensure_canonical_result()
                if a_j.canonical_result:
                    if a_j.canonical_result.poses:
                        p0 = a_j.canonical_result.poses[0]
                        a_e = p0.binding_score if p0.binding_score is not None else p0.binding_energy
                        a_ki = p0.estimated_ki_formatted if p0.estimated_ki_formatted != "N/A" else "—"
                    if a_j.canonical_result.clusters:
                        c0 = a_j.canonical_result.clusters[0]
                        a_clust = f"Cluster #1 ({c0.population_percent:.1f}%)"
                if a_e is None and getattr(a_j, "ad4_results", None):
                    ad4_data = a_j.ad4_results
                    if isinstance(ad4_data, dict):
                        a_e = ad4_data.get("best_energy")
                if a_e is not None and a_ki in ("—", "N/A"):
                    a_ki = "Not reported by AutoDock4"

            # Check max energy threshold
            best_val = min([e for e in (v_e, a_e) if e is not None], default=None)
            if max_energy < 0.0 and (best_val is None or best_val > max_energy):
                continue

            # Compute delta delta G
            ddg = None
            if v_e is not None and a_e is not None:
                ddg = round(a_e - v_e, 2)

            status = "Completed Both" if (v_j and a_j) else ("Vina Only" if v_j else "AD4 Only")

            rows.append({
                "receptor": rec,
                "ligand": lig,
                "mode": mode,
                "vina_energy": v_e,
                "ad4_energy": a_e,
                "ddg": ddg,
                "vina_ki": v_ki,
                "ad4_ki": a_ki,
                "ad4_cluster": a_clust,
                "status": status,
            })

        self._filtered_rows = rows
        self._render_table(rows)

        # Update KPIs
        self.kpi_total_pairs.setText(str(len(rows)))
        v_energies = [r["vina_energy"] for r in rows if r["vina_energy"] is not None]
        a_energies = [r["ad4_energy"] for r in rows if r["ad4_energy"] is not None]
        ddg_vals = [r["ddg"] for r in rows if r["ddg"] is not None]

        self.kpi_best_vina.setText(f"{min(v_energies):.2f} kcal/mol" if v_energies else "—")
        self.kpi_best_ad4.setText(f"{min(a_energies):.2f} kcal/mol" if a_energies else "—")
        mean_ddg = sum(ddg_vals) / len(ddg_vals) if ddg_vals else None
        self.kpi_mean_delta.setText(f"{mean_ddg:+.2f} kcal/mol" if mean_ddg is not None else "—")

    def _render_table(self, rows: List[Dict[str, Any]]) -> None:
        """Render rows into QTableWidget."""
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(r["receptor"]))
            self.table.setItem(i, 1, QTableWidgetItem(r["ligand"]))
            self.table.setItem(i, 2, QTableWidgetItem(r["mode"]))

            v_item = QTableWidgetItem(f"{r['vina_energy']:.2f}" if r["vina_energy"] is not None else "—")
            v_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(i, 3, v_item)

            a_item = QTableWidgetItem(f"{r['ad4_energy']:.2f}" if r["ad4_energy"] is not None else "—")
            a_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(i, 4, a_item)

            ddg_item = QTableWidgetItem(f"{r['ddg']:+.2f}" if r["ddg"] is not None else "—")
            ddg_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(i, 5, ddg_item)

            self.table.setItem(i, 6, QTableWidgetItem(str(r["vina_ki"])))
            self.table.setItem(i, 7, QTableWidgetItem(str(r["ad4_ki"])))
            self.table.setItem(i, 8, QTableWidgetItem(str(r["ad4_cluster"])))
            self.table.setItem(i, 9, QTableWidgetItem(str(r["status"])))

    def _export_csv(self) -> None:
        """Export current filtered comparison rows to CSV."""
        if not self._filtered_rows:
            QMessageBox.information(self, "No Data", "No comparison rows available to export.")
            return

        f, _ = QFileDialog.getSaveFileName(
            self, "Export Comparison CSV",
            str(Path.cwd() / "docking_cross_comparison.csv"),
            "CSV Files (*.csv)"
        )
        if not f:
            return

        try:
            with open(f, "w", newline="", encoding="utf-8") as fp:
                writer = csv.DictWriter(fp, fieldnames=[
                    "receptor", "ligand", "mode", "vina_energy", "ad4_energy",
                    "ddg", "vina_ki", "ad4_ki", "ad4_cluster", "status"
                ])
                writer.writeheader()
                writer.writerows(self._filtered_rows)
            QMessageBox.information(self, "Export Complete", f"Comparison table successfully exported to:\n{f}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to write CSV: {e}")
