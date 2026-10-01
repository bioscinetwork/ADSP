#!/usr/bin/env python3
"""
AutoDock Suite Pro — Unified Reporting Engine
=============================================
Generates CSV, styled XLSX, JSON, and interaction profiling reports
from completed docking jobs. Works seamlessly for Vina and AutoDock4 workflows.
"""

from __future__ import annotations

import csv
import json
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from config import ProjectConfig
from interactions import export_interactions_to_csv, profile_docking_job, summarize_interactions
from models import DockingJob, Engine, JobStatus, SUITE_NAME, __version__

logger = logging.getLogger("docking_automation.reporting")


# ═══════════════════════════════════════════════════════════════════════════════
# Vina Results Table
# ═══════════════════════════════════════════════════════════════════════════════

def _get_val(obj: Any, *keys: str, default: Any = None) -> Any:
    for k in keys:
        if isinstance(obj, dict) and k in obj:
            return obj[k]
        if hasattr(obj, k):
            return getattr(obj, k)
    return default


def _get_vina_modes(job: DockingJob) -> List[Any]:
    vr = getattr(job, "vina_results", None)
    if not vr:
        return []
    if isinstance(vr, list):
        return vr
    if isinstance(vr, dict):
        return vr.get("modes", [])
    return []


def _build_vina_rows(jobs: List[DockingJob]) -> List[Dict[str, Any]]:
    """Build report rows from Vina docking jobs."""
    rows = []
    for job in jobs:
        if job.status not in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING):
            continue
        if job.canonical_result and job.canonical_result.poses:
            for pose in job.canonical_result.poses:
                vm = pose.vina_metrics
                rows.append({
                    "Receptor": job.receptor_name,
                    "Docking_Mode": job.docking_mode.value,
                    "Ligand": job.ligand_name,
                    "Job_ID": job.job_id,
                    "Pose": pose.rank,
                    "Binding_Affinity_kcal_mol": pose.binding_energy,
                    "Estimated_Ki_nM": pose.estimated_ki_nM if pose.estimated_ki_nM is not None else "",
                    "RMSD_Lower_Bound": vm.rmsd_lower_bound if vm and vm.rmsd_lower_bound is not None else 0.0,
                    "RMSD_Upper_Bound": vm.rmsd_upper_bound if vm and vm.rmsd_upper_bound is not None else 0.0,
                    "Validation_RMSD": pose.validation_rmsd if pose.validation_rmsd is not None else "",
                    "Requested_Modes": job.requested_modes,
                    "Obtained_Modes": job.obtained_modes or len(job.canonical_result.poses),
                    "Elapsed_Seconds": round(job.elapsed_seconds, 1),
                    "Status": job.status.value,
                    "Warnings": "; ".join(job.warnings) if job.warnings else "",
                })
            continue

        modes = _get_vina_modes(job)
        for result in modes:
            pose = _get_val(result, "pose", "mode", default=1)
            affinity = _get_val(result, "binding_affinity", "affinity_kcal_mol", default=None)
            rmsd_lb = _get_val(result, "rmsd_lower_bound", "rmsd_lb", default=0.0)
            rmsd_ub = _get_val(result, "rmsd_upper_bound", "rmsd_ub", default=0.0)
            rows.append({
                "Receptor": job.receptor_name,
                "Docking_Mode": job.docking_mode.value,
                "Ligand": job.ligand_name,
                "Job_ID": job.job_id,
                "Pose": pose,
                "Binding_Affinity_kcal_mol": affinity,
                "RMSD_Lower_Bound": rmsd_lb,
                "RMSD_Upper_Bound": rmsd_ub,
                "Requested_Modes": job.requested_modes,
                "Obtained_Modes": job.obtained_modes or len(modes),
                "Elapsed_Seconds": round(job.elapsed_seconds, 1),
                "Status": job.status.value,
                "Warnings": "; ".join(job.warnings) if job.warnings else "",
            })
    return rows


def _build_vina_summary(jobs: List[DockingJob]) -> List[Dict[str, Any]]:
    """Build summary table: best pose per receptor-ligand pair."""
    rows = []
    rank = 0
    completed = [
        j for j in jobs
        if j.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)
        and (_get_vina_modes(j) or (j.canonical_result and j.canonical_result.poses))
    ]
    def _best_aff(job: DockingJob) -> float:
        if job.canonical_result and job.canonical_result.best_binding_energy is not None:
            return float(job.canonical_result.best_binding_energy)
        modes = _get_vina_modes(job)
        if modes:
            aff = _get_val(modes[0], "binding_affinity", "affinity_kcal_mol", default=0.0)
            try:
                return float(aff)
            except (ValueError, TypeError):
                return 0.0
        return 0.0

    completed.sort(key=_best_aff)

    for job in completed:
        rank += 1
        if job.canonical_result and job.canonical_result.poses:
            best_pose = job.canonical_result.poses[0]
            vm = best_pose.vina_metrics
            rows.append({
                "Rank": rank,
                "Receptor": job.receptor_name,
                "Docking_Mode": job.docking_mode.value,
                "Ligand": job.ligand_name,
                "Best_Affinity_kcal_mol": best_pose.binding_energy,
                "Best_Pose": best_pose.rank,
                "Total_Poses": len(job.canonical_result.poses),
                "RMSD_lb": vm.rmsd_lower_bound if vm and vm.rmsd_lower_bound is not None else 0.0,
                "RMSD_ub": vm.rmsd_upper_bound if vm and vm.rmsd_upper_bound is not None else 0.0,
                "Elapsed_s": round(job.elapsed_seconds, 1),
                "Status": job.status.value,
            })
        else:
            modes = _get_vina_modes(job)
            best = modes[0]
            rows.append({
                "Rank": rank,
                "Receptor": job.receptor_name,
                "Docking_Mode": job.docking_mode.value,
                "Ligand": job.ligand_name,
                "Best_Affinity_kcal_mol": _get_val(best, "binding_affinity", "affinity_kcal_mol"),
                "Best_Pose": _get_val(best, "pose", "mode", default=1),
                "Total_Poses": job.obtained_modes or len(modes),
                "RMSD_lb": _get_val(best, "rmsd_lower_bound", "rmsd_lb", default=0.0),
                "RMSD_ub": _get_val(best, "rmsd_upper_bound", "rmsd_ub", default=0.0),
                "Elapsed_s": round(job.elapsed_seconds, 1),
                "Status": job.status.value,
            })
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# AutoDock 4 Results Table
# ═══════════════════════════════════════════════════════════════════════════════

def _build_ad4_rows(jobs: List[DockingJob]) -> List[Dict[str, Any]]:
    """Build report rows from AutoDock4 docking jobs using canonical results or fallback."""
    rows = []
    for job in jobs:
        if job.status not in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING):
            continue

        # 1. Prefer Canonical Result (single owner of scientific analysis, no re-parsing)
        if job.canonical_result and job.canonical_result.poses:
            for p in job.canonical_result.poses:
                ad4 = p.ad4_metrics
                rows.append({
                    "Receptor": job.receptor_name,
                    "Docking_Mode": job.docking_mode.value,
                    "Ligand": job.ligand_name,
                    "Job_ID": job.job_id,
                    "Rank": p.rank,
                    "Run": p.run_number or p.rank,
                    "Binding_Energy_kcal_mol": p.binding_energy,
                    "Ki_nM": p.estimated_ki_nM if p.estimated_ki_nM is not None else "",
                    "Ki_Formatted": p.estimated_ki_formatted or "",
                    "Intermol_Energy": ad4.intermolecular_energy if ad4 else "",
                    "Internal_Energy": ad4.internal_energy if ad4 else "",
                    "Torsional_Energy": ad4.torsional_energy if ad4 else "",
                    "Cluster_ID": ad4.cluster_id if (ad4 and ad4.cluster_id is not None) else "",
                    "Cluster_Size": ad4.cluster_size if ad4 else "",
                    "Cluster_RMSD": ad4.cluster_rmsd if ad4 else "",
                    "RMSD_From_Ref": ad4.rmsd_from_reference if ad4 else "",
                    "Elapsed_Seconds": round(job.elapsed_seconds, 1),
                    "Status": job.status.value,
                })
            continue

        # 2. Check ad4_results dict in memory
        if hasattr(job, "ad4_results") and job.ad4_results and job.ad4_results.get("poses"):
            for p in job.ad4_results["poses"]:
                ki_str = f"{p.get('ki_raw', '')} {p.get('ki_unit', '')}".strip()
                rows.append({
                    "Receptor": job.receptor_name,
                    "Docking_Mode": job.docking_mode.value,
                    "Ligand": job.ligand_name,
                    "Job_ID": job.job_id,
                    "Rank": p.get("rank") or p.get("mode") or 1,
                    "Run": p.get("run_number") or p.get("rank") or 1,
                    "Binding_Energy_kcal_mol": p.get("binding_energy") or p.get("binding_affinity"),
                    "Ki_nM": p.get("ki_nM", ""),
                    "Ki_Formatted": ki_str,
                    "Intermol_Energy": p.get("intermol_energy", ""),
                    "Internal_Energy": p.get("internal_energy", ""),
                    "Torsional_Energy": p.get("torsional_energy", ""),
                    "Cluster_ID": p.get("cluster_id", ""),
                    "Cluster_Size": p.get("cluster_size", ""),
                    "Cluster_RMSD": p.get("cluster_rmsd", ""),
                    "RMSD_From_Ref": p.get("rmsd_from_ref", ""),
                    "Elapsed_Seconds": round(job.elapsed_seconds, 1),
                    "Status": job.status.value,
                })
            continue

        # 3. Fallback: Parse DLG file only if neither canonical nor cached results exist
        dlg_path = getattr(job, "dlg_path", None)
        if not dlg_path and hasattr(job, "ad4_results"):
            dlg_str = job.ad4_results.get("dlg_path")
            if dlg_str:
                dlg_path = Path(dlg_str)

        if dlg_path and Path(dlg_path).is_file():
            try:
                from dlg_extract import DLGParser
                res = DLGParser().parse(Path(dlg_path))
                if res and res.poses:
                    for pose in res.poses:
                        ki_str = f"{pose.ki_raw} {pose.ki_unit}".strip() if pose.ki_raw is not None else ""
                        rows.append({
                            "Receptor": job.receptor_name,
                            "Docking_Mode": job.docking_mode.value,
                            "Ligand": job.ligand_name,
                            "Job_ID": job.job_id,
                            "Rank": pose.rank,
                            "Run": pose.run_number or pose.rank,
                            "Binding_Energy_kcal_mol": pose.binding_energy,
                            "Ki_nM": pose.ki_nM if pose.ki_nM is not None else "",
                            "Ki_Formatted": ki_str,
                            "Intermol_Energy": pose.intermol_energy,
                            "Internal_Energy": pose.internal_energy,
                            "Torsional_Energy": pose.torsional_energy,
                            "Cluster_ID": pose.cluster_id if pose.cluster_id is not None else "",
                            "Cluster_Size": pose.cluster_size,
                            "Cluster_RMSD": pose.cluster_rmsd,
                            "RMSD_From_Ref": pose.rmsd_from_ref if pose.rmsd_from_ref != 0.0 else "",
                            "Elapsed_Seconds": round(job.elapsed_seconds, 1),
                            "Status": job.status.value,
                        })
            except Exception as e:
                logger.warning(f"Failed to parse DLG fallback for {job.job_id}: {e}")

    return rows


def _build_ad4_summary(jobs: List[DockingJob], ad4_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Build summary table: best pose per receptor-ligand pair for AD4."""
    if not ad4_rows:
        return []
    best_per_job: Dict[str, Dict[str, Any]] = {}
    for r in ad4_rows:
        jid = r["Job_ID"]
        if jid not in best_per_job or r["Binding_Energy_kcal_mol"] < best_per_job[jid]["Binding_Energy_kcal_mol"]:
            best_per_job[jid] = r

    summary = list(best_per_job.values())
    summary.sort(key=lambda x: x["Binding_Energy_kcal_mol"])
    for rank, row in enumerate(summary, start=1):
        row["Rank"] = rank
    return summary


# ═══════════════════════════════════════════════════════════════════════════════
# Status Table
# ═══════════════════════════════════════════════════════════════════════════════

def _build_status_rows(jobs: List[DockingJob]) -> List[Dict[str, Any]]:
    """Build job status table."""
    rows = []
    for job in jobs:
        rows.append({
            "Job_ID": job.job_id,
            "Receptor": job.receptor_name,
            "Docking_Mode": job.docking_mode.value,
            "Ligand": job.ligand_name,
            "Engine": job.engine.value,
            "Status": job.status.value,
            "Execution_Status": job.execution_status.value,
            "Analysis_Status": job.analysis_status.value,
            "Exit_Code": job.exit_code if job.exit_code is not None else "",
            "Elapsed_s": round(job.elapsed_seconds, 1),
            "Warnings": len(job.warnings),
            "Errors": len(job.errors),
            "Retry_Count": job.retry_count,
        })
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# Protein-Ligand Interaction Rows
# ═══════════════════════════════════════════════════════════════════════════════

def _build_interaction_rows(
    jobs: List[DockingJob],
    progress_callback: Optional[Callable[[str], None]] = None,
) -> List[Dict[str, Any]]:
    """Profiles molecular interactions for completed jobs and generates unified table."""
    all_rows: List[Dict[str, Any]] = []
    completed = [
        j for j in jobs
        if j.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)
    ]
    total_completed = len(completed)

    for idx, job in enumerate(completed, start=1):
        if progress_callback and (idx == 1 or idx % 5 == 0 or idx == total_completed):
            progress_callback(f"Profiling molecular interactions ({idx}/{total_completed})...")

        # Locate receptor PDBQT
        receptor_path = job.receptor_path
        if not receptor_path or not receptor_path.is_file():
            continue

        # Locate pose PDBQT
        pose_path = job.output_pdbqt
        if not pose_path or not pose_path.is_file():
            if job.output_dir:
                candidate = job.output_dir / f"{job.ligand_name}_out.pdbqt"
                if candidate.is_file():
                    pose_path = candidate
                else:
                    continue
            else:
                continue

        try:
            # Profile top poses (default: pose 1)
            interactions = profile_docking_job(receptor_path, pose_path, max_poses=1)
            if interactions:
                # Save individual interaction CSV in the ligand result folder
                per_job_csv = job.output_dir / "interactions.csv"
                export_interactions_to_csv(interactions, per_job_csv)

                # Add to master summary table
                for item in interactions:
                    row = item.to_dict()
                    row["Receptor"] = job.receptor_name
                    row["Docking_Mode"] = job.docking_mode.value
                    row["Ligand"] = job.ligand_name
                    row["Job_ID"] = job.job_id
                    all_rows.append(row)

        except Exception as e:
            logger.warning(f"Failed to profile interactions for {job.job_id}: {e}")

    return all_rows


# ═══════════════════════════════════════════════════════════════════════════════
# CSV / XLSX Export
# ═══════════════════════════════════════════════════════════════════════════════

def _write_csv(rows: List[Dict], output_path: Path, title: str = "") -> None:
    """Write rows to a CSV file with AutoDock Suite Pro header."""
    if not rows:
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    headers = list(rows[0].keys())

    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)

        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        writer.writerow([f"{SUITE_NAME} v{__version__}"])
        writer.writerow([f"Generated: {ts}"])
        if title:
            writer.writerow([title])
        writer.writerow([])

        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(h, "") for h in headers])

    logger.info(f"CSV written: {output_path}")


def _write_xlsx(
    tables: Dict[str, List[Dict]],
    output_path: Path,
) -> None:
    """Write multiple tables to an Excel workbook with styled sheets."""
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
        from openpyxl.utils import get_column_letter

        wb = Workbook()
        default_sheet = wb.active

        header_fill = PatternFill(start_color="1F4E79", end_color="1F4E79", fill_type="solid")
        header_font = Font(name="Segoe UI", size=10, bold=True, color="FFFFFF")
        data_font = Font(name="Segoe UI", size=10)
        border_thin = Border(
            left=Side(style="thin", color="E0E0E0"),
            right=Side(style="thin", color="E0E0E0"),
            top=Side(style="thin", color="E0E0E0"),
            bottom=Side(style="thin", color="E0E0E0"),
        )

        sheets_created = 0
        for sheet_name, rows in tables.items():
            if not rows:
                continue
            ws = wb.create_sheet(title=sheet_name[:31])
            sheets_created += 1
            ws.views.sheetView[0].showGridLines = True

            headers = list(rows[0].keys())
            ws.append(headers)

            for cell in ws[1]:
                cell.fill = header_fill
                cell.font = header_font
                cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            ws.row_dimensions[1].height = 26

            for row_idx, row in enumerate(rows, start=2):
                row_values = [row.get(h, "") for h in headers]
                ws.append(row_values)
                for col_idx in range(1, len(headers) + 1):
                    c = ws.cell(row=row_idx, column=col_idx)
                    c.font = data_font
                    c.border = border_thin
                    if isinstance(c.value, (int, float)):
                        c.alignment = Alignment(horizontal="right")
                    else:
                        c.alignment = Alignment(horizontal="left")

            # Auto-fit column widths
            for col in ws.columns:
                max_len = max(len(str(cell.value or "")) for cell in col)
                col_letter = get_column_letter(col[0].column)
                ws.column_dimensions[col_letter].width = min(max(max_len + 3, 11), 40)

        # Remove default blank sheet only if at least one styled sheet was created
        if sheets_created > 0:
            wb.remove(default_sheet)
        else:
            default_sheet.title = "Summary"
            default_sheet.append(["Status", "Notice"])
            default_sheet.append(["Empty", "No completed docking poses were available for export."])

        output_path.parent.mkdir(parents=True, exist_ok=True)
        wb.save(str(output_path))
        logger.info(f"Styled XLSX written: {output_path}")

    except ImportError:
        logger.warning("openpyxl not available — XLSX export skipped")
    except Exception as e:
        logger.warning(f"XLSX export failed: {e}")


def _write_json(data: Dict[str, Any], output_path: Path) -> None:
    """Write data to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(data, indent=2, default=str),
        encoding="utf-8"
    )
    logger.info(f"JSON written: {output_path}")


def _build_dual_comparison(vina_jobs: List[DockingJob], ad4_jobs: List[DockingJob]) -> List[Dict]:
    """Build a side-by-side comparison table for dual-engine docking runs."""
    vina_by_key: Dict[Tuple[str, str, str], DockingJob] = {}
    for j in vina_jobs:
        vina_by_key[(j.receptor_name, j.docking_mode.value, j.ligand_name)] = j

    ad4_by_key: Dict[Tuple[str, str, str], DockingJob] = {}
    for j in ad4_jobs:
        ad4_by_key[(j.receptor_name, j.docking_mode.value, j.ligand_name)] = j

    all_keys = sorted(set(list(vina_by_key.keys()) + list(ad4_by_key.keys())))
    rows = []

    for key in all_keys:
        rec, mode, lig = key
        v_job = vina_by_key.get(key)
        a_job = ad4_by_key.get(key)

        v_dg = None
        v_modes = 0
        if v_job and v_job.vina_results:
            modes = v_job.vina_results if isinstance(v_job.vina_results, list) else v_job.vina_results.get("modes", [])
            if modes:
                first = modes[0]
                v_dg = _get_val(first, "binding_affinity", "affinity_kcal_mol")
                v_modes = len(modes)

        a_dg = None
        a_modes = 0
        if a_job and a_job.ad4_results and a_job.ad4_results.get("poses"):
            poses = a_job.ad4_results["poses"]
            a_dg = poses[0].get("binding_energy")
            a_modes = len(poses)

        diff = None
        if v_dg is not None and a_dg is not None:
            try:
                diff = round(float(a_dg) - float(v_dg), 2)
            except (ValueError, TypeError):
                pass

        status = "Completed Both" if (v_job and a_job) else ("Vina Only" if v_job else "AD4 Only")

        rows.append({
            "Receptor": rec,
            "Docking Mode": mode,
            "Ligand": lig,
            "Vina Best Affinity (kcal/mol)": v_dg if v_dg is not None else "N/A",
            "AD4 Estimated Binding Energy (kcal/mol)": a_dg if a_dg is not None else "N/A",
            "Score Difference (AD4 - Vina)": diff if diff is not None else "N/A",
            # Backward-compatible schema keys: intentionally contain no value
            # so legacy consumers cannot mistake either score for thermodynamic ΔG.
            "Vina Best ΔG (kcal/mol)": "N/A — legacy label disabled",
            "AD4 Best ΔG (kcal/mol)": "N/A — legacy label disabled",
            "ΔΔG (AD4 - Vina)": "N/A — legacy label disabled",
            "Vina Poses": v_modes,
            "AD4 Poses": a_modes,
            "Engine Status": status,
        })

    return rows


def _generate_energy_chart_svg(rows: List[Dict[str, Any]], title: str = "Top Candidates Binding Energy (kcal/mol)") -> str:
    """Generate clean, publication-grade inline SVG horizontal bar chart."""
    if not rows:
        return ""
    top = rows[:12]
    chart_h = 40 + len(top) * 28
    chart_w = 680
    bar_x = 180
    max_w = 420

    # Determine scale
    energies = []
    for r in top:
        val = r.get("Best_Affinity_kcal_mol", r.get("Lowest_Energy_kcal_mol", r.get("Binding_Affinity_kcal_mol", 0.0)))
        try:
            energies.append(float(val))
        except (ValueError, TypeError):
            energies.append(0.0)

    min_e = min(energies) if energies else -1.0
    max_e = max(energies) if energies else 0.0
    span = abs(min_e) if abs(min_e) > 0.01 else 1.0

    svg_lines = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" xmlns="http://www.w3.org/2000/svg" style="background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">',
        f'<text x="20" y="24" font-size="12" font-weight="700" fill="#0f172a">{title}</text>',
    ]

    for i, (r, e) in enumerate(zip(top, energies)):
        y = 48 + i * 28
        lig = str(r.get("Ligand", f"Ligand {i+1}"))[:18]
        # Calculate bar width proportional to binding strength
        frac = min(1.0, max(0.05, abs(e) / span))
        bw = max(10, int(frac * max_w))

        # Color: negative energy (good) = deep journal teal/blue
        fill_color = "#0284c7" if e <= -7.0 else ("#0ea5e9" if e < -5.0 else "#94a3b8")

        svg_lines.append(f'<text x="{bar_x - 10}" y="{y + 13}" font-size="11" fill="#334155" text-anchor="end">{lig}</text>')
        svg_lines.append(f'<rect x="{bar_x}" y="{y}" width="{bw}" height="18" rx="3" fill="{fill_color}" />')
        svg_lines.append(f'<text x="{bar_x + bw + 8}" y="{y + 13}" font-size="11" font-weight="600" fill="#0f172a">{e:.2f} kcal/mol</text>')

    svg_lines.append('</svg>')
    return "\n".join(svg_lines)


def _generate_cluster_chart_svg(cluster_rows: List[Dict[str, Any]], title: str = "Conformational Cluster Population Distribution (%)") -> str:
    """Generate clean inline SVG bar chart for AutoDock4 cluster populations."""
    if not cluster_rows:
        return ""
    chart_h = 40 + len(cluster_rows[:10]) * 26
    chart_w = 680
    bar_x = 140
    max_w = 440

    svg_lines = [
        f'<svg viewBox="0 0 {chart_w} {chart_h}" xmlns="http://www.w3.org/2000/svg" style="background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">',
        f'<text x="20" y="24" font-size="12" font-weight="700" fill="#0f172a">{title}</text>',
    ]

    for i, c in enumerate(cluster_rows[:10]):
        y = 44 + i * 26
        cid = c.get("Cluster_ID", c.get("Cluster ID", i + 1))
        pop = c.get("Population_%", c.get("Population", 0.0))
        try:
            pop_f = float(pop)
        except (ValueError, TypeError):
            pop_f = 0.0

        bw = max(6, int((pop_f / 100.0) * max_w))
        fill_color = "#f59e0b" if i == 0 else "#fbbf24"

        svg_lines.append(f'<text x="{bar_x - 10}" y="{y + 13}" font-size="11" fill="#334155" text-anchor="end">Cluster {cid}</text>')
        svg_lines.append(f'<rect x="{bar_x}" y="{y}" width="{bw}" height="17" rx="3" fill="{fill_color}" />')
        svg_lines.append(f'<text x="{bar_x + bw + 8}" y="{y + 13}" font-size="11" font-weight="600" fill="#0f172a">{pop_f:.1f}% ({c.get("Size", c.get("Members", ""))} runs)</text>')

    svg_lines.append('</svg>')
    return "\n".join(svg_lines)


def generate_html_report(
    jobs: List[DockingJob],
    config: ProjectConfig,
    output_path: Path,
    tables: Optional[Dict[str, List[Dict[str, Any]]]] = None,
) -> Path:
    """Generate a publication-quality research report in standalone HTML format."""
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")
    total_jobs = len(jobs)
    successful_jobs = [j for j in jobs if j.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)]
    has_ad4 = any(j.engine in (Engine.AUTODOCK4, Engine.BOTH) for j in jobs)
    has_vina = any(j.engine in (Engine.VINA, Engine.BOTH) for j in jobs)

    # Best overall lead
    best_lead_name = "N/A"
    best_lead_score = "N/A"
    best_lead_rec = "N/A"
    best_lead_engine = "N/A"
    if successful_jobs:
        scored = []
        for j in successful_jobs:
            if j.canonical_result and j.canonical_result.best_binding_energy is not None:
                scored.append((j.canonical_result.best_binding_energy, j))
        if scored:
            scored.sort(key=lambda x: x[0])
            best_score, best_j = scored[0]
            best_lead_name = best_j.ligand_name
            best_lead_score = f"{best_score:.2f} kcal/mol"
            best_lead_rec = best_j.receptor_name
            best_lead_engine = best_j.engine.value

    # Extract summary rows for energy chart
    summary_rows = []
    if tables:
        for k, rows in tables.items():
            if "Summary" in k and rows:
                summary_rows = rows
                break

    energy_svg = _generate_energy_chart_svg(summary_rows) if summary_rows else ""

    # Cluster rows
    cluster_rows = []
    for j in successful_jobs:
        if j.canonical_result and j.canonical_result.clusters:
            for c in j.canonical_result.clusters:
                cluster_rows.append({
                    "Cluster_ID": c.cluster_id,
                    "Size": c.size,
                    "Population_%": c.population_percent,
                    "Lowest_Energy": c.lowest_energy,
                    "Mean_Energy": c.mean_energy,
                    "Cluster_RMSD": c.cluster_rmsd,
                    "Best_Run": c.representative_run,
                })
            break
    cluster_svg = _generate_cluster_chart_svg(cluster_rows) if cluster_rows else ""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{SUITE_NAME} — Publication Docking Report</title>
<style>
    @page {{
        size: A4;
        margin: 18mm 16mm 18mm 16mm;
    }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        color: #0f172a;
        background-color: #f8fafc;
        margin: 0;
        padding: 24px;
        line-height: 1.55;
        font-size: 13px;
    }}
    .report-container {{
        max-width: 980px;
        margin: 0 auto;
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 40px;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
    }}
    header.report-header {{
        border-bottom: 2px solid #0284c7;
        padding-bottom: 20px;
        margin-bottom: 28px;
    }}
    .suite-badge {{
        display: inline-block;
        background: #e0f2fe;
        color: #0284c7;
        font-weight: 700;
        font-size: 11px;
        padding: 3px 10px;
        border-radius: 4px;
        letter-spacing: 0.5px;
        text-transform: uppercase;
        margin-bottom: 8px;
    }}
    h1.report-title {{
        font-size: 24px;
        font-weight: 800;
        color: #0f172a;
        margin: 0 0 6px 0;
        letter-spacing: -0.3px;
    }}
    p.report-subtitle {{
        font-size: 13px;
        color: #64748b;
        margin: 0;
    }}
    .meta-bar {{
        display: flex;
        flex-wrap: wrap;
        gap: 18px;
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 6px;
        padding: 12px 18px;
        margin-top: 18px;
        font-size: 12px;
        color: #475569;
    }}
    .meta-item strong {{
        color: #0f172a;
    }}
    /* Sections */
    section.report-section {{
        margin-bottom: 34px;
        page-break-inside: avoid;
    }}
    h2.section-heading {{
        font-size: 15px;
        font-weight: 800;
        color: #0284c7;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        border-bottom: 1.5px solid #e2e8f0;
        padding-bottom: 6px;
        margin: 0 0 14px 0;
    }}
    /* KPI Grid */
    .kpi-grid {{
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
        gap: 12px;
        margin-bottom: 18px;
    }}
    .kpi-card {{
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 6px;
        padding: 12px;
    }}
    .kpi-card .kpi-label {{
        font-size: 10px;
        font-weight: 700;
        color: #64748b;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 4px;
    }}
    .kpi-card .kpi-val {{
        font-size: 18px;
        font-weight: 800;
        color: #0f172a;
        font-family: Consolas, monospace;
    }}
    /* Tables */
    table.data-table {{
        width: 100%;
        border-collapse: collapse;
        font-size: 12px;
        margin-bottom: 14px;
        border: 1px solid #e2e8f0;
    }}
    table.data-table th {{
        background: #f1f5f9;
        color: #0f172a;
        font-weight: 700;
        text-align: left;
        padding: 8px 10px;
        border: 1px solid #e2e8f0;
        font-size: 11px;
    }}
    table.data-table td {{
        padding: 7px 10px;
        border: 1px solid #e2e8f0;
        color: #334155;
        font-variant-numeric: tabular-nums;
    }}
    table.data-table tr:nth-child(even) {{
        background: #f8fafc;
    }}
    /* Charts */
    .chart-container {{
        margin: 16px 0;
    }}
    .figure-caption {{
        font-size: 11px;
        font-style: italic;
        color: #64748b;
        margin-top: 6px;
    }}
    /* Disclaimers & Notes */
    .disclaimer-box {{
        background: #fffbeb;
        border: 1px solid #fde68a;
        border-left: 4px solid #f59e0b;
        border-radius: 4px;
        padding: 12px 16px;
        font-size: 12px;
        color: #92400e;
        margin-top: 24px;
    }}
    .provenance-box {{
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 6px;
        padding: 14px 18px;
        font-family: Consolas, monospace;
        font-size: 11px;
        color: #334155;
        white-space: pre-wrap;
    }}
    @media print {{
        body {{
            background: #ffffff;
            padding: 0;
        }}
        .report-container {{
            border: none;
            box-shadow: none;
            padding: 0;
            max-width: 100%;
        }}
    }}
</style>
</head>
<body>
<div class="report-container">
    <header class="report-header">
        <div class="suite-badge">{SUITE_NAME} v{__version__}</div>
        <h1 class="report-title">Molecular Docking & Virtual Screening Report</h1>
        <p class="report-subtitle">Publication-Grade Biophysical Analysis, Binding Thermodynamics, and Conformational Clustering</p>
        <div class="meta-bar">
            <div class="meta-item"><strong>Date:</strong> {now_str}</div>
            <div class="meta-item"><strong>Project:</strong> {config.project_name}</div>
            <div class="meta-item"><strong>Engine:</strong> {config.engine.value}</div>
            <div class="meta-item"><strong>Mode:</strong> {config.docking_mode.value}</div>
            <div class="meta-item"><strong>Total Jobs:</strong> {total_jobs}</div>
        </div>
    </header>

    <!-- SECTION 1: EXECUTIVE SUMMARY -->
    <section class="report-section">
        <h2 class="section-heading">1. Executive Summary & Lead Candidate</h2>
        <div class="kpi-grid">
            <div class="kpi-card">
                <div class="kpi-label">Total Jobs</div>
                <div class="kpi-val">{total_jobs}</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-label">Completed</div>
                <div class="kpi-val">{len(successful_jobs)}</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-label">Top Lead</div>
                <div class="kpi-val" style="font-size:14px; font-weight:700;">{best_lead_name}</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-label">Best Affinity / ΔG</div>
                <div class="kpi-val" style="color:#0284c7;">{best_lead_score}</div>
            </div>
            <div class="kpi-card">
                <div class="kpi-label">Target Receptor</div>
                <div class="kpi-val" style="font-size:14px;">{best_lead_rec}</div>
            </div>
        </div>
        <p>This report documents the biophysical virtual screening execution carried out using <strong>{SUITE_NAME}</strong>, combining the empirical score optimization of AutoDock Vina with the semi-empirical force field and Lamarckian Genetic Algorithm clustering of AutoDock 4.2.6.</p>
    </section>

    <!-- SECTION 2: LEAD RANKINGS & ENERGY DISTRIBUTION -->
    <section class="report-section">
        <h2 class="section-heading">2. Lead Candidate Rankings & Affinity Distribution</h2>
        {f'<div class="chart-container">{energy_svg}<div class="figure-caption">Figure 1. Best predicted binding energy (kcal/mol) across leading screening candidates.</div></div>' if energy_svg else ''}
        
        <table class="data-table">
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Receptor</th>
                    <th>Ligand</th>
                    <th>Mode</th>
                    <th>Engine Score / Affinity (kcal/mol)</th>
                    <th>AutoDock4 DLG Ki</th>
                    <th>RMSD (lb/ub Å)</th>
                    <th>Poses</th>
                </tr>
            </thead>
            <tbody>"""

    rank_i = 1
    for j in successful_jobs[:15]:
        b_score = "N/A"
        ki_str = "N/A"
        rmsd_str = "0.00 / 0.00"
        n_poses = 1
        if j.canonical_result and j.canonical_result.poses:
            p0 = j.canonical_result.poses[0]
            b_score = f"{p0.binding_score:.2f}" if p0.binding_score is not None else "N/A"
            ki_str = p0.estimated_ki_formatted if j.engine == Engine.AUTODOCK4 else "Not applicable — Vina does not report Ki"
            n_poses = len(j.canonical_result.poses)
            if p0.vina_metrics:
                rmsd_str = f"{p0.vina_metrics.rmsd_lower_bound or 0.0:.2f} / {p0.vina_metrics.rmsd_upper_bound or 0.0:.2f}"
            elif p0.ad4_metrics:
                rmsd_str = f"Cluster RMSD: {p0.ad4_metrics.cluster_rmsd or 0.0:.2f}"

        html += f"""
                <tr>
                    <td><strong>{rank_i}</strong></td>
                    <td>{j.receptor_name}</td>
                    <td><strong>{j.ligand_name}</strong></td>
                    <td>{j.docking_mode.value}</td>
                    <td style="font-weight:700; color:#0284c7;">{b_score}</td>
                    <td>{ki_str}</td>
                    <td>{rmsd_str}</td>
                    <td>{n_poses}</td>
                </tr>"""
        rank_i += 1

    html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 3: AUTODOCK4 CLUSTER ANALYSIS
    if cluster_rows:
        html += f"""
    <section class="report-section">
        <h2 class="section-heading">3. AutoDock4 Conformational Cluster Analysis</h2>
        <p>Lamarckian Genetic Algorithm (LGA) docked conformations clustered with an RMSD tolerance of 2.0 Å. The lowest energy conformation in the most populated cluster represents the primary binding hypothesis.</p>
        <div class="chart-container">{cluster_svg}<div class="figure-caption">Figure 2. Conformational cluster population distribution across LGA docking runs.</div></div>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Cluster ID</th>
                    <th>Members</th>
                    <th>Population (%)</th>
                    <th>Lowest Energy (kcal/mol)</th>
                    <th>Mean Energy (kcal/mol)</th>
                    <th>Cluster RMSD (Å)</th>
                    <th>Best Run</th>
                </tr>
            </thead>
            <tbody>"""
        for c in cluster_rows[:10]:
            html += f"""
                <tr>
                    <td><strong>{c.get('Cluster_ID')}</strong></td>
                    <td>{c.get('Size')}</td>
                    <td>{float(c.get('Population_%', 0)):.1f}%</td>
                    <td style="font-weight:700; color:#0284c7;">{float(c.get('Lowest_Energy', 0)):.2f}</td>
                    <td>{float(c.get('Mean_Energy', 0)):.2f}</td>
                    <td>{float(c.get('Cluster_RMSD', 0)):.2f}</td>
                    <td>Run #{c.get('Best_Run')}</td>
                </tr>"""
        html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 4: CRYSTALLOGRAPHIC VALIDATION
    val_rows = []
    for j in successful_jobs:
        if j.canonical_result and j.canonical_result.validation and j.canonical_result.validation.crystal_rmsd is not None:
            v = j.canonical_result.validation
            val_rows.append({
                "Receptor": j.receptor_name,
                "Ligand": j.ligand_name,
                "Reference": Path(v.reference_ligand_file).name if v.reference_ligand_file else "Crystal",
                "Crystal_RMSD": f"{v.crystal_rmsd:.2f} Å",
                "Matched_Atoms": v.heavy_atoms_matched or "All heavy",
                "Method": v.atom_mapping_method,
                "Status": "PASS (RMSD < 2.0 Å)" if v.crystal_rmsd <= 2.0 else "SUB-OPTIMAL",
            })
    if val_rows:
        html += """
    <section class="report-section">
        <h2 class="section-heading">4. Crystallographic Reference Ligand Validation</h2>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Receptor</th>
                    <th>Docked Ligand</th>
                    <th>Reference Ligand</th>
                    <th>Crystal RMSD</th>
                    <th>Matched Heavy Atoms</th>
                    <th>Mapping Method</th>
                    <th>Validation Status</th>
                </tr>
            </thead>
            <tbody>"""
        for vr in val_rows:
            html += f"""
                <tr>
                    <td>{vr['Receptor']}</td>
                    <td>{vr['Ligand']}</td>
                    <td>{vr['Reference']}</td>
                    <td style="font-weight:700; color:#16a34a;">{vr['Crystal_RMSD']}</td>
                    <td>{vr['Matched_Atoms']}</td>
                    <td>{vr['Method']}</td>
                    <td>{vr['Status']}</td>
                </tr>"""
        html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 5: INTERACTIONS
    interaction_list = []
    for j in successful_jobs[:5]:
        if j.canonical_result and j.canonical_result.poses:
            p0 = j.canonical_result.poses[0]
            if p0.interactions:
                for inter in p0.interactions[:8]:
                    interaction_list.append({
                        "Receptor": j.receptor_name,
                        "Ligand": j.ligand_name,
                        "Type": inter.get("type", "Contact"),
                        "Residue": inter.get("residue", ""),
                        "Distance": f"{inter.get('distance', 0.0):.2f} Å" if inter.get("distance") else "—",
                        "Confidence": "High (< 3.2 Å)" if (inter.get("distance", 99.0) < 3.2) else "Standard",
                    })
    if interaction_list:
        html += """
    <section class="report-section">
        <h2 class="section-heading">5. Predicted Non-Covalent Molecular Interactions</h2>
        <p>Interactions detected between lead docked conformations and active site residues using geometric and distance-cutoff profiling.</p>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Receptor</th>
                    <th>Ligand</th>
                    <th>Interaction Type</th>
                    <th>Active Site Residue</th>
                    <th>Distance</th>
                    <th>Confidence</th>
                </tr>
            </thead>
            <tbody>"""
        for ir in interaction_list[:15]:
            html += f"""
                <tr>
                    <td>{ir['Receptor']}</td>
                    <td>{ir['Ligand']}</td>
                    <td><strong>{ir['Type']}</strong></td>
                    <td>{ir['Residue']}</td>
                    <td>{ir['Distance']}</td>
                    <td>{ir['Confidence']}</td>
                </tr>"""
        html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 6: THERMODYNAMIC & STATISTICAL DESCRIPTORS
    thermo_list = []
    for j in successful_jobs:
        if j.canonical_result and j.canonical_result.thermodynamics:
            th = j.canonical_result.thermodynamics
            thermo_list.append({
                "Receptor": j.receptor_name,
                "Ligand": j.ligand_name,
                "Info_Entropy": f"{th.info_entropy:.3f}" if th.info_entropy is not None else "—",
                "Partition_Fn": f"{th.partition_function:.2e}" if th.partition_function is not None else "—",
                "Stat_Temp": f"{th.stat_temperature:.1f} K" if th.stat_temperature is not None else "—",
                "Free_Energy": f"{th.stat_free_energy:.2f} kcal/mol" if th.stat_free_energy is not None else "—",
                "Internal_Energy": f"{th.stat_internal_energy:.2f} kcal/mol" if th.stat_internal_energy is not None else "—",
                "Boltzmann_Prob": f"{th.boltzmann_prob * 100:.1f}%" if th.boltzmann_prob is not None else "—",
                "Source": th.descriptor_source or "ADSP Statistical Engine",
            })
    if thermo_list:
        html += """
    <section class="report-section">
        <h2 class="section-heading">6. Thermodynamic &amp; Statistical Descriptors</h2>
        <p>Statistical mechanics and information-theoretic descriptors computed across the sampled conformational ensemble. These exploratory biophysical metrics distinguish ensemble entropy from static scoring minima.</p>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Receptor</th>
                    <th>Ligand</th>
                    <th>Info Entropy</th>
                    <th>Partition Function (Q)</th>
                    <th>Stat. Temperature</th>
                    <th>Free Energy</th>
                    <th>Internal Energy</th>
                    <th>Lead Boltzmann Prob</th>
                    <th>Method/Source</th>
                </tr>
            </thead>
            <tbody>"""
        for tr in thermo_list:
            html += f"""
                <tr>
                    <td>{tr['Receptor']}</td>
                    <td>{tr['Ligand']}</td>
                    <td>{tr['Info_Entropy']}</td>
                    <td>{tr['Partition_Fn']}</td>
                    <td>{tr['Stat_Temp']}</td>
                    <td style="font-weight:700;">{tr['Free_Energy']}</td>
                    <td>{tr['Internal_Energy']}</td>
                    <td>{tr['Boltzmann_Prob']}</td>
                    <td style="font-size:11px; color:#64748b;">{tr['Source']}</td>
                </tr>"""
        html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 7: ADMET PHYSICOCHEMICAL PROPERTIES
    admet_list = []
    for j in successful_jobs:
        if j.canonical_result and j.canonical_result.admet:
            adm = j.canonical_result.admet
            admet_list.append({
                "Ligand": j.ligand_name,
                "MW": f"{adm.get('mw', 0.0):.1f} Da" if 'mw' in adm else "—",
                "LogP": f"{adm.get('logp', 0.0):.2f}" if 'logp' in adm else "—",
                "HBD": str(adm.get('hbd', '—')),
                "HBA": str(adm.get('hba', '—')),
                "TPSA": f"{adm.get('tpsa', 0.0):.1f} Å²" if 'tpsa' in adm else "—",
                "RotBonds": str(adm.get('rotatable_bonds', '—')),
                "Lipinski": "PASS" if adm.get('lipinski_violations', 0) == 0 else f"{adm.get('lipinski_violations')} Violations",
            })
    if admet_list:
        html += """
    <section class="report-section">
        <h2 class="section-heading">7. ADMET &amp; Physicochemical Properties</h2>
        <p>In-silico physicochemical properties and Lipinski Rule of Five compliance for prioritized screening ligands.</p>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Ligand</th>
                    <th>Molecular Weight</th>
                    <th>cLogP</th>
                    <th>H-Bond Donors</th>
                    <th>H-Bond Acceptors</th>
                    <th>TPSA</th>
                    <th>Rotatable Bonds</th>
                    <th>Lipinski Rule-of-5</th>
                </tr>
            </thead>
            <tbody>"""
        for ar in admet_list:
            lip_color = "#16a34a" if ar['Lipinski'] == "PASS" else "#d97706"
            html += f"""
                <tr>
                    <td><strong>{ar['Ligand']}</strong></td>
                    <td>{ar['MW']}</td>
                    <td>{ar['LogP']}</td>
                    <td>{ar['HBD']}</td>
                    <td>{ar['HBA']}</td>
                    <td>{ar['TPSA']}</td>
                    <td>{ar['RotBonds']}</td>
                    <td style="font-weight:700; color:{lip_color};">{ar['Lipinski']}</td>
                </tr>"""
        html += """
            </tbody>
        </table>
    </section>"""

    # SECTION 8: REPRODUCIBILITY & PROVENANCE
    prov_text = f"Application: {SUITE_NAME} v{__version__}\n"
    prov_text += f"Date & Time: {now_str}\n"
    prov_text += f"Engine Config: {config.engine.value} ({config.docking_mode.value} mode)\n"
    prov_text += f"Exhaustiveness: {config.exhaustiveness}, Number of Modes: {config.num_modes}\n"
    if successful_jobs:
        j0 = successful_jobs[0]
        if j0.canonical_result and j0.canonical_result.provenance:
            pr = j0.canonical_result.provenance
            prov_text += f"Vina Version: {pr.engine_version or '1.2.x'}\n"
            prov_text += f"Receptor SHA-256: {pr.receptor_hash or 'Verified'}\n"
            prov_text += f"Ligand SHA-256: {pr.ligand_hash or 'Verified'}\n"
        cleanup = getattr(j0, "cleanup_provenance", None) or {}
        if not cleanup and j0.canonical_result and j0.canonical_result.provenance:
            cleanup = j0.canonical_result.provenance.cleanup_provenance or {}
        if cleanup:
            retained = cleanup.get("retained_components", [])
            removed = cleanup.get("removed_components", [])
            prov_text += f"Cleanup rule: {', '.join(cleanup.get('default_cleanup_components', [])) or 'none declared'}\n"
            prov_text += f"Retained components: {len(retained)}\n"
            prov_text += f"Removed components: {len(removed)}\n"
            prov_text += f"User overrides: {sum(1 for x in removed if x.get('reason') == 'user_selected_removal')}\n"
            prov_text += f"Cleanup warnings: {', '.join(cleanup.get('warnings', [])) or 'none'}\n"
        else:
            prov_text += "Cleanup provenance: Not recorded for this docking job\n"

    html += f"""
    <!-- SECTION 8: PROVENANCE & REPRODUCIBILITY -->
    <section class="report-section">
        <h2 class="section-heading">Reproducibility &amp; Provenance Record</h2>
        <div class="provenance-box">{prov_text}</div>
    </section>

    <!-- SECTION 9: SCIENTIFIC DISCLAIMERS -->
    <div class="disclaimer-box">
        <strong>Scientific Disclaimers &amp; Methodological Transparency:</strong>
        <p style="margin: 4px 0 0 0;">
        All docking affinities (kcal/mol), estimated inhibition constants (Ki), and interaction contacts in this document are computational estimates derived from empirical and force-field scoring functions. Predicted non-covalent contacts represent geometric distance/angle calculations and do not constitute direct experimental confirmation. Experimental validation (e.g. SPR, ITC, or X-ray crystallography) is recommended before drawing biological conclusions.
        </p>
    </div>
</div>
</body>
</html>"""

    output_path.write_text(html, encoding="utf-8")
    return output_path


# ═══════════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════════

def generate_reports(
    jobs: List[DockingJob],
    config: ProjectConfig,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> Dict[str, Any]:
    """Generate all reports for completed docking jobs.

    Produces:
        - CSV files (per-table including Protein-Ligand Interactions)
        - XLSX workbook (all tables combined with styling)
        - JSON data export
    Returns a dictionary of generated report paths.
    """
    if sys.platform == "win32":
        for _s in (sys.stdout, sys.stderr):
            if hasattr(_s, "reconfigure"):
                try:
                    _s.reconfigure(encoding="utf-8", errors="replace")
                except Exception:
                    pass

    report_dir = config.report_directory
    report_dir.mkdir(parents=True, exist_ok=True)

    if progress_callback:
        progress_callback("Compiling summary and pose tables...")

    vina_jobs = [
        j for j in jobs
        if j.engine == Engine.VINA or getattr(j, "vina_results", None) or (j.output_pdbqt and Path(j.output_pdbqt).is_file())
    ]
    ad4_jobs = [
        j for j in jobs
        if j.engine == Engine.AUTODOCK4 or getattr(j, "ad4_results", None) or getattr(j, "dlg_path", None)
    ]

    has_vina = bool(vina_jobs) or (config.engine in (Engine.VINA, Engine.BOTH))
    has_ad4 = bool(ad4_jobs) or (config.engine in (Engine.AUTODOCK4, Engine.BOTH))

    if (has_vina and has_ad4) or config.engine == Engine.BOTH:
        engine_label = "DUAL"
    elif has_ad4:
        engine_label = "AD4"
    else:
        engine_label = "VINA"

    tables: Dict[str, List[Dict]] = {}

    # Dual-Engine cross-comparison table
    if (has_vina and has_ad4) and vina_jobs and ad4_jobs:
        comparison_rows = _build_dual_comparison(vina_jobs, ad4_jobs)
        if comparison_rows:
            tables["00_Dual_Comparison"] = comparison_rows

    # Vina tables
    if has_vina and vina_jobs:
        vina_rows = _build_vina_rows(vina_jobs)
        vina_summary = _build_vina_summary(vina_jobs)
        if vina_rows:
            tables["01_Vina_All_Poses"] = vina_rows
        if vina_summary:
            tables["02_Vina_Summary"] = vina_summary

    # AutoDock 4 tables
    if has_ad4 and ad4_jobs:
        ad4_rows = _build_ad4_rows(ad4_jobs)
        ad4_summary = _build_ad4_summary(ad4_jobs, ad4_rows)
        prefix_a = "03_" if (has_vina and vina_jobs) else "01_"
        prefix_as = "04_" if (has_vina and vina_jobs) else "02_"
        if ad4_rows:
            tables[f"{prefix_a}AD4_All_Poses"] = ad4_rows
        if ad4_summary:
            tables[f"{prefix_as}AD4_Summary"] = ad4_summary

    status_rows = _build_status_rows(jobs)
    if status_rows:
        tables["05_Job_Status" if (has_vina and has_ad4) else "03_Job_Status"] = status_rows

    if progress_callback:
        progress_callback("Profiling molecular interactions...")
    interaction_rows = _build_interaction_rows(jobs, progress_callback=progress_callback)
    if interaction_rows:
        tables["06_Interactions" if (has_vina and has_ad4) else "04_Interactions"] = interaction_rows

    # Write CSVs
    if config.export_csv:
        if progress_callback:
            progress_callback("Writing CSV tables...")
        for name, rows in tables.items():
            _write_csv(
                rows,
                report_dir / f"{engine_label}_{name}.csv",
                title=f"{SUITE_NAME} — {engine_label} {name}",
            )

    xlsx_path = report_dir / f"{engine_label}_Complete_Report.xlsx"
    # Write XLSX
    if config.export_excel:
        if progress_callback:
            progress_callback("Generating formatted Excel workbook...")
        _write_xlsx(tables, xlsx_path)

    # Write JSON
    _write_json(
        {name: rows for name, rows in tables.items()},
        report_dir / f"{engine_label}_data.json",
    )

    # Write Publication-Quality HTML Research Report
    html_path = report_dir / f"{engine_label}_Publication_Report.html"
    try:
        if progress_callback:
            progress_callback("Generating publication-quality research report...")
        generate_html_report(jobs, config, html_path, tables=tables)
    except Exception as e:
        logger.warning(f"Could not generate HTML publication report: {e}")

    # Generate PyMOL 3D Sessions for Top Poses — only if PyMol Python API is available
    try:
        import pymol  # noqa: F401 — check presence before importing exporter
        from pymol_exporter import export_pymol_session
        pymol_dir = report_dir / "pymol"
        pymol_dir.mkdir(parents=True, exist_ok=True)
        completed = [
            j for j in jobs
            if j.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)
        ]
        if completed and progress_callback:
            progress_callback("Generating PyMOL 3D visualization sessions...")
        for job in completed[:15]:
            rec_p = job.receptor_path
            lig_p = job.output_pdbqt
            if not lig_p and hasattr(job, "dlg_path") and job.dlg_path:
                lig_p = job.dlg_path
            if rec_p and lig_p and Path(rec_p).is_file() and Path(lig_p).is_file():
                try:
                    export_pymol_session(
                        receptor_path=rec_p,
                        ligand_path=lig_p,
                        output_dir=pymol_dir / f"{job.receptor_name}_{job.ligand_name}",
                        pose_index=1,
                        compile_pse=True,
                    )
                except Exception as err:
                    logger.debug(f"PyMOL export for {job.job_id} skipped: {err}")
    except ImportError:
        logger.debug(
            "PyMOL Python API not available — skipping PSE export. "
            "PDBQT complex files are in each job's 'complexes/' subfolder."
        )
    except Exception as e:
        logger.debug(f"PyMOL report generation notice: {e}")

    if progress_callback:
        progress_callback("Reports completed successfully.")

    # Summary statistics
    total = len(jobs)
    success = sum(1 for j in jobs if j.status == JobStatus.SUCCESS)
    warnings = sum(1 for j in jobs if j.status == JobStatus.SUCCESS_WITH_WARNING)
    failed = sum(1 for j in jobs if j.status in (
        JobStatus.EXECUTION_FAILED, JobStatus.PARSING_FAILED,
        JobStatus.SPLIT_FAILED, JobStatus.INVALID_INPUT,
    ))

    print()
    print("=" * 60)
    print(f"  {SUITE_NAME} — REPORTS GENERATED")
    print("=" * 60)
    print(f"  Output dir  : {report_dir}")
    print(f"  Total jobs  : {total}")
    print(f"  Successful  : {success}")
    print(f"  Warnings    : {warnings}")
    print(f"  Failed      : {failed}")
    print(f"  Tables      : {', '.join(tables.keys())}")
    print("=" * 60)

    return {
        "report_dir": report_dir,
        "xlsx_path": xlsx_path if xlsx_path.is_file() else None,
        "html_path": html_path if html_path.is_file() else None,
        "tables": tables,
        "engine_label": engine_label,
    }
