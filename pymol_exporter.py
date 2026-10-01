#!/usr/bin/env python3
"""
AutoDock Suite Pro — PyMOL 3D Session & Visualization Exporter
==============================================================
Generates publication-ready PyMOL sessions (.pse) and reproducible PyMOL
scripts (.pml) for docked receptor-ligand complexes.

Features:
  - Automatic detection of PyMOL installations (Windows, macOS, Linux, Conda)
  - Receptor cartoon with secondary-structure styling and optional transparent surface
  - Docked ligand poses with ball-and-stick styling and multi-pose state flipping
  - Active site / pocket residues highlighted as sticks with residue labels
  - Full non-covalent interaction profiling (Hydrogen bonds, Salt bridges,
    Pi-stacking, Halogen bonds) rendered as color-coded 3D dashed lines with
    measured Ångström distances
  - 3D Wireframe Grid Box visualization based on docking parameters
  - 1-Click PyMOL desktop launch from GUI or CLI
  - Headless compilation to native binary .pse session files
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from interactions import (
    Interaction, PDBQTAtom, analyze_pose_interactions,
    dist3d, parse_docked_poses, parse_receptor_pdbqt,
)

logger = logging.getLogger("docking_automation.pymol_exporter")


# ==============================================================================
# PYMOL INSTALLATION DISCOVERY
# ==============================================================================

COMMON_WINDOWS_PYMOL_PATHS = [
    # Schrodinger AppData (current user standard)
    Path(os.path.expandvars(r"%LOCALAPPDATA%\Schrodinger\PyMOL\Scripts\pymol.exe")),
    Path(os.path.expandvars(r"%LOCALAPPDATA%\Schrodinger\PyMOL\PyMOLWin.exe")),
    Path(os.path.expandvars(r"%LOCALAPPDATA%\Schrodinger\PyMOL\python.exe")),
    # Standard 64-bit Program Files
    Path(r"C:\Program Files\PyMOL\PyMOLWin.exe"),
    Path(r"C:\Program Files\PyMOL\PyMOL.exe"),
    Path(r"C:\Program Files\Schrodinger\PyMOL\PyMOL.exe"),
    Path(r"C:\Program Files\Schrodinger\PyMOL\PyMOLWin.exe"),
    # 32-bit Program Files
    Path(r"C:\Program Files (x86)\PyMOL\PyMOL.exe"),
    Path(r"C:\Program Files (x86)\PyMOL\PyMOLWin.exe"),
    # Direct Root Installations
    Path(r"C:\PyMOL\PyMOL.exe"),
    Path(r"C:\PyMOL\PyMOLWin.exe"),
    # Anaconda / Miniconda base & envs
    Path(os.path.expanduser(r"~\anaconda3\Scripts\pymol.exe")),
    Path(os.path.expanduser(r"~\miniconda3\Scripts\pymol.exe")),
    Path(os.path.expanduser(r"~\anaconda3\python.exe")),
    Path(os.path.expanduser(r"~\miniconda3\python.exe")),
]

COMMON_UNIX_PYMOL_PATHS = [
    Path("/usr/bin/pymol"),
    Path("/usr/local/bin/pymol"),
    Path("/opt/pymol/bin/pymol"),
    Path("/Applications/PyMOL.app/Contents/MacOS/PyMOL"),
    Path(os.path.expanduser("~/miniconda3/bin/pymol")),
    Path(os.path.expanduser("~/anaconda3/bin/pymol")),
]


def find_pymol_executable(custom_path: Optional[Path | str] = None) -> Optional[Path]:
    """Auto-detects the PyMOL GUI executable on the system.

    Searches:
      1. Explicit user-provided custom_path
      2. Environment variable PYMOL_PATH
      3. System PATH (via shutil.which)
      4. Known platform-specific installation directories

    Returns:
        Path to PyMOL executable if found, or None.
    """
    # 1. Custom path
    if custom_path:
        p = Path(custom_path).resolve()
        if p.is_file():
            return p
        if p.is_dir():
            for name in ("pymol.exe", "PyMOLWin.exe", "pymol"):
                candidate = p / name
                if candidate.is_file():
                    return candidate

    # 2. PYMOL_PATH environment variable
    env_path = os.environ.get("PYMOL_PATH")
    if env_path:
        p = Path(env_path).resolve()
        if p.is_file():
            return p
        if p.is_dir():
            for name in ("pymol.exe", "PyMOLWin.exe", "pymol"):
                candidate = p / name
                if candidate.is_file():
                    return candidate

    # 3. System PATH
    for name in ("pymol", "pymol.exe", "PyMOLWin.exe", "pymolwin"):
        found = shutil.which(name)
        if found:
            p = Path(found).resolve()
            if p.is_file():
                return p

    # 4. Standard directories
    candidate_list = COMMON_WINDOWS_PYMOL_PATHS if sys.platform == "win32" else COMMON_UNIX_PYMOL_PATHS
    for candidate in candidate_list:
        try:
            if candidate.is_file():
                return candidate.resolve()
        except OSError:
            continue

    return None


def find_pymol_python(custom_path: Optional[Path | str] = None) -> Optional[Path]:
    """Finds a Python interpreter that has the `pymol` module installed.

    Used for headless execution and programmatic .pse session generation.
    """
    # First check current interpreter
    try:
        import pymol  # type: ignore # noqa: F401
        return Path(sys.executable).resolve()
    except ImportError:
        pass

    # Next check next to detected pymol executable
    pymol_exe = find_pymol_executable(custom_path)
    if pymol_exe:
        parent = pymol_exe.parent
        # Candidate 1: parent/python.exe
        candidates = [
            parent / "python.exe",
            parent.parent / "python.exe",
            parent / "python",
            parent.parent / "python",
        ]
        for c in candidates:
            if c.is_file():
                return c.resolve()

    return None


def is_pymol_available(custom_path: Optional[Path | str] = None) -> bool:
    """Returns True if PyMOL can be executed or invoked programmatically."""
    return find_pymol_executable(custom_path) is not None or find_pymol_python(custom_path) is not None


# ==============================================================================
# PDB COORDINATE CONVERSION UTILITIES
# ==============================================================================

def pdbqt_to_clean_pdb_lines(atoms: List[PDBQTAtom], res_name_override: Optional[str] = None) -> List[str]:
    """Converts a list of PDBQTAtom instances into standard PDB format text lines."""
    lines: List[str] = []
    for a in atoms:
        record_type = "HETATM" if a.is_hetero else "ATOM  "
        res_name = (res_name_override or a.res_name)[:3].ljust(3)
        atom_name = a.name.strip()
        # Standard PDB atom name formatting: 4 chars
        if len(atom_name) < 4 and not atom_name[0].isdigit():
            fmt_atom_name = f" {atom_name:<3s}"
        else:
            fmt_atom_name = f"{atom_name:<4s}"

        chain = (a.chain or "A")[:1]
        res_seq = a.res_seq if a.res_seq is not None else 1
        elem = a.element[:2].rjust(2)

        line = (
            f"{record_type}{a.serial:5d} {fmt_atom_name} {res_name} {chain}{res_seq:4d}    "
            f"{a.x:8.3f}{a.y:8.3f}{a.z:8.3f}{1.00:6.2f}{20.00:6.2f}          {elem}"
        )
        lines.append(line)
    return lines


# ==============================================================================
# PYMOL SCRIPT (.PML) GENERATOR
# ==============================================================================

def generate_pymol_script(
    receptor_rel_or_abs: str,
    ligand_rel_or_abs: str,
    receptor_atoms: List[PDBQTAtom],
    poses: List[Tuple[int, List[PDBQTAtom], List[PDBQTAtom]]],
    active_pose_index: int = 1,
    interactions: Optional[List[Interaction]] = None,
    active_lig_atoms: Optional[List[PDBQTAtom]] = None,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    title: str = "AutoDock Suite Pro Docking Complex",
    session_pse_name: Optional[str] = None,
) -> str:
    """Generate a publication-quality PyMOL .pml script.

    Targets journal figures at Nature / JACS / JCIM standard:
      • White background, ray_trace_mode 1, antialias 2, proper ambient/specular
      • Secondary-structure coloring: helices=salmon, sheets=paleyellow, loops=white
      • Ligand ball-and-stick, gold carbons + full CPK heteroatom colors
      • Per-type interaction lines with calibrated widths, gaps, and radii:
          H-bond=yellow dashed, Salt bridge=magenta solid (dash_gap=0),
          Pi=cyan dashed, Halogen=orange dashed,
          Hydrophobic=slate dashes (label_size 0, disabled by default)
      • Binding pocket surface (wheat, 55 % transparent, disabled by default)
      • Compact Cα residue labels (resn+resi) within 4.5 Å of ligand
      • Wireframe grid box if grid_center/grid_size provided (steel blue)
      • 4 auto-rendered 300 DPI PNGs: front, side, back, top (1600×1200 px)
      • Clean grouped object hierarchy for publication figure editing
    """
    from datetime import datetime

    rec_posix = Path(receptor_rel_or_abs).as_posix()
    lig_posix = Path(ligand_rel_or_abs).as_posix()
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    stem = (
        Path(receptor_rel_or_abs).stem
        + "_"
        + Path(ligand_rel_or_abs).stem.replace("_poses", "")
    )

    rec_atom_map = {(a.chain or "A", a.res_seq, a.name): a for a in receptor_atoms}
    if active_lig_atoms is None and poses:
        selected_pose = [(m, l, f) for m, l, f in poses if m == active_pose_index]
        active_lig_atoms = selected_pose[0][1] if selected_pose else poses[0][1]

    num_poses = len(poses) if poses else 1
    active_lig_sel = f"ligand_pose_{active_pose_index}"

    pml: List[str] = []

    # ── Header ────────────────────────────────────────────────────────────────
    pml += [
        f"# {'=' * 76}",
        f"# {title}",
        f"# Generated by AutoDock Suite Pro  |  {now}",
        "# Publication-quality PyMOL session — journal figure ready",
        f"# {'=' * 76}",
        "",
    ]

    # ── Global initialisation ─────────────────────────────────────────────────
    pml += [
        "# ── Global Initialisation ──────────────────────────────────────────────────",
        "reinitialize",
        "",
        "# Background & viewport  (journal standard: white)",
        "bg_color white",
        "set orthoscopic,         0",
        "viewport 1600, 1200",
        "",
        "# ── Ray-trace quality (mode 1 = production ray-trace) ──────────────────────",
        "set ray_trace_mode,         1",
        "set ray_trace_gain,         0.15",
        "set ray_trace_color,        black",
        "set ray_shadows,            1",
        "set ray_shadow_decay_factor, 0.10",
        "set ray_shadow_decay_range,  1.8",
        "set ray_interior_color,     grey80",
        "set antialias,              2",
        "set ambient,                0.30",
        "set specular,               0.55",
        "set specular_intensity,     0.55",
        "set shininess,              50",
        "set reflect,                0.35",
        "set reflect_power,          1",
        "set direct,                 0.45",
        "set depth_cue,              1",
        "set fog_start,              0.45",
        "",
        "# ── Global bond/atom display ────────────────────────────────────────────────",
        "set stick_radius,           0.18",
        "set stick_ball,             on",
        "set stick_ball_ratio,       1.5",
        "set sphere_scale,           0.22",
        "set sphere_quality,         4",
        "set stick_quality,          15",
        "",
        "# ── Interaction-line defaults ────────────────────────────────────────────────",
        "set dash_gap,               0.25",
        "set dash_width,             3.5",
        "set dash_radius,            0.055",
        "set dash_length,            0.40",
        "set label_size,             13",
        "set label_font_id,          7",
        "set label_color,            black",
        "set label_outline_color,    white",
        "set label_position,         (0, 0.8, 0)",
        "",
    ]

    # ── Receptor ──────────────────────────────────────────────────────────────
    pml += [
        "# ── Receptor (Protein) ─────────────────────────────────────────────────────",
        f'load "{rec_posix}", receptor',
        "hide everything, receptor",
        "",
        "# Cartoon with secondary-structure spectrum coloring",
        "show cartoon, receptor",
        "color white,      receptor",
        "color paleyellow, receptor and ss s",
        "color salmon,     receptor and ss h",
        "set cartoon_fancy_helices,     1",
        "set cartoon_highlight_color,   grey90",
        "set cartoon_side_chain_helper, 1",
        "set cartoon_helix_radius,      1.0",
        "set cartoon_loop_radius,       0.35",
        "set cartoon_tube_radius,       0.40",
        "",
        "# Transparent molecular surface (disabled by default; enable for surface figure panels)",
        "create receptor_surface, receptor",
        "hide everything, receptor_surface",
        "show surface,    receptor_surface",
        "color grey90,    receptor_surface",
        "set transparency,    0.55, receptor_surface",
        "set surface_quality, 2",
        "set surface_type,    0",
        "disable receptor_surface",
        "",
    ]

    # ── Ligand poses ──────────────────────────────────────────────────────────
    pml += [
        "# ── Docked Ligand Poses ────────────────────────────────────────────────────",
        f'load "{lig_posix}", ligand_poses',
        "hide everything, ligand_poses",
        "",
        f"# Active pose: Pose {active_pose_index}  (GOLD carbons — stands out on white/salmon background)",
        f"create {active_lig_sel}, ligand_poses, {active_pose_index}, 1",
        f"show sticks,  {active_lig_sel}",
        f"show spheres, {active_lig_sel}",
        f"set sphere_scale, 0.22, {active_lig_sel}",
        f"color gold,       {active_lig_sel} and elem C",
        f"util.cnc          {active_lig_sel}",
        f"# Full CPK heteroatom coloring",
        f"color blue,     {active_lig_sel} and elem N",
        f"color red,      {active_lig_sel} and elem O",
        f"color yellow,   {active_lig_sel} and elem S",
        f"color orange,   {active_lig_sel} and elem P",
        f"color tv_green, {active_lig_sel} and elem F",
        "",
    ]

    if num_poses > 1:
        other_pose_names: List[str] = []
        for p_idx in range(1, num_poses + 1):
            if p_idx == active_pose_index:
                continue
            obj = f"ligand_pose_{p_idx}"
            pml += [
                f"create {obj}, ligand_poses, {p_idx}, 1",
                f"show sticks,  {obj}",
                f"color marine, {obj} and elem C",
                f"util.cnc {obj}",
                f"set stick_radius, 0.12, {obj}",
                f"disable {obj}",
            ]
            other_pose_names.append(obj)
        if other_pose_names:
            pml.append(f"group Alternative_Poses, {' '.join(other_pose_names)}")
        pml += ["disable ligand_poses", ""]

    # ── Pocket residues & surface ─────────────────────────────────────────────
    pml += [
        "# ── Binding Pocket Residues ─────────────────────────────────────────────────",
        f"select pocket_residues, byres (receptor within 4.5 of {active_lig_sel})",
        "show sticks, pocket_residues",
        "color grey80, pocket_residues and elem C",
        "util.cnc pocket_residues",
        "",
        "# Compact Cα residue labels (resn+resi e.g. ASP102)",
        "label pocket_residues and name CA, '%-3s%s' % (resn, resi)",
        "set label_size,         12",
        "set label_font_id,       7",
        "set label_color,        black",
        "set label_outline_color, white",
        "",
        "# Pocket volume surface (wheat, disabled by default)",
        f"select pocket_volume, byres (receptor within 6.0 of {active_lig_sel})",
        "create pocket_surface_obj, pocket_volume",
        "hide everything, pocket_surface_obj",
        "show surface,    pocket_surface_obj",
        "color wheat,     pocket_surface_obj",
        "set transparency,    0.55, pocket_surface_obj",
        "set surface_quality, 2",
        "disable pocket_surface_obj",
        "",
    ]

    # ── Non-covalent interactions ─────────────────────────────────────────────
    pml.append("# ── Non-Covalent Interactions ──────────────────────────────────────────────")
    hbond_objs:   List[str] = []
    salt_objs:    List[str] = []
    pi_objs:      List[str] = []
    halogen_objs: List[str] = []
    hydro_objs:   List[str] = []

    if interactions:
        for idx, item in enumerate(interactions, start=1):
            rec_atom  = item.receptor_atom
            rec_seq   = item.receptor_res_seq
            rec_chain = item.receptor_chain or "A"
            rec_resn  = item.receptor_res_name
            lig_atom  = item.ligand_atom
            dist      = item.distance_angstrom
            itype     = item.interaction_type

            sel_rec = (
                f"receptor and chain {rec_chain} and resi {rec_seq} and name {rec_atom}"
            )

            # Resolve exact ligand atom serial by proximity matching
            lig_serial = None
            if active_lig_atoms:
                rec_a = rec_atom_map.get((rec_chain, rec_seq, rec_atom))
                if rec_a:
                    best_la = min(
                        active_lig_atoms,
                        key=lambda la: abs(dist3d(rec_a.coords, la.coords) - dist),
                    )
                    if abs(dist3d(rec_a.coords, best_la.coords) - dist) < 0.2:
                        lig_serial = best_la.serial

            sel_lig = (
                f"{active_lig_sel} and id {lig_serial}"
                if lig_serial is not None
                else f"{active_lig_sel} and name {lig_atom}"
            )

            safe_resn = re.sub(r"[^A-Za-z0-9_]", "_", rec_resn)
            safe_lat  = re.sub(r"[^A-Za-z0-9_]", "_", lig_atom)

            if itype == "Hydrogen Bond":
                obj = f"hb_{idx}_{safe_resn}{rec_seq}_{safe_lat}"
                pml += [
                    f"distance {obj}, {sel_rec}, {sel_lig}",
                    f"color yellow,    {obj}",
                    f"set dash_width,  3.5,   {obj}",
                    f"set dash_gap,    0.20,  {obj}",
                    f"set dash_radius, 0.060, {obj}",
                ]
                hbond_objs.append(obj)

            elif itype == "Salt Bridge":
                obj = f"sb_{idx}_{safe_resn}{rec_seq}_{safe_lat}"
                pml += [
                    f"distance {obj}, {sel_rec}, {sel_lig}",
                    f"color magenta,   {obj}",
                    f"set dash_width,  4.0,   {obj}",
                    f"set dash_gap,    0.0,   {obj}",
                    f"set dash_radius, 0.075, {obj}",
                ]
                salt_objs.append(obj)

            elif "Pi" in itype or "pi" in itype:
                obj = f"pi_{idx}_{safe_resn}{rec_seq}_{safe_lat}"
                pml += [
                    f"distance {obj}, {sel_rec}, {sel_lig}",
                    f"color cyan,      {obj}",
                    f"set dash_width,  3.0,   {obj}",
                    f"set dash_gap,    0.30,  {obj}",
                    f"set dash_radius, 0.055, {obj}",
                ]
                pi_objs.append(obj)

            elif itype == "Halogen Bond":
                obj = f"xb_{idx}_{safe_resn}{rec_seq}_{safe_lat}"
                pml += [
                    f"distance {obj}, {sel_rec}, {sel_lig}",
                    f"color orange,    {obj}",
                    f"set dash_width,  3.5,   {obj}",
                    f"set dash_gap,    0.25,  {obj}",
                    f"set dash_radius, 0.060, {obj}",
                ]
                halogen_objs.append(obj)

            elif "Hydrophobic" in itype or "hydrophobic" in itype:
                obj = f"hp_{idx}_{safe_resn}{rec_seq}_{safe_lat}"
                pml += [
                    f"distance {obj}, {sel_rec}, {sel_lig}",
                    f"color slate,     {obj}",
                    f"set dash_width,  2.0,   {obj}",
                    f"set dash_gap,    0.40,  {obj}",
                    f"set dash_radius, 0.040, {obj}",
                    f"set label_size,  0,     {obj}",
                ]
                hydro_objs.append(obj)

        pml.append("")
        all_groups: List[str] = []
        if hbond_objs:
            pml.append(f"group Hydrogen_Bonds,      {' '.join(hbond_objs)}")
            all_groups.append("Hydrogen_Bonds")
        if salt_objs:
            pml.append(f"group Salt_Bridges,        {' '.join(salt_objs)}")
            all_groups.append("Salt_Bridges")
        if pi_objs:
            pml.append(f"group Pi_Interactions,     {' '.join(pi_objs)}")
            all_groups.append("Pi_Interactions")
        if halogen_objs:
            pml.append(f"group Halogen_Bonds,       {' '.join(halogen_objs)}")
            all_groups.append("Halogen_Bonds")
        if hydro_objs:
            pml.append(f"group Hydrophobic_Contacts, {' '.join(hydro_objs)}")
            pml.append("disable Hydrophobic_Contacts")
            all_groups.append("Hydrophobic_Contacts")
        if all_groups:
            pml.append(f"group Interactions,        {' '.join(all_groups)}")
        pml.append("")

    # ── 3D Wireframe Grid Box ─────────────────────────────────────────────────
    if grid_center and grid_size:
        cx, cy, cz = grid_center
        sx, sy, sz = grid_size
        mn_x, mx_x = cx - sx / 2, cx + sx / 2
        mn_y, mx_y = cy - sy / 2, cy + sy / 2
        mn_z, mx_z = cz - sz / 2, cz + sz / 2

        edges = [
            (mn_x, mn_y, mn_z, mx_x, mn_y, mn_z),
            (mx_x, mn_y, mn_z, mx_x, mx_y, mn_z),
            (mx_x, mx_y, mn_z, mn_x, mx_y, mn_z),
            (mn_x, mx_y, mn_z, mn_x, mn_y, mn_z),
            (mn_x, mn_y, mx_z, mx_x, mn_y, mx_z),
            (mx_x, mn_y, mx_z, mx_x, mx_y, mx_z),
            (mx_x, mx_y, mx_z, mn_x, mx_y, mx_z),
            (mn_x, mx_y, mx_z, mn_x, mn_y, mx_z),
            (mn_x, mn_y, mn_z, mn_x, mn_y, mx_z),
            (mx_x, mn_y, mn_z, mx_x, mn_y, mx_z),
            (mx_x, mx_y, mn_z, mx_x, mx_y, mx_z),
            (mn_x, mx_y, mn_z, mn_x, mx_y, mx_z),
        ]
        vertex_lines = "\n".join(
            f"    VERTEX, {x1:.3f}, {y1:.3f}, {z1:.3f},"
            f" VERTEX, {x2:.3f}, {y2:.3f}, {z2:.3f},"
            for x1, y1, z1, x2, y2, z2 in edges
        )

        pml += [
            "# ── Docking Search Space Grid Box ───────────────────────────────────────────",
            "python",
            "from pymol.cgo import BEGIN, LINES, COLOR, VERTEX, END",
            "from pymol import cmd",
            "grid_box = [",
            "    BEGIN, LINES,",
            "    COLOR, 0.10, 0.65, 0.95,   # steel blue",
            vertex_lines,
            "    END,",
            "]",
            "cmd.load_cgo(grid_box, 'Docking_Grid_Box')",
            "python end",
            "set cgo_line_width, 1.8, Docking_Grid_Box",
            "",
        ]

    # ── Camera & Focus ────────────────────────────────────────────────────────
    pml += [
        "# ── Camera Orientation — centred on binding pocket ──────────────────────────",
        f"center {active_lig_sel}",
        f"orient {active_lig_sel} or pocket_residues",
        f"zoom   {active_lig_sel} or pocket_residues, buffer=5.0",
        "",
    ]

    # ── Clean object hierarchy ────────────────────────────────────────────────
    pml += [
        "# ── Object Hierarchy ───────────────────────────────────────────────────────",
        "group Receptor, receptor receptor_surface",
        f"group Ligand,   {active_lig_sel} ligand_poses",
        "group Pocket,   pocket_residues pocket_surface_obj",
        "",
        "deselect",
        "set selection_width, 0",
        "",
    ]

    # ── Publication PNG renders (4 views, 300 DPI) ─────────────────────────────
    pml += [
        "# ── Auto Multi-Angle Publication PNG Renders (300 DPI, 1600 × 1200 px) ─────",
        "# Front view",
        "ray 1600, 1200",
        f'png "{stem}_front.png", dpi=300',
        "",
        "# 90 ° rotation — side view",
        "turn y, 90",
        "ray 1600, 1200",
        f'png "{stem}_side.png", dpi=300',
        "",
        "# 180 ° total — back view",
        "turn y, 90",
        "ray 1600, 1200",
        f'png "{stem}_back.png", dpi=300',
        "",
        "# Return to front",
        "turn y, -180",
        "",
        "# Top-down overhead view of binding pocket",
        "turn x, 90",
        f"zoom {active_lig_sel}, buffer=6.0",
        "ray 1600, 1200",
        f'png "{stem}_top.png", dpi=300',
        "turn x, -90",
        f"zoom {active_lig_sel} or pocket_residues, buffer=5.0",
        "",
    ]

    # ── Save PSE ──────────────────────────────────────────────────────────────
    if session_pse_name:
        pml += [
            "# ── Save Binary Session (.pse) ──────────────────────────────────────────────",
            f'save "{session_pse_name}"',
            "",
        ]

    # ── Figure tips ───────────────────────────────────────────────────────────
    pml += [
        "# ── Journal Figure Tips ─────────────────────────────────────────────────────",
        "# • Enable 'receptor_surface' or 'pocket_surface_obj' for surface figure panels.",
        "# • Disable 'Hydrophobic_Contacts' group to reduce line clutter.",
        "# • Use  set ray_trace_mode, 3  for a bold outline (black-on-white) look.",
        "# • Manually adjust zoom/orient then re-run:  ray 1600, 1200  +  png 'fig.png', dpi=300",
        "# • For dark-background figures: bg_color black  (change label_color to white first).",
        f"# • 4 PNG figures exported: {stem}_front / _side / _back / _top  (300 DPI each)",
    ]

    return "\n".join(pml)


# ==============================================================================
# PYMOL SESSION EXPORTER & COMPILER
# ==============================================================================

def export_pymol_session(
    receptor_path: Path | str,
    ligand_path: Path | str,
    output_dir: Optional[Path | str] = None,
    pose_index: int = 1,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    compile_pse: bool = True,
    custom_pymol_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    """Generates a complete PyMOL export bundle (.pml script, PDB files, and .pse session).

    Args:
        receptor_path: Path to receptor PDB or PDBQT.
        ligand_path: Path to docked ligand PDBQT.
        output_dir: Directory where files will be created (defaults to ./pymol_session).
        pose_index: Active pose number (1-based).
        grid_center: (cx, cy, cz) tuple.
        grid_size: (sx, sy, sz) tuple.
        compile_pse: If True, attempts to compile native .pse via local PyMOL.
        custom_pymol_path: Optional path to pymol executable.

    Returns:
        Dict with keys:
          - 'dir': Path to session folder
          - 'pml': Path to .pml script
          - 'pse': Path to .pse file (or None if compilation failed or skipped)
          - 'receptor_pdb': Path to clean receptor PDB
          - 'ligand_pdb': Path to clean ligand PDB
          - 'interactions_count': Total interactions identified
    """
    rec_p = Path(receptor_path).resolve()
    lig_p = Path(ligand_path).resolve()

    if not rec_p.is_file():
        raise FileNotFoundError(f"Receptor file not found: {rec_p}")
    if not lig_p.is_file():
        raise FileNotFoundError(f"Ligand file not found: {lig_p}")

    rec_stem = rec_p.stem
    lig_stem = lig_p.stem

    if output_dir:
        out_dir = Path(output_dir).resolve()
    else:
        out_dir = rec_p.parent / f"pymol_{rec_stem}_{lig_stem}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parse coordinates
    rec_atoms = parse_receptor_pdbqt(rec_p)
    poses = parse_docked_poses(lig_p)

    # 2. Export clean PDB files into output_dir for maximum PyMOL compatibility
    clean_rec_pdb = out_dir / f"{rec_stem}.pdb"
    rec_pdb_lines = pdbqt_to_clean_pdb_lines(rec_atoms)
    clean_rec_pdb.write_text("\n".join(rec_pdb_lines) + "\nEND\n", encoding="utf-8")

    clean_lig_pdb = out_dir / f"{lig_stem}_poses.pdb"
    all_pose_lines: List[str] = []
    for m_idx, lig_atoms, flex_atoms in poses:
        all_pose_lines.append(f"MODEL        {m_idx:4d}")
        all_pose_lines.extend(pdbqt_to_clean_pdb_lines(lig_atoms + flex_atoms, res_name_override="LIG"))
        all_pose_lines.append("ENDMDL")
    if not all_pose_lines:
        # Fallback if poses empty
        all_pose_lines = ["MODEL        1", "ENDMDL"]
    clean_lig_pdb.write_text("\n".join(all_pose_lines) + "\nEND\n", encoding="utf-8")

    # 3. Analyze interactions for requested pose
    selected_pose = [(m, l, f) for m, l, f in poses if m == pose_index]
    active_lig_atoms = selected_pose[0][1] if selected_pose else (poses[0][1] if poses else [])
    active_flex_atoms = selected_pose[0][2] if selected_pose else (poses[0][2] if poses else [])

    interactions = analyze_pose_interactions(
        receptor_atoms=rec_atoms,
        ligand_atoms=active_lig_atoms,
        pose_index=pose_index,
        flex_atoms=active_flex_atoms if active_flex_atoms else None,
    )

    # 4. Generate .pml script
    pse_file = out_dir / f"{rec_stem}_{lig_stem}_docking.pse"
    pml_file = out_dir / f"{rec_stem}_{lig_stem}_docking.pml"

    pml_text = generate_pymol_script(
        receptor_rel_or_abs=clean_rec_pdb.name,
        ligand_rel_or_abs=clean_lig_pdb.name,
        receptor_atoms=rec_atoms,
        poses=poses,
        active_pose_index=pose_index,
        interactions=interactions,
        active_lig_atoms=active_lig_atoms,
        grid_center=grid_center,
        grid_size=grid_size,
        title=f"{rec_stem} :: {lig_stem} (Pose {pose_index})",
        session_pse_name=pse_file.name,
    )
    pml_file.write_text(pml_text, encoding="utf-8")
    logger.info(f"Generated PyMOL script: {pml_file}")

    # 5. Compile .pse session file if requested and PyMOL is available
    compiled_pse: Optional[Path] = None
    if compile_pse:
        compiled_pse = compile_pymol_session(
            pml_path=pml_file,
            pse_output_path=pse_file,
            work_dir=out_dir,
            custom_pymol_path=custom_pymol_path,
        )

    return {
        "dir": out_dir,
        "pml": pml_file,
        "pse": compiled_pse,
        "receptor_pdb": clean_rec_pdb,
        "ligand_pdb": clean_lig_pdb,
        "interactions_count": len(interactions),
    }


def compile_pymol_session(
    pml_path: Path | str,
    pse_output_path: Path | str,
    work_dir: Optional[Path | str] = None,
    custom_pymol_path: Optional[Path | str] = None,
) -> Optional[Path]:
    """Compiles a .pml script into a native binary .pse session file headlessly."""
    pml = Path(pml_path).resolve()
    pse = Path(pse_output_path).resolve()
    cwd = Path(work_dir).resolve() if work_dir else pml.parent

    # Remove any stale .pse so the freshness check below is reliable across
    # repeated exports to the same output directory.
    if pse.is_file():
        try:
            pse.unlink()
        except OSError:
            pass

    # Strategy 1: Run PyMOL CLI in batch mode (-c -q) — fastest and most reliable
    pymol_exe = find_pymol_executable(custom_pymol_path)
    if pymol_exe:
        try:
            logger.info(f"Compiling .pse via PyMOL CLI: {pymol_exe}")
            res = subprocess.run(
                [str(pymol_exe), "-c", "-q", str(pml), "-d", f'save "{pse.as_posix()}"; quit'],
                cwd=str(cwd),
                capture_output=True,
                text=True,
                timeout=25,
            )
            if pse.is_file() and pse.stat().st_size > 0:
                logger.info(f"Compiled binary PyMOL session (.pse) via CLI: {pse} ({pse.stat().st_size} bytes)")
                return pse
        except Exception as e:
            logger.warning(f"PyMOL CLI compilation failed: {e}")

    # Strategy 2: Python interpreter with pymol module fallback
    pymol_python = find_pymol_python(custom_pymol_path)
    if pymol_python:
        try:
            logger.info(f"Compiling .pse via PyMOL Python interpreter: {pymol_python}")
            py_code = (
                "import os, sys\n"
                "os.environ['PYOPENGL_PLATFORM'] = 'osmesa'\n"
                "import pymol\n"
                "pymol.pymol_argv = ['pymol', '-cq']\n"
                "pymol.finish_launching()\n"
                "from pymol import cmd\n"
                f"cmd.load(r'{pml.as_posix()}')\n"
                f"cmd.save(r'{pse.as_posix()}')\n"
                "cmd.quit()\n"
            )
            res = subprocess.run(
                [str(pymol_python), "-c", py_code],
                cwd=str(cwd),
                capture_output=True,
                text=True,
                timeout=25,
            )
            if pse.is_file() and pse.stat().st_size > 0:
                logger.info(f"Compiled binary PyMOL session (.pse) via Python: {pse} ({pse.stat().st_size} bytes)")
                return pse
        except Exception as e:
            logger.warning(f"PyMOL Python compilation failed: {e}")

    logger.warning("Could not compile .pse session file; .pml script remains available.")
    return None


# ==============================================================================
# PYMOL INTERACTIVE GUI LAUNCHER
# ==============================================================================

def launch_pymol(
    target_file: Path | str,
    custom_pymol_path: Optional[Path | str] = None,
) -> bool:
    """Launches PyMOL desktop application with the specified session or script loaded.

    Args:
        target_file: Path to .pse or .pml file.
        custom_pymol_path: Optional path to PyMOL executable.

    Returns:
        True if PyMOL process was successfully launched, False otherwise.
    """
    target = Path(target_file).resolve()
    if not target.is_file():
        raise FileNotFoundError(f"File to open in PyMOL not found: {target}")

    pymol_exe = find_pymol_executable(custom_pymol_path)
    if not pymol_exe:
        logger.error("PyMOL executable not found on this system.")
        return False

    cwd = target.parent
    logger.info(f"Launching PyMOL ({pymol_exe}) with target: {target}")

    try:
        # Spawn detached process so GUI stays responsive
        if sys.platform == "win32":
            subprocess.Popen(
                [str(pymol_exe), str(target)],
                cwd=str(cwd),
                creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP,
                close_fds=True,
            )
        else:
            subprocess.Popen(
                [str(pymol_exe), str(target)],
                cwd=str(cwd),
                start_new_session=True,
                close_fds=True,
            )
        return True
    except Exception as e:
        logger.error(f"Failed to launch PyMOL: {e}")
        return False


# ==============================================================================
# CLI EXECUTION (STANDALONE TESTING)
# ==============================================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="AutoDock Suite Pro — PyMOL 3D Session & Visualization Exporter"
    )
    parser.add_argument("-r", "--receptor", required=True, help="Path to receptor PDB or PDBQT")
    parser.add_argument("-l", "--ligand", required=True, help="Path to docked ligand PDBQT")
    parser.add_argument("-o", "--output-dir", default=None, help="Output directory")
    parser.add_argument("-p", "--pose", type=int, default=1, help="Pose index (1-based)")
    parser.add_argument("--cx", type=float, default=None, help="Grid box center X")
    parser.add_argument("--cy", type=float, default=None, help="Grid box center Y")
    parser.add_argument("--cz", type=float, default=None, help="Grid box center Z")
    parser.add_argument("--sx", type=float, default=None, help="Grid box size X")
    parser.add_argument("--sy", type=float, default=None, help="Grid box size Y")
    parser.add_argument("--sz", type=float, default=None, help="Grid box size Z")
    parser.add_argument("--pymol-exe", default=None, help="Custom path to PyMOL executable")
    parser.add_argument("--launch", action="store_true", help="Launch PyMOL GUI after export")
    args = parser.parse_args()

    center = (args.cx, args.cy, args.cz) if args.cx is not None and args.cy is not None and args.cz is not None else None
    size = (args.sx, args.sy, args.sz) if args.sx is not None and args.sy is not None and args.sz is not None else None

    print(f"Exporting PyMOL session:\n  Receptor: {args.receptor}\n  Ligand:   {args.ligand}\n  Pose:     {args.pose}")
    res = export_pymol_session(
        receptor_path=args.receptor,
        ligand_path=args.ligand,
        output_dir=args.output_dir,
        pose_index=args.pose,
        grid_center=center,
        grid_size=size,
        custom_pymol_path=args.pymol_exe,
    )

    print("\nExport completed successfully:")
    print(f"  Directory:    {res['dir']}")
    print(f"  Script:       {res['pml']}")
    print(f"  Session PSE:  {res['pse']}")
    print(f"  Interactions: {res['interactions_count']}")

    if args.launch:
        target = res['pse'] or res['pml']
        print(f"\nLaunching PyMOL with {target} ...")
        launch_pymol(target, custom_pymol_path=args.pymol_exe)
