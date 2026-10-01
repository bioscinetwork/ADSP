"""
interaction_diagram.py
======================
Publication-grade, LigPlot+-inspired 2D protein–ligand interaction diagram
rendered as a fully self-contained offline HTML+SVG document.

Architecture
------------
The diagram is built in five passes:
  1. Molecule loading  – RDKit or OpenBabel fallback
  2. Residue grouping  – group all contacts per residue, pick dominant type
  3. Orbital layout    – evenly distribute residues on an ellipse, then
                         iteratively push overlapping labels apart (greedy
                         angular collision avoidance)
  4. SVG build         – lines first, badges second (painter's algorithm)
  5. HTML wrap         – responsive, 100 % offline, dark/light themes

Interaction visual vocabulary
------------------------------
  Hydrogen Bond         – animated green dashed line, distance badge, atom labels
  Hydrophobic           – classic LigPlot+ eyelash arc (brick-red/amber)
  Salt Bridge           – magenta thick solid line, +/– charge symbols
  Pi-Stacking           – purple dashed line, ⬡ ring symbols on residue badge
  Pi-Cation             – cyan dashed line, ⬡ on ligand side
  Halogen Bond          – orange dashed line, X marker
  (fallback)            – slate dashed line

Public API
----------
  render_interaction_diagram_html(ligand_path, interactions, obabel_exe,
                                   width, height, theme) -> Optional[str]
  get_interaction_color(interaction_type, theme) -> str
  is_available() -> bool
"""

from __future__ import annotations

import html
import math
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class InteractionDiagramResult:
    """Standardized result containing true SVG, wrapped HTML, pose ID and metadata."""
    svg: str
    html: str
    pose_index: int
    metadata: Dict[str, Any] = field(default_factory=dict)
    contact_count: int = 0
    residue_count: int = 0
    interaction_classes: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Interaction-type priority, colors, and dash patterns
# ---------------------------------------------------------------------------

#: Priority list – first match wins when a residue has multiple contact types
_INTERACTION_PRIORITY: List[str] = [
    "Hydrogen Bond",
    "Salt Bridge",
    "Halogen Bond",
    "Pi-Stacking (Face-to-face)",
    "Pi-Stacking (T-shaped)",
    "Pi-Cation",
    "Hydrophobic",          # canonical short name (also matches "Hydrophobic Contact")
    "Hydrophobic Contact",
]

# Light-theme palette (journal / LigPlot Paper style)
_COLOUR_MAP: Dict[str, str] = {
    "hydrogen bond":              "#16a34a",   # Classic LigPlot green
    "salt bridge":                "#be185d",   # Deep magenta
    "halogen bond":               "#c2410c",   # Burnt orange
    "pi-stacking (face-to-face)": "#6d28d9",   # Deep violet
    "pi-stacking (t-shaped)":     "#6d28d9",
    "pi-cation":                  "#0369a1",   # Ocean blue
    "hydrophobic":                "#b91c1c",   # LigPlot brick red
    "hydrophobic contact":        "#b91c1c",
}

# Dark-theme palette (fluorescent / neon)
_COLOUR_MAP_DARK: Dict[str, str] = {
    "hydrogen bond":              "#22c55e",   # Neon emerald
    "salt bridge":                "#f43f5e",   # Hot rose / magenta
    "halogen bond":               "#fb923c",   # Bright amber-orange
    "pi-stacking (face-to-face)": "#a855f7",   # Electric violet
    "pi-stacking (t-shaped)":     "#a855f7",
    "pi-cation":                  "#38bdf8",   # Sky cyan
    "hydrophobic":                "#f59e0b",   # Glowing amber
    "hydrophobic contact":        "#f59e0b",
}

# Stroke-dasharray values per interaction type (empty = solid)
_DASH_MAP: Dict[str, str] = {
    "hydrogen bond":              "7,4",
    "salt bridge":                "",
    "halogen bond":               "5,3",
    "pi-stacking (face-to-face)": "5,3",
    "pi-stacking (t-shaped)":     "5,3",
    "pi-cation":                  "4,3",
    "hydrophobic":                "",
    "hydrophobic contact":        "",
}


def get_interaction_color(interaction_type: str, theme: str = "dark") -> str:
    """Return the canonical hex color for an interaction type label."""
    key = str(interaction_type).strip().lower()
    palette = _COLOUR_MAP_DARK if theme.lower() == "dark" else _COLOUR_MAP
    return palette.get(key, "#94a3b8")


def is_available() -> bool:
    """Return True when RDKit is installed and importable."""
    try:
        from rdkit import Chem  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Molecule loading helpers
# ---------------------------------------------------------------------------

def _load_mol_rdkit(ligand_path: Path):
    """Load a ligand file into an RDKit Mol with 2-D coordinates."""
def load_pose_molecule(
    ligand_path: str | Path,
    pose_index: int = 1,
    run_number: Optional[int] = None,
    original_ligand_path: Optional[str | Path] = None,
    obabel_exe: Optional[str | Path] = None,
) -> Optional[Any]:
    """
    Loads an RDKit Mol object corresponding strictly to the specified pose.

    1. If original_ligand_path is provided (SDF/MOL/MOL2/SMILES), loads that
       for high-fidelity chemical bond layout and aromaticity.
    2. Otherwise extracts the exact single-pose PDBQT coordinate block for pose_index.
    3. Converts to RDKit Mol and computes 2D coordinates.
    Never loads Pose 1 when Pose 2, 3, etc. is requested.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None

    # Step 1: Chemically rich original structure if provided
    if original_ligand_path:
        orig_p = Path(original_ligand_path)
        if orig_p.is_file():
            suf = orig_p.suffix.lower()
            mol = None
            try:
                if suf in (".sdf", ".mol"):
                    mol = Chem.MolFromMolFile(str(orig_p), sanitize=True, removeHs=True)
                elif suf == ".mol2":
                    mol = Chem.MolFromMol2File(str(orig_p), sanitize=True, removeHs=True)
                elif suf in (".smi", ".smiles"):
                    smi_text = orig_p.read_text(encoding="utf-8", errors="replace").strip().split()[0]
                    mol = Chem.MolFromSmiles(smi_text)
                if mol is not None:
                    AllChem.Compute2DCoords(mol)
                    return mol
            except Exception:
                pass

    p = Path(ligand_path)
    if not p.is_file():
        return None

    suffix = p.suffix.lower()

    # Step 2: If ligand_path is an SDF/MOL/MOL2, read it directly
    if suffix in (".sdf", ".mol"):
        try:
            mol = Chem.MolFromMolFile(str(p), sanitize=True, removeHs=True)
            if mol is not None:
                AllChem.Compute2DCoords(mol)
                return mol
        except Exception:
            pass
    elif suffix == ".mol2":
        try:
            mol = Chem.MolFromMol2File(str(p), sanitize=True, removeHs=True)
            if mol is not None:
                AllChem.Compute2DCoords(mol)
                return mol
        except Exception:
            pass

    # Step 3: Extract exact single-pose PDBQT block for pose_index
    from interactions import extract_pose_pdbqt_block
    pose_block = extract_pose_pdbqt_block(p, pose_index=pose_index, run_number=run_number)
    if pose_block is None:
        try:
            pose_block = p.read_text(encoding="utf-8", errors="replace")
        except Exception:
            return None

    # Check for REMARK SMILES inside this specific pose block
    for line in pose_block.splitlines():
        m = re.match(r"REMARK\s+SMILES\s+(\S+)", line, re.IGNORECASE)
        if m:
            try:
                mol = Chem.MolFromSmiles(m.group(1).strip())
                if mol is not None:
                    AllChem.Compute2DCoords(mol)
                    return mol
            except Exception:
                pass

    # Parse ATOM / HETATM lines from this single pose block as a clean PDB block
    pdb_lines: List[str] = []
    for line in pose_block.splitlines():
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            pdb_lines.append(line[:66].ljust(66) + "\n")
    if pdb_lines:
        pdb_block = "".join(pdb_lines) + "END\n"
        try:
            mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    pass
                AllChem.Compute2DCoords(mol)
                return mol
        except Exception:
            pass

    # Step 4: Fallback to OpenBabel for this specific single-pose block if provided
    if obabel_exe:
        try:
            res = subprocess.run(
                [str(obabel_exe), "-ipdbqt", "-osdf", "--gen2D"],
                input=pose_block.encode("utf-8"),
                capture_output=True, timeout=5,
            )
            sdf_block = res.stdout.decode(errors="replace")
            if sdf_block.strip():
                mol = Chem.MolFromMolBlock(sdf_block, sanitize=True, removeHs=True)
                if mol is not None:
                    AllChem.Compute2DCoords(mol)
                    return mol
        except Exception:
            pass

    return None


def _load_mol_rdkit(
    ligand_path: Path,
    pose_index: int = 1,
    run_number: Optional[int] = None,
    original_ligand_path: Optional[Path] = None,
    obabel_exe: Optional[str] = None,
):
    """Backwards-compatible wrapper delegating to load_pose_molecule."""
    return load_pose_molecule(
        ligand_path,
        pose_index=pose_index,
        run_number=run_number,
        original_ligand_path=original_ligand_path,
        obabel_exe=obabel_exe,
    )


# ---------------------------------------------------------------------------
# RDKit SVG depiction
# ---------------------------------------------------------------------------

def _rdkit_mol_svg(mol, lig_w: int, lig_h: int, theme: str) -> Optional[str]:
    """Render molecule with RDKit; returns raw SVG string or None."""
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(lig_w, lig_h)
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = True
        opts.bondLineWidth = 2.2
        opts.clearBackground = False
        # Atom-color overrides for dark theme readability
        if theme.lower() == "dark":
            try:
                from rdkit.Chem.Draw.rdMolDraw2D import SetDarkMode
                SetDarkMode(opts)
            except Exception:
                pass
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return None


def _extract_svg_inner(svg_doc: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract viewBox attribute and inner body from an RDKit SVG string."""
    svg_doc = re.sub(r"<\?xml[^?]*\?>", "", svg_doc).strip()
    vb_m = re.search(r'<svg[^>]*\bviewBox=["\']([\d\s.]+)["\']', svg_doc)
    viewbox = vb_m.group(1).strip() if vb_m else None
    inner_m = re.search(r"<svg[^>]*>(.*)</svg>", svg_doc, re.DOTALL)
    inner = inner_m.group(1).strip() if inner_m else None
    return viewbox, inner


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _ellipse_point(cx: float, cy: float, rx: float, ry: float,
                   angle_deg: float) -> Tuple[float, float]:
    """Point on an axis-aligned ellipse at *angle_deg* (0 = top, clockwise)."""
    rad = math.radians(angle_deg - 90.0)
    return cx + rx * math.cos(rad), cy + ry * math.sin(rad)


def _rect_edge_pt(cx: float, cy: float, hw: float, hh: float,
                  tx: float, ty: float) -> Tuple[float, float]:
    """
    Intersection of the ray from (cx, cy) toward (tx, ty) with the rectangle
    of half-width hw, half-height hh centered at (cx, cy).
    """
    dx, dy = tx - cx, ty - cy
    if dx == 0 and dy == 0:
        return cx, cy
    sx = abs(hw / dx) if dx != 0 else 1e9
    sy = abs(hh / dy) if dy != 0 else 1e9
    s = min(sx, sy)
    return cx + dx * s, cy + dy * s


def _angular_dist(a: float, b: float) -> float:
    """Signed angular distance b – a wrapped to (−180, 180]."""
    d = (b - a) % 360.0
    if d > 180.0:
        d -= 360.0
    return d


def _collision_free_angles(n: int, label_span: float = 26.0) -> List[float]:
    """
    Generate *n* evenly-spaced initial angles around a circle starting at the
    top (−90° offset already baked in via _ellipse_point), then push apart any
    that are too close using a greedy iterative relaxation.

    *label_span* is the angular width (degrees) a badge is assumed to occupy
    at the orbit radius; used as the minimum separation between two badge
    centres.
    """
    if n == 0:
        return []
    # Start evenly spaced
    angles = [(360.0 / n) * i for i in range(n)]

    # Iterative relaxation (max 60 passes)
    for _ in range(60):
        moved = False
        for i in range(n):
            j = (i + 1) % n
            diff = _angular_dist(angles[i], angles[j])
            if abs(diff) < label_span:
                push = (label_span - abs(diff)) / 2.0 + 0.5
                # push i back, j forward
                angles[i] = (angles[i] - push) % 360.0
                angles[j] = (angles[j] + push) % 360.0
                moved = True
        if not moved:
            break
    return angles


# ---------------------------------------------------------------------------
# LigPlot+ eyelash arc
# ---------------------------------------------------------------------------

def _eyelash_arc_svg(
    cx: float, cy: float,
    radius: float,
    angle_deg: float,
    color: str,
    span_deg: float = 50.0,
    num_spokes: int = 6,
    spoke_len: float = 8.0,
    stroke_w: float = 2.0,
) -> str:
    """
    Generate SVG markup for LigPlot+'s signature hydrophobic eyelash arc.

    The arc is centred at (cx, cy) at the given radius; spokes radiate
    outward (away from the ligand centre).
    """
    start_deg = angle_deg - span_deg / 2.0
    end_deg   = angle_deg + span_deg / 2.0

    r1 = math.radians(start_deg - 90.0)
    r2 = math.radians(end_deg   - 90.0)

    x1 = cx + radius * math.cos(r1)
    y1 = cy + radius * math.sin(r1)
    x2 = cx + radius * math.cos(r2)
    y2 = cy + radius * math.sin(r2)

    # large-arc flag: 0 when span < 180°
    large_arc = 1 if span_deg >= 180.0 else 0

    parts = [
        f'<path d="M {x1:.2f},{y1:.2f} A {radius:.2f},{radius:.2f} 0 {large_arc},1 '
        f'{x2:.2f},{y2:.2f}" fill="none" stroke="{color}" '
        f'stroke-width="{stroke_w:.1f}" stroke-linecap="round"/>'
    ]

    for i in range(num_spokes):
        frac = i / (num_spokes - 1) if num_spokes > 1 else 0.5
        a_rad = math.radians((start_deg + frac * span_deg) - 90.0)
        cos_a, sin_a = math.cos(a_rad), math.sin(a_rad)
        sx1 = cx + radius * cos_a
        sy1 = cy + radius * sin_a
        sx2 = cx + (radius + spoke_len) * cos_a
        sy2 = cy + (radius + spoke_len) * sin_a
        parts.append(
            f'<line x1="{sx1:.2f}" y1="{sy1:.2f}" x2="{sx2:.2f}" y2="{sy2:.2f}" '
            f'stroke="{color}" stroke-width="{stroke_w * 0.75:.1f}" stroke-linecap="round"/>'
        )

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Badge & line builders
# ---------------------------------------------------------------------------

def _pill_badge(
    bx: float, by: float,
    line1: str, line2: str,
    color: str, badge_bg: str, text_col: str,
    extra_count: int = 0,
    tooltip: str = "",
    anim_id: str = "",
) -> str:
    """
    Build a two-line pill badge SVG group centred at (bx, by).

    line1  – large text (e.g. 'ASP 102')
    line2  – small text (e.g. 'OD1')
    extra_count – number of additional interactions ('+ N more' tag)
    tooltip     – SVG <title> text for hover info
    """
    char_w1 = 7.2
    char_w2 = 6.0
    w1 = len(line1) * char_w1
    w2 = len(line2) * char_w2 if line2 else 0.0
    bw = max(52.0, max(w1, w2) + 18.0)
    bh = 30.0 if line2 else 20.0
    bh_half = bh / 2.0
    rx = bh_half  # full pill rounding

    parts = ["<g>"]
    if tooltip:
        parts.append(f"  <title>{html.escape(tooltip)}</title>")

    # Shadow glow
    parts.append(
        f'  <rect x="{bx - bw/2 - 1:.1f}" y="{by - bh_half - 1:.1f}" '
        f'width="{bw + 2:.1f}" height="{bh + 2:.1f}" rx="{rx:.1f}" '
        f'fill="{color}" opacity="0.18"/>'
    )
    # Main body
    parts.append(
        f'  <rect x="{bx - bw/2:.1f}" y="{by - bh_half:.1f}" '
        f'width="{bw:.1f}" height="{bh:.1f}" rx="{rx:.1f}" '
        f'fill="{badge_bg}" stroke="{color}" stroke-width="1.6"/>'
    )

    if line2:
        # Two-line layout
        parts.append(
            f'  <text x="{bx:.1f}" y="{by - 2:.1f}" text-anchor="middle" '
            f'font-size="10" font-weight="700" fill="{text_col}" '
            f'font-family="Segoe UI,system-ui,sans-serif">{html.escape(line1)}</text>'
        )
        parts.append(
            f'  <text x="{bx:.1f}" y="{by + 10:.1f}" text-anchor="middle" '
            f'font-size="8.5" font-weight="400" fill="{color}" '
            f'font-family="monospace">{html.escape(line2)}</text>'
        )
    else:
        parts.append(
            f'  <text x="{bx:.1f}" y="{by + 4:.1f}" text-anchor="middle" '
            f'font-size="10" font-weight="700" fill="{text_col}" '
            f'font-family="Segoe UI,system-ui,sans-serif">{html.escape(line1)}</text>'
        )

    # "+ N more" count badge
    if extra_count > 0:
        cbx = bx + bw / 2 + 7.0
        cby = by - bh_half + 7.0
        parts.append(
            f'  <circle cx="{cbx:.1f}" cy="{cby:.1f}" r="7" fill="{color}" opacity="0.9"/>'
        )
        parts.append(
            f'  <text x="{cbx:.1f}" y="{cby + 4:.1f}" text-anchor="middle" '
            f'font-size="7.5" font-weight="700" fill="white">+{extra_count}</text>'
        )

    parts.append("</g>")
    return "\n".join(parts)


def _interaction_line_svg(
    x1: float, y1: float, x2: float, y2: float,
    itype_key: str,
    color: str,
    dist: float = 0.0,
    card_bg: str = "#111827",
    anim_css_class: str = "",
) -> str:
    """
    Build the interaction line SVG between ligand boundary and badge anchor.

    For hydrogen bonds an animated dashed line and a distance badge are added.
    For salt bridges a thick solid line is drawn.
    """
    dash = _DASH_MAP.get(itype_key, "5,3")
    sw = 2.4 if "salt" in itype_key else (2.2 if "hydrogen" in itype_key else 1.8)

    anim_attr = f' class="{anim_css_class}"' if anim_css_class else ""
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""

    parts = [
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
        f'stroke="{color}" stroke-width="{sw:.1f}"{dash_attr} '
        f'stroke-linecap="round"{anim_attr}/>'
    ]

    # Salt bridge charge symbols at ¼ and ¾ positions
    if "salt" in itype_key:
        q1x, q1y = x1 + (x2 - x1) * 0.30, y1 + (y2 - y1) * 0.30
        q2x, q2y = x1 + (x2 - x1) * 0.70, y1 + (y2 - y1) * 0.70
        for qx, qy, sym in [(q1x, q1y, "−"), (q2x, q2y, "+")]:
            parts.append(
                f'<circle cx="{qx:.1f}" cy="{qy:.1f}" r="5.5" '
                f'fill="{color}" opacity="0.85"/>'
            )
            parts.append(
                f'<text x="{qx:.1f}" y="{qy + 3.5:.1f}" text-anchor="middle" '
                f'font-size="9" font-weight="900" fill="white">{sym}</text>'
            )

    # Pi-stacking ring symbols at midpoint
    if "pi-stacking" in itype_key or "pi-cation" in itype_key:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        parts.append(
            f'<text x="{mx:.1f}" y="{my + 5:.1f}" text-anchor="middle" '
            f'font-size="13" fill="{color}" opacity="0.9">⬡</text>'
        )

    # Halogen bond X marker
    if "halogen" in itype_key:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        parts.append(
            f'<text x="{mx:.1f}" y="{my + 4:.1f}" text-anchor="middle" '
            f'font-size="11" font-weight="900" fill="{color}" opacity="0.9">✕</text>'
        )

    # H-bond distance badge
    if "hydrogen" in itype_key and dist > 0:
        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        dist_txt = f"{dist:.2f}\u00c5"  # Å character
        tw = len(dist_txt) * 6.2 + 8
        parts.append(
            f'<rect x="{mx - tw/2:.1f}" y="{my - 9:.1f}" width="{tw:.1f}" height="14" '
            f'rx="5" fill="{card_bg}" stroke="{color}" stroke-width="0.9" opacity="0.95"/>'
        )
        parts.append(
            f'<text x="{mx:.1f}" y="{my + 2:.1f}" text-anchor="middle" '
            f'font-size="8.5" font-family="monospace" font-weight="700" fill="{color}">'
            f'{html.escape(dist_txt)}</text>'
        )

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Ligand-atom indicator dot
# ---------------------------------------------------------------------------

def _ligand_atom_dot(cx: float, cy: float, lig_hw: float, lig_hh: float,
                     angle_deg: float, color: str) -> str:
    """
    Draw a small glowing dot on the ligand rectangle boundary to indicate
    which face the interaction exits from.  We project a point on the rect
    edge in the direction of the interacting residue.
    """
    rad = math.radians(angle_deg - 90.0)
    # Approx edge point inside the ligand box (80 % of half-extent)
    ix = cx + lig_hw * 0.55 * math.cos(rad)
    iy = cy + lig_hh * 0.55 * math.sin(rad)
    return (
        f'<circle cx="{ix:.1f}" cy="{iy:.1f}" r="4" fill="{color}" opacity="0.85"/>'
        f'<circle cx="{ix:.1f}" cy="{iy:.1f}" r="7" fill="{color}" opacity="0.18"/>'
    )


# ---------------------------------------------------------------------------
# Legend builder
# ---------------------------------------------------------------------------

def _build_legend(
    y: float, W: float,
    present_types: List[str],
    theme: str,
    subtext_color: str,
    card_bg: str,
) -> str:
    """Build a compact SVG legend row showing only the interaction types that
    actually appear in this diagram."""

    # canonical ordering for display
    display_order = [
        ("Hydrogen Bond",              "H-Bond",       "hbond"),
        ("Salt Bridge",                "Salt Bridge",  "salt"),
        ("Halogen Bond",               "Halogen",      "halogen"),
        ("Pi-Stacking (Face-to-face)", "π-Stack",      "pi"),
        ("Pi-Stacking (T-shaped)",     "π-Stack T",    "pi"),
        ("Pi-Cation",                  "π-Cation",     "pi"),
        ("Hydrophobic",                "Hydrophobic",  "hydro"),
        ("Hydrophobic Contact",        "Hydrophobic",  "hydro"),
    ]

    present_lower = {t.lower() for t in present_types}
    items = []
    seen_labels = set()
    for full, short, shape in display_order:
        if full.lower() in present_lower and short not in seen_labels:
            c = get_interaction_color(full, theme)
            items.append((short, c, shape))
            seen_labels.add(short)

    if not items:
        return ""

    # Each item occupies approximately 90 px
    item_w = 90.0
    total_w = len(items) * item_w
    start_x = max(10.0, (W - total_w) / 2.0)
    lx = start_x
    parts: List[str] = []

    for label, c, shape in items:
        icon_cx = lx + 8.0
        icon_cy = y - 5.0

        if shape == "hbond":
            parts.append(
                f'<line x1="{lx:.0f}" y1="{icon_cy:.0f}" x2="{lx + 16:.0f}" y2="{icon_cy:.0f}" '
                f'stroke="{c}" stroke-width="2.2" stroke-dasharray="5,3" stroke-linecap="round"/>'
            )
        elif shape == "salt":
            parts.append(
                f'<line x1="{lx:.0f}" y1="{icon_cy:.0f}" x2="{lx + 16:.0f}" y2="{icon_cy:.0f}" '
                f'stroke="{c}" stroke-width="3" stroke-linecap="round"/>'
            )
            parts.append(
                f'<circle cx="{lx + 4:.0f}" cy="{icon_cy:.0f}" r="3.5" fill="{c}" opacity="0.85"/>'
                f'<text x="{lx + 4:.0f}" y="{icon_cy + 2.5:.0f}" text-anchor="middle" '
                f'font-size="6" fill="white">−</text>'
            )
            parts.append(
                f'<circle cx="{lx + 12:.0f}" cy="{icon_cy:.0f}" r="3.5" fill="{c}" opacity="0.85"/>'
                f'<text x="{lx + 12:.0f}" y="{icon_cy + 2.5:.0f}" text-anchor="middle" '
                f'font-size="6" fill="white">+</text>'
            )
        elif shape == "pi":
            parts.append(
                f'<line x1="{lx:.0f}" y1="{icon_cy:.0f}" x2="{lx + 16:.0f}" y2="{icon_cy:.0f}" '
                f'stroke="{c}" stroke-width="2" stroke-dasharray="4,2.5"/>'
            )
            parts.append(
                f'<text x="{lx + 8:.0f}" y="{icon_cy + 5:.0f}" text-anchor="middle" '
                f'font-size="10" fill="{c}">⬡</text>'
            )
        elif shape == "halogen":
            parts.append(
                f'<line x1="{lx:.0f}" y1="{icon_cy:.0f}" x2="{lx + 16:.0f}" y2="{icon_cy:.0f}" '
                f'stroke="{c}" stroke-width="2" stroke-dasharray="4,2.5"/>'
            )
            parts.append(
                f'<text x="{lx + 8:.0f}" y="{icon_cy + 4:.0f}" text-anchor="middle" '
                f'font-size="9" font-weight="900" fill="{c}">✕</text>'
            )
        elif shape == "hydro":
            # Mini eyelash arc
            r_mini = 5.0
            parts.append(
                f'<path d="M {lx:.0f},{icon_cy:.0f} A {r_mini:.0f},{r_mini:.0f} 0 0,1 '
                f'{lx + 16:.0f},{icon_cy:.0f}" fill="none" stroke="{c}" stroke-width="2" '
                f'stroke-linecap="round"/>'
            )
            for k in range(4):
                sx = lx + 2 + k * 4
                parts.append(
                    f'<line x1="{sx:.0f}" y1="{icon_cy:.0f}" x2="{sx:.0f}" y2="{icon_cy - 5:.0f}" '
                    f'stroke="{c}" stroke-width="1.3" stroke-linecap="round"/>'
                )

        # Label
        parts.append(
            f'<text x="{lx + 20:.0f}" y="{y:.0f}" font-size="8.5" '
            f'font-family="Segoe UI,system-ui,sans-serif" fill="{subtext_color}">'
            f'{html.escape(label)}</text>'
        )
        lx += item_w

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# No-interaction placeholder
# ---------------------------------------------------------------------------

def _no_interaction_placeholder(
    cx: float, cy: float, diag_h: float, W: float,
    subtext_color: str, card_bg: str, border_color: str,
) -> str:
    """Return SVG markup for when there are no interactions to display."""
    ew, eh = W * 0.7, diag_h * 0.35
    ex, ey = cx - ew / 2, cy - eh / 2
    parts = [
        f'<rect x="{ex:.1f}" y="{ey:.1f}" width="{ew:.1f}" height="{eh:.1f}" '
        f'rx="12" fill="{card_bg}" stroke="{border_color}" stroke-width="1" stroke-dasharray="6,3"/>',
        f'<text x="{cx:.1f}" y="{cy - 8:.1f}" text-anchor="middle" font-size="14" '
        f'fill="{subtext_color}" font-family="Segoe UI,system-ui,sans-serif">🔬</text>',
        f'<text x="{cx:.1f}" y="{cy + 10:.1f}" text-anchor="middle" font-size="10" '
        f'fill="{subtext_color}" font-family="Segoe UI,system-ui,sans-serif">'
        f'No interactions detected</text>',
    ]
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Main public renderer
# ---------------------------------------------------------------------------

def render_interaction_diagram(
    ligand_path,
    interactions: List[Any],
    pose_index: int = 1,
    run_number: Optional[int] = None,
    original_ligand_path: Optional[Any] = None,
    pose_metadata: Optional[Dict[str, Any]] = None,
    obabel_exe=None,
    width: int = 480,
    height: int = 380,
    theme: str = "dark",
    publication_mode: bool = False,
) -> InteractionDiagramResult:
    """
    Generate a self-contained, publication-quality 2D protein-ligand
    interaction diagram for a specific docking pose.

    Returns an InteractionDiagramResult containing both the exact outer SVG
    and the self-contained HTML representation.
    """
    lig_p = Path(ligand_path) if ligand_path else Path("ligand.pdbqt")
    is_pub = bool(publication_mode)
    is_dark = (theme.lower() == "dark") and not is_pub

    # ── Theme tokens ──────────────────────────────────────────────────────
    if is_pub:
        bg_color      = "#ffffff"
        card_bg       = "#ffffff"
        panel_bg      = "#f8fafc"
        text_color    = "#0f172a"
        subtext_color = "#475569"
        border_color  = "#94a3b8"
        badge_bg      = "#f1f5f9"
        grid_color    = "#ffffff"
    elif is_dark:
        bg_color      = "#0b0f1a"
        card_bg       = "#111827"
        panel_bg      = "#0f1724"
        text_color    = "#f1f5f9"
        subtext_color = "#94a3b8"
        border_color  = "#1e293b"
        badge_bg      = "#1c2535"
        grid_color    = "#1a2233"
    else:
        bg_color      = "#f8fafc"
        card_bg       = "#ffffff"
        panel_bg      = "#f0f4f8"
        text_color    = "#0f172a"
        subtext_color = "#64748b"
        border_color  = "#cbd5e1"
        badge_bg      = "#f1f5f9"
        grid_color    = "#e8eef4"

    # ── Layout constants ──────────────────────────────────────────────────
    W  = float(width)
    H  = float(height)

    header_h = 36.0
    legend_h = 38.0
    diag_h   = H - header_h - legend_h        # working area height
    cx       = W / 2.0
    cy       = header_h + diag_h / 2.0        # ligand centre

    # Ligand panel size (square)
    lig_size = min(W * 0.40, diag_h * 0.48)
    lig_hw   = lig_size / 2.0
    lig_hh   = lig_size / 2.0

    # Orbital radii
    orbit_rx = W * 0.43
    orbit_ry = diag_h * 0.46

    # Badge half-dimensions
    badge_half_w = 46.0
    badge_half_h = 18.0

    # ── Load molecule for EXACT selected pose ──────────────────────────────
    mol = None
    if is_available() and ligand_path:
        mol = load_pose_molecule(
            ligand_path=lig_p,
            pose_index=pose_index,
            run_number=run_number,
            original_ligand_path=original_ligand_path,
            obabel_exe=obabel_exe,
        )

    # ── 2-D SVG depiction ─────────────────────────────────────────────────
    lig_svg_inner: Optional[str] = None
    lig_vb:        Optional[str] = None
    mol_theme = "light" if is_pub else theme
    if mol is not None:
        raw_svg = _rdkit_mol_svg(mol, int(lig_size), int(lig_size), mol_theme)
        if raw_svg:
            lig_vb, lig_svg_inner = _extract_svg_inner(raw_svg)

    # ── Group & prioritise interactions per residue ───────────────────────
    inter_list = getattr(interactions, "interactions", interactions) or []
    residue_data: Dict[str, Dict[str, Any]] = {}

    for it in inter_list:
        res       = getattr(it, "receptor_residue",  "UNK").strip()
        res_name  = getattr(it, "receptor_res_name", "").strip()
        res_seq   = getattr(it, "receptor_res_seq",  0)
        chain     = getattr(it, "receptor_chain",    "").strip()
        itype     = getattr(it, "interaction_type",  "Hydrophobic").strip()
        dist      = float(getattr(it, "distance_angstrom", 0.0) or 0.0)
        rec_atom  = getattr(it, "receptor_atom",     "").strip()
        lig_atom  = getattr(it, "ligand_atom",       "").strip()
        angle     = getattr(it, "angle_deg",         None)
        basis     = getattr(it, "detection_basis",   "").strip()

        # Build a display label for the residue: "ASP 102" style
        if res_name and res_seq:
            display_label = f"{res_name} {res_seq}"
        elif res_name:
            display_label = res_name
        else:
            display_label = res

        angle_str = f" · {angle:.1f}°" if angle is not None else ""
        info_line = f"{itype}: {rec_atom}→{lig_atom} {dist:.2f}Å{angle_str}" if dist else f"{itype}"
        if basis:
            info_line += f" ({basis})"

        if res not in residue_data:
            residue_data[res] = {
                "residue":        res,
                "display_label":  display_label,
                "chain":          chain,
                "dominant_type":  itype,
                "min_dist":       dist,
                "rec_atom":       rec_atom,
                "lig_atom":       lig_atom,
                "angle":          angle,
                "all_types":      [itype],
                "all_details":    [info_line],
            }
        else:
            residue_data[res]["all_types"].append(itype)
            residue_data[res]["all_details"].append(info_line)
            if dist > 0 and (residue_data[res]["min_dist"] == 0 or dist < residue_data[res]["min_dist"]):
                residue_data[res]["min_dist"] = dist
                residue_data[res]["rec_atom"] = rec_atom
                residue_data[res]["lig_atom"] = lig_atom
                residue_data[res]["angle"] = angle

    # Determine dominant interaction type by priority
    priority_index = {t: i for i, t in enumerate(_INTERACTION_PRIORITY)}
    for rinfo in residue_data.values():
        best = min(
            rinfo["all_types"],
            key=lambda t: priority_index.get(t, 99),
        )
        rinfo["dominant_type"] = best

    # Sort and limit to top 16 residues
    MAX_BADGES = 16
    sorted_residues = sorted(
        residue_data.values(),
        key=lambda x: (
            priority_index.get(x["dominant_type"], 99),
            x["min_dist"] if x["min_dist"] > 0 else 99.0,
        ),
    )[:MAX_BADGES]

    n_res = len(sorted_residues)
    lig_title = lig_p.stem.replace("_out", "").replace("_", " ")

    # Compute collision-free angles
    avg_r = (orbit_rx + orbit_ry) / 2.0
    label_span_deg = math.degrees(2 * math.atan2(badge_half_w, avg_r)) + 4.0
    angles = _collision_free_angles(n_res, label_span=max(label_span_deg, 22.0))

    # ── SVG construction ──────────────────────────────────────────────────
    svg: List[str] = []

    # Root SVG element
    svg.append(
        f'<svg width="100%" height="100%" viewBox="0 0 {W:.1f} {H:.1f}" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:adsp="https://autodocksuite.pro/schema" '
        f'style="background:{bg_color};font-family:Segoe UI,-apple-system,system-ui,sans-serif;">'
    )

    # ── ADSP Structural Metadata ──────────────────────────────────────────
    pmeta = pose_metadata or {}
    engine_name = html.escape(str(pmeta.get("engine", "AutoDock / Vina")))
    mode_name   = html.escape(str(pmeta.get("docking_mode", "Rigid / Flexible")))
    job_id_str  = html.escape(str(pmeta.get("job_id", "")))
    svg.append(f"""<metadata id="adsp-metadata">
  <adsp:provenance>
    <adsp:generator>AutoDockSuite Pro</adsp:generator>
    <adsp:pose_index>{pose_index}</adsp:pose_index>
    <adsp:run_number>{run_number if run_number is not None else ""}</adsp:run_number>
    <adsp:ligand>{html.escape(lig_title)}</adsp:ligand>
    <adsp:engine>{engine_name}</adsp:engine>
    <adsp:docking_mode>{mode_name}</adsp:docking_mode>
    <adsp:job_id>{job_id_str}</adsp:job_id>
    <adsp:contacts>{len(inter_list)}</adsp:contacts>
    <adsp:residues>{len(residue_data)}</adsp:residues>
    <adsp:publication_mode>{str(is_pub).lower()}</adsp:publication_mode>
  </adsp:provenance>
</metadata>""")

    # ── CSS styles ────────────────────────────────────────────────────────
    if is_pub:
        # Static publication styles (vector-friendly, no animation or glow)
        svg.append("""<defs>
  <style>
    .hbond-line { stroke-dasharray: 6 3; }
    .badge-group text { font-family: Segoe UI, system-ui, sans-serif; }
  </style>
</defs>""")
    else:
        # Interactive dark/light studio styles
        svg.append("""<defs>
  <style>
    .hbond-line {
      stroke-dasharray: 7 4;
      animation: marchDash 1.4s linear infinite;
    }
    @keyframes marchDash {
      to { stroke-dashoffset: -22; }
    }
    .badge-group:hover rect { filter: brightness(1.18); cursor: pointer; }
  </style>
  <filter id="glow" x="-30%" y="-30%" width="160%" height="160%">
    <feGaussianBlur stdDeviation="2.5" result="blur"/>
    <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>
  </filter>
</defs>""")

    # ── Grid background (Screen mode only) ────────────────────────────────
    if not is_pub:
        svg.append(
            f'<pattern id="dotgrid" x="0" y="0" width="20" height="20" patternUnits="userSpaceOnUse">'
            f'<circle cx="10" cy="10" r="0.8" fill="{grid_color}"/>'
            f'</pattern>'
            f'<rect width="{W}" height="{H}" fill="url(#dotgrid)"/>'
        )

    # ── Header bar ────────────────────────────────────────────────────
    svg.append(
        f'<rect x="0" y="0" width="{W:.1f}" height="{header_h:.1f}" '
        f'fill="{panel_bg}" opacity="{"1.0" if is_pub else "0.9"}"/>'
    )
    # Header title with explicit POSE identity
    header_subtitle = f" · Pose {pose_index}"
    if run_number is not None and run_number != pose_index:
        header_subtitle += f" (Run {run_number})"

    svg.append(
        f'<text x="14" y="23" font-size="12.5" '
        f'font-weight="700" fill="{text_color}" letter-spacing="0.3">'
        f'2D Interaction Map'
        f'<tspan font-weight="700" fill="{get_interaction_color("Hydrogen Bond", "light" if is_pub else theme)}">{html.escape(header_subtitle)}</tspan>'
        f'<tspan font-weight="400" fill="{subtext_color}"> · {html.escape(lig_title)}</tspan>'
        f'</text>'
    )

    # Header summary metrics
    if inter_list:
        n_total = len(inter_list)
        n_res_total = len(residue_data)
        # Class counts
        type_counts: Dict[str, int] = {}
        for it in inter_list:
            t = getattr(it, "interaction_type", "Contact")
            type_counts[t] = type_counts.get(t, 0) + 1

        top_classes = [f"{t[:4]} {c}" for t, c in sorted(type_counts.items(), key=lambda kv: -kv[1])[:3]]
        summary_str = f"{n_total} contacts · {n_res_total} res · " + " · ".join(top_classes)
        svg.append(
            f'<text x="{W - 14:.1f}" y="23" text-anchor="end" font-size="9" '
            f'fill="{subtext_color}">{html.escape(summary_str)}</text>'
        )

    # ── No-interaction placeholder ─────────────────────────────────────
    if n_res == 0:
        svg.append(_no_interaction_placeholder(
            cx, cy, diag_h, W, subtext_color, card_bg, border_color,
        ))
        svg.append("</svg>")
        full_svg = "\n".join(svg)
        doc_html = _wrap_html(full_svg, bg_color, lig_title, pose_index)
        return InteractionDiagramResult(
            svg=full_svg,
            html=doc_html,
            pose_index=pose_index,
            metadata={"pose_index": pose_index, "contacts": 0, "residues": 0},
            contact_count=0,
            residue_count=0,
            interaction_classes=[],
        )

    # ── Ligand panel ──────────────────────────────────────────────────
    panel_x = cx - lig_hw
    panel_y = cy - lig_hh
    if not is_pub:
        # Outer soft border in studio mode
        svg.append(
            f'<rect x="{panel_x - 3:.1f}" y="{panel_y - 3:.1f}" '
            f'width="{lig_size + 6:.1f}" height="{lig_size + 6:.1f}" '
            f'rx="14" fill="{border_color}" opacity="0.3"/>'
        )
    # Main panel box
    svg.append(
        f'<rect x="{panel_x:.1f}" y="{panel_y:.1f}" '
        f'width="{lig_size:.1f}" height="{lig_size:.1f}" '
        f'rx="10" fill="{card_bg}" stroke="{border_color}" stroke-width="{"1.2" if is_pub else "1.0"}"/>'
    )
    # Pose watermark in panel
    svg.append(
        f'<text x="{cx:.1f}" y="{panel_y + 11:.1f}" text-anchor="middle" '
        f'font-size="7.5" font-weight="600" fill="{subtext_color}" opacity="0.8">LIGAND · POSE {pose_index}</text>'
    )

    # Embed molecule SVG or placeholder
    if lig_svg_inner and lig_vb:
        vb_parts = lig_vb.split()
        vbw = float(vb_parts[2]) if len(vb_parts) == 4 else lig_size
        vbh = float(vb_parts[3]) if len(vb_parts) == 4 else lig_size
        svg.append(
            f'<svg x="{panel_x + 4:.1f}" y="{panel_y + 10:.1f}" '
            f'width="{lig_size - 8:.1f}" height="{lig_size - 14:.1f}" '
            f'viewBox="0 0 {vbw:.1f} {vbh:.1f}">'
        )
        svg.append(lig_svg_inner)
        svg.append("</svg>")
    else:
        # Graceful text placeholder
        svg.append(
            f'<rect x="{panel_x + 8:.1f}" y="{panel_y + 12:.1f}" '
            f'width="{lig_size - 16:.1f}" height="{lig_size - 20:.1f}" '
            f'rx="8" fill="none" stroke="{subtext_color}" '
            f'stroke-width="0.9" stroke-dasharray="4,4"/>'
        )
        svg.append(
            f'<text x="{cx:.1f}" y="{cy:.1f}" text-anchor="middle" '
            f'font-size="10" font-weight="600" fill="{subtext_color}">{html.escape(lig_title)}</text>'
        )
        svg.append(
            f'<text x="{cx:.1f}" y="{cy + 14:.1f}" text-anchor="middle" '
            f'font-size="8" fill="{subtext_color}" opacity="0.8">Pose {pose_index}</text>'
        )

    # ── Pass 1: Interaction lines & indicators ─────────────────────────
    atom_dot_svgs: List[str] = []
    badge_data: List[Tuple[float, float, Dict[str, Any]]] = []

    for idx, rinfo in enumerate(sorted_residues):
        angle      = angles[idx]
        itype      = rinfo["dominant_type"]
        itype_key  = itype.lower()
        color      = get_interaction_color(itype, "light" if is_pub else theme)
        dist_val   = rinfo["min_dist"]
        rec_atom   = rinfo["rec_atom"]
        lig_atom   = rinfo["lig_atom"]

        # Badge anchor point on the orbit ellipse
        bx, by = _ellipse_point(cx, cy, orbit_rx, orbit_ry, angle)

        # Clamp to keep inside the viewBox with margin
        margin_x = badge_half_w + 4.0
        margin_y = badge_half_h + 6.0
        bx = max(margin_x, min(W - margin_x, bx))
        by = max(header_h + margin_y, min(H - legend_h - margin_y, by))

        # Interaction line origin: edge of the ligand panel
        lx, ly = _rect_edge_pt(cx, cy, lig_hw * 2, lig_hh * 2, bx, by)

        # Hydrophobic: eyelash arc
        if "hydrophobic" in itype_key:
            arc_r = lig_hw * 0.28
            arc_angle = math.degrees(math.atan2(by - cy, bx - cx)) + 90.0
            svg.append(_eyelash_arc_svg(
                lx, ly, arc_r, arc_angle + 180.0, color,
                span_deg=54.0, num_spokes=6, spoke_len=9.0, stroke_w=2.0,
            ))
        else:
            anim_cls = ("" if is_pub else "hbond-line") if "hydrogen" in itype_key else ""
            svg.append(_interaction_line_svg(
                lx, ly, bx, by, itype_key, color, dist_val,
                card_bg=card_bg, anim_css_class=anim_cls,
            ))

            # Atom indicator dot
            if "hydrogen" in itype_key or "salt" in itype_key or "halogen" in itype_key:
                atom_dot_svgs.append(
                    _ligand_atom_dot(cx, cy, lig_hw, lig_hh, angle, color)
                )

        badge_data.append((bx, by, rinfo))

    # Atom dots
    for dot_svg in atom_dot_svgs:
        svg.append(dot_svg)

    # ── Pass 2: Residue badges (rendered on top of lines) ─────────────
    for bx, by, rinfo in badge_data:
        itype      = rinfo["dominant_type"]
        color      = get_interaction_color(itype, "light" if is_pub else theme)
        display_lb = rinfo["display_label"]
        rec_atom   = rinfo["rec_atom"]
        all_types  = rinfo["all_types"]
        all_details = rinfo["all_details"]

        line1 = display_lb
        line2 = rec_atom if rec_atom else ""

        itype_key = itype.lower()
        if "pi-stacking" in itype_key or "pi-cation" in itype_key:
            line1 = "⬡ " + line1
        elif "halogen" in itype_key:
            line1 = "✕ " + line1

        tooltip_lines = [
            f"Residue: {rinfo['residue']}",
            f"Pose: {pose_index}",
        ] + all_details
        tooltip = "\n".join(tooltip_lines)

        unique_types = list(dict.fromkeys(all_types))
        extra = len(unique_types) - 1

        svg.append(_pill_badge(
            bx, by, line1, line2, color, badge_bg, text_color,
            extra_count=extra,
            tooltip=tooltip,
        ))

    # ── Legend (Active classes only) ──────────────────────────────────
    leg_y = H - 10.0
    present_types = list({rinfo["dominant_type"] for _, _, rinfo in badge_data})
    svg.append(
        f'<rect x="0" y="{H - legend_h:.1f}" width="{W:.1f}" height="{legend_h:.1f}" '
        f'fill="{panel_bg}" opacity="{"1.0" if is_pub else "0.85"}"/>'
    )
    legend_svg = _build_legend(
        leg_y - 6, W, present_types, "light" if is_pub else theme, subtext_color, card_bg,
    )
    svg.append(legend_svg)

    svg.append("</svg>")
    full_svg = "\n".join(svg)
    doc_html = _wrap_html(full_svg, bg_color, lig_title, pose_index)

    metadata_dict = {
        "pose_index": pose_index,
        "run_number": run_number,
        "ligand": lig_title,
        "contacts": len(inter_list),
        "residues": len(residue_data),
        "interaction_classes": present_types,
        "publication_mode": is_pub,
        "theme": theme,
    }
    if pose_metadata:
        metadata_dict.update(pose_metadata)

    return InteractionDiagramResult(
        svg=full_svg,
        html=doc_html,
        pose_index=pose_index,
        metadata=metadata_dict,
        contact_count=len(inter_list),
        residue_count=len(residue_data),
        interaction_classes=present_types,
    )


# ---------------------------------------------------------------------------
# Compatibility and convenience wrappers
# ---------------------------------------------------------------------------

def render_interaction_diagram_html(
    ligand_path,
    interactions: List[Any],
    obabel_exe=None,
    width: int = 480,
    height: int = 380,
    theme: str = "dark",
    pose_index: int = 1,
    run_number: Optional[int] = None,
    original_ligand_path: Optional[Any] = None,
    pose_metadata: Optional[Dict[str, Any]] = None,
    publication_mode: bool = False,
) -> Optional[str]:
    """
    Generate a self-contained HTML 2D interaction diagram.
    """
    res = render_interaction_diagram(
        ligand_path=ligand_path,
        interactions=interactions,
        pose_index=pose_index,
        run_number=run_number,
        original_ligand_path=original_ligand_path,
        pose_metadata=pose_metadata,
        obabel_exe=obabel_exe,
        width=width,
        height=height,
        theme=theme,
        publication_mode=publication_mode,
    )
    return res.html if res else None


def render_interaction_diagram_svg(
    ligand_path,
    interactions: List[Any],
    obabel_exe=None,
    width: int = 480,
    height: int = 380,
    theme: str = "dark",
    pose_index: int = 1,
    run_number: Optional[int] = None,
    original_ligand_path: Optional[Any] = None,
    pose_metadata: Optional[Dict[str, Any]] = None,
    publication_mode: bool = False,
) -> Optional[str]:
    """
    Generate an unadulterated, valid outer SVG 2D interaction diagram.
    """
    res = render_interaction_diagram(
        ligand_path=ligand_path,
        interactions=interactions,
        pose_index=pose_index,
        run_number=run_number,
        original_ligand_path=original_ligand_path,
        pose_metadata=pose_metadata,
        obabel_exe=obabel_exe,
        width=width,
        height=height,
        theme=theme,
        publication_mode=publication_mode,
    )
    return res.svg if res else None


# ---------------------------------------------------------------------------
# HTML wrapper
# ---------------------------------------------------------------------------

def _wrap_html(full_svg: str, bg_color: str, title: str, pose_index: int = 1) -> str:
    """Wrap the SVG in a minimal, self-contained, responsive HTML document."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1.0"/>
<title>2D Interaction Map · {html.escape(title)} · Pose {pose_index}</title>
<style>
  *{{box-sizing:border-box;margin:0;padding:0}}
  html,body{{
    width:100%;height:100%;
    overflow:hidden;
    background:{bg_color};
    display:flex;
    align-items:center;
    justify-content:center;
  }}
  svg{{
    max-width:100%;
    max-height:100%;
    display:block;
    user-select:none;
    -webkit-user-select:none;
  }}
</style>
</head>
<body>
{full_svg}
</body>
</html>"""

