#!/usr/bin/env python3
"""
AutoDock Suite Pro — Interactive 3D Molecular WebGL Visualizer
=============================================================
Generates standalone interactive 3D WebGL visualizations using 3Dmol.js.
Visualizes:
  - Receptor cartoon & transparent molecular surface
  - Active site target residues (rendered as sticks with labels)
  - Biophysical 3D wireframe Grid Box overlaying the binding pocket
  - Docked ligand poses with pose switcher (Mode 1 to N)
  - Real-time hydrogen bond dashed lines and interaction distance labels
  - 1-Click browser launch from Python / GUI
"""

from __future__ import annotations

import json
import logging
import tempfile
import webbrowser
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from interactions import (
    PDBQTAtom, detect_hydrogen_bonds, parse_docked_poses, parse_receptor_pdbqt,
)

logger = logging.getLogger("docking_automation.viewer_3d")


def _read_file_safe(file_path: Path | str) -> str:
    p = Path(file_path)
    if not p.is_file():
        return ""
    try:
        return p.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""


def generate_3d_viewer_html(
    receptor_path: Path | str,
    ligand_path: Optional[Path | str] = None,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    target_residues: Optional[List[str]] = None,
    title: str = "AutoDock Suite Pro — 3D Molecular Inspector",
) -> str:
    """Builds a complete, standalone HTML string embedding 3Dmol.js with interactive controls.

    Args:
        receptor_path: Path to receptor PDB or PDBQT.
        ligand_path: Optional path to multi-pose docked ligand PDBQT or single pose.
        grid_center: (cx, cy, cz) coordinates.
        grid_size: (sx, sy, sz) dimensions.
        target_residues: List of residue identifiers (e.g. ['GLN:206', 'ASP:102']).
        title: Page title.

    Returns:
        HTML string.
    """
    receptor_content = _read_file_safe(receptor_path)
    receptor_ext = Path(receptor_path).suffix.lstrip(".").lower()
    if receptor_ext not in ("pdb", "pdbqt"):
        receptor_ext = "pdbqt"

    # Process ligand poses if available
    poses_json: List[Dict[str, Any]] = []
    hbonds_per_pose: List[List[Dict[str, Any]]] = []

    if ligand_path and Path(ligand_path).is_file():
        rec_atoms = parse_receptor_pdbqt(receptor_path)
        rec_atom_map = {(a.res_name, a.res_seq, a.name): a for a in rec_atoms}
        parsed_poses = parse_docked_poses(ligand_path)

        for pose_idx, lig_atoms, flex_atoms in parsed_poses:
            all_pose_atoms = lig_atoms + flex_atoms
            lig_atom_map = {a.name: a for a in lig_atoms}
            hbonds = detect_hydrogen_bonds(rec_atoms, lig_atoms, pose_idx=pose_idx)
            pose_text_lines = []
            for a in all_pose_atoms:
                # Format simple PDB ATOM record for 3Dmol
                line = (
                    f"ATOM  {a.serial:5d} {a.name:<4s} "
                    f"LIG A   1    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
                    f"  1.00 20.00          {a.element:>2s}"
                )
                pose_text_lines.append(line)
            pose_pdb_text = "\n".join(pose_text_lines) + "\nEND\n"

            poses_json.append({
                "pose": pose_idx,
                "pdb": pose_pdb_text,
                "num_atoms": len(lig_atoms),
            })

            hb_list = []
            for hb in hbonds:
                rec_a = rec_atom_map.get((hb.receptor_res_name, hb.receptor_res_seq, hb.receptor_atom))
                lig_a = lig_atom_map.get(hb.ligand_atom)
                if rec_a and lig_a:
                    hb_list.append({
                        "res": f"{hb.receptor_res_name} {hb.receptor_res_seq}",
                        "dist": round(hb.distance_angstrom, 2),
                        "rec_x": rec_a.x,
                        "rec_y": rec_a.y,
                        "rec_z": rec_a.z,
                        "lig_x": lig_a.x,
                        "lig_y": lig_a.y,
                        "lig_z": lig_a.z,
                        "lig_atom": hb.ligand_atom,
                        "rec_atom": hb.receptor_atom,
                    })
            hbonds_per_pose.append(hb_list)

    # Escape contents for embedded JavaScript
    safe_rec_content = json.dumps(receptor_content)
    safe_poses_json = json.dumps(poses_json)
    safe_hbonds_json = json.dumps(hbonds_per_pose)
    safe_target_res = json.dumps(target_residues or [])

    center_obj = json.dumps(
        {"x": grid_center[0], "y": grid_center[1], "z": grid_center[2]}
        if grid_center else None
    )
    size_obj = json.dumps(
        {"x": grid_size[0], "y": grid_size[1], "z": grid_size[2]}
        if grid_size else None
    )

    rec_name = Path(receptor_path).stem
    lig_name = Path(ligand_path).stem if ligand_path else "None"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root {{
      --bg-main: #0d1117;
      --bg-panel: #161b22;
      --border-color: #30363d;
      --accent: #58a6ff;
      --accent-green: #3fb950;
      --accent-orange: #d29922;
      --text: #c9d1d9;
      --text-bright: #f0f6fc;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg-main);
      color: var(--text);
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }}
    header {{
      background: var(--bg-panel);
      border-bottom: 1px solid var(--border-color);
      padding: 10px 20px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 15px;
    }}
    header h1 {{
      font-size: 16px;
      font-weight: 700;
      color: var(--accent);
      letter-spacing: 0.5px;
      display: flex;
      align-items: center;
      gap: 10px;
    }}
    .badge {{
      font-size: 11px;
      padding: 2px 8px;
      border-radius: 12px;
      background: #1f6feb33;
      color: var(--accent);
      border: 1px solid #1f6feb;
    }}
    .toolbar {{
      display: flex;
      align-items: center;
      gap: 12px;
      flex-wrap: wrap;
    }}
    .btn {{
      background: #21262d;
      border: 1px solid var(--border-color);
      color: var(--text-bright);
      padding: 6px 14px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 600;
      cursor: pointer;
      transition: all 0.15s ease;
      display: inline-flex;
      align-items: center;
      gap: 6px;
    }}
    .btn:hover {{
      background: #30363d;
      border-color: #8b949e;
    }}
    .btn.active {{
      background: #1f6feb;
      border-color: #388bfd;
      color: #fff;
    }}
    select {{
      background: #21262d;
      border: 1px solid var(--border-color);
      color: var(--text-bright);
      padding: 6px 10px;
      border-radius: 6px;
      font-size: 12px;
    }}
    #viewer-container {{
      flex: 1;
      position: relative;
      width: 100%;
      height: 100%;
    }}
    #gldiv {{
      width: 100%;
      height: 100%;
    }}
    #legend-panel {{
      position: absolute;
      top: 15px;
      left: 15px;
      background: rgba(22, 27, 34, 0.85);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 12px 16px;
      font-size: 12px;
      line-height: 1.6;
      max-width: 280px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.4);
      z-index: 10;
    }}
    #legend-panel h3 {{
      font-size: 13px;
      color: var(--text-bright);
      margin-bottom: 6px;
      border-bottom: 1px solid var(--border-color);
      padding-bottom: 4px;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin-top: 4px;
    }}
    .color-swatch {{
      width: 12px;
      height: 12px;
      border-radius: 3px;
      display: inline-block;
    }}
    #interaction-card {{
      position: absolute;
      bottom: 15px;
      right: 15px;
      background: rgba(22, 27, 34, 0.88);
      backdrop-filter: blur(8px);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 12px 16px;
      font-size: 12px;
      max-height: 220px;
      width: 320px;
      overflow-y: auto;
      box-shadow: 0 4px 12px rgba(0,0,0,0.4);
      z-index: 10;
    }}
    #interaction-card h4 {{
      font-size: 12px;
      color: var(--accent-green);
      margin-bottom: 6px;
    }}
    .hb-row {{
      display: flex;
      justify-content: space-between;
      padding: 3px 0;
      border-bottom: 1px dashed #21262d;
    }}
  </style>
</head>
<body>
  <header>
    <h1>
      <span>AutoDock Suite Pro</span>
      <span class="badge">3D WebGL Viewer</span>
    </h1>
    <div class="toolbar">
      <label for="pose-select" style="font-size: 12px; color: var(--text);">Docked Pose:</label>
      <select id="pose-select" onchange="onPoseChange(this.value)">
        <!-- Poses added dynamically -->
      </select>
      <button class="btn" id="btn-surface" onclick="toggleSurface()">🌐 Toggle Surface</button>
      <button class="btn" id="btn-box" onclick="toggleGridBox()">📦 Toggle Grid Box</button>
      <button class="btn" id="btn-hbonds" onclick="toggleHBonds()">⚡ Toggle H-Bonds</button>
      <button class="btn" onclick="resetCamera()">🎯 Reset View</button>
      <button class="btn" onclick="downloadPMLScript()" style="border-color:#3fb950; color:#3fb950;">🔬 Export PyMOL (.pml)</button>
    </div>
  </header>

  <div id="viewer-container">
    <div id="gldiv"></div>

    <div id="legend-panel">
      <h3>Target & Ligand Info</h3>
      <div><strong>Receptor:</strong> {rec_name}</div>
      <div><strong>Ligand:</strong> {lig_name}</div>
      <div class="legend-item"><span class="color-swatch" style="background:#58a6ff;"></span> Receptor Cartoon</div>
      <div class="legend-item"><span class="color-swatch" style="background:#00d26a;"></span> Ligand Carbons</div>
      <div class="legend-item"><span class="color-swatch" style="background:#ff9f43;"></span> Active-Site Residues</div>
      <div class="legend-item"><span class="color-swatch" style="background:#00ffff;"></span> Grid Box Bounding Volume</div>
      <div class="legend-item"><span class="color-swatch" style="background:#ffd32a;"></span> H-Bond Contacts (d ≤ 3.5 Å)</div>
    </div>

    <div id="interaction-card">
      <h4>Identified Pose Interactions</h4>
      <div id="interaction-list">Loading interactions...</div>
    </div>
  </div>

  <script>
    const receptorData = {safe_rec_content};
    const receptorExt = "{receptor_ext}";
    const posesData = {safe_poses_json};
    const hbondsData = {safe_hbonds_json};
    const targetResidues = {safe_target_res};
    const gridCenter = {center_obj};
    const gridSize = {size_obj};

    let glviewer = null;
    let recModel = null;
    let ligModel = null;
    let surfMesh = null;
    let boxShapes = [];
    let hbondLines = [];
    let showSurface = false;
    let showGrid = true;
    let showHBonds = true;
    let currentPoseIndex = 0;

    window.addEventListener('DOMContentLoaded', () => {{
      const element = document.getElementById('gldiv');
      const config = {{ backgroundColor: '#0d1117' }};
      glviewer = $3Dmol.createViewer(element, config);

      // Load Receptor
      if (receptorData && receptorData.trim().length > 0) {{
        recModel = glviewer.addModel(receptorData, receptorExt);
        recModel.setStyle({{}}, {{
          cartoon: {{ color: 'spectrum', opacity: 0.9 }},
        }});

        // Highlight target residues as sticks
        if (targetResidues && targetResidues.length > 0) {{
          targetResidues.forEach(resStr => {{
            const parts = resStr.split(':');
            const resName = parts[0];
            const resNum = parts.length > 1 ? parseInt(parts[1]) : null;
            const sel = resNum ? {{ resi: resNum }} : {{ resn: resName }};
            recModel.setStyle(sel, {{
              stick: {{ colorscheme: 'orangeCarbon', radius: 0.25 }},
              cartoon: {{ color: '#ff9f43' }}
            }});
          }});
        }}
      }}

      // Populate Poses Dropdown
      const sel = document.getElementById('pose-select');
      if (posesData.length > 0) {{
        posesData.forEach((p, idx) => {{
          const opt = document.createElement('option');
          opt.value = idx;
          opt.textContent = `Pose ${{p.pose}} (${{p.num_atoms}} atoms)`;
          sel.appendChild(opt);
        }});
        loadPose(0);
      }} else {{
        const opt = document.createElement('option');
        opt.textContent = "No ligand loaded";
        sel.appendChild(opt);
        document.getElementById('interaction-list').innerHTML = "<em>No docked ligand loaded.</em>";
      }}

      // Draw Grid Box Wireframe
      if (gridCenter && gridSize) {{
        drawGridBox(gridCenter, gridSize);
      }}

      glviewer.zoomTo();
      glviewer.render();
    }});

    function drawGridBox(center, size) {{
      const minX = center.x - size.x / 2.0;
      const maxX = center.x + size.x / 2.0;
      const minY = center.y - size.y / 2.0;
      const maxY = center.y + size.y / 2.0;
      const minZ = center.z - size.z / 2.0;
      const maxZ = center.z + size.z / 2.0;

      const corners = [
        {{x: minX, y: minY, z: minZ}},
        {{x: maxX, y: minY, z: minZ}},
        {{x: maxX, y: maxY, z: minZ}},
        {{x: minX, y: maxY, z: minZ}},
        {{x: minX, y: minY, z: maxZ}},
        {{x: maxX, y: minY, z: maxZ}},
        {{x: maxX, y: maxY, z: maxZ}},
        {{x: minX, y: maxY, z: maxZ}},
      ];

      const edges = [
        [0, 1], [1, 2], [2, 3], [3, 0], // bottom face
        [4, 5], [5, 6], [6, 7], [7, 4], // top face
        [0, 4], [1, 5], [2, 6], [3, 7]  // pillars
      ];

      edges.forEach(([i, j]) => {{
        const cyl = glviewer.addCylinder({{
          start: corners[i],
          end: corners[j],
          radius: 0.12,
          color: '#00ffff',
          opacity: 0.85
        }});
        boxShapes.push(cyl);
      }});
    }}

    function loadPose(index) {{
      currentPoseIndex = index;
      if (ligModel) {{
        glviewer.removeModel(ligModel);
        ligModel = null;
      }}
      clearHBonds();

      if (!posesData[index]) return;

      const pose = posesData[index];
      ligModel = glviewer.addModel(pose.pdb, 'pdb');
      ligModel.setStyle({{}}, {{
        stick: {{ colorscheme: 'greenCarbon', radius: 0.35 }}
      }});

      // Draw H-Bonds
      const hbs = hbondsData[index] || [];
      const listDiv = document.getElementById('interaction-list');
      listDiv.innerHTML = "";

      if (hbs.length === 0) {{
        listDiv.innerHTML = "<em>No hydrogen bonds detected within 3.5 Å.</em>";
      }} else {{
        hbs.forEach(hb => {{
          const row = document.createElement('div');
          row.className = 'hb-row';
          row.innerHTML = `<span><strong>${{hb.res}}</strong> (${{hb.rec_atom}} ··· ${{hb.lig_atom}})</span><span>${{hb.dist}} Å</span>`;
          listDiv.appendChild(row);

          if (showHBonds) {{
            const line = glviewer.addLine({{
              start: {{ x: hb.rec_x, y: hb.rec_y, z: hb.rec_z }},
              end: {{ x: hb.lig_x, y: hb.lig_y, z: hb.lig_z }},
              color: '#ffd32a',
              dashed: true,
              dashLength: 0.3,
              gapLength: 0.2
            }});
            hbondLines.push(line);
          }}
        }});
      }}

      glviewer.render();
    }}

    function clearHBonds() {{
      hbondLines.forEach(l => glviewer.removeShape(l));
      hbondLines = [];
    }}

    function onPoseChange(idx) {{
      loadPose(parseInt(idx));
    }}

    function toggleSurface() {{
      showSurface = !showSurface;
      const btn = document.getElementById('btn-surface');
      if (showSurface) {{
        btn.classList.add('active');
        if (!surfMesh && recModel) {{
          surfMesh = glviewer.addSurface($3Dmol.SurfaceType.VDW, {{
            opacity: 0.35,
            color: '#58a6ff'
          }}, {{}});
        }}
      }} else {{
        btn.classList.remove('active');
        if (surfMesh) {{
          glviewer.removeSurface(surfMesh);
          surfMesh = null;
        }}
      }}
      glviewer.render();
    }}

    function toggleGridBox() {{
      showGrid = !showGrid;
      const btn = document.getElementById('btn-box');
      boxShapes.forEach(shape => shape.hidden = !showGrid);
      if (showGrid) btn.classList.remove('active');
      else btn.classList.add('active');
      glviewer.render();
    }}

    function toggleHBonds() {{
      showHBonds = !showHBonds;
      const btn = document.getElementById('btn-hbonds');
      if (showHBonds) {{
        btn.classList.remove('active');
        loadPose(currentPoseIndex);
      }} else {{
        btn.classList.add('active');
        clearHBonds();
        glviewer.render();
      }}
    }}

    function resetCamera() {{
      glviewer.zoomTo();
      glviewer.render();
    }}

    function downloadPMLScript() {{
      let pml = `# PyMOL Script for {rec_name} :: {lig_name}\\n`;
      pml += `reinitialize\\n`;
      pml += `set ray_shadows, 0\\nset cartoon_fancy_helices, 1\\nset cartoon_side_chain_helper, 1\\n`;
      pml += `set stick_radius, 0.22\\nset stick_ball, on\\nset stick_ball_ratio, 1.35\\n`;
      pml += `set dash_gap, 0.22\\nset dash_width, 3.0\\nset label_size, 14\\nbg_color black\\n\\n`;
      pml += `# Load structures\\n`;
      pml += `load "{rec_name}.{receptor_ext}", receptor\\n`;
      pml += `show cartoon, receptor\\ncolor slate, receptor and elem C\\nutil.cbc receptor\\n\\n`;
      pml += `load "{lig_name}.pdbqt", ligand_poses\\n`;
      pml += `create ligand_pose_1, ligand_poses, 1, 1\\n`;
      pml += `show sticks, ligand_pose_1\\nshow spheres, ligand_pose_1\\nset sphere_scale, 0.24, ligand_pose_1\\n`;
      pml += `color limegreen, ligand_pose_1 and elem C\\nutil.cnc ligand_pose_1\\n\\n`;
      pml += `select pocket_residues, byres (receptor within 4.5 of ligand_pose_1)\\n`;
      pml += `show sticks, pocket_residues\\ncolor grey80, pocket_residues and elem C\\nutil.cnc pocket_residues\\n`;
      pml += `label pocket_residues and name CA, '  %s %s' % (resn, resi)\\n\\n`;
      pml += `center ligand_pose_1\\norient ligand_pose_1 or pocket_residues\\nzoom ligand_pose_1 or pocket_residues, buffer=4.0\\n`;

      const blob = new Blob([pml], {{ type: 'text/plain' }});
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = `{rec_name}_{lig_name}_docking.pml`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(a.href);
    }}
  </script>
</body>
</html>
"""
    return html


def launch_3d_viewer(
    receptor_path: Path | str,
    ligand_path: Optional[Path | str] = None,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    target_residues: Optional[List[str]] = None,
    output_html: Optional[Path] = None,
) -> Path:
    """Generates the interactive 3D WebGL page and launches it in the user's default browser."""
    html = generate_3d_viewer_html(
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        grid_center=grid_center,
        grid_size=grid_size,
        target_residues=target_residues,
    )

    if not output_html:
        temp_dir = Path(tempfile.gettempdir()) / "autodock_3d_viewer"
        temp_dir.mkdir(parents=True, exist_ok=True)
        rec_stem = Path(receptor_path).stem
        lig_stem = Path(ligand_path).stem if ligand_path else "grid"
        output_html = temp_dir / f"viewer_{rec_stem}_{lig_stem}.html"

    output_html.write_text(html, encoding="utf-8")
    logger.info(f"Generated 3D WebGL viewer at {output_html}")
    webbrowser.open(f"file:///{output_html.resolve().as_posix()}")
    return output_html


def launch_pymol_viewer(
    receptor_path: Path | str,
    ligand_path: Path | str,
    pose_index: int = 1,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    output_dir: Optional[Path | str] = None,
    custom_pymol_path: Optional[Path | str] = None,
) -> bool:
    """Exports a complete session and launches the native PyMOL desktop application."""
    from pymol_exporter import export_pymol_session, launch_pymol
    bundle = export_pymol_session(
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        output_dir=output_dir,
        pose_index=pose_index,
        grid_center=grid_center,
        grid_size=grid_size,
        compile_pse=True,
        custom_pymol_path=custom_pymol_path,
    )
    target = bundle.get("pse") or bundle.get("pml")
    return launch_pymol(target, custom_pymol_path=custom_pymol_path)


def export_pymol_bundle(
    receptor_path: Path | str,
    ligand_path: Path | str,
    output_dir: Optional[Path | str] = None,
    pose_index: int = 1,
    grid_center: Optional[Tuple[float, float, float]] = None,
    grid_size: Optional[Tuple[float, float, float]] = None,
    compile_pse: bool = True,
    custom_pymol_path: Optional[Path | str] = None,
) -> Dict[str, Any]:
    """Generates a complete PyMOL .pml script, .pse session, and clean coordinates bundle."""
    from pymol_exporter import export_pymol_session
    return export_pymol_session(
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        output_dir=output_dir,
        pose_index=pose_index,
        grid_center=grid_center,
        grid_size=grid_size,
        compile_pse=compile_pse,
        custom_pymol_path=custom_pymol_path,
    )
