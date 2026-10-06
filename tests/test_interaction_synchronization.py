"""
Regression Test Suite for AutoDockSuite Pro Interaction Subsystem
Proves:
1. Exact pose selection (Pose 1 vs Pose 2)
2. No silent fallback to Pose 1 on invalid/missing pose
3. True pose coordinate extraction
4. Pose isolation in 2D diagram generation
5. Complete outer SVG export without regex truncation
6. Self-contained offline HTML export
7. Stale asynchronous worker result rejection
8. Sequential pose switching (Pose 1 -> 2 -> 3 -> 1)
"""

import os
import sys
import tempfile
from pathlib import Path
import pytest

import interactions
from interactions import (
    PoseLookupStatus,
    PoseNotFoundError,
    locate_docked_pose,
    extract_pose_pdbqt_block,
    profile_docking_pose,
    profile_docking_job,
    Interaction,
)
from gui.interaction_diagram import (
    render_interaction_diagram,
    render_interaction_diagram_html,
    render_interaction_diagram_svg,
    load_pose_molecule,
    InteractionDiagramResult,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

RECEPTOR_PDBQT = """ATOM      1  N   ASP A 102      10.000  10.000   8.000  1.00 20.00    -0.200 N 
ATOM      2  CA  ASP A 102      10.000  10.000   9.000  1.00 20.00     0.100 C 
ATOM      3  CB  ASP A 102      10.000  10.000   9.500  1.00 20.00     0.000 C 
ATOM      4  CG  ASP A 102      10.000  10.000   9.800  1.00 20.00     0.400 C 
ATOM      5  OD1 ASP A 102      10.000  10.000  10.000  1.00 20.00    -0.600 OA
ATOM      6  N   SER A 195      30.000  30.000  28.000  1.00 20.00    -0.200 N 
ATOM      7  CA  SER A 195      30.000  30.000  29.000  1.00 20.00     0.100 C 
ATOM      8  CB  SER A 195      30.000  30.000  29.500  1.00 20.00     0.000 C 
ATOM      9  OG  SER A 195      30.000  30.000  30.000  1.00 20.00    -0.600 OA
"""

# Multi-pose ligand:
# Model 1 has atom O1 close to ASP 102 OD1 (distance ~2.8 A)
# Model 2 has atom O1 close to SER 195 OG (distance ~2.8 A), far from ASP 102
# Model 3 has atom O1 at (50, 50, 50), far from everything
MULTI_POSE_LIGAND_PDBQT = """MODEL 1
REMARK VINA RESULT:      -8.5      0.000      0.000
ATOM      1  O1  LIG L   1      10.000  10.000  12.800  1.00 20.00    -0.400 OA
ATOM      2  C2  LIG L   1      10.000  10.000  14.200  1.00 20.00     0.100 C 
ATOM      3  C3  LIG L   1      10.000  10.000  15.500  1.00 20.00     0.000 C 
ENDMDL
MODEL 2
REMARK VINA RESULT:      -7.2      2.100      3.400
ATOM      1  O1  LIG L   1      30.000  30.000  32.800  1.00 20.00    -0.400 OA
ATOM      2  C2  LIG L   1      30.000  30.000  34.200  1.00 20.00     0.100 C 
ATOM      3  C3  LIG L   1      30.000  30.000  35.500  1.00 20.00     0.000 C 
ENDMDL
MODEL 3
REMARK VINA RESULT:      -6.1      4.500      6.200
ATOM      1  O1  LIG L   1      50.000  50.000  52.800  1.00 20.00    -0.400 OA
ATOM      2  C2  LIG L   1      50.000  50.000  54.200  1.00 20.00     0.100 C 
ATOM      3  C3  LIG L   1      50.000  50.000  55.500  1.00 20.00     0.000 C 
ENDMDL
"""


@pytest.fixture
def docking_fixture(tmp_path):
    rec_file = tmp_path / "receptor.pdbqt"
    rec_file.write_text(RECEPTOR_PDBQT, encoding="utf-8")
    lig_file = tmp_path / "ligand_out.pdbqt"
    lig_file.write_text(MULTI_POSE_LIGAND_PDBQT, encoding="utf-8")
    return rec_file, lig_file


# ---------------------------------------------------------------------------
# Test 1: Load multi-pose output. Request Pose 1. Verify pose_index == 1.
# ---------------------------------------------------------------------------

def test_1_pose_1_selection(docking_fixture):
    rec_path, lig_path = docking_fixture
    inters = profile_docking_pose(rec_path, lig_path, pose_index=1)

    assert len(inters) > 0, "Pose 1 should find contacts"
    assert all(it.pose_index == 1 for it in inters), "All interactions must have pose_index == 1"
    
    # Pose 1 is located near ASP 102
    asp_contacts = [it for it in inters if it.receptor_res_seq == 102 and it.receptor_res_name == "ASP"]
    assert len(asp_contacts) > 0, "Pose 1 must form contact with ASP 102"
    assert pytest.approx(asp_contacts[0].distance_angstrom, abs=0.1) == 2.80


# ---------------------------------------------------------------------------
# Test 2: Request Pose 2. Verify pose_index == 2 and not same coords as Pose 1.
# ---------------------------------------------------------------------------

def test_2_pose_2_selection_different_coordinates(docking_fixture):
    rec_path, lig_path = docking_fixture
    inters_2 = profile_docking_pose(rec_path, lig_path, pose_index=2)

    assert len(inters_2) > 0, "Pose 2 should find contacts"
    assert all(it.pose_index == 2 for it in inters_2), "All interactions must have pose_index == 2"

    # Pose 2 is near SER 195, NOT ASP 102
    ser_contacts = [it for it in inters_2 if it.receptor_res_seq == 195 and it.receptor_res_name == "SER"]
    asp_contacts = [it for it in inters_2 if it.receptor_res_seq == 102 and it.receptor_res_name == "ASP"]

    assert len(ser_contacts) > 0, "Pose 2 must form contact with SER 195"
    assert len(asp_contacts) == 0, "Pose 2 must NOT form contact with ASP 102"

    # Verify physical coordinates extracted for Pose 1 and Pose 2 are distinct
    block_1 = extract_pose_pdbqt_block(lig_path, pose_index=1)
    block_2 = extract_pose_pdbqt_block(lig_path, pose_index=2)
    assert block_1 != block_2
    assert "10.000  10.000  12.800" in block_1
    assert "30.000  30.000  32.800" in block_2


# ---------------------------------------------------------------------------
# Test 3: Request nonexistent Pose 999. Verify POSE NOT FOUND, NOT Pose 1.
# ---------------------------------------------------------------------------

def test_3_nonexistent_pose_fails_explicitly(docking_fixture):
    rec_path, lig_path = docking_fixture

    # locate_docked_pose must report POSE_NOT_FOUND
    poses = interactions.parse_docked_poses(lig_path)
    status, match = locate_docked_pose(poses, pose_index=999)
    assert status == PoseLookupStatus.POSE_NOT_FOUND
    assert match is None

    # extract_pose_pdbqt_block with raise_on_missing must raise PoseNotFoundError
    with pytest.raises(PoseNotFoundError) as exc_info:
        extract_pose_pdbqt_block(lig_path, pose_index=999, raise_on_missing=True)
    assert exc_info.value.status == PoseLookupStatus.POSE_NOT_FOUND

    # extract_pose_pdbqt_block without raise_on_missing returns None (NEVER Pose 1 block)
    block_999 = extract_pose_pdbqt_block(lig_path, pose_index=999, raise_on_missing=False)
    assert block_999 is None

    # profile_docking_pose with raise_on_missing must raise PoseNotFoundError
    with pytest.raises(PoseNotFoundError):
        profile_docking_pose(rec_path, lig_path, pose_index=999, raise_on_missing=True)

    # profile_docking_pose without raise_on_missing returns diagnostic profile with POSE_NOT_FOUND
    prof_999 = profile_docking_pose(rec_path, lig_path, pose_index=999, raise_on_missing=False)
    assert prof_999.status == PoseLookupStatus.POSE_NOT_FOUND
    assert len(prof_999.interactions) == 0

    # profile_docking_job must NOT silently return Pose 1's interactions
    with pytest.raises(PoseNotFoundError):
        profile_docking_job(rec_path, lig_path, pose_index=999, raise_on_missing=True)

    # When raise_on_missing=False on profile_docking_job, it returns [] (NEVER Pose 1)
    res_999 = profile_docking_job(rec_path, lig_path, pose_index=999, raise_on_missing=False)
    assert res_999 == []


# ---------------------------------------------------------------------------
# Test 4: Render Pose 1 & 2. Verify diagrams contain different pose metadata
#         and different interaction content.
# ---------------------------------------------------------------------------

def test_4_diagram_rendering_pose_isolation(docking_fixture):
    rec_path, lig_path = docking_fixture
    inters_1 = profile_docking_pose(rec_path, lig_path, pose_index=1)
    inters_2 = profile_docking_pose(rec_path, lig_path, pose_index=2)

    diag_1 = render_interaction_diagram(lig_path, inters_1, pose_index=1)
    diag_2 = render_interaction_diagram(lig_path, inters_2, pose_index=2)

    # Pose index verification
    assert diag_1.pose_index == 1
    assert diag_2.pose_index == 2
    assert diag_1.metadata["pose_index"] == 1
    assert diag_2.metadata["pose_index"] == 2

    # Structural SVG metadata tags
    assert "<adsp:pose_index>1</adsp:pose_index>" in diag_1.svg
    assert "<adsp:pose_index>2</adsp:pose_index>" in diag_2.svg

    # HTML headers
    assert "Pose 1" in diag_1.html
    assert "Pose 2" in diag_2.html

    # Residue contents
    assert "ASP" in diag_1.svg
    assert "SER" in diag_2.svg
    assert diag_1.svg != diag_2.svg


# ---------------------------------------------------------------------------
# Test 5: Export Pose 2 as SVG. Verify valid, complete, and contains Pose 2.
# ---------------------------------------------------------------------------

def test_5_export_pose_2_as_svg(docking_fixture, tmp_path):
    rec_path, lig_path = docking_fixture
    inters_2 = profile_docking_pose(rec_path, lig_path, pose_index=2)
    diag_2 = render_interaction_diagram(lig_path, inters_2, pose_index=2)

    export_path = tmp_path / "test_export_pose2.svg"
    export_path.write_text(diag_2.svg, encoding="utf-8")

    assert export_path.exists()
    assert export_path.stat().st_size > 0

    content = export_path.read_text(encoding="utf-8").strip()
    assert content.startswith("<svg")
    assert content.endswith("</svg>")
    assert "<adsp:pose_index>2</adsp:pose_index>" in content
    assert "SER" in content
    assert "Pose 2" in content


# ---------------------------------------------------------------------------
# Test 6: Export Pose 2 as HTML. Verify self-contained offline document.
# ---------------------------------------------------------------------------

def test_6_export_pose_2_as_html(docking_fixture, tmp_path):
    rec_path, lig_path = docking_fixture
    inters_2 = profile_docking_pose(rec_path, lig_path, pose_index=2)
    diag_2 = render_interaction_diagram(lig_path, inters_2, pose_index=2)

    export_path = tmp_path / "test_export_pose2.html"
    export_path.write_text(diag_2.html, encoding="utf-8")

    assert export_path.exists()
    assert export_path.stat().st_size > 0

    content = export_path.read_text(encoding="utf-8").strip()
    assert "<!DOCTYPE html>" in content
    assert "Pose 2" in content
    assert "<svg" in content
    # Offline check: no remote scripts or external style sheets
    offline_check_content = content.replace("http://www.w3.org/2000/svg", "").replace("https://autodocksuite.pro/schema", "")
    assert "http://" not in offline_check_content
    assert "https://" not in offline_check_content


# ---------------------------------------------------------------------------
# Test 7: Stale asynchronous worker results cannot overwrite newer pose.
# ---------------------------------------------------------------------------

def test_7_stale_asynchronous_worker_rejection(docking_fixture):
    """
    Simulate the exact race condition:
    1. User initiates Pose 1 (generation 1).
    2. User immediately clicks Pose 2 (generation 2).
    3. Pose 1 worker finishes late and attempts to deliver Pose 1.
    4. Validate that Pose 1 is discarded and only Pose 2 is accepted.
    """
    from PySide6.QtWidgets import QApplication
    from config import ProjectConfig
    from models import DockingJob, Engine, DockingMode
    from gui_qt.results_tab import ResultsTab

    app = QApplication.instance() or QApplication(sys.argv)

    rec_path, lig_path = docking_fixture
    cfg = ProjectConfig()
    results_tab = ResultsTab(cfg, on_config_change=lambda: None)

    job = DockingJob(
        job_id="test_job_001",
        receptor_name="receptor.pdbqt",
        ligand_name="ligand_out.pdbqt",
        engine=Engine.VINA,
        docking_mode=DockingMode.RIGID,
        output_dir=str(lig_path.parent),
    )
    job.ensure_canonical_result()
    results_tab._current_job = job

    # User selects Pose 1
    results_tab._current_pose = 1
    results_tab._render_generation = 1
    gen_1 = 1

    # User immediately clicks Pose 2
    results_tab._current_pose = 2
    results_tab._render_generation = 2
    gen_2 = 2

    # Pose 1 worker finishes late and tries to emit results with gen_1
    inters_pose_1 = profile_docking_pose(rec_path, lig_path, pose_index=1)
    diag_pose_1 = render_interaction_diagram(lig_path, inters_pose_1, pose_index=1)

    # Deliver stale results for Pose 1
    results_tab._on_interactions_ready(inters_pose_1, job.job_id, pose_idx=1, gen=gen_1)
    results_tab._on_diagram_ready(diag_pose_1, job.job_id, pose_idx=1, gen=gen_1)

    # Verify stale results were rejected!
    assert results_tab._current_pose == 2
    assert results_tab._current_interactions != inters_pose_1
    assert results_tab._current_diagram_res != diag_pose_1

    # Now Pose 2 worker delivers result with gen_2
    inters_pose_2 = profile_docking_pose(rec_path, lig_path, pose_index=2)
    diag_pose_2 = render_interaction_diagram(lig_path, inters_pose_2, pose_index=2)

    results_tab._on_interactions_ready(inters_pose_2, job.job_id, pose_idx=2, gen=gen_2)
    results_tab._on_diagram_ready(diag_pose_2, job.job_id, pose_idx=2, gen=gen_2)

    # Verify Pose 2 was accepted!
    assert results_tab._current_pose == 2
    assert results_tab._current_interactions == inters_pose_2
    assert results_tab._current_diagram_res == diag_pose_2
    assert "<adsp:pose_index>2</adsp:pose_index>" in results_tab._current_svg


# ---------------------------------------------------------------------------
# Test 8: Sequential pose switching: Pose 1 -> Pose 2 -> Pose 3 -> Pose 1.
# ---------------------------------------------------------------------------

def test_8_sequential_pose_switching(docking_fixture):
    rec_path, lig_path = docking_fixture
    
    # 1. Pose 1
    p1 = profile_docking_pose(rec_path, lig_path, pose_index=1)
    d1 = render_interaction_diagram(lig_path, p1, pose_index=1)
    assert d1.pose_index == 1
    assert "ASP" in d1.svg

    # 2. Pose 2
    p2 = profile_docking_pose(rec_path, lig_path, pose_index=2)
    d2 = render_interaction_diagram(lig_path, p2, pose_index=2)
    assert d2.pose_index == 2
    assert "SER" in d2.svg

    # 3. Pose 3 (no close contacts)
    p3 = profile_docking_pose(rec_path, lig_path, pose_index=3)
    d3 = render_interaction_diagram(lig_path, p3, pose_index=3)
    assert d3.pose_index == 3
    assert d3.contact_count == 0

    # 4. Switch back to Pose 1
    p1_again = profile_docking_pose(rec_path, lig_path, pose_index=1)
    d1_again = render_interaction_diagram(lig_path, p1_again, pose_index=1)
    assert d1_again.pose_index == 1
    assert "ASP" in d1_again.svg
    assert d1_again.svg == d1.svg
