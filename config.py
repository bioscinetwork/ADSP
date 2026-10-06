#!/usr/bin/env python3
"""
AutoDock Suite Pro — Configuration Module
==========================================
Loads TOML configuration, validates it, and generates template configs.

Folder Convention (Option A):
    receptors/
        RECEPTOR_NAME/
            rigid/          ← Rigid docking PDBQT goes here
            flex/           ← Flexible docking PDBQTs go here (optional)
            config.txt      ← Vina grid config (or auto-generated)
    ligands/
        LIGAND_NAME.pdbqt   ← Prepared ligand PDBQT files
"""

from __future__ import annotations

import os
import shutil
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from executables import resolve_executable
from models import DockingMode, Engine


# ═══════════════════════════════════════════════════════════════════════════════
# Configuration Dataclass
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ProjectConfig:
    """Parsed and validated project configuration."""

    # Project metadata
    project_name: str = ""
    project_root: Path = field(default_factory=Path)

    # Executables
    vina_executable: Path = field(default_factory=Path)
    vina_split_executable: Path = field(default_factory=Path)
    autogrid4_executable: Path = field(default_factory=Path)
    autodock4_executable: Path = field(default_factory=Path)
    obabel_executable: Path = field(default_factory=Path)

    # Input directories
    receptor_directory: Path = field(default_factory=Path)
    ligand_directory: Path = field(default_factory=Path)

    # Output directories
    result_directory: Path = field(default_factory=Path)
    log_directory: Path = field(default_factory=Path)
    report_directory: Path = field(default_factory=Path)

    # Engine and docking mode
    engine: Engine = Engine.VINA
    docking_mode: DockingMode = DockingMode.RIGID   # RIGID or FLEXIBLE

    # Receptor/ligand selection (empty = all)
    selected_receptors: List[str] = field(default_factory=list)
    selected_ligands: List[str] = field(default_factory=list)

    # Prepare options
    auto_generate_config: bool = False  # Auto-create config.txt if missing

    # Vina parameters
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 6.0
    cpu: str = "AUTO"
    seed: str = "AUTO"

    # Retry policy (Vina only)
    retry_enabled: bool = False
    retry_exhaustiveness: List[int] = field(default_factory=list)

    # Execution flags
    sequential: bool = True
    max_workers: int = 1
    resume: bool = True
    dry_run: bool = False

    # Reporting
    export_csv: bool = True
    export_excel: bool = True

    # AutoDock4 specific
    gpf_template: Optional[Path] = None
    dpf_template: Optional[Path] = None
    ad4_parameter_profile: str = "ad4_standard_4.2"
    reuse_existing_maps: bool = False
    ga_run: int = 100
    ga_pop_size: int = 150
    ga_num_evals: int = 2500000
    ga_num_generations: int = 27000
    ad4_num_modes: int = 100
    rmstol: float = 2.0
    ad4_seed: str = "AUTO"
    ad4_algorithm: str = "LGA"

    # AutoDock4 GA fine-grained parameters (exposed in GUI)
    ga_elitism: int = 1
    ga_mutation_rate: float = 0.02
    ga_crossover_rate: float = 0.80
    ga_window_size: int = 10
    ga_cauchy_alpha: float = 0.0
    ga_cauchy_beta: float = 1.0

    # AutoDock4 Solis-Wets local search parameters (LGA/LS)
    sw_max_its: int = 300
    sw_max_succ: int = 4
    sw_max_fail: int = 4
    sw_rho: float = 1.0
    sw_lb_rho: float = 0.01
    ls_search_freq: float = 0.06
    outlev: int = 1
    unbound_model: str = "bound"

    # Vina advanced options
    vina_scoring: str = "vina"   # "vina", "vinardo", "ad4"
    vina_min_rmsd: float = 1.0

    # Cocrystallized ligand
    cocrystallized_ligand_name: str = ""
    cocrystallized_ligand_action: str = ""  # REMOVE or KEEP

    # UI Theme
    ui_theme: str = "Noir"

    def __post_init__(self) -> None:
        """Ensure all executables are auto-resolved to real paths instead of empty/dot."""
        if not self.vina_executable or str(self.vina_executable).strip() in (".", ""):
            self.vina_executable = resolve_executable(None, "vina.exe")
        if not self.vina_split_executable or str(self.vina_split_executable).strip() in (".", ""):
            self.vina_split_executable = resolve_executable(None, "vina_split.exe")
        if not self.autogrid4_executable or str(self.autogrid4_executable).strip() in (".", ""):
            self.autogrid4_executable = resolve_executable(None, "autogrid4.exe")
        if not self.autodock4_executable or str(self.autodock4_executable).strip() in (".", ""):
            self.autodock4_executable = resolve_executable(None, "autodock4.exe")
        if not self.obabel_executable or str(self.obabel_executable).strip() in (".", ""):
            self.obabel_executable = resolve_executable(None, "obabel.exe")

    @property
    def vina_path(self) -> Path:
        return self.vina_executable

    @vina_path.setter
    def vina_path(self, val: Path | str) -> None:
        self.vina_executable = Path(val) if val else Path()

    @property
    def vina_split_path(self) -> Path:
        return self.vina_split_executable

    @vina_split_path.setter
    def vina_split_path(self, val: Path | str) -> None:
        self.vina_split_executable = Path(val) if val else Path()

    @property
    def autodock4_path(self) -> Path:
        return self.autodock4_executable

    @autodock4_path.setter
    def autodock4_path(self, val: Path | str) -> None:
        self.autodock4_executable = Path(val) if val else Path()

    @property
    def autogrid4_path(self) -> Path:
        return self.autogrid4_executable

    @autogrid4_path.setter
    def autogrid4_path(self, val: Path | str) -> None:
        self.autogrid4_executable = Path(val) if val else Path()

    @property
    def obabel_path(self) -> Path:
        return self.obabel_executable

    @obabel_path.setter
    def obabel_path(self, val: Path | str) -> None:
        self.obabel_executable = Path(val) if val else Path()


# ═══════════════════════════════════════════════════════════════════════════════
# Config Loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_config(config_path: Path | str) -> ProjectConfig:
    """Load a TOML configuration file and return a ProjectConfig."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, "rb") as f:
        data = tomllib.load(f)

    config = ProjectConfig()

    # Project section
    proj = data.get("project", {})
    config.project_name = proj.get("name", config_path.stem)
    root_str = proj.get("root", "")
    if not root_str or root_str == ".":
        config.project_root = config_path.parent.resolve()
    else:
        p_root = Path(root_str)
        if not p_root.is_absolute():
            config.project_root = (config_path.parent / p_root).resolve()
        else:
            config.project_root = p_root.resolve()

    # Executables section (prioritizes local workspace/suite bin/ before system paths)
    exe = data.get("executables", {})

    def _resolve_exe(configured_val: Optional[str], default_name: str) -> Path:
        val = configured_val or f"bin/{default_name}"
        if val and not Path(val).is_absolute():
            cand = config.project_root / val
            if cand.is_file():
                return cand.resolve()
        return resolve_executable(val, default_name)

    config.vina_executable = _resolve_exe(exe.get("vina"), "vina.exe")
    config.vina_split_executable = _resolve_exe(exe.get("vina_split"), "vina_split.exe")
    config.autogrid4_executable = _resolve_exe(exe.get("autogrid4"), "autogrid4.exe")
    config.autodock4_executable = _resolve_exe(exe.get("autodock4"), "autodock4.exe")
    config.obabel_executable = _resolve_exe(exe.get("obabel"), "obabel.exe")

    # Input directories
    inputs = data.get("inputs", {})
    config.receptor_directory = _resolve_input_dir(
        inputs.get("receptor_directory", "receptors"),
        config.project_root,
        ["data/receptors", "data/Macromolecules", "Macromolecules", "macromolecules", "receptors", "Receptors"],
    )
    config.ligand_directory = _resolve_input_dir(
        inputs.get("ligand_directory", "ligands"),
        config.project_root,
        ["data/ligands", "data/compounds", "Ligands", "ligands", "compounds", "Compounds"],
    )

    # Receptor/ligand selection
    config.selected_receptors = inputs.get("selected_receptors", [])
    config.selected_ligands = inputs.get("selected_ligands", [])

    # Output directories
    outputs = data.get("outputs", {})
    config.result_directory = _resolve_path(
        outputs.get("result_directory", "results"), config.project_root
    )
    config.log_directory = _resolve_path(
        outputs.get("log_directory", "logs"), config.project_root
    )
    config.report_directory = _resolve_path(
        outputs.get("report_directory", "reports"), config.project_root
    )

    # Engine
    engine_str = data.get("engine", {}).get("type", "VINA").upper()
    config.engine = Engine(engine_str)

    # Docking mode
    docking = data.get("docking", {})
    mode_str = docking.get("mode", "AUTO").upper()
    try:
        config.docking_mode = DockingMode(mode_str)
    except ValueError:
        config.docking_mode = DockingMode.AUTO

    # Prepare options
    prepare = data.get("prepare", {})
    config.auto_generate_config = prepare.get("auto_generate_config", False)

    # Vina parameters
    vina = data.get("vina", {})
    config.exhaustiveness = vina.get("exhaustiveness", 8)
    config.num_modes = vina.get("num_modes", 9)
    config.energy_range = vina.get("energy_range", 6.0)
    config.cpu = str(vina.get("cpu", "AUTO")).upper()
    config.seed = str(vina.get("seed", "AUTO")).upper()
    config.vina_scoring = str(vina.get("scoring", "vina")).lower()
    config.vina_min_rmsd = float(vina.get("min_rmsd", 1.0))

    # Retry
    retry = data.get("retry", {})
    config.retry_enabled = retry.get("enabled", False)
    config.retry_exhaustiveness = retry.get("exhaustiveness_values", [])

    # Execution
    execution = data.get("execution", {})
    config.sequential = execution.get("sequential", True)
    mw = execution.get("max_workers", 1)
    if str(mw).upper() == "AUTO":
        import os
        config.max_workers = max(1, (os.cpu_count() or 2) // 2)
    else:
        try:
            config.max_workers = max(1, int(mw))
        except ValueError:
            config.max_workers = 1
    if config.max_workers > 1:
        config.sequential = False
    config.resume = execution.get("resume", True)
    config.dry_run = execution.get("dry_run", False)

    # Reporting
    reporting = data.get("reporting", {})
    config.export_csv = reporting.get("csv", True)
    config.export_excel = reporting.get("excel", True)

    # AutoDock4 parameters
    ad4 = data.get("autodock4", {})
    config.ad4_parameter_profile = str(
        ad4.get("parameter_profile", "ad4_standard_4.2")
    ).strip()
    if ad4.get("gpf_template"):
        config.gpf_template = _resolve_path(ad4["gpf_template"], config.project_root)
    if ad4.get("dpf_template"):
        config.dpf_template = _resolve_path(ad4["dpf_template"], config.project_root)
    config.reuse_existing_maps = ad4.get("reuse_existing_maps", False)
    config.ga_run = ad4.get("ga_run", 100)
    config.ga_pop_size = ad4.get("ga_pop_size", 150)
    config.ga_num_evals = ad4.get("ga_num_evals", 2500000)
    config.ga_num_generations = ad4.get("ga_num_generations", 27000)
    config.ad4_num_modes = ad4.get("num_modes", 100)
    config.rmstol = ad4.get("rmstol", 2.0)
    config.ad4_seed = str(ad4.get("seed", "AUTO")).upper()
    config.ad4_algorithm = str(ad4.get("algorithm", "LGA")).upper()
    config.ga_elitism = int(ad4.get("ga_elitism", 1))
    config.ga_mutation_rate = float(ad4.get("ga_mutation_rate", 0.02))
    config.ga_crossover_rate = float(ad4.get("ga_crossover_rate", 0.80))
    config.ga_window_size = int(ad4.get("ga_window_size", 10))
    config.ga_cauchy_alpha = float(ad4.get("ga_cauchy_alpha", 0.0))
    config.ga_cauchy_beta = float(ad4.get("ga_cauchy_beta", 1.0))
    config.sw_max_its = int(ad4.get("sw_max_its", 300))
    config.sw_max_succ = int(ad4.get("sw_max_succ", 4))
    config.sw_max_fail = int(ad4.get("sw_max_fail", 4))
    config.sw_rho = float(ad4.get("sw_rho", 1.0))
    config.sw_lb_rho = float(ad4.get("sw_lb_rho", 0.01))
    config.ls_search_freq = float(ad4.get("ls_search_freq", 0.06))
    config.outlev = int(ad4.get("outlev", 1))
    config.unbound_model = str(ad4.get("unbound_model", "bound")).lower()

    # Cocrystallized ligand
    cocryst = data.get("cocrystallized_ligand", {})
    config.cocrystallized_ligand_name = cocryst.get("name", "")
    config.cocrystallized_ligand_action = cocryst.get("action", "").upper()

    # UI Theme
    ui_sec = data.get("ui", {})
    config.ui_theme = str(ui_sec.get("theme", data.get("theme", "Noir")))

    return config


def _resolve_path(path_str: str, root: Path) -> Path:
    """Resolve a path relative to the project root if not absolute."""
    p = Path(path_str)
    if p.is_absolute():
        return p.resolve()
    return (root.resolve() / p).resolve()


def set_workspace_directory(config: ProjectConfig, workspace_dir: Path | str) -> None:
    """Set and activate a new workspace directory for the project.

    1. Resolves workspace_dir to an absolute Path.
    2. If project_config.toml exists in workspace_dir, loads it and applies its settings to config.
    3. Updates project_root, project_name, receptor_directory, ligand_directory,
       result_directory, log_directory, and report_directory to be inside workspace_dir.
    4. Creates standard subdirectories (receptors, ligands, results, results/DLG, results/DPF,
       results/VINA, results/BSNDVP_RESULTS, DLG, BSNDVP_RESULTS, logs, reports).
    5. Copies dlg_extract.py into the workspace if not already present.
    6. Switches Python's current working directory (os.chdir) to workspace_dir.
    """
    ws = Path(workspace_dir).resolve()
    ws.mkdir(parents=True, exist_ok=True)

    toml_path = ws / "project_config.toml"
    if toml_path.is_file():
        try:
            loaded = load_config(toml_path)
            for f_name in loaded.__dataclass_fields__:
                val = getattr(loaded, f_name)
                if f_name.endswith("_executable") and (not val or str(val).strip() in (".", "")):
                    continue
                setattr(config, f_name, val)
        except Exception:
            pass

    config.project_root = ws
    config.project_name = ws.name

    config.receptor_directory = _resolve_input_dir(
        "receptors", ws, ["data/receptors", "data/Macromolecules", "Macromolecules", "macromolecules", "receptors", "Receptors"]
    )
    config.ligand_directory = _resolve_input_dir(
        "ligands", ws, ["data/ligands", "data/compounds", "Ligands", "ligands", "compounds", "Compounds"]
    )

    config.result_directory = (ws / "results").resolve()
    config.log_directory = (ws / "logs").resolve()
    config.report_directory = (ws / "reports").resolve()

    for d in [
        config.receptor_directory,
        config.ligand_directory,
        config.result_directory,
        config.result_directory / "DLG",
        config.result_directory / "DPF",
        config.result_directory / "VINA",
        config.result_directory / "BSNDVP_RESULTS",
        ws / "DLG",
        ws / "BSNDVP_RESULTS",
        config.log_directory,
        config.report_directory,
    ]:
        try:
            d.mkdir(parents=True, exist_ok=True)
        except Exception:
            pass

    # Copy dlg_extract.py into workspace if missing
    suite_dlg_extract = Path(__file__).resolve().parent / "dlg_extract.py"
    target_dlg_extract = ws / "dlg_extract.py"
    if suite_dlg_extract.is_file() and not target_dlg_extract.is_file() and ws != suite_dlg_extract.parent.resolve():
        try:
            shutil.copy2(suite_dlg_extract, target_dlg_extract)
        except Exception:
            pass

    # Switch process working directory
    try:
        os.chdir(ws)
    except Exception:
        pass


def _resolve_input_dir(path_str: str, root: Path, candidates: List[str]) -> Path:
    """Resolve input directory with multi-location auto-discovery.

    Checks:
    1. Absolute path if specified
    2. root / path_str
    3. root / .. / path_str (parent workspace)
    4. Alternative names (e.g. Macromolecules, Ligands) in root and parent workspace
    """
    p = Path(path_str)
    if p.is_absolute():
        try:
            if p.exists():
                return p
        except OSError:
            pass

    # Check root / path_str
    try:
        if (root / p).exists():
            return root / p
    except OSError:
        pass

    # Check parent of root
    try:
        parent_cand = root / ".." / p
        if parent_cand.exists():
            return parent_cand.resolve()
    except OSError:
        pass

    # Check candidates in root and parent
    for cand in candidates:
        try:
            cand_path = root / cand
            if cand_path.exists():
                return cand_path
        except OSError:
            pass
        try:
            cand_parent = root / ".." / cand
            if cand_parent.exists():
                return cand_parent.resolve()
        except OSError:
            pass

    return root / p


# ═══════════════════════════════════════════════════════════════════════════════
# Config Validation
# ═══════════════════════════════════════════════════════════════════════════════

def validate_config(config: ProjectConfig) -> List[str]:
    """Validate configuration and return a list of error messages.

    Returns an empty list if configuration is valid.
    """
    errors: List[str] = []

    # Validate project root
    try:
        if not config.project_root.exists():
            errors.append(f"Project root does not exist: {config.project_root}")
    except OSError as e:
        errors.append(f"Project root is an invalid path: {config.project_root} ({e})")

    # Validate input directories
    try:
        if not config.receptor_directory.exists():
            errors.append(f"Receptor directory does not exist: {config.receptor_directory}")
    except OSError as e:
        errors.append(f"Receptor directory path is invalid: {config.receptor_directory} ({e})")

    try:
        if not config.ligand_directory.exists():
            errors.append(f"Ligand directory does not exist: {config.ligand_directory}")
    except OSError as e:
        errors.append(f"Ligand directory path is invalid: {config.ligand_directory} ({e})")

    # Validate engine-specific executables
    if config.engine in (Engine.VINA, Engine.BOTH):
        if not config.vina_executable or not config.vina_executable.exists():
            errors.append(
                f"Vina executable not found: {config.vina_executable}\n"
                f"  → Set 'vina' in [executables] to the correct path."
            )
        if not config.vina_split_executable or not config.vina_split_executable.exists():
            errors.append(
                f"vina_split executable not found: {config.vina_split_executable}\n"
                f"  → Set 'vina_split' in [executables] to the correct path."
            )
    if config.engine in (Engine.AUTODOCK4, Engine.BOTH):
        from ad4_compatibility import PARAMETER_PROFILES, validate_profile_id
        if config.ad4_parameter_profile not in PARAMETER_PROFILES:
            errors.append(
                f"Unknown AutoDock4 parameter profile: {config.ad4_parameter_profile}. "
                f"Available profiles: {', '.join(sorted(PARAMETER_PROFILES))}"
            )
        else:
            profile_result = validate_profile_id(config.ad4_parameter_profile)
            if not profile_result.valid:
                errors.append(
                    f"Invalid AutoDock4 parameter profile '{config.ad4_parameter_profile}': "
                    + "; ".join(profile_result.errors)
                )
        if not config.autogrid4_executable or not config.autogrid4_executable.exists():
            errors.append(
                f"AutoGrid4 executable not found: {config.autogrid4_executable}\n"
                f"  → Set 'autogrid4' in [executables] to the correct path."
            )
        if not config.autodock4_executable or not config.autodock4_executable.exists():
            errors.append(
                f"AutoDock4 executable not found: {config.autodock4_executable}\n"
                f"  → Set 'autodock4' in [executables] to the correct path."
            )

    # Validate Vina parameters
    if config.engine in (Engine.VINA, Engine.BOTH):
        if config.exhaustiveness < 1:
            errors.append(f"exhaustiveness must be >= 1, got {config.exhaustiveness}")
        if config.num_modes < 1:
            errors.append(f"num_modes must be >= 1, got {config.num_modes}")
        if config.energy_range <= 0:
            errors.append(f"energy_range must be > 0, got {config.energy_range}")

    # Validate retry settings
    if config.retry_enabled and not config.retry_exhaustiveness:
        errors.append("Retry is enabled but no retry exhaustiveness values specified.")

    return errors


# ═══════════════════════════════════════════════════════════════════════════════
# Template Config Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_template_config(output_path: Path) -> None:
    """Write an example TOML configuration file with the new folder structure."""
    template = '''# ============================================================
# AUTODOCK SUITE PRO — PROJECT CONFIGURATION v0.2.0
# ============================================================
# Edit this file to match your project setup.
# Lines starting with # are comments.
# Paths must use forward slashes or escaped backslashes.
#
# FOLDER STRUCTURE EXPECTED:
# ─────────────────────────────────────────────────────────
#  receptors/
#    RECEPTOR_NAME/
#      rigid/                    ← Place your rigid-receptor PDBQT here
#        receptor.pdbqt
#      flex/                     ← (Optional) Flexible docking files
#        receptor_rigid.pdbqt    ← Rigid part (protein backbone)
#        receptor_flex.pdbqt     ← Flexible part (moveable sidechains)
#      config.txt                ← Vina grid config (or auto-generated)
#  ligands/
#    compound_001.pdbqt          ← Prepared ligand PDBQTs (any names)
#    compound_002.pdbqt
#
# Run "python main.py --scaffold" to create this structure automatically.
# Run "python main.py --prepare" to auto-generate config.txt / GPF files.
# ─────────────────────────────────────────────────────────

[project]
name = "Docking_Project"
root = "."

# ── Executables ──────────────────────────────────────────────
[executables]
# AutoDock Vina & Vina Split
vina       = "bin/vina.exe"
vina_split = "bin/vina_split.exe"

# AutoDock 4 & AutoGrid 4
autogrid4  = "bin/autogrid4.exe"
autodock4  = "bin/autodock4.exe"
obabel     = "bin/obabel.exe"

# ── Input Directories ────────────────────────────────────────
[inputs]
# Directories can be absolute or relative to project root
receptor_directory = "receptors"
ligand_directory   = "ligands"

# Optional: select specific receptors/ligands (empty = all)
# selected_receptors = ["receptor_01", "receptor_02"]
# selected_ligands   = ["ligand_001", "ligand_005"]

# ── Output Directories ───────────────────────────────────────
[outputs]
result_directory = "results"
log_directory    = "logs"
report_directory = "reports"

# ── Engine Selection ─────────────────────────────────────────
[engine]
# Options: "VINA" or "AUTODOCK4"
type = "VINA"

# ── Docking Mode ─────────────────────────────────────────────
[docking]
# Options: "AUTO", "RIGID", "FLEXIBLE", "BOTH"
#
# AUTO     → Auto-detects from folders:
#            - rigid/ only        → RIGID docking
#            - flex/ only         → FLEXIBLE docking
#            - BOTH folders exist → Runs BOTH rigid and flexible docking!
# RIGID    → Force rigid docking only (ignores flex/ folder)
# FLEXIBLE → Force flexible docking only (ignores rigid/ folder)
# BOTH     → Explicitly queue both rigid and flexible docking
mode = "AUTO"

# ── Preparation Options ──────────────────────────────────────
[prepare]
# If true, a missing config.txt is auto-generated using grid_center
# and grid_size values below (Vina only). Requires center/size set.
auto_generate_config = false

# Default grid center (Angstroms) — used when auto_generate_config = true
# grid_center_x = 0.0
# grid_center_y = 0.0
# grid_center_z = 0.0

# Default grid box size (Angstroms)
# grid_size_x = 25.0
# grid_size_y = 25.0
# grid_size_z = 25.0

# ── Vina Parameters ──────────────────────────────────────────
[vina]
exhaustiveness = 8
num_modes      = 9        # Maximum requested poses (actual may be fewer)
energy_range   = 6.0      # kcal/mol
cpu            = "AUTO"   # Number of CPUs or "AUTO"
seed           = "AUTO"   # Random seed or "AUTO"

# ── Retry Policy (Vina only) ─────────────────────────────────
[retry]
enabled = false
# exhaustiveness_values = [16, 32]

# ── AutoDock4 Parameters ─────────────────────────────────────
[autodock4]
# gpf_template = "config/template.gpf"
# dpf_template = "config/template.dpf"
parameter_profile = "ad4_standard_4.2"  # Standard atom-type compatibility profile; no specialized metal workflow
reuse_existing_maps  = false
ga_run               = 100
ga_pop_size          = 150
ga_num_evals         = 2500000
ga_num_generations   = 27000
num_modes            = 100
rmstol               = 2.0
seed                 = "AUTO"

# ── Cocrystallized Ligand ────────────────────────────────────
[cocrystallized_ligand]
# name   = "<LIGAND_RESNAME>"
# action = "REMOVE"   # or "KEEP"

# ── Execution Settings ───────────────────────────────────────
[execution]
sequential  = true    # Process jobs one at a time (set false for parallel)
max_workers = 1       # Number of concurrent workers (1 = sequential, 2-16 or "AUTO")
resume      = true    # Skip previously completed jobs
dry_run     = false   # Validate without executing

# ── Reporting ─────────────────────────────────────────────────
[reporting]
csv   = true
excel = true
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(template, encoding="utf-8")
