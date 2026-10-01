"""
AutoDock Suite Pro — Config Sync / Save / Load (gui/settings.py)
================================================================
Decouples all config persistence logic from the GUI tabs.
"""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict

from config import ProjectConfig, load_config
from models import DockingMode, Engine


# ---------------------------------------------------------------------------
# Sync UI → Config
# ---------------------------------------------------------------------------

def sync_config_from_ui(config: ProjectConfig, ui_refs: Dict[str, Any]) -> None:
    """Pull every user-editable field from the UI into the ProjectConfig object.

    Args:
        config:   The ProjectConfig instance to update in-place.
        ui_refs:  Dict of {field_name: ctk_widget} mapping.
    """
    def _get(name: str, default: str = "") -> str:
        w = ui_refs.get(name)
        if w is None:
            return default
        try:
            return w.get().strip()
        except Exception:
            return default

    def _path(name: str) -> Path | None:
        v = _get(name)
        return Path(v) if v else None

    # Executables
    if p := _path("vina"):
        config.vina_executable = p
    if p := _path("vina_split"):
        config.vina_split_executable = p
    if p := _path("autodock4"):
        config.autodock4_executable = p
    if p := _path("autogrid4"):
        config.autogrid4_executable = p
    if p := _path("obabel"):
        config.obabel_executable = p

    # Directories
    if p := _path("receptor_dir"):
        config.receptor_directory = p
    if p := _path("ligand_dir"):
        config.ligand_directory = p
    if p := _path("result_dir"):
        config.result_directory = p
    if p := _path("log_dir"):
        config.log_directory = p
    if p := _path("report_dir"):
        config.report_directory = p

    # Engine
    engine_val = _get("engine", "VINA")
    try:
        config.engine = Engine(engine_val.upper())
    except ValueError:
        pass

    # Docking mode
    mode_val = _get("docking_mode", "AUTO")
    try:
        config.docking_mode = DockingMode(mode_val.upper())
    except ValueError:
        pass

    # Vina parameters
    try:
        config.exhaustiveness = int(_get("exhaustiveness", "8"))
    except ValueError:
        pass
    try:
        config.num_modes = int(_get("num_modes", "9"))
    except ValueError:
        pass
    try:
        config.energy_range = float(_get("energy_range", "3.0"))
    except ValueError:
        pass
    config.seed = _get("seed", "AUTO")
    config.cpu = _get("cpu", "AUTO")
    config.vina_scoring = _get("vina_scoring", getattr(config, "vina_scoring", "vina"))
    try:
        config.vina_min_rmsd = float(_get("vina_min_rmsd", str(getattr(config, "vina_min_rmsd", 1.0))))
    except ValueError:
        pass

    # AutoDock4 parameters
    ad4_alg = _get("ad4_algorithm", getattr(config, "ad4_algorithm", "LGA"))
    config.ad4_algorithm = ad4_alg.upper()
    try:
        config.ga_run = int(_get("ga_run", str(getattr(config, "ga_run", 100))))
    except ValueError:
        pass
    try:
        config.ga_pop_size = int(_get("ga_pop_size", str(getattr(config, "ga_pop_size", 150))))
    except ValueError:
        pass
    try:
        config.ga_num_evals = int(_get("ga_num_evals", str(getattr(config, "ga_num_evals", 2500000))))
    except ValueError:
        pass
    try:
        config.ga_num_generations = int(_get("ga_num_generations", str(getattr(config, "ga_num_generations", 27000))))
    except ValueError:
        pass
    try:
        config.ga_elitism = int(_get("ga_elitism", str(getattr(config, "ga_elitism", 1))))
    except ValueError:
        pass
    try:
        config.ga_mutation_rate = float(_get("ga_mutation_rate", str(getattr(config, "ga_mutation_rate", 0.02))))
    except ValueError:
        pass
    try:
        config.ga_crossover_rate = float(_get("ga_crossover_rate", str(getattr(config, "ga_crossover_rate", 0.80))))
    except ValueError:
        pass
    try:
        config.sw_max_its = int(_get("sw_max_its", str(getattr(config, "sw_max_its", 300))))
    except ValueError:
        pass
    try:
        config.sw_max_succ = int(_get("sw_max_succ", str(getattr(config, "sw_max_succ", 4))))
    except ValueError:
        pass
    try:
        config.sw_max_fail = int(_get("sw_max_fail", str(getattr(config, "sw_max_fail", 4))))
    except ValueError:
        pass
    try:
        config.sw_rho = float(_get("sw_rho", str(getattr(config, "sw_rho", 1.0))))
    except ValueError:
        pass
    try:
        config.sw_lb_rho = float(_get("sw_lb_rho", str(getattr(config, "sw_lb_rho", 0.01))))
    except ValueError:
        pass
    try:
        config.ls_search_freq = float(_get("ls_search_freq", str(getattr(config, "ls_search_freq", 0.06))))
    except ValueError:
        pass
    try:
        config.rmstol = float(_get("rmstol", str(getattr(config, "rmstol", 2.0))))
    except ValueError:
        pass
    try:
        config.outlev = int(_get("outlev", str(getattr(config, "outlev", 1))))
    except ValueError:
        pass
    config.unbound_model = _get("unbound_model", getattr(config, "unbound_model", "bound"))
    config.ad4_seed = _get("ad4_seed", getattr(config, "ad4_seed", "AUTO"))


# ---------------------------------------------------------------------------
# Save Config → TOML
# ---------------------------------------------------------------------------

def save_config_to_toml(config: ProjectConfig, path: Path) -> None:
    """Serialises the current ProjectConfig back to a project_config.toml file.

    Only writes the fields that matter for reproducibility. Comments are
    regenerated fresh (not round-tripped) to keep the file clean.
    """

    def _exe(p: Path | None) -> str:
        return str(p).replace("\\", "/") if p else ""

    lines = [
        "# ============================================================\n",
        "# AUTODOCK SUITE PRO — PROJECT CONFIGURATION v0.3.0\n",
        "# ============================================================\n",
        "\n",
        "[project]\n",
        f'name = "{config.project_name}"\n',
        f'root = "."\n',
        "\n",
        "[executables]\n",
        f'vina       = "{_exe(config.vina_executable)}"\n',
        f'vina_split = "{_exe(config.vina_split_executable)}"\n',
        f'autogrid4  = "{_exe(config.autogrid4_executable)}"\n',
        f'autodock4  = "{_exe(config.autodock4_executable)}"\n',
        f'obabel     = "{_exe(config.obabel_executable)}"\n',
        "\n",
        "[inputs]\n",
        f'receptor_directory = "{str(config.receptor_directory).replace(chr(92), "/")}"\n',
        f'ligand_directory   = "{str(config.ligand_directory).replace(chr(92), "/")}"\n',
        "\n",
        "[outputs]\n",
        f'result_directory = "{str(config.result_directory).replace(chr(92), "/")}"\n',
        f'log_directory    = "{str(config.log_directory).replace(chr(92), "/")}"\n',
        f'report_directory = "{str(config.report_directory).replace(chr(92), "/")}"\n',
        "\n",
        "[engine]\n",
        f'type = "{config.engine.value}"\n',
        "\n",
        "[docking]\n",
        f'mode = "{config.docking_mode.value}"\n',
        "\n",
        "[vina]\n",
        f"exhaustiveness = {config.exhaustiveness}\n",
        f"num_modes      = {config.num_modes}\n",
        f"energy_range   = {config.energy_range:.1f}\n",
        f'cpu            = "{config.cpu}"\n',
        f'seed           = "{config.seed}"\n',
        f'scoring        = "{getattr(config, "vina_scoring", "vina")}"\n',
        f"min_rmsd       = {getattr(config, 'vina_min_rmsd', 1.0):.1f}\n",
        "\n",
        "[retry]\n",
        f"enabled = {str(config.retry_enabled).lower()}\n",
        "\n",
        "[autodock4]\n",
        f'algorithm          = "{config.ad4_algorithm}"\n',
        f"ga_run             = {config.ga_run}\n",
        f"ga_pop_size        = {config.ga_pop_size}\n",
        f"ga_num_evals       = {config.ga_num_evals}\n",
        f"ga_num_generations = {getattr(config, 'ga_num_generations', 27000)}\n",
        f"ga_elitism         = {getattr(config, 'ga_elitism', 1)}\n",
        f"ga_mutation_rate   = {getattr(config, 'ga_mutation_rate', 0.02)}\n",
        f"ga_crossover_rate  = {getattr(config, 'ga_crossover_rate', 0.80)}\n",
        f"sw_max_its         = {getattr(config, 'sw_max_its', 300)}\n",
        f"sw_max_succ        = {getattr(config, 'sw_max_succ', 4)}\n",
        f"sw_max_fail        = {getattr(config, 'sw_max_fail', 4)}\n",
        f"sw_rho             = {getattr(config, 'sw_rho', 1.0)}\n",
        f"sw_lb_rho          = {getattr(config, 'sw_lb_rho', 0.01)}\n",
        f"ls_search_freq     = {getattr(config, 'ls_search_freq', 0.06)}\n",
        f"rmstol             = {config.rmstol:.1f}\n",
        f"outlev             = {getattr(config, 'outlev', 1)}\n",
        f'unbound_model      = "{getattr(config, "unbound_model", "bound")}"\n',
        f'seed               = "{getattr(config, "ad4_seed", "AUTO")}"\n',
        f"reuse_existing_maps = {str(getattr(config, 'reuse_existing_maps', False)).lower()}\n",
        "\n",
        "[execution]\n",
        f"sequential  = {str(config.sequential).lower()}\n",
        f"max_workers = {config.max_workers}\n",
        f"resume      = {str(config.resume).lower()}\n",
        f"dry_run     = {str(config.dry_run).lower()}\n",
        "\n",
        "[reporting]\n",
        f"csv   = true\n",
        f"excel = true\n",
        "\n",
        "[ui]\n",
        f'theme = "{getattr(config, "ui_theme", "Noir")}"\n',
    ]

    path.write_text("".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Load Config
# ---------------------------------------------------------------------------

def load_config_safe(toml_path: Path) -> ProjectConfig:
    """Loads a ProjectConfig from TOML, returning a default on any error."""
    try:
        return load_config(toml_path)
    except Exception:
        return ProjectConfig()
