#!/usr/bin/env python3
"""
AutoDock Suite Pro — Core Data Models
======================================
Enumerations, dataclasses, and type definitions shared across all modules.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ═══════════════════════════════════════════════════════════════════════════════
# Suite Version
# ═══════════════════════════════════════════════════════════════════════════════

__version__ = "0.3.0"
SUITE_NAME = "AutoDock Suite Pro"
SUITE_BANNER = f"{SUITE_NAME} v{__version__}"


# ═══════════════════════════════════════════════════════════════════════════════
# Enumerations
# ═══════════════════════════════════════════════════════════════════════════════

class Engine(str, Enum):
    """Docking engine selection."""
    VINA = "VINA"
    AUTODOCK4 = "AUTODOCK4"
    BOTH = "BOTH"          # Run Vina then AutoDock4 sequentially in one pass


class DockingMode(str, Enum):
    """Rigid or flexible docking mode."""
    AUTO = "AUTO"
    RIGID = "RIGID"
    FLEXIBLE = "FLEXIBLE"
    BOTH = "BOTH"


class ExecutionStatus(str, Enum):
    """Status of the docking execution step."""
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    CANCEL_REQUESTED = "CANCEL_REQUESTED"
    CANCELLED = "CANCELLED"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"
    TIMEOUT = "TIMEOUT"
    SKIPPED = "SKIPPED"


class AnalysisStatus(str, Enum):
    """Status of post-docking analysis (parsing, splitting, etc.)."""
    PENDING = "PENDING"
    SUCCESS = "SUCCESS"
    WARNING = "WARNING"
    FAILED = "FAILED"


class JobStatus(str, Enum):
    """Final composite status for a docking job."""
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    CANCEL_REQUESTED = "CANCEL_REQUESTED"
    CANCELLED = "CANCELLED"
    SUCCESS = "SUCCESS"
    SUCCESS_WITH_WARNING = "SUCCESS_WITH_WARNING"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    PARSING_FAILED = "PARSING_FAILED"
    SPLIT_FAILED = "SPLIT_FAILED"
    POST_PROCESSING_FAILED = "POST_PROCESSING_FAILED"
    REPORT_GENERATION_FAILED = "REPORT_GENERATION_FAILED"
    INVALID_INPUT = "INVALID_INPUT"
    SKIPPED = "SKIPPED"


class ResumeMode(str, Enum):
    """How to handle previously completed jobs."""
    RESUME = "RESUME"             # Skip completed, run remaining
    RERUN_FAILED = "RERUN_FAILED" # Rerun failed + run remaining
    FORCE = "FORCE"               # Rerun everything


# ═══════════════════════════════════════════════════════════════════════════════
# Scientific Formatting and Missing Value Semantics
# ═══════════════════════════════════════════════════════════════════════════════

def format_metric(
    val: Optional[float | int],
    unit: str = "",
    precision: int = 2,
    na_str: str = "N/A — not reported by this engine",
) -> str:
    """Format scientific metrics without silently converting missing data to zero."""
    if val is None:
        return na_str
    if isinstance(val, int):
        return f"{val} {unit}".strip()
    if math.isnan(val) or math.isinf(val):
        return na_str
    formatted = f"{val:.{precision}f}"
    return f"{formatted} {unit}".strip() if unit else formatted


def calculate_inhibition_constant(
    delta_g: Optional[float], temp_k: Optional[float] = None
) -> Tuple[Optional[float], str]:
    """Convert an energy-like value to a concentration proxy (legacy only).

    Returns:
        (ki_molar, formatted_string_with_unit).
        If delta_g is None, returns (None, "N/A — not calculated").
    """
    if delta_g is None or math.isnan(delta_g) or math.isinf(delta_g):
        return (None, "N/A — not calculated")

    # A temperature must be explicit. This legacy helper must never be used to
    # derive Ki from Vina affinity; only AutoDock4's own DLG Ki is authoritative.
    if temp_k is None or not math.isfinite(temp_k) or temp_k <= 0:
        return (None, "N/A — temperature not supplied")

    # R = 1.98720425864083e-3 kcal/(mol·K)
    rt = 0.00198720425864083 * temp_k
    try:
        ki = math.exp(delta_g / rt)
    except OverflowError:
        return (None, "> 1 M (very low affinity)")

    if ki < 1e-15:
        return (ki, f"{ki * 1e18:.2f} aM")
    elif ki < 1e-12:
        return (ki, f"{ki * 1e15:.2f} fM")
    elif ki < 1e-9:
        return (ki, f"{ki * 1e12:.2f} pM")
    elif ki < 1e-6:
        return (ki, f"{ki * 1e9:.2f} nM")
    elif ki < 1e-3:
        return (ki, f"{ki * 1e6:.2f} μM")
    elif ki < 1.0:
        return (ki, f"{ki * 1e3:.2f} mM")
    return (ki, f"{ki:.2f} M")


# ═══════════════════════════════════════════════════════════════════════════════
# Engine-Specific Scientific Metric Structures
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class VinaMetrics:
    """AutoDock Vina-specific pose metrics."""
    rmsd_lower_bound: Optional[float] = None  # RMSD l.b. relative to best pose
    rmsd_upper_bound: Optional[float] = None  # RMSD u.b. relative to best pose

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rmsd_lower_bound": self.rmsd_lower_bound,
            "rmsd_upper_bound": self.rmsd_upper_bound,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "VinaMetrics":
        return cls(
            rmsd_lower_bound=data.get("rmsd_lower_bound"),
            rmsd_upper_bound=data.get("rmsd_upper_bound"),
        )


@dataclass
class AutoDock4Metrics:
    """AutoDock 4-specific pose metrics."""
    run_number: Optional[int] = None
    intermolecular_energy: Optional[float] = None
    internal_energy: Optional[float] = None
    torsional_energy: Optional[float] = None
    unbound_energy: Optional[float] = None
    cluster_id: Optional[int] = None
    cluster_size: Optional[int] = None
    cluster_population_percent: Optional[float] = None
    cluster_population: Optional[float] = None
    cluster_rmsd: Optional[float] = None
    rmsd_from_reference: Optional[float] = None
    receptor_strain: Optional[float] = None
    cluster_rank: Optional[int] = None
    dlg_source: Optional[str] = None
    vdW_hbond_desolvation_energy: Optional[float] = None
    electrostatic_energy: Optional[float] = None

    def __post_init__(self):
        if self.cluster_population is not None and self.cluster_population_percent is None:
            self.cluster_population_percent = self.cluster_population
        elif self.cluster_population_percent is not None and self.cluster_population is None:
            self.cluster_population = self.cluster_population_percent

    def to_dict(self) -> Dict[str, Any]:
        return {
            "run_number": self.run_number,
            "intermolecular_energy": self.intermolecular_energy,
            "internal_energy": self.internal_energy,
            "torsional_energy": self.torsional_energy,
            "unbound_energy": self.unbound_energy,
            "cluster_id": self.cluster_id,
            "cluster_size": self.cluster_size,
            "cluster_population_percent": self.cluster_population_percent,
            "cluster_population": self.cluster_population,
            "cluster_rmsd": self.cluster_rmsd,
            "rmsd_from_reference": self.rmsd_from_reference,
            "receptor_strain": self.receptor_strain,
            "cluster_rank": self.cluster_rank,
            "dlg_source": self.dlg_source,
            "vdW_hbond_desolvation_energy": self.vdW_hbond_desolvation_energy,
            "electrostatic_energy": self.electrostatic_energy,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "AutoDock4Metrics":
        pop = data.get("cluster_population", data.get("cluster_population_percent"))
        return cls(
            run_number=data.get("run_number"),
            intermolecular_energy=data.get("intermolecular_energy"),
            internal_energy=data.get("internal_energy"),
            torsional_energy=data.get("torsional_energy"),
            unbound_energy=data.get("unbound_energy"),
            cluster_id=data.get("cluster_id"),
            cluster_size=data.get("cluster_size"),
            cluster_population_percent=pop,
            cluster_population=pop,
            cluster_rmsd=data.get("cluster_rmsd"),
            rmsd_from_reference=data.get("rmsd_from_reference"),
            receptor_strain=data.get("receptor_strain"),
            cluster_rank=data.get("cluster_rank"),
            dlg_source=data.get("dlg_source"),
            vdW_hbond_desolvation_energy=data.get("vdW_hbond_desolvation_energy"),
            electrostatic_energy=data.get("electrostatic_energy"),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Canonical Pose & Cluster Models
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CanonicalPose:
    """Canonical representation of a single docked conformation."""
    rank: int = 1
    binding_score: Optional[float] = None  # legacy sorting alias; see engine-specific metrics
    run_number: Optional[int] = None
    estimated_ki: Optional[float] = None
    estimated_ki_formatted: str = "N/A"
    vina_metrics: Optional[VinaMetrics] = None
    ad4_metrics: Optional[AutoDock4Metrics] = None
    validation_rmsd: Optional[float] = None
    interactions: List[Dict[str, Any]] = field(default_factory=list)
    pose_file: Optional[Path] = None
    raw_atoms: List[Dict[str, Any]] = field(default_factory=list)
    binding_energy: Optional[float] = None
    estimated_ki_nM: Optional[float] = None
    source_model: Optional[int] = None
    source_run: Optional[int] = None
    source_file: Optional[str] = None

    def __post_init__(self):
        if self.binding_energy is not None and self.binding_score is None:
            self.binding_score = self.binding_energy
        elif self.binding_score is not None and self.binding_energy is None:
            self.binding_energy = self.binding_score

        if self.estimated_ki_nM is not None and self.estimated_ki is None:
            self.estimated_ki = self.estimated_ki_nM * 1e-9
        elif self.estimated_ki is not None and self.estimated_ki_nM is None:
            self.estimated_ki_nM = self.estimated_ki * 1e9

    @property
    def pose_id(self) -> str:
        """Canonical pose identifier, e.g. 'pose_1'. Consistent with Interaction.pose_id."""
        return f"pose_{self.rank}"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rank": self.rank,
            "binding_score": self.binding_score,
            "binding_energy": self.binding_score,
            "run_number": self.run_number,
            "estimated_ki": self.estimated_ki,
            "estimated_ki_nM": self.estimated_ki_nM,
            "estimated_ki_formatted": self.estimated_ki_formatted,
            "vina_metrics": self.vina_metrics.to_dict() if self.vina_metrics else None,
            "ad4_metrics": self.ad4_metrics.to_dict() if self.ad4_metrics else None,
            "validation_rmsd": self.validation_rmsd,
            "interactions": self.interactions,
            "pose_file": str(self.pose_file) if self.pose_file else None,
            "source_model": self.source_model,
            "source_run": self.source_run,
            "source_file": self.source_file,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "CanonicalPose":
        v_m = VinaMetrics.from_dict(data["vina_metrics"]) if data.get("vina_metrics") else None
        ad_m = AutoDock4Metrics.from_dict(data["ad4_metrics"]) if data.get("ad4_metrics") else None
        b_score = data.get("binding_score", data.get("binding_energy"))
        ki_val = data.get("estimated_ki")
        ki_nm = data.get("estimated_ki_nM")
        if ki_val is None and ki_nm is not None:
            ki_val = ki_nm * 1e-9
        elif ki_val is not None and ki_nm is None:
            ki_nm = ki_val * 1e9
        return cls(
            rank=data.get("rank", 1),
            binding_score=b_score,
            binding_energy=b_score,
            run_number=data.get("run_number"),
            estimated_ki=ki_val,
            estimated_ki_nM=ki_nm,
            estimated_ki_formatted=data.get("estimated_ki_formatted", "N/A"),
            vina_metrics=v_m,
            ad4_metrics=ad_m,
            validation_rmsd=data.get("validation_rmsd"),
            interactions=data.get("interactions", []),
            pose_file=Path(data["pose_file"]) if data.get("pose_file") else None,
            source_model=data.get("source_model"),
            source_run=data.get("source_run"),
            source_file=data.get("source_file"),
        )


@dataclass
class ClusterInfo:
    """AutoDock4 conformational cluster statistics."""
    cluster_id: int = 0
    size: int = 0
    population_percent: float = 0.0
    lowest_energy: float = 0.0
    mean_energy: float = 0.0
    cluster_rmsd: float = 0.0
    representative_run: int = 0
    representative_pose: int = 0
    runs: List[int] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cluster_id": self.cluster_id,
            "size": self.size,
            "population_percent": self.population_percent,
            "lowest_energy": self.lowest_energy,
            "mean_energy": self.mean_energy,
            "cluster_rmsd": self.cluster_rmsd,
            "representative_run": self.representative_run,
            "representative_pose": self.representative_pose,
            "runs": self.runs,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ClusterInfo":
        return cls(
            cluster_id=data.get("cluster_id", 0),
            size=data.get("size", 0),
            population_percent=data.get("population_percent", 0.0),
            lowest_energy=data.get("lowest_energy", 0.0),
            mean_energy=data.get("mean_energy", 0.0),
            cluster_rmsd=data.get("cluster_rmsd", 0.0),
            representative_run=data.get("representative_run", 0),
            representative_pose=data.get("representative_pose", 0),
            runs=data.get("runs", []),
        )


@dataclass
class ThermodynamicAnalysis:
    """Thermodynamic and statistical mechanics analysis."""
    info_entropy: Optional[float] = None
    info_entropy_rmstol: Optional[float] = None
    partition_function: Optional[float] = None
    stat_temperature: Optional[float] = None
    stat_free_energy: Optional[float] = None
    stat_internal_energy: Optional[float] = None
    stat_entropy: Optional[float] = None
    boltzmann_prob: Optional[float] = None
    descriptor_source: str = "ADSP Statistical Engine (dlg_extract)"
    temperature_reported: bool = False
    source_section: str = "AutoDock4 DLG"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "info_entropy": self.info_entropy,
            "info_entropy_rmstol": self.info_entropy_rmstol,
            "partition_function": self.partition_function,
            "stat_temperature": self.stat_temperature,
            "stat_free_energy": self.stat_free_energy,
            "stat_internal_energy": self.stat_internal_energy,
            "stat_entropy": self.stat_entropy,
            "boltzmann_prob": self.boltzmann_prob,
            "descriptor_source": self.descriptor_source,
            "temperature_reported": self.temperature_reported,
            "source_section": self.source_section,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ThermodynamicAnalysis":
        return cls(
            info_entropy=data.get("info_entropy"),
            info_entropy_rmstol=data.get("info_entropy_rmstol"),
            partition_function=data.get("partition_function"),
            stat_temperature=data.get("stat_temperature"),
            stat_free_energy=data.get("stat_free_energy"),
            stat_internal_energy=data.get("stat_internal_energy"),
            stat_entropy=data.get("stat_entropy"),
            boltzmann_prob=data.get("boltzmann_prob"),
            descriptor_source=data.get("descriptor_source", "ADSP Statistical Engine (dlg_extract)"),
            temperature_reported=data.get("temperature_reported", data.get("stat_temperature") is not None),
            source_section=data.get("source_section", "AutoDock4 DLG"),
        )


@dataclass
class ValidationMetrics:
    """Rigorous crystallographic / reference ligand validation."""
    reference_ligand: str = ""
    reference_ligand_file: str = ""
    docked_ligand: str = ""
    rmsd: Optional[float] = None
    crystal_rmsd: Optional[float] = None
    docking_search_rmsd: Optional[float] = None
    heavy_atoms_matched: Optional[int] = None
    matched_atoms: Optional[int] = None
    mapping_method: str = "Kabsch / MCS"
    atom_mapping_method: str = "Kabsch / MCS"
    status: str = "Validation not performed"
    validation_performed: bool = False
    details: str = ""

    def __post_init__(self):
        if self.reference_ligand_file and not self.reference_ligand:
            self.reference_ligand = self.reference_ligand_file
        elif self.reference_ligand and not self.reference_ligand_file:
            self.reference_ligand_file = self.reference_ligand

        if self.crystal_rmsd is not None and self.rmsd is None:
            self.rmsd = self.crystal_rmsd
        elif self.rmsd is not None and self.crystal_rmsd is None:
            self.crystal_rmsd = self.rmsd

        if self.matched_atoms is not None and self.heavy_atoms_matched is None:
            self.heavy_atoms_matched = self.matched_atoms
        elif self.heavy_atoms_matched is not None and self.matched_atoms is None:
            self.matched_atoms = self.heavy_atoms_matched

        if self.atom_mapping_method and self.mapping_method == "Kabsch / MCS":
            self.mapping_method = self.atom_mapping_method
        elif self.mapping_method and self.atom_mapping_method == "Kabsch / MCS":
            self.atom_mapping_method = self.mapping_method

        if self.validation_performed and self.status == "Validation not performed":
            self.status = "SUCCESS"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reference_ligand": self.reference_ligand,
            "reference_ligand_file": self.reference_ligand_file,
            "docked_ligand": self.docked_ligand,
            "rmsd": self.rmsd,
            "crystal_rmsd": self.crystal_rmsd,
            "docking_search_rmsd": self.docking_search_rmsd,
            "heavy_atoms_matched": self.heavy_atoms_matched,
            "matched_atoms": self.matched_atoms,
            "mapping_method": self.mapping_method,
            "atom_mapping_method": self.atom_mapping_method,
            "status": self.status,
            "validation_performed": self.validation_performed,
            "details": self.details,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ValidationMetrics":
        ref = data.get("reference_ligand", data.get("reference_ligand_file", ""))
        c_rmsd = data.get("crystal_rmsd", data.get("rmsd"))
        m_atoms = data.get("matched_atoms", data.get("heavy_atoms_matched"))
        method = data.get("atom_mapping_method", data.get("mapping_method", "Kabsch / MCS"))
        return cls(
            reference_ligand=ref,
            reference_ligand_file=ref,
            docked_ligand=data.get("docked_ligand", ""),
            rmsd=c_rmsd,
            crystal_rmsd=c_rmsd,
            docking_search_rmsd=data.get("docking_search_rmsd"),
            heavy_atoms_matched=m_atoms,
            matched_atoms=m_atoms,
            mapping_method=method,
            atom_mapping_method=method,
            status=data.get("status", "Validation not performed"),
            validation_performed=data.get("validation_performed", False),
            details=data.get("details", ""),
        )


@dataclass
class ProvenanceRecord:
    """Complete scientific provenance and execution reproducibility manifest."""
    application_version: str = __version__
    app_version: Optional[str] = None
    schema_version: str = "2.0.0"
    engine: str = ""
    engine_version: str = ""
    docking_mode: str = ""
    receptor_name: str = ""
    receptor_path: Optional[str] = None
    receptor_hash: str = ""
    ligand_name: str = ""
    ligand_path: Optional[str] = None
    ligand_hash: str = ""
    config_path: Optional[str] = None
    config_hash: str = ""
    executable_path: Optional[str] = None
    executable_hash: str = ""
    start_time: str = ""
    completion_time: str = ""
    elapsed_seconds: float = 0.0
    output_files: List[str] = field(default_factory=list)
    output_status: str = ""
    parser_version: str = "5.0.0"
    analysis_status: str = ""
    reference_ligand_path: Optional[str] = None
    reference_ligand_hash: Optional[str] = None
    cleanup_provenance: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.app_version is not None and (self.application_version == __version__ or not self.application_version):
            self.application_version = self.app_version
        elif self.application_version and self.app_version is None:
            self.app_version = self.application_version

    def to_dict(self) -> Dict[str, Any]:
        return {
            "application_version": self.application_version,
            "app_version": self.application_version,
            "schema_version": self.schema_version,
            "engine": self.engine,
            "engine_version": self.engine_version,
            "docking_mode": self.docking_mode,
            "receptor_name": self.receptor_name,
            "receptor_path": self.receptor_path,
            "receptor_hash": self.receptor_hash,
            "ligand_name": self.ligand_name,
            "ligand_path": self.ligand_path,
            "ligand_hash": self.ligand_hash,
            "config_path": self.config_path,
            "config_hash": self.config_hash,
            "executable_path": self.executable_path,
            "executable_hash": self.executable_hash,
            "start_time": self.start_time,
            "completion_time": self.completion_time,
            "elapsed_seconds": self.elapsed_seconds,
            "output_files": self.output_files,
            "output_status": self.output_status,
            "parser_version": self.parser_version,
            "analysis_status": self.analysis_status,
            "reference_ligand_path": self.reference_ligand_path,
            "reference_ligand_hash": self.reference_ligand_hash,
            "cleanup_provenance": self.cleanup_provenance,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ProvenanceRecord":
        app_ver = data.get("application_version", data.get("app_version", __version__))
        return cls(
            application_version=app_ver,
            app_version=app_ver,
            schema_version=data.get("schema_version", "2.0.0"),
            engine=data.get("engine", ""),
            engine_version=data.get("engine_version", ""),
            docking_mode=data.get("docking_mode", ""),
            receptor_name=data.get("receptor_name", ""),
            receptor_path=data.get("receptor_path"),
            receptor_hash=data.get("receptor_hash", ""),
            ligand_name=data.get("ligand_name", ""),
            ligand_path=data.get("ligand_path"),
            ligand_hash=data.get("ligand_hash", ""),
            config_path=data.get("config_path"),
            config_hash=data.get("config_hash", ""),
            executable_path=data.get("executable_path"),
            executable_hash=data.get("executable_hash", ""),
            start_time=data.get("start_time", ""),
            completion_time=data.get("completion_time", ""),
            elapsed_seconds=data.get("elapsed_seconds", 0.0),
            output_files=data.get("output_files", []),
            output_status=data.get("output_status", ""),
            parser_version=data.get("parser_version", "5.0.0"),
            analysis_status=data.get("analysis_status", ""),
            reference_ligand_path=data.get("reference_ligand_path"),
            reference_ligand_hash=data.get("reference_ligand_hash"),
            cleanup_provenance=data.get("cleanup_provenance", {}),
        )


@dataclass
class DockingResult:
    """Unified, engine-aware canonical scientific docking result."""
    job_id: str = ""
    engine: Engine = Engine.VINA
    docking_mode: DockingMode = DockingMode.RIGID
    summary: Dict[str, Any] = field(default_factory=dict)
    poses: List[CanonicalPose] = field(default_factory=list)
    clusters: List[ClusterInfo] = field(default_factory=list)
    thermodynamics: Optional[ThermodynamicAnalysis] = None
    validation: Optional[ValidationMetrics] = None
    admet: Dict[str, Any] = field(default_factory=dict)
    provenance: Optional[ProvenanceRecord] = None
    success: bool = True
    error_message: str = ""
    schema_version: str = "2.0.0"
    best_binding_energy: Optional[float] = None
    best_binding_score: Optional[float] = None

    def __post_init__(self):
        if self.best_binding_energy is not None and self.best_binding_score is None:
            self.best_binding_score = self.best_binding_energy
        elif self.best_binding_score is not None and self.best_binding_energy is None:
            self.best_binding_energy = self.best_binding_score
        elif self.best_binding_energy is None and self.best_binding_score is None:
            valid = [p.binding_score for p in self.poses if p.binding_score is not None]
            if valid:
                val = min(valid)
                self.best_binding_score = val
                self.best_binding_energy = val

    @property
    def best_energy(self) -> Optional[float]:
        return self.best_binding_energy

    @property
    def best_pose(self) -> Optional[CanonicalPose]:
        if not self.poses:
            return None
        valid = [p for p in self.poses if p.binding_score is not None]
        return min(valid, key=lambda p: p.binding_score) if valid else self.poses[0]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "job_id": self.job_id,
            "engine": self.engine.value if hasattr(self.engine, "value") else str(self.engine),
            "docking_mode": self.docking_mode.value if hasattr(self.docking_mode, "value") else str(self.docking_mode),
            "best_binding_energy": self.best_binding_energy,
            "best_binding_score": self.best_binding_score,
            "summary": self.summary,
            "poses": [p.to_dict() for p in self.poses],
            "clusters": [c.to_dict() for c in self.clusters],
            "thermodynamics": self.thermodynamics.to_dict() if self.thermodynamics else None,
            "validation": self.validation.to_dict() if self.validation else None,
            "admet": self.admet,
            "provenance": self.provenance.to_dict() if self.provenance else None,
            "success": self.success,
            "error_message": self.error_message,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DockingResult":
        engine_str = data.get("engine", "VINA")
        mode_str = data.get("docking_mode", "RIGID")
        b_e = data.get("best_binding_energy", data.get("best_binding_score"))
        return cls(
            job_id=data.get("job_id", ""),
            engine=Engine(engine_str) if engine_str in Engine._value2member_map_ else Engine.VINA,
            docking_mode=DockingMode(mode_str) if mode_str in DockingMode._value2member_map_ else DockingMode.RIGID,
            summary=data.get("summary", {}),
            poses=[CanonicalPose.from_dict(p) for p in data.get("poses", [])],
            clusters=[ClusterInfo.from_dict(c) for c in data.get("clusters", [])],
            thermodynamics=ThermodynamicAnalysis.from_dict(data["thermodynamics"]) if data.get("thermodynamics") else None,
            validation=ValidationMetrics.from_dict(data["validation"]) if data.get("validation") else None,
            admet=data.get("admet", {}),
            provenance=ProvenanceRecord.from_dict(data["provenance"]) if data.get("provenance") else None,
            success=data.get("success", True),
            error_message=data.get("error_message", ""),
            schema_version=data.get("schema_version", "2.0.0"),
            best_binding_energy=b_e,
            best_binding_score=b_e,
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Legacy Vina Result Dataclass (Preserved for 100% Backward Compatibility)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class VinaResult:
    """A single Vina docking pose result."""
    pose: int = 0
    binding_affinity: float = 0.0
    rmsd_lower_bound: float = 0.0
    rmsd_upper_bound: float = 0.0

    @property
    def mode(self) -> int:
        return self.pose

    @property
    def rmsd_lb(self) -> float:
        return self.rmsd_lower_bound

    @property
    def rmsd_ub(self) -> float:
        return self.rmsd_upper_bound


# ═══════════════════════════════════════════════════════════════════════════════
# Docking Job Dataclass
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class DockingJob:
    """Represents a single receptor–ligand docking job."""
    job_id: str = ""
    receptor_name: str = ""
    ligand_name: str = ""
    receptor_path: Optional[Path] = None        # Rigid receptor PDBQT
    ligand_path: Optional[Path] = None
    config_path: Optional[Path] = None          # Vina config.txt
    flex_receptor_path: Optional[Path] = None   # Flexible sidechain PDBQT (flex mode only)
    output_dir: Optional[Path] = None
    engine: Engine = Engine.VINA
    docking_mode: DockingMode = DockingMode.RIGID

    # Status tracking
    status: JobStatus = JobStatus.PENDING
    execution_status: ExecutionStatus = ExecutionStatus.PENDING
    analysis_status: AnalysisStatus = AnalysisStatus.PENDING

    # Canonical Result Model (v2.0)
    canonical_result: Optional[DockingResult] = None

    # Vina-specific legacy results (backward-compatible)
    requested_modes: int = 9
    obtained_modes: int = 0
    vina_results: List[VinaResult] = field(default_factory=list)

    # AD4-specific legacy results (backward-compatible)
    ad4_results: Dict[str, Any] = field(default_factory=dict)

    # Execution details
    command_executed: str = ""
    exit_code: Optional[int] = None
    stdout: str = ""
    stderr: str = ""
    elapsed_seconds: float = 0.0

    # Warnings and errors
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    cleanup_provenance: Dict[str, Any] = field(default_factory=dict)

    # Output file paths (Vina)
    output_pdbqt: Optional[Path] = None
    output_log: Optional[Path] = None
    split_poses: List[Path] = field(default_factory=list)

    # Output file paths (AutoDock4)
    gpf_path: Optional[Path] = None    # Grid Parameter File
    glg_path: Optional[Path] = None    # Grid Log File (autogrid4 output)
    dpf_path: Optional[Path] = None    # Docking Parameter File
    dlg_path: Optional[Path] = None    # Docking Log File (autodock4 output)

    # Retry tracking
    retry_count: int = 0
    retry_history: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def best_binding_energy(self) -> Optional[float]:
        """Return the best binding energy/score regardless of engine without fabricating zeros."""
        if self.canonical_result and self.canonical_result.best_binding_score is not None:
            return self.canonical_result.best_binding_score
        if self.vina_results:
            return float(self.vina_results[0].binding_affinity)
        if self.ad4_results:
            if "best_energy" in self.ad4_results:
                return float(self.ad4_results["best_energy"])
            poses = self.ad4_results.get("poses", [])
            if poses:
                p0 = poses[0]
                val = p0.get("binding_energy", p0.get("binding_affinity", p0.get("energy")))
                if val is not None:
                    return float(val)
        return None

    def ensure_canonical_result(self) -> DockingResult:
        """Synthesize or return canonical_result from legacy fields if absent."""
        if self.canonical_result is not None:
            return self.canonical_result

        # Build from existing engine results
        result = DockingResult(
            job_id=self.job_id,
            engine=self.engine,
            docking_mode=self.docking_mode,
            success=(self.status in (JobStatus.SUCCESS, JobStatus.SUCCESS_WITH_WARNING)),
            error_message="; ".join(self.errors) if self.errors else "",
        )

        if self.engine == Engine.VINA and self.vina_results:
            for r in self.vina_results:
                pose = CanonicalPose(
                    rank=r.pose,
                    run_number=None,
                    binding_score=r.binding_affinity,
                    estimated_ki=None,
                    estimated_ki_formatted="Not applicable — Vina does not report Ki",
                    vina_metrics=VinaMetrics(
                        rmsd_lower_bound=r.rmsd_lower_bound,
                        rmsd_upper_bound=r.rmsd_upper_bound,
                    ),
                    ad4_metrics=None,
                )
                result.poses.append(pose)
            result.summary = {
                "engine": "AutoDock Vina",
                "modes": len(self.vina_results),
                "best_affinity": self.vina_results[0].binding_affinity if self.vina_results else None,
            }
        elif self.engine == Engine.AUTODOCK4 and self.ad4_results:
            raw_poses = self.ad4_results.get("poses", [])
            for r in raw_poses:
                dg = r.get("binding_energy", r.get("binding_affinity"))
                ki_val = r.get("ki_nM")
                ki_str = f"{r.get('ki_raw')} {r.get('ki_unit')}".strip() if r.get("ki_raw") is not None else ""
                pose = CanonicalPose(
                    rank=r.get("rank", 1),
                    run_number=r.get("run_number"),
                    binding_score=dg,
                    estimated_ki=ki_val,
                    estimated_ki_formatted=ki_str if ki_str else "Not reported by AutoDock4",
                    vina_metrics=None,
                    ad4_metrics=AutoDock4Metrics(
                        run_number=r.get("run_number"),
                        intermolecular_energy=r.get("intermol_energy"),
                        internal_energy=r.get("internal_energy"),
                        torsional_energy=r.get("torsional_energy"),
                        unbound_energy=r.get("unbound_energy"),
                        cluster_id=r.get("cluster_id"),
                        cluster_size=r.get("cluster_size"),
                        cluster_population_percent=r.get("cluster_population_percent"),
                        cluster_rmsd=r.get("cluster_rmsd"),
                        rmsd_from_reference=r.get("rmsd_from_ref"),
                        receptor_strain=r.get("receptor_strain"),
                    ),
                )
                result.poses.append(pose)
            clusters_raw = self.ad4_results.get("clusters", [])
            for c in clusters_raw:
                result.clusters.append(ClusterInfo.from_dict(c) if isinstance(c, dict) else c)

            result.summary = {
                "engine": "AutoDock4",
                "ga_runs": self.ad4_results.get("num_runs", len(raw_poses)),
                "docked_poses": len(raw_poses),
                "clusters": len(result.clusters),
                "best_binding_energy": self.ad4_results.get("best_energy"),
            }

        self.canonical_result = result
        return result

    def to_dict(self) -> Dict[str, Any]:
        """Serialize job for JSON status persistence."""
        # Ensure canonical result is up to date
        if self.canonical_result is None and (self.vina_results or self.ad4_results):
            self.ensure_canonical_result()

        return {
            "schema_version": "2.0.0",
            "job_id": self.job_id,
            "receptor_name": self.receptor_name,
            "ligand_name": self.ligand_name,
            "receptor_path": str(self.receptor_path) if self.receptor_path else None,
            "ligand_path": str(self.ligand_path) if self.ligand_path else None,
            "config_path": str(self.config_path) if self.config_path else None,
            "flex_receptor_path": str(self.flex_receptor_path) if self.flex_receptor_path else None,
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "engine": self.engine.value,
            "docking_mode": self.docking_mode.value,
            "status": self.status.value,
            "execution_status": self.execution_status.value,
            "analysis_status": self.analysis_status.value,
            "requested_modes": self.requested_modes,
            "obtained_modes": self.obtained_modes,
            "exit_code": self.exit_code,
            "elapsed_seconds": self.elapsed_seconds,
            "warnings": self.warnings,
            "errors": self.errors,
            "cleanup_provenance": self.cleanup_provenance,
            "output_pdbqt": str(self.output_pdbqt) if self.output_pdbqt else None,
            "output_log": str(self.output_log) if self.output_log else None,
            "split_poses": [str(p) for p in self.split_poses],
            "vina_results": [
                {
                    "pose": r.pose,
                    "binding_affinity": r.binding_affinity,
                    "rmsd_lower_bound": r.rmsd_lower_bound,
                    "rmsd_upper_bound": r.rmsd_upper_bound,
                }
                for r in self.vina_results
            ],
            "gpf_path": str(self.gpf_path) if self.gpf_path else None,
            "glg_path": str(self.glg_path) if self.glg_path else None,
            "dpf_path": str(self.dpf_path) if self.dpf_path else None,
            "dlg_path": str(self.dlg_path) if self.dlg_path else None,
            "ad4_results": self.ad4_results,
            "canonical_result": self.canonical_result.to_dict() if self.canonical_result else None,
            "retry_count": self.retry_count,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DockingJob":
        """Deserialize job from JSON status with full backward compatibility."""
        job = cls()
        job.job_id = data.get("job_id", "")
        job.receptor_name = data.get("receptor_name", "")
        job.ligand_name = data.get("ligand_name", "")
        job.receptor_path = Path(data["receptor_path"]) if data.get("receptor_path") else None
        job.ligand_path = Path(data["ligand_path"]) if data.get("ligand_path") else None
        job.config_path = Path(data["config_path"]) if data.get("config_path") else None
        job.flex_receptor_path = Path(data["flex_receptor_path"]) if data.get("flex_receptor_path") else None
        job.output_dir = Path(data["output_dir"]) if data.get("output_dir") else None
        job.engine = Engine(data.get("engine", "VINA"))
        job.docking_mode = DockingMode(data.get("docking_mode", "RIGID"))
        job.status = JobStatus(data.get("status", "PENDING"))
        job.execution_status = ExecutionStatus(data.get("execution_status", "PENDING"))
        job.analysis_status = AnalysisStatus(data.get("analysis_status", "PENDING"))
        job.requested_modes = data.get("requested_modes", 9)
        job.obtained_modes = data.get("obtained_modes", 0)
        job.exit_code = data.get("exit_code")
        job.elapsed_seconds = data.get("elapsed_seconds", 0.0)
        job.warnings = data.get("warnings", [])
        job.errors = data.get("errors", [])
        job.cleanup_provenance = data.get("cleanup_provenance", {})
        job.output_pdbqt = Path(data["output_pdbqt"]) if data.get("output_pdbqt") else None
        job.output_log = Path(data["output_log"]) if data.get("output_log") else None
        job.split_poses = [Path(p) for p in data.get("split_poses", [])]
        job.vina_results = [
            VinaResult(
                pose=r.get("pose", 0),
                binding_affinity=r.get("binding_affinity", 0.0),
                rmsd_lower_bound=r.get("rmsd_lower_bound", 0.0),
                rmsd_upper_bound=r.get("rmsd_upper_bound", 0.0),
            )
            for r in data.get("vina_results", [])
        ]
        job.gpf_path = Path(data["gpf_path"]) if data.get("gpf_path") else None
        job.glg_path = Path(data["glg_path"]) if data.get("glg_path") else None
        job.dpf_path = Path(data["dpf_path"]) if data.get("dpf_path") else None
        job.dlg_path = Path(data["dlg_path"]) if data.get("dlg_path") else None
        job.ad4_results = data.get("ad4_results", {})
        job.retry_count = data.get("retry_count", 0)

        # Deserialize canonical result if available, or synthesize
        if data.get("canonical_result"):
            try:
                job.canonical_result = DockingResult.from_dict(data["canonical_result"])
            except Exception:
                job.ensure_canonical_result()
        else:
            job.ensure_canonical_result()

        return job


# ═══════════════════════════════════════════════════════════════════════════════
# Job ID Generation
# ═══════════════════════════════════════════════════════════════════════════════

def generate_job_id(receptor_name: str, ligand_name: str,
                    docking_mode: DockingMode = DockingMode.RIGID,
                    engine: Optional["Engine"] = None) -> str:
    """Generate a deterministic, filesystem-safe job ID.

    AutoDock4 jobs receive an ``ad4_`` prefix so they never collide with Vina
    job IDs when both engines are run on the same receptor/ligand pair and their
    results are stored in a shared ``job_status.json``.

    Examples::

        generate_job_id("2v5z", "compound_001")
            → "2v5z__compound_001"
        generate_job_id("2v5z", "compound_001", engine=Engine.AUTODOCK4)
            → "ad4_2v5z__compound_001"
        generate_job_id("2v5z", "compound_001", DockingMode.FLEXIBLE)
            → "2v5z_flex__compound_001"
        generate_job_id("2v5z", "compound_001", DockingMode.FLEXIBLE, Engine.AUTODOCK4)
            → "ad4_2v5z_flex__compound_001"
    """
    rec = _sanitize_name(receptor_name)
    lig = _sanitize_name(ligand_name)
    if docking_mode == DockingMode.FLEXIBLE and not rec.lower().endswith("_flex"):
        rec = f"{rec}_flex"
    # Namespace AutoDock4 jobs so they are distinct from Vina jobs at the same
    # receptor/ligand in the shared status file.
    if engine is not None and engine == Engine.AUTODOCK4:
        return f"ad4_{rec}__{lig}"
    return f"{rec}__{lig}"


def _sanitize_name(name: str) -> str:
    """Sanitize a name for use in filenames while preserving biological identifiers."""
    # Remove file extension if present
    name = Path(name).stem
    # Replace problematic characters with underscores
    name = re.sub(r'[<>:"/\\|?*\s]+', '_', name)
    # Remove leading/trailing underscores
    name = name.strip('_')
    return name or "unknown"
