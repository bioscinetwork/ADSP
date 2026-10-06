# Docking Automation Suite — Developer & Architecture Guide

**Version**: 1.0.0  
**Target Audience**: Computational Chemists, Scientific Software Engineers, and Pipeline Developers  

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Module Map & Component Responsibilities](#2-module-map--component-responsibilities)
3. [Data Flow & Lifecycle](#3-data-flow--lifecycle)
4. [Core Data Models (`models.py`)](#4-core-data-models-modelspy)
5. [Configuration Engine (`config.py`)](#5-configuration-engine-configpy)
6. [Job Management & State Persistence (`job_manager.py`)](#6-job-management--state-persistence-job_managerpy)
7. [Vina Workflow Implementation](#7-vina-workflow-implementation)
   - [Parser (`vina_parser.py`)](#parser-vina_parserpy)
   - [Splitter (`vina_splitter.py`)](#splitter-vina_splitterpy)
   - [Workflow Orchestrator (`vina_workflow.py`)](#workflow-orchestrator-vina_workflowpy)
8. [AutoDock 4 Workflow & BSNDVP™ Integration (`autodock4_workflow.py`)](#8-autodock-4-workflow--bsndvp-integration-autodock4_workflowpy)
9. [Reporting Engine (`reporting.py`)](#9-reporting-engine-reportingpy)
10. [Error Handling & Resiliency Strategy](#10-error-handling--resiliency-strategy)
11. [Extending the Suite](#11-extending-the-suite)
    - [Adding a New Docking Engine](#adding-a-new-docking-engine)
    - [Adding a Custom Exporter](#adding-a-custom-exporter)
12. [Testing & Verification](#12-testing--verification)

---

## 1. Architecture Overview

The **Docking Automation Suite** is engineered with a modular, decoupled architecture adhering to the following design principles:

- **Single Responsibility Principle**: Distinct modules handle validation, execution, parsing, splitting, reporting, and state management.
- **Fail-Safe Execution**: Every docking job runs independently within a sandboxed directory structure. Failure in one ligand–receptor pair never aborts the overall pipeline.
- **Deterministic State Persistence**: State is updated to `job_status.json` after every individual docking calculation, allowing instant resume from crashes or power interruptions.
- **Engine Agnostic Abstraction**: The core queue, progress tracker, and reporting engine interact with uniform `DockingJob` models regardless of whether Vina or AutoDock4 executes the simulation.

```mermaid
graph TD
    CLI[main.py / CLI] --> Config[config.py]
    CLI --> Val[validators.py / executables.py]
    CLI --> JM[job_manager.py]
    JM --> VW[vina_workflow.py]
    JM --> AW[autodock4_workflow.py]
    
    VW --> VP[vina_parser.py]
    VW --> VS[vina_splitter.py]
    
    AW --> DLG[dlg_extract.py (BSNDVP™)]
    
    VW --> REP[reporting.py]
    AW --> REP
    
    REP --> CSV[CSV Reports]
    REP --> XLSX[Excel Spreadsheets]
    REP --> JSON[JSON Metadata]
```

---

## 2. Module Map & Component Responsibilities

| Module | Lines | Primary Purpose |
|---|---|---|
| `models.py` | ~250 | Central dataclasses (`DockingJob`, `VinaResult`, `AD4Result`) and enums (`Engine`, `JobStatus`, `ResumeMode`). |
| `config.py` | ~385 | TOML loader, schema validator, default parameter manager, and template generator. |
| `executables.py` | ~110 | Cross-platform binary verification, version sniffing (`--version`), and CPU core detection. |
| `validators.py` | ~230 | PDBQT structure validation, Vina config parsing, grid parameter validation, and post-docking sanity checks. |
| `logging_utils.py` | ~90 | Hierarchical logging: global execution logs + isolated per-job log files with timestamps. |
| `job_manager.py` | ~400 | Pairwise matrix builder, resume filtering, atomic state serialization, and terminal progress indicators. |
| `vina_parser.py` | ~120 | Regex-based parser for multi-model Vina PDBQT files and affinity score extraction. |
| `vina_splitter.py` | ~180 | Robust `vina_split` wrapper with automated fallback manual PDBQT model chunker. |
| `vina_workflow.py` | ~580 | End-to-end AutoDock Vina orchestrator: argument generation, retries, pose splitting, and metadata capture. |
| `autodock4_workflow.py` | ~610 | AutoDock 4 pipeline: GPF/DPF parameter generation, AutoGrid 4 mapping, AutoDock 4 GA search, and BSNDVP™ handoff. |
| `reporting.py` | ~250 | Multi-sheet Excel workbook generator (via `openpyxl`), standardized CSV tables, and summary statistics. |
| `main.py` | ~230 | Unified command-line interface with `--init`, `--validate`, `--dry-run`, `--resume`, and `--report-only`. |
| `dlg_extract.py` | ~3,300 | Pre-existing BSNDVP™ AutoDock 4 comprehensive DLG parser and clustering analysis engine. |
| `viewer_3d.py` | ~630 | Interactive standalone WebGL 3D molecular visualizer embedding `3Dmol.js`. |
| `pymol_exporter.py` | ~750 | PyMOL session (.pse) and script (.pml) generator with 3D non-covalent interaction rendering. |

---

## 3. Data Flow & Lifecycle

1. **Initialization Phase**:
   - `main.py` parses CLI arguments.
   - `config.py` reads `project_config.toml`, validates required keys, types, and path existences.
   - `executables.py` tests binary execution and logs detected versions.
2. **Queue Generation**:
   - `job_manager.build_job_queue()` traverses `receptor_directory` and `ligand_directory`.
   - Generates all $N \times M$ `DockingJob` instances.
   - `filter_jobs()` inspects existing `job_status.json` and filters out completed tasks if `--resume` is active.
3. **Execution Phase**:
   - Jobs are processed sequentially.
   - Each job directory is populated with specific configs and log sinks.
   - External processes (`vina.exe` or `autogrid4.exe` + `autodock4.exe`) are monitored with strict timeout enforcement.
4. **Post-Processing Phase**:
   - Output PDBQT/DLG is parsed and validated for structural integrity.
   - Poses are split into individual numbered PDBQT files.
   - `job_status.json` is updated atomically.
5. **Reporting Phase**:
   - `reporting.generate_reports()` aggregates all `DockingJob` results.
   - Generates Excel workbooks with ranked leaderboards and detailed pose metrics.

---

## 4. Core Data Models (`models.py`)

### `DockingJob`
The central entity passed across the pipeline:
```python
@dataclass
class DockingJob:
    job_id: str
    receptor_name: str
    ligand_name: str
    receptor_pdbqt: Path
    ligand_pdbqt: Path
    engine: Engine
    
    # Execution states
    status: JobStatus = JobStatus.PENDING
    execution_status: ExecutionStatus = ExecutionStatus.NOT_STARTED
    analysis_status: AnalysisStatus = AnalysisStatus.NOT_STARTED
    
    # Results
    vina_results: List[VinaResult] = field(default_factory=list)
    ad4_results: Dict[str, Any] = field(default_factory=dict)
    
    # Diagnostics
    elapsed_seconds: float = 0.0
    exit_code: Optional[int] = None
    warnings: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
```

### `VinaResult`
Encapsulates an individual pose from Vina multi-model output:
```python
@dataclass
class VinaResult:
    pose: int
    binding_affinity: float  # kcal/mol
    rmsd_lower_bound: float
    rmsd_upper_bound: float
    pdbqt_path: Optional[Path] = None
```

---

## 5. Configuration Engine (`config.py`)

Configuration is stored in human-readable TOML:
- Loads via standard Python 3.11+ `tomllib`.
- Validates that paths exist, parameters fall within physical ranges (e.g., $exhaustiveness \ge 1$), and required executables are reachable.
- Includes `generate_template_config(output_path)` to emit a fully annotated starter configuration with standard paths (`C:/Program Files (x86)/MGLTools-1.5.7/vina.exe`).
- Features `resolve_executable()` in `executables.py` with multi-PC automatic discovery across standard installation folders (`MGLTools`, `ADT3`, `Vina`, `PATH`), ensuring compatibility on different PCs without manual config adjustments.

---

## 6. Job Management & State Persistence (`job_manager.py`)

### State File Schema (`results/job_status.json`)
The status of every job is preserved as a JSON array of records:
```json
[
  {
    "job_id": "mpro_dimer_vs_ligand_01",
    "receptor_name": "mpro_dimer",
    "ligand_name": "ligand_01",
    "engine": "VINA",
    "status": "SUCCESS",
    "requested_modes": 9,
    "obtained_modes": 9,
    "elapsed_seconds": 14.2,
    "vina_results": [
      {
        "pose": 1,
        "binding_affinity": -8.4,
        "rmsd_lower_bound": 0.0,
        "rmsd_upper_bound": 0.0
      }
    ]
  }
]
```

### Atomic File Writes
State writes employ temporary files (`.tmp`) followed by atomic rename operations to prevent file corruption in case of unexpected process termination.

---

## 7. Vina Workflow Implementation

### Parser (`vina_parser.py`)
Uses compiled regular expressions to parse Vina PDBQT models and affinity remarks:
```python
RE_MODEL_START = re.compile(r"^MODEL\s+(\d+)", re.MULTILINE)
RE_AFFINITY = re.compile(
    r"^REMARK\s+VINA\s+RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
    re.MULTILINE
)
```

### Splitter (`vina_splitter.py`)
Splits multi-model output using the official `vina_split.exe`. If `vina_split` is missing or fails, an internal fallback parser splits by `MODEL ... ENDMDL` records and preserves headers, ensuring zero pose loss.

### Workflow Orchestrator (`vina_workflow.py`)
Orchestrates:
1. Dynamic config parameter synthesis (`input_config.txt`).
2. Subprocess execution with timeout capture.
3. Post-docking model validation.
4. Retry escalation with higher exhaustiveness if poses obtained < requested modes.

---

## 8. AutoDock 4 Workflow & BSNDVP™ Integration (`autodock4_workflow.py`)

### GPF & DPF Generation
When explicit `.gpf` or `.dpf` template files are not provided, the module dynamically synthesizes valid parameter files:
- Determines receptor grid center and bounding dimensions from `<receptor>.txt`.
- Infers atom types directly from receptor and ligand PDBQT atom records.
- Configures Lamarckian Genetic Algorithm parameters (`ga_run`, `ga_num_evals`, `rmstol`).

### Integration with `dlg_extract.py`
Completed `.dlg` files are collected and passed directly to the pre-existing BSNDVP™ extraction engine:
```python
spec = importlib.util.spec_from_file_location("dlg_extract", str(dlg_extract_path))
dlg_extract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dlg_extract)

pipeline = dlg_extract.DLGPipeline(dlg_dir, output_dir, dlg_analysis_dir)
pipeline.run()
```

---

## 9. Reporting Engine (`reporting.py`)

Produces three output formats:
1. **Spreadsheets (`.xlsx`)**: Formatted using `openpyxl` with column headers, auto-adjusted widths, and numeric styling.
2. **Delimited Tables (`.csv`)**: Standard RFC 4180 CSV tables with headers.
3. **Structured Data (`.json`)**: Machine-readable output for programmatic pipelines.

---

## 10. Error Handling & Resiliency Strategy

| Failure Mode | Mitigation Strategy |
|---|---|
| Missing executable | Detected during pre-flight validation; informative error printed before docking starts. |
| Invalid PDBQT format | Validated before docking; job flagged as `INVALID_INPUT` and logged without crashing batch. |
| Timeout / Hung Process | Subprocesses monitored with timeout limits; process killed and marked `EXECUTION_FAILED`. |
| Corrupt output PDBQT | Validated by `validate_vina_output()`; flagged as `PARSING_FAILED`. |
| Power failure / Crash | Status saved after every job; resume with `python main.py --resume`. |

---

## 11. Extending the Suite

### Adding a New Docking Engine
1. Add an entry to the `Engine` enum in `models.py` (e.g., `GNINA = "GNINA"`).
2. Create `<engine>_workflow.py` implementing `run_<engine>_workflow(config, resume_mode)`.
3. Add engine-specific parameters to `ProjectConfig` in `config.py`.
4. Register the new engine in `main.py`.

### Adding a Custom Exporter
1. Define your generator function in `reporting.py`:
   ```python
   def export_custom_format(jobs: List[DockingJob], output_path: Path) -> None:
       ...
   ```
2. Call your exporter within `generate_reports()`.

---

## 12. Testing & Verification

Run the verification suite to ensure all components function correctly:
```bash
# Verify environment and binaries
python main.py --version

# Run validation against template configuration
python main.py --validate --config project_config.toml

# Execute dry-run simulation
python main.py --dry-run --config project_config.toml
```
