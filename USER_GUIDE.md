# AutoDock Suite Pro — User Guide

**Version**: 0.3.0  
**Engines Supported**: AutoDock Vina 1.2.7+ and AutoDock 4.2.6  

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Prerequisites & System Setup](#2-prerequisites--system-setup)
3. [How to Use This Tool — Where It Lives & How to Run It](#3-how-to-use-this-tool--where-it-lives--how-to-run-it)
4. [Project Setup Workflow](#4-project-setup-workflow)
   - [Step 1: Create Your Project Folder](#step-1-create-your-project-folder)
   - [Step 2: Generate the Config Template](#step-2-generate-the-config-template)
   - [Step 3: Scaffold the Folder Structure](#step-3-scaffold-the-folder-structure)
   - [Step 4: Prepare Your Receptor PDBQT Files](#step-4-prepare-your-receptor-pdbqt-files)
   - [Step 5: Prepare Your Ligand PDBQT Files](#step-5-prepare-your-ligand-pdbqt-files)
   - [Step 6: Set Up the Grid Config](#step-6-set-up-the-grid-config)
   - [Step 7: Edit project_config.toml](#step-7-edit-project_configtoml)
   - [Step 8: Validate Before Docking](#step-8-validate-before-docking)
   - [Step 9: Run Docking](#step-9-run-docking)
5. [Folder Structure Reference](#5-folder-structure-reference)
   - [Rigid Docking Layout](#rigid-docking-layout)
   - [Flexible Docking Layout](#flexible-docking-layout)
6. [Flexible Docking — Full Setup Guide](#6-flexible-docking--full-setup-guide)
7. [Configuration Reference (`project_config.toml`)](#7-configuration-reference-project_configtoml)
8. [Running Docking — All Commands](#8-running-docking--all-commands)
9. [Understanding Results & Output Structure](#9-understanding-results--output-structure)
10. [Troubleshooting & FAQs](#10-troubleshooting--faqs)
11. [Modern Desktop Graphical User Interface (GUI)](#11-modern-desktop-graphical-user-interface-gui)
12. [Protein–Ligand Interaction Profiler](#12-proteinligand-interaction-profiler)
13. [Multi-Core Parallel Docking](#13-multi-core-parallel-docking)
14. [Standalone Portable Executable (`.exe`) Compilation](#14-standalone-portable-executable-exe-compilation)

---

## 1. Introduction

**AutoDock Suite Pro** is a cross-engine molecular docking automation platform for high-throughput computational drug discovery. It does NOT prepare your receptor or ligand files — that is done externally using MGLTools scripts or any compatible tool. What it does is:

- Automatically dock **all ligands against all receptors** (Cartesian product: N × M jobs)
- Fully supports **AutoDock Vina 1.2.7** and **AutoDock 4.2.6**
- Handles **rigid** and **flexible** receptor docking
- Automatically generates **GPF** (Grid Parameter File), **GLG** (Grid Log), **DPF** (Docking Parameter File), and **DLG** (Docking Log) for AutoDock4
- Automatically auto-generates `config.txt` files for Vina when they are missing
- Crash-resilient: resume exactly where you left off if interrupted
- Retry policy: automatically retries with higher exhaustiveness if insufficient poses
- Generates **CSV**, **Excel**, and **JSON** summary reports
- Full **BSNDVP™** DLG extraction for AutoDock4 results

---

## 2. Prerequisites & System Setup

### 2.1 Software Requirements
- **Python** 3.10+ (tested on Python 3.14, Windows 64-bit)
- **MGLTools 1.5.7** — installs Vina, Vina Split, AutoGrid4, AutoDock4, and the preparation scripts:
  - Default path: `C:/Program Files (x86)/MGLTools-1.5.7/`
  - Key executables: `vina.exe`, `vina_split.exe`, `autogrid4.exe`, `autodock4.exe`
  - Preparation scripts: `prepare_receptor4.py`, `prepare_ligand4.py`, `prepare_flexreceptor4.py`

> [!NOTE]
> **Auto-Detection**: If your MGLTools is installed in a non-standard path, set it manually in `project_config.toml` under `[executables]`. The suite scans common paths and system PATH automatically.

### 2.2 Python Dependencies
```bash
pip install openpyxl
```
All other libraries (`tomllib`, `csv`, `json`, `subprocess`, `pathlib`) are built into Python 3.10+.

### 2.3 Verify Installed Tools
```bash
python main.py --version
```
Expected output:
```
======================================================================
  AutoDock Suite Pro v0.2.0
  Cross-Engine Molecular Docking Automation Platform
======================================================================
  Python           : 3.14.x (64bit)
  Operating System : Windows-11
  CPU Logical Cores: 4
----------------------------------------------------------------------
  AutoDock Vina    : v1.2.7
  Vina Split       : v1.2.7
  AutoGrid 4       : 4.2.6
  AutoDock 4       : 4.2.6
======================================================================
```

---

## 3. How to Use This Tool — Where It Lives & How to Run It

### Where does AutoDock Suite Pro live?
The suite is a collection of Python scripts. You are **free to place the tool folder into any working directory** (e.g. `C:\AutoDockSuitePro\` or in your user profile) — it is fully self-contained.

### Three ways to use it

### Option A — Copy the tool into your working directory *(recommended — simplest)*

Copy the entire `Docking Automation` folder directly into your project folder:

```
C:\MyResearch\DrugProject1\
  ├── main.py                  ← all tool scripts live here alongside your data
  ├── vina_workflow.py
  ├── autodock4_workflow.py
  ├── config.py
  ├── job_manager.py
  ├── prepare.py
  ├── ... (all other .py scripts)
  ├── project_config.toml
  ├── receptors\
  ├── ligands\
  └── results\
```

Then from your project folder:
```bash
cd C:\MyResearch\DrugProject1
python main.py --scaffold
python main.py --validate
python main.py
```

In `project_config.toml`:
```toml
[project]
root = "."          # current folder — receptors\ and ligands\ are right here
```

> [!TIP]
> This is the **cleanest and most portable** setup. Everything is in one place. You can zip the entire project folder and move it to another PC — just copy it and run.

### Option B — Run the tool from your project folder (without copying)
```bash
cd C:\MyResearch\DrugProject1\
python "C:\AutoDockSuitePro\main.py" --validate
python "C:\AutoDockSuitePro\main.py"
```
The suite finds `project_config.toml` automatically in the current folder.

### Option C — Run from the tool folder, point to your project
```bash
cd "C:\AutoDockSuitePro\"
python main.py --config "C:\MyResearch\DrugProject1\project_config.toml"
```

> [!TIP]
> **If you have multiple projects**, use Option C — one copy of the tool, many projects each with their own `project_config.toml`.

### What if I have multiple projects?
Each project gets its own folder with its own `project_config.toml`. Example:
```
C:\MyResearch\
  Project_Mpro\
    project_config.toml
    receptors\
    ligands\
  Project_EGFR\
    project_config.toml
    receptors\
    ligands\
```
Switch between them by navigating to each project folder before running.

---

## 4. Project Setup Workflow

### Step 1: Create Your Project Folder
Create a folder anywhere on your PC for your docking project:
```
C:\MyResearch\DrugProject1\
```

### Step 2: Generate the Config Template
Navigate into your project folder and run:
```bash
cd C:\MyResearch\DrugProject1
python "C:\...\Docking Automation\main.py" --init
```
This creates `project_config.toml` in your project folder with all options commented and explained.

### Step 3: Scaffold the Folder Structure
```bash
python "C:\...\Docking Automation\main.py" --scaffold
```
This creates the full folder structure with `README.txt` guides telling you exactly what to put where:
```
DrugProject1\
  receptors\
    README.txt                  ← Instructions
    example_receptor\
      rigid\                    ← Put your rigid PDBQT here
      flex\                     ← Put flexible PDBQTs here (optional)
      config.txt                ← Grid config stub
  ligands\
    README.txt                  ← Instructions
```

### Step 4: Prepare Your Receptor PDBQT Files

AutoDock Suite Pro does NOT prepare PDBQT files — this is done externally. Use **MGLTools scripts** from the command line:

**Using prepare_receptor4.py (MGLTools):**
```bash
cd C:\MyResearch\DrugProject1\

# Prepare your receptor
python "C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py" ^
  -r your_protein.pdb ^
  -o receptors\protein1\rigid\protein1.pdbqt ^
  -A checkhydrogens
```
Common flags:
- `-A checkhydrogens` — adds hydrogens if missing
- `-U nphs` — merge non-polar hydrogens
- `-e` — delete waters

**Alternatively with OpenBabel:**
```bash
obabel your_protein.pdb -O receptors\protein1\rigid\protein1.pdbqt
```

> [!IMPORTANT]
> The prepared PDBQT must go into `receptors\RECEPTOR_NAME\rigid\` — this is where the suite will look for it.

### Step 5: Prepare Your Ligand PDBQT Files

**Using prepare_ligand4.py (MGLTools):**
```bash
python "C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py" ^
  -l compound_001.mol2 ^
  -o ligands\compound_001.pdbqt
```

**Batch preparation (for many ligands from an SDF file with OpenBabel):**
```bash
obabel compounds.sdf -O ligands\compound_.pdbqt --split
```

> [!IMPORTANT]
> All ligand PDBQTs go directly into the `ligands\` folder — NOT in subfolders.

### Step 6: Set Up the Grid Config

For **Vina**, each receptor needs a `config.txt` specifying the grid box center and size.

**Option A — Use the `--prepare` command (recommended):**
```bash
python "C:\...\Docking Automation\main.py" --prepare
```
The suite will detect which receptors are missing `config.txt` and interactively ask you for the grid center (x, y, z) and size. You can get these coordinates from:
- **UCSF Chimera** → Tools → Surface/Binding → Grid Box
- **PyMOL** → selection center + box tool
- **AutoDockTools GUI** → Grid → Grid Box

**Option B — Create config.txt manually:**
In `receptors\protein1\config.txt`:
```ini
# Grid Box Center (Angstroms)
center_x = 12.540
center_y = -3.210
center_z = 24.890

# Grid Box Size (Angstroms)
size_x = 25.0
size_y = 25.0
size_z = 25.0

# Search parameters
exhaustiveness = 8
num_modes = 9
energy_range = 3.0
```

### Step 7: Edit project_config.toml

Open the generated `project_config.toml` and at minimum verify:
1. **Executable paths** in `[executables]` — update if MGLTools is not in the default path
2. **Engine** in `[engine]` — set `type = "VINA"` or `"AUTODOCK4"`
3. **Docking mode** in `[docking]` — set `mode = "RIGID"` or `"FLEXIBLE"`
4. **Vina parameters** in `[vina]` — adjust `exhaustiveness`, `num_modes` as needed

### Step 8: Validate Before Docking
```bash
python "C:\...\Docking Automation\main.py" --validate
```
✅ Expected output:
```
[OK] All configuration checks passed.
[OK] Receptors found : 3
[OK] Ligands found   : 20 (*.pdbqt)
```

Also run a dry-run to see the planned job matrix:
```bash
python "C:\...\Docking Automation\main.py" --dry-run
```

### Step 9: Run Docking
```bash
python "C:\...\Docking Automation\main.py"
```
That's it. The suite will dock all ligands against all receptors and generate reports.

---

## 5. Folder Structure Reference

### Rigid Docking Layout

```
MyProject\
├── project_config.toml
├── receptors\
│   ├── protein1\
│   │   ├── rigid\
│   │   │   └── protein1.pdbqt     ← Prepared rigid receptor PDBQT
│   │   └── config.txt             ← Grid box config for Vina
│   └── protein2\
│       ├── rigid\
│       │   └── protein2.pdbqt
│       └── config.txt
├── ligands\
│   ├── compound_001.pdbqt
│   ├── compound_002.pdbqt
│   └── ...
├── results\                       ← Auto-created during docking
├── logs\                          ← Auto-created during docking
└── reports\                       ← Auto-created after docking
```

### Flexible Docking Layout

```
MyProject\
├── project_config.toml
├── receptors\
│   └── protein1\
│       ├── flex\
│       │   ├── protein1_rigid.pdbqt        ← Backbone (passed to --receptor)
│       │   └── protein1_flex.pdbqt         ← Sidechains (passed to --flex)
│       └── config.txt
└── ligands\
    └── ...
```

### Dual-Mode Layout (Both Rigid and Flexible)

If you place **both** a `rigid\` folder and a `flex\` folder inside the same receptor directory:

```
receptors\
└── 2V5Z\
    ├── rigid\
    │   └── 2v5z_clean.pdbqt                ← Rigid receptor
    ├── flex\
    │   ├── 2v5z_clean_rigid.pdbqt          ← Flexible docking backbone
    │   └── 2v5z_clean_flex.pdbqt           ← Flexible docking sidechains
    └── config.txt                          ← Grid box parameters
```

> [!TIP]
> **Automatic Dual-Mode Execution**: When both `rigid\` and `flex\` folders exist, AutoDock Suite Pro automatically executes **both** rigid docking and flexible docking!
> - Rigid outputs go to: `results\VINA\2V5Z\<ligand>\`
> - Flexible outputs go to: `results\VINA\2V5Z_flex\<ligand>\`
> - Both sets of scores appear side-by-side in your CSV and Excel reports with a `Docking_Mode` column (`RIGID` vs `FLEXIBLE`).
>
> To force only one mode, use `--mode RIGID` or `--mode FLEXIBLE`.

---

## 6. Flexible Docking — Full Setup Guide

### What you need
For flexible docking, MGLTools splits your receptor into two PDBQT files:
- **Rigid part** (`*_rigid.pdbqt`): the protein backbone (static)
- **Flex part** (`*_flex.pdbqt`): the selected flexible sidechain atoms (moveable)

### How to split the receptor with MGLTools
```bash
python "C:\Program Files (x86)\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_flexreceptor4.py" ^
  -r protein1.pdbqt ^
  -s A:315_A:316_A:317
```
- `-r` = your prepared receptor PDBQT
- `-s` = residues to make flexible in format `Chain:ResNum` separated by `_`

This produces:
- `protein1_rigid.pdbqt`
- `protein1_flex.pdbqt`

### Place them in the correct subfolder
```
receptors\protein1\
  flex\
    protein1_rigid.pdbqt    ← matches *rigid*.pdbqt
    protein1_flex.pdbqt     ← matches *flex*.pdbqt
  config.txt
```

### What Vina command the suite builds
```
vina.exe
  --receptor protein1_rigid.pdbqt    ← backbone
  --flex     protein1_flex.pdbqt     ← flexible sidechains
  --ligand   compound_001.pdbqt
  --config   config.txt
  --out      compound_001_out.pdbqt
```

### Naming flexibility
The file finder accepts any names as long as they follow one of these patterns:
- `*rigid*` or `*backbone*` → treated as rigid backbone
- `*flex*` or `*sidechain*` or `*flexible*` → treated as flex sidechains

---

## 7. Configuration Reference (`project_config.toml`)

```toml
[project]
name = "SARS-CoV-2_Main_Protease"
root = "."                          # "." = project_config.toml's own folder

[executables]
vina       = "C:/Program Files (x86)/MGLTools-1.5.7/vina.exe"
vina_split = "C:/Program Files (x86)/MGLTools-1.5.7/vina_split.exe"
autogrid4  = "C:/Program Files (x86)/MGLTools-1.5.7/autogrid4.exe"
autodock4  = "C:/Program Files (x86)/MGLTools-1.5.7/autodock4.exe"

[inputs]
receptor_directory = "receptors"    # Relative to project root
ligand_directory   = "ligands"

# Dock only specific receptors or ligands (leave commented for all):
# selected_receptors = ["protein1", "protein2"]
# selected_ligands   = ["compound_001", "compound_005"]

[outputs]
result_directory = "results"
log_directory    = "logs"
report_directory = "reports"

[engine]
type = "VINA"     # "VINA" or "AUTODOCK4"

[docking]
mode = "RIGID"    # "RIGID" or "FLEXIBLE"
                  # NOTE: also auto-detected per receptor via flex/ subfolder

[prepare]
# Set true to auto-generate config.txt for receptors that are missing one.
# This requires grid_center and grid_size to be set (or use --prepare to be prompted).
auto_generate_config = false

[vina]
exhaustiveness = 8        # Search depth: 8=standard, 16=thorough, 32=very thorough
num_modes      = 9        # Max poses to generate (actual may be fewer)
energy_range   = 3.0      # Max energy above top pose to include (kcal/mol)
cpu            = "AUTO"   # "AUTO" = use all CPU cores
seed           = "AUTO"   # Integer for reproducible results, "AUTO" for random

[retry]
enabled = false
# If Vina returns fewer poses than requested, retry with higher exhaustiveness:
# exhaustiveness_values = [16, 32]

[autodock4]
reuse_existing_maps  = false    # Skip AutoGrid4 if .map files already exist
ga_run               = 100      # Lamarckian GA runs (more = more thorough but slower)
ga_pop_size          = 150      # GA population size
ga_num_evals         = 2500000  # Max energy evaluations per run
ga_num_generations   = 27000    # Max generations
num_modes            = 100      # Max poses to cluster
rmstol               = 2.0      # Cluster RMSD tolerance (Angstroms)
seed                 = "AUTO"

[execution]
sequential = true    # Always recommended (one job at a time)
resume     = true    # Skip completed jobs if run is interrupted
dry_run    = false   # Set true to simulate without actually docking

[reporting]
csv   = true
excel = true
```

---

## 8. Running Docking — All Commands

### First-time setup
```bash
# 1. Check environment
python main.py --version

# 2. Create config template
python main.py --init

# 3. Create folder structure
python main.py --scaffold

# 4. Auto-generate missing config.txt files (Vina)
python main.py --prepare

# 5. Validate everything before docking
python main.py --validate

# 6. Dry run — see what would be docked
python main.py --dry-run
```

### Running docking
```bash
# Run with default engine (from project_config.toml)
python main.py

# Run with Vina explicitly
python main.py --engine VINA

# Run with AutoDock4 explicitly
python main.py --engine AUTODOCK4

# Point to a different config file (project in another folder)
python main.py --config "C:\MyResearch\AnotherProject\project_config.toml"
```

### Docking Mode Selection (Rigid, Flexible, or Both)
```bash
# Auto mode (default) — checks folders; if both rigid/ and flex/ exist, runs BOTH!
python main.py
python main.py --mode AUTO

# Force rigid docking only (ignores flex/ folder even if present)
python main.py --mode RIGID

# Force flexible docking only (ignores rigid/ folder even if present)
python main.py --mode FLEXIBLE

# Explicitly queue both rigid and flexible docking
python main.py --mode BOTH
```

### Resume / restart controls
```bash
# Resume — skip already-completed jobs (default behavior)
python main.py --resume

# Rerun only failed jobs
python main.py --rerun-failed

# Force rerun everything from scratch
python main.py --force
```

### Filtering — dock only specific pairs
```bash
# One receptor, all ligands
python main.py --receptor protein1

# All receptors, one ligand
python main.py --ligand compound_001

# One receptor vs one ligand
python main.py --receptor protein1 --ligand compound_001
```

### Reports & utilities
```bash
# Re-generate reports without re-running docking
python main.py --report-only

# Verbose logging (shows detailed debug output)
python main.py -v
```

---

## 9. Understanding Results & Output Structure

### Directory Layout After Docking

```
MyProject\
├── results\
│   ├── job_status.json              ← Job state database (resume support)
│   ├── run_metadata.json            ← Version + environment record
│   ├── VINA\
│   │   └── protein1\
│   │       └── compound_001\
│   │           ├── input\           ← Copy of inputs used (receptor, ligand, config)
│   │           ├── output\
│   │           │   ├── compound_001_out.pdbqt   ← All poses (multi-model)
│   │           │   └── compound_001_log.txt     ← Vina stdout log
│   │           └── poses\
│   │               ├── compound_001_out_ligand_01.pdbqt   ← Best pose
│   │               ├── compound_001_out_ligand_02.pdbqt
│   │               └── ...
│   └── AD4\
│       ├── DLG_ANALYSIS\
│       │   └── BSNDVP_RESULTS\      ← AutoDock4 full analysis reports
│       └── protein1\
│           └── compound_001\
│               ├── protein1.gpf     ← Grid Parameter File
│               ├── protein1.glg     ← Grid Log File (AutoGrid4 output)
│               ├── protein1_compound_001.dpf    ← Docking Parameter File
│               └── protein1_compound_001.dlg    ← Docking Log (AutoDock4 output)
├── logs\
│   └── job_VINA_protein1__compound_001.log   ← Per-job detailed log
└── reports\
    ├── VINA_Complete_Report.xlsx    ← Excel workbook
    └── VINA_Complete_Report.csv     ← CSV export
```

### Vina Results (Excel Report)
`reports\VINA_Complete_Report.xlsx` contains:
- **01_Vina_All_Poses** — every pose for every job with affinity (kcal/mol), RMSD_LB, RMSD_UB
- **02_Vina_Summary** — ranked table of best pose per receptor–ligand pair
- **03_Job_Status** — execution health, elapsed time, warnings, errors

### AutoDock4 Files Explained
| File | What it is |
|---|---|
| `.gpf` | Grid Parameter File — tells AutoGrid4 where and how to compute grids |
| `.glg` | Grid Log File — AutoGrid4's output log (check for errors here) |
| `.dpf` | Docking Parameter File — tells AutoDock4 the GA parameters |
| `.dlg` | Docking Log File — full AutoDock4 results including all poses and cluster analysis |

---

## 10. Troubleshooting & FAQs

### Q: I ran `python main.py` but got "No receptors found"
**Cause**: Receptors are not in the expected subfolder layout.  
**Fix**: Each receptor must be in its own subfolder with a `rigid\` subdirectory:
```
receptors\
  protein1\             ← subfolder named after your receptor
    rigid\
      protein1.pdbqt    ← PDBQT file inside rigid\
    config.txt
```
Run `python main.py --scaffold` to see the expected structure.

### Q: I got "No config.txt found for receptor"
**Cause**: Vina needs a grid box configuration for each receptor.  
**Fix**: Run `python main.py --prepare` and enter the grid center/size when prompted. Or manually create `receptors\RECEPTOR_NAME\config.txt`.

### Q: Flexible docking is not running — it's doing rigid instead
**Cause**: The `flex\` subfolder is missing or empty, or the PDBQT files are not named correctly.  
**Fix**:
1. Verify `receptors\protein1\flex\` exists and contains `.pdbqt` files
2. One file must have `rigid` or `backbone` in its name → passed to `--receptor`
3. One file must have `flex`, `sidechain`, or `flexible` in its name → passed to `--flex`

Run `python main.py --validate` — it will show which mode was detected for each receptor.

### Q: Why do log files say 9 dockings/modes requested, but output PDBQT has fewer models (e.g., 4 or 6)?
**Explanation**: In AutoDock Vina, this is **completely normal and mathematically expected**. Here is how Vina's pose selection algorithm works:

1. **`num_modes` (default: 9)** is an **upper limit**, not a guarantee. Vina will never output *more* than 9 models, but it will only output models that meet its scientific quality thresholds.
2. **`energy_range` (default: 3.0 kcal/mol)**: Vina strictly discards any pose whose binding affinity is more than `energy_range` kcal/mol worse than the top pose (Pose 1).
   - *Example from your run*: If Pose 1 is **-9.9 kcal/mol**, then any pose with energy worse than **-6.9 kcal/mol** (-9.9 + 3.0) is discarded. If only 6 poses were found between -9.9 and -7.0, you will get exactly 6 models.
3. **RMSD Clustering**: Vina clusters similar conformations (usually within 1.0 Å RMSD) and only keeps the best energy representative for each cluster. Duplicate poses are removed.

**How to get more models**:
If you want Vina to retain poses with higher energy differences, increase `energy_range` in `project_config.toml`:
```toml
[vina]
num_modes    = 9
energy_range = 4.5    # Increase from 3.0 to 4.5 or 5.0 kcal/mol
```

### Q: How do I run AutoDock 4 instead of Vina? How should I structure my folders?
**Folder Structure**: Exactly the same as Vina! Both rigid and flexible docking are fully supported in AutoDock 4:
```
MyProject\
├── project_config.toml
├── receptors\
│   └── 2V5Z\
│       ├── rigid\
│       │   └── 2v5z_clean.pdbqt       ← Rigid receptor PDBQT
│       ├── flex\                      ← (Optional) AutoDock 4 Flexible Docking
│       │   ├── 2v5z_clean_rigid.pdbqt ← Rigid backbone
│       │   └── 2v5z_clean_flex.pdbqt  ← Flexible sidechains (flexres)
│       └── config.txt                 ← Grid box center & size (used to generate GPF)
└── ligands\
    ├── 2v5z_d_sag.pdbqt
    └── ...
```

**How to Run**:
You can run AutoDock 4 in either of two ways:
1. **Via CLI override (quickest)**:
   ```bash
   python main.py --engine AUTODOCK4
   ```
2. **Via `project_config.toml`**:
   Change line 59 in `project_config.toml`:
   ```toml
   [engine]
   type = "AUTODOCK4"
   ```
   Then simply run:
   ```bash
   python main.py
   ```

**Search Algorithm Options (in `project_config.toml`)**:
Under `[autodock4]`, you can select your search algorithm:
```toml
[autodock4]
algorithm = "LGA"   # Options: "LGA" (Lamarckian GA), "GA" (Traditional GA), or "LS" (Local Search)
ga_run    = 100     # Number of docking runs
ga_evals  = 2500000 # Number of energy evaluations
```

**What AutoDock Suite Pro Does Automatically for AutoDock 4**:
- **Automatic GPF Generation**: Gathers moving atom types (including flexible sidechain atom types in flex mode), computes grid points from `config.txt`, and generates `receptor.gpf`.
- **AutoGrid 4 Execution**: Runs `autogrid4.exe` to compute potential energy maps and produces the **GLG** (`receptor.glg`).
- **Automatic DPF Generation**: Generates `receptor_ligand.dpf` with chosen algorithm (LGA / GA / LS) and includes the `flexres` directive automatically for flexible docking.
- **AutoDock 4 Execution**: Runs `autodock4.exe` to perform molecular docking and writes the **DLG** (`receptor_ligand.dlg`).
- **BSNDVP™ Post-Docking Extraction**: Automatically parses the DLG, extracts docked ligand poses and moved sidechains into SDF and PDBQT, clusters conformations, classifies binding modes, and produces publication-ready Excel & CSV reports in `results\AD4\DLG_ANALYSIS\BSNDVP_RESULTS\`!

### Q: Vina exits with a non-zero error code
**Solution**: Check the job log in `logs\job_VINA_*.log` — the exact error is logged there. Common causes:
- Receptor PDBQT has atom type errors (re-prepare with `prepare_receptor4.py`)
- Grid box center is outside the receptor (check `config.txt` coordinates)
- Missing or corrupt ligand PDBQT

### Q: AutoGrid4 fails with "cannot find maps"
**Cause**: AutoGrid4 requires the receptor PDBQT to have valid AutoDock atom types.  
**Fix**: Re-prepare the receptor with `prepare_receptor4.py -A hydrogens` and check that no atoms have type `X` or blank types.

### Q: Permission denied when writing the Excel report
**Fix**: Close `VINA_Complete_Report.xlsx` in Excel before running. Windows locks open files.

### Q: My project is on a different drive / network path
**Fix**: Use an absolute path in `project_config.toml`:
```toml
[project]
root = "D:/MyResearch/DrugProject1"

[inputs]
receptor_directory = "D:/MyResearch/DrugProject1/receptors"
ligand_directory   = "D:/MyResearch/DrugProject1/ligands"
```

### Q: How do I dock one compound at a time to test?
```bash
python main.py --receptor protein1 --ligand compound_001 --dry-run
python main.py --receptor protein1 --ligand compound_001
```

### Q: Can I use this on multiple PCs with different software paths?
Yes. The suite auto-detects common MGLTools installation paths. If it can't find the executables, set them explicitly:
```toml
[executables]
vina = "D:/Software/MGLTools/vina.exe"
```
Each PC can have its own `project_config.toml` pointing to the same shared receptor/ligand folder.

---

## 11. Modern Desktop Graphical User Interface (GUI)

AutoDock Suite Pro includes a modern, dark/light-themed graphical interface powered by **CustomTkinter**.

### 11.1 Launching the GUI
You can launch the GUI in three simple ways:
1. **Double-click** `run_gui.bat` in the project root
2. **Double-click** `main.py` directly from Windows Explorer (automatically starts GUI on desktop systems)
3. Run from PowerShell / Command Prompt:
   ```bash
   python main.py --gui
   ```

### 11.2 Level 3 Scientific Platform & Desktop Architecture
AutoDock Suite Pro features a state-of-the-art computational chemistry workstation interface engineered in PySide6/Qt:

#### 1. Centralized Scientific Design System & 5 Thematic Palettes
Every interface element follows a centralized design token scale (`gui_qt/styles.py`) governing typography, spacing, border radii, and accessible high-contrast color tokens:
- **Noir** (Default): Minimalist, distraction-free monochrome workstation with maximum contrast.
- **Dark Studio**: Deep slate and obsidian workstation aesthetic optimized for long computational sessions.
- **Scientific Light**: Publication-grade, journal-ready crisp laboratory aesthetic with dark navy typography.
- **Molecular Plasma**: Sophisticated deep indigo and vibrant violet palette tailored for structural biology.
- **Bio-Neutral**: Restrained, professional emerald and forest-teal aesthetic for biochemical profiling.
- **Live Theme Swatch Preview**: Interactive preview widget in the header and settings displaying background, surface, accent, and text swatches before applying themes.

#### 2. Modern Sidebar Navigation & Header Architecture
- **Sidebar Workflow Navigation**: Intuitive vertical navigation panel organizing the 8 scientific workspaces with accessible icons and hover/selection feedback.
- **Global Application Header**:
  - Embedded branding and subtitle: *Molecular Docking & Virtual Screening Workstation*.
  - **Workspace Indicator Pill**: Shows the active anchored workspace directory.
  - **Global Status Badge**: Real-time semantic status pills (`● READY`, `● RUNNING`, `● COMPLETED`, `● CANCELLED`, `● FAILED`).
  - **Engine Status Indicators**: Rapid visual confirmation of detected Vina and AutoDock4 binaries.
  - **Non-Blocking In-App Notifications**: Subtle transient alert banner for job completions, export events, and warnings without intrusive blocking popups.

#### 3. Laboratory Workspaces
1. **🗂 Workspace Studio**:
   - Laboratory KPI summary cards displaying live counts for `RECEPTORS`, `LIGANDS`, `COMPLETED JOBS`, and `ACTIVE / READY`.
   - Receptor and ligand inventory tables with inline deletion, directory browsing, and configuration anchoring.
2. **🧪 Molecular Preparation Studio**:
   - **Visual Progression Stepper**: High-visibility pipeline indicator showing progression stages:
     `INPUT → INSPECTION → CLEANING → CHARGES / H → FLEXIBLE RESIDUES → VALIDATION → READY`
   - Active-site residue selection and automated rigid/flexible sidechain splitting.
3. **🎯 Grid Box Configuration**:
   - Precise 3D coordinate bounding box editor with calculated volume in Å³.
   - Automated GPF and Vina `config.txt` generation and synchronization.
4. **📚 Screening Library & ADMET Filter**:
   - In-silico physicochemical descriptor calculator (MW, LogP, TPSA, HBD, HBA, RotBonds).
   - Drug-likeness evaluation against Lipinski Rule of 5 and Veber bioavailability criteria.
5. **▶ Docking Execution Console & Timeline**:
   - **Experiment Configuration & Monitor**: Live summary of active Engine, Docking Mode, Exhaustiveness, Modes, and Output Directory.
   - **Dual-Engine Execution Status**: When running in `BOTH` mode, separate badges track `Vina` and `AutoDock4` stages independently (`● Running`, `✓ Complete`, `○ Pending`).
   - **Scientific Execution Timeline**: Live stage stepper tracking `Preparation → Grid Verification → Docking Search → DLG Parsing → Biophysical Analysis → Report Generation`.
6. **📊 Results & Biophysical Analytics Workstation**:
   - **Overview Tab**: Key metrics summary card, lead pose information, and embedded **Vector SVG Conformation Energy Spectrum** plot.
   - **Poses Tab**: Sortable data grid with engine badges, binding energies, estimated $K_i$, RMSD, and one-click PyMOL/WebGL inspection.
   - **Clusters Tab (AD4)**: Cluster population table, representative pose energy, and embedded **Vector SVG Cluster Population Distribution** bar chart.
   - **Validation Tab**: Crystallographic reference ligand superposition, matched atom counts, and reference RMSD assessment.
   - **Interactions Tab**: 2D interaction diagram and non-covalent contact tables (H-bonds, salt bridges, $\pi$-stacking, halogen bonds).
   - **Thermodynamics Tab**: Statistical mechanics ensemble metrics (Shannon information entropy, partition function $Q$, statistical temperature, internal energy, Boltzmann lead probability).
   - **ADMET Tab**: Structured physicochemical report and drug-likeness compliance.
   - **Provenance Tab**: Full cryptographic SHA-256 manifest and engine provenance record with one-click clipboard copy.
7. **🔬 Experiment Comparison Studio**:
   - Cross-experiment analysis workspace comparing docking jobs across different receptors, ligands, engines, and docking modes.
   - **Scoring Delta Analysis ($\Delta\Delta G = \text{AD4} - \text{Vina}$)** to contrast force-field vs. empirical scoring functions.
   - Real-time filtering by ligand search, target receptor, engine, and binding energy cutoff with CSV export.
8. **⚙ Modern Grouped Settings**:
   - Structured card layout for Appearance, Docking Engines (with binary resolution and live version verification), Biophysical Analysis, Reporting, and Reproducibility.

#### 4. Publication-Quality Research Reporting
- **Standalone Scientific HTML Report** (`generate_html_report()`):
  - Formatted as a peer-reviewed research document with clean typography, numbered sections, and executive summary metadata.
  - Features inline resolution-independent vector SVG figures with professional scientific captions:
    * *Figure 1. Binding affinity distribution across top-ranked docked poses.*
    * *Figure 2. Conformational cluster populations and representative binding energies.*
  - Comprehensive coverage: Target & Experiment Summary, Top Docking Conformations, Conformational Cluster Analysis, Reference Ligand Validation, Predicted Non-Covalent Interactions, Thermodynamic & Statistical Mechanics, ADMET Profiling, Reproducibility & Provenance, and Scientific Disclaimers.
  - Dedicated print and PDF stylesheets (`@media print`) for 1-click archiving and document generation.
- **Multiformat Data Exports**: XLSX multi-sheet workbooks, per-analysis CSV files, and canonical JSON records.



---

## 12. Protein–Ligand Interaction Profiler

AutoDock Suite Pro features a **Pure Python 3 Interaction Profiler** (`interactions.py`) that analyzes 3D coordinates directly from PDBQT files without requiring PyMOL, Discovery Studio, or external licenses.

### 12.1 Supported Molecular Interactions
| Interaction Type | Geometrical Criteria | Key Residues / Atoms |
| :--- | :--- | :--- |
| **Hydrogen Bonds** | Donor-Acceptor $d \le 3.5\text{ \AA}$, polar H distance $d \le 2.5\text{ \AA}$ | ASN, GLN, SER, THR, TYR, TRP, LYS, ARG, backbone N/O |
| **Hydrophobic Contacts** | Non-polar Carbon-to-Carbon $d \le 4.0\text{ \AA}$ | ALA, VAL, LEU, ILE, PRO, PHE, MET, TRP, TYR |
| **Salt Bridges** | Cation to Anion distance $d \le 4.0\text{ \AA}$ | LYS (NZ), ARG (NH1, NH2), protonated HIS with carboxylate/phosphate |
| **$\pi$-Stacking** | Aromatic Ring Centroid-to-Centroid $d \le 5.0\text{ \AA}$ | PHE, TYR, TRP, HIS rings with ligand aromatic rings |
| **$\pi$-Cation** | Positive Group to Aromatic Centroid $d \le 4.5\text{ \AA}$ | LYS, ARG, or ligand positive amine to aromatic ring |
| **Halogen Bonds** | Ligand Halogen (Cl, Br, I) to Lewis base (O, N, S) $d \le 3.8\text{ \AA}$ | Halogenated ligands interacting with backbone or sidechain |

### 12.2 Flexible Residue Awareness
When flexible docking is used (Vina or AutoDock 4), the profiler **automatically substitutes the flexible residue coordinates** output by the docking engine for those residues, ensuring 100% true physical contacts are calculated!

### 12.3 Reports & Output Files
- **Master Interaction Table**: Saved to `reports/VINA_04_Interactions.csv` (or `reports/AD4_02_Interactions.csv`) and embedded as a dedicated worksheet in `Complete_Report.xlsx`.
- **Per-Ligand CSV**: Saved directly in each result folder as `results/<ENGINE>/<RECEPTOR>/<LIGAND>/interactions.csv`.
- **Standalone CLI Usage**:
  ```bash
  python interactions.py -r receptors/2V5Z/rigid/2v5z_clean.pdbqt -l results/VINA/2V5Z_flex/matrine/output/matrine_out.pdbqt -o interactions.csv
  ```

---

## 13. Multi-Core Parallel Docking

Speed up large batch virtual screening runs by running jobs concurrently across multiple CPU cores.

### 13.1 Enabling Parallel Execution
1. **Via Command Line**:
   ```bash
   # Run with 4 concurrent docking workers
   python main.py --parallel 4

   # Run with AUTO workers (half of logical cores)
   python main.py --parallel
   ```
2. **Via `project_config.toml`**:
   ```toml
   [execution]
   sequential  = false
   max_workers = 4     # or "AUTO"
   ```
3. **Via GUI**:
   Adjust the **Parallel Workers** slider on the Dashboard.

### 13.2 Thread-Safe Resume & Progress Tracking
The worker pool uses a synchronized lock for saving `results/job_status.json` and console progress updates. If interrupted at any time, run with `--resume` to continue immediately without losing completed jobs.

---

## 14. Standalone Portable Executable (`.exe`) Compilation

You can package the entire AutoDock Suite Pro application (Python runtime + CustomTkinter + binaries + scripts) into a standalone Windows folder with `AutoDockSuitePro.exe` that runs on any Windows PC with zero installation.

### 14.1 Building the Executable
1. Simply double-click `build.bat`, or run from terminal:
   ```bash
   python build_exe.py
   ```
2. The compilation will:
   - Install PyInstaller if missing
   - Discover local engine binaries (`vina.exe`, `vina_split.exe`, `autodock4.exe`, `autogrid4.exe`) and bundle them directly into the distribution (`bin/` and `_internal/bin/`)
   - Bundle all scientific dependencies: `meeko`, `rdkit`, `openbabel`, `scipy`, `numpy`, `openpyxl`, and `customtkinter`
   - Bundle `project_config.toml`, `receptors/`, and `ligands/`
   - Produce a 100% self-contained portable distribution folder at:
     ```
     dist\AutoDockSuitePro\
     ```
3. To distribute, simply copy or zip the `dist\AutoDockSuitePro\` folder. Double-clicking `AutoDockSuitePro.exe` launches the full modern GUI with all docking engines, molecular preparation tools, and 3D viewers functioning out-of-the-box!

---

## 15. Molecular Preparation Studio (AutoDockTools & PyRx Python 3 Upgrade)

AutoDock Suite Pro includes a native **Molecular Preparation Studio** powered by Scripps' modern Python 3 **Meeko 0.7+**, **OpenBabel**, and **RDKit**. You no longer need legacy Python 2 MGLTools (`prepare_receptor4.py`, `prepare_ligand4.py`) or external command-line tools.

### 15.1 Prepare Macromolecule (Receptor PDB/CIF → PDBQT)
- **What it does**:
  1. Automatically strips crystallographic waters (`HOH`, `WAT`) and buffer co-solvents (`EDO`, `GOL`, `SO4`, `PO4`, `DMS`).
  2. Adds polar hydrogens at physiological pH (7.4).
  3. Computes Gasteiger partial atomic charges.
  4. Assigns AutoDock 4 atom types (`C`, `A`, `OA`, `SA`, `HD`, `NA`, etc.).
  5. Outputs a clean, rigid `.pdbqt` directly into `receptors/<NAME>/rigid/<NAME>.pdbqt`.
- **In the GUI**:
  - Open the **Molecular Preparation** tab.
  - Under *Prepare Macromolecule*, browse to your `.pdb` or `.cif` file.
  - Click **Convert to Receptor PDBQT**.

### 15.2 Prepare Ligands (SDF / MOL2 / PDB / SMI → PDBQT)
- **What it does**:
  1. Reads chemical structures in SDF, MOL2, PDB, or SMILES format.
  2. Generates 3D conformers with MMFF energy minimization.
  3. Automatically computes Gasteiger charges and merges non-polar hydrogens.
  4. Detects rotatable bonds and constructs the hierarchical torsion tree (`ROOT`, `BRANCH`, `TORSDOF`) using Scripps Meeko.
  5. Saves directly into `ligands/<NAME>.pdbqt` and calculates physicochemical descriptors.
- **In the GUI**:
  - Open the **Molecular Preparation** tab.
  - Under *Prepare Ligand*, browse to your compound file.
  - Click **Convert to Flexible Ligand PDBQT**.

### 15.3 Flexible Residues Splitter ("Make Residues Flexible")
- **What it does**:
  1. Inspects the receptor and extracts all amino acid residues.
  2. Selects target active-site residues (e.g. `ILE 199`, `TYR 326`, `GLN 206`).
  3. Splits the receptor into:
     - `receptors/<NAME>/flex/<NAME>_rigid.pdbqt` (retains backbone N, CA, C, O)
     - `receptors/<NAME>/flex/<NAME>_flex.pdbqt` (sidechains with `BEGIN_RES`, `ROOT`, `BRANCH`, `END_RES`)
- **In the GUI**:
  - Open the **Molecular Preparation** tab.
  - Under *Flexible Receptor Splitter*, select your receptor from the dropdown.
  - Enter the target residues to make flexible (comma-separated).
  - Click **Generate Rigid & Flexible PDBQT Pair**.

---

## 16. Screening Library & Chemical Descriptors Explorer

AutoDock Suite Pro provides a dedicated cheminformatics explorer (**Screening Library** tab) showing:
- **Molecular Weight (MW)**: In g/mol
- **Calculated LogP (SlogP)**: Octanol-water partition coefficient
- **H-Bond Donors (HBD) & Acceptors (HBA)**
- **Rotatable Torsions (TORSDOF)**: Flexibility metric
- **Heavy Atom Count**: Non-hydrogen atoms (used for Ligand Efficiency calculation)
- **Lipinski's Rule of 5 Compliance**: Highlights lead-like drug candidates (`Pass` or `Fail`)

---

## 17. Interactive 3D Visualization: WebGL & Native PyMOL Integration

AutoDock Suite Pro includes comprehensive 3D molecular visualization capabilities through both an in-browser WebGL engine and full native **PyMOL** desktop integration (`pymol_exporter.py`, `viewer_3d.py`).

### 17.1. In-Browser 3D WebGL Inspector (3Dmol.js)

Click **🌐 WebGL** on the Results tab:
- **Zero-Setup Inspection**: Launches directly in your default web browser without needing third-party software.
- **Docked Pose Switcher**: Toggle through all docked poses (Modes 1 to $N$) with instant re-centering.
- **Pocket Residues**: Active site target residues highlighted as orange sticks with residue sequence labels.
- **Real-Time H-Bond Dashes**: Yellow dashed lines with distance measurements for hydrogen bond contacts ($d \le 3.5\text{ \AA}$).
- **Docking Grid Box**: 3D wireframe bounding box overlaying the binding pocket cavity.
- **Surface Display**: Toggle semi-transparent molecular surface.
- **Export PyMOL Script**: 1-click button in the WebGL toolbar to download a PyMOL `.pml` script.

---

### 17.2. Native PyMOL 3D Session & Script Integration

AutoDock Suite Pro connects directly with **PyMOL** (Open-Source and Incentive editions by Schrödinger):

#### 1. One-Click Desktop Launch ("🔬 Open in PyMOL")
- Automatically detects installed PyMOL on Windows (`%LOCALAPPDATA%\Schrodinger\PyMOL`, `Program Files`, Anaconda/Miniconda) as well as macOS and Linux.
- Compiles a native binary PyMOL session (`.pse`) and opens the PyMOL desktop application with your entire docking result pre-configured and oriented.
- If PyMOL is located in a custom path, you can select `PyMOL.exe` once and the suite remembers it.

#### 2. Complete Session Export ("💾 Export PyMOL")
- Exports a complete session bundle into your chosen destination:
  * `<receptor>_<ligand>_docking.pml`: Universal, reproducible PyMOL macro script.
  * `<receptor>_<ligand>_docking.pse`: Compiled native binary PyMOL session file.
  * Standardized clean PDB coordinates for both receptor and multi-model ligand poses.

#### 3. Complete Result Representation in PyMOL
- **Receptor Architecture**: Secondary-structure cartoon coloring with alpha helices and beta sheets.
- **Toggleable Pocket Cavity Surface**: Semi-transparent VDW molecular surface (`receptor_surface`, 65% transparency) toggleable in the PyMOL object tree.
- **Docked Ligand Poses**:
  * Active pose rendered in vibrant lime-green carbons with ball-and-stick sphere overlays.
  * Alternative poses organized in the `Alternative_Poses` PyMOL group (hidden by default, easily enabled).
- **Active Site Residues**: All binding pocket residues within $4.5\text{ \AA}$ rendered as clean gray/white sticks with Helvetica residue labels (`ARG 100`, `GLU 84`).
- **Comprehensive Non-Covalent Interactions** (with measured $\text{\AA}$ distances):
  * **Hydrogen Bonds**: Yellow dashed lines ($3.2\text{ pt}$ width) with measured distances.
  * **Salt Bridges**: Magenta dashed lines ($3.5\text{ pt}$ width).
  * **$\pi$-Stacking & $\pi$-Cation**: Cyan dashed lines ($2.8\text{ pt}$ width).
  * **Halogen Bonds**: Bright orange dashed lines ($3.2\text{ pt}$ width).
  * **Hydrophobic Contacts**: Slate dashed lines organized in the `Hydrophobic_Contacts` group.
  * All interactions organized in a master `Interactions` group in the PyMOL panel for 1-click toggling.
- **3D Wireframe Docking Grid Box**: Compiled Graphics Object (`Docking_Grid_Box`) matching the exact center $(x, y, z)$ and dimensions $(s_x, s_y, s_z)$ of your search space.
- **Publication Presets**: Pre-configured depth-cueing, ray-tracing shadows disabled (`ray_shadows, 0`), dark background (switch to white for journal papers with `bg_color white`), and optimal camera focus.

---

### 17.3. Command-Line PyMOL Export

You can also generate PyMOL sessions directly from the terminal or in automated scripts:

```bash
python pymol_exporter.py -r receptors/2V5Z/rigid/2V5Z.pdbqt -l results/VINA/2V5Z/berberine/output/berberine_out.pdbqt -p 1 --cx 51.886 --cy 156.453 --cz 28.559 --sx 25.0 --sy 25.0 --sz 25.0 -o pymol_export --launch
```

Parameters:
- `-r / --receptor`: Path to receptor PDB or PDBQT.
- `-l / --ligand`: Path to docked ligand pose PDBQT.
- `-p / --pose`: Pose index (default: 1).
- `--cx, --cy, --cz`: Grid box center coordinates.
- `--sx, --sy, --sz`: Grid box dimensions in Ångströms.
- `-o / --output-dir`: Output folder for `.pml`, `.pse`, and coordinate files.
- `--launch`: Automatically open the resulting session in PyMOL.


