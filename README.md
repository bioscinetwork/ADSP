# AutoDock Suite Platform (ADSP) v0.3.0

AutoDock Suite Platform (ADSP) is a research-oriented molecular docking workflow platform providing unified automation, preparation, docking execution, pose parsing, interaction analysis, complex building, and auditable reporting for **AutoDock4** and **AutoDock Vina**.

ADSP is designed for reproducible computational chemistry and structural biology research, featuring hardware-accelerated PySide6 desktop GUI and comprehensive CLI automation.

---

## Key Capabilities & Architecture

- **Dual-Engine Architecture**: Independent, isolated execution pathways for **AutoDock Vina** (empirical scoring, Monte Carlo global optimization) and **AutoDock4** (semi-empirical force field, Lamarckian Genetic Algorithm).
- **Strict Thermodynamic Isolation**:
  - **AutoDock Vina**: Binding affinity strictly reported in kcal/mol. No synthetic thermodynamic inference ($K_i$, $\Delta H$, $\Delta S$ are never fabricated).
  - **AutoDock4**: Binding free energy and inhibition constant ($K_i$) are strictly DLG-driven from AutoDock4's internal thermodynamic analysis.
- **PDBQT & Structure Validation**: Fixed-column format validation enforcing atom names, coordinates, partial charges, and atom type fields with Open Babel round-trip compatibility.
- **Graph-Constrained RMSD**: Hungarian algorithm automorphism symmetry matching and Maximum Common Substructure (MCS) via RDKit/obrms. Zero coordinate-order fallback; structural mismatch raises explicit failures rather than reporting misleading 0.0 Å.
- **Auditable Reporting**: Automated multi-format reports (CSV, XLSX, JSON, TXT, HTML) with complete provenance tracking (engine, hashes, parameters, timestamps).

---

## Repository Structure

```
D:\Projects\ADSP
├── main.py                     # Unified CLI entry point
├── config.py                   # Project configuration and directory resolution
├── models.py                   # Canonical data schemas and provenance models
├── prepare.py                  # Macromolecule and ligand preparation engine
├── validators.py               # Fixed-column PDBQT and geometry validation
├── dlg_extract.py              # Canonical AutoDock4 DLG extraction engine
├── vina_workflow.py            # Vina docking workflow orchestration
├── autodock4_workflow.py       # AutoDock4 & AutoGrid4 workflow orchestration
├── reporting.py                # Multi-format report generation (CSV/XLSX/JSON/TXT)
├── interactions.py             # Protein-ligand interaction profiler (HB, vdW, pi)
├── complex_builder.py          # Receptor-ligand complex generation & validation
├── executables.py              # Portable runtime binary & toolchain resolution
│
├── gui_qt/                     # Canonical PySide6 (Qt6) hardware-accelerated GUI
├── gui/                        # Legacy/fallback desktop GUI and icon assets
│
├── tests/                      # Automated test suite (124 pytest tests)
│   ├── test_suite.py
│   ├── test_packaging.py
│   ├── test_ad4_compatibility.py
│   └── ...
│
├── benchmarks/                 # Curated crystallographic benchmark fixtures
│   ├── 1CA2/                   # Human Carbonic Anhydrase II (Zn²⁺)
│   ├── 1MBN/                   # Sperm Whale Myoglobin (Fe²⁺/HEM)
│   └── 2NV6/                   # M. tuberculosis InhA / INH-NAD adduct
│
├── data/                       # Structural data directories
│   ├── receptors/              # Target macromolecule files
│   ├── ligands/                # Compound libraries
│   └── reference/              # Reference ligand conformers (ADP, HRM, Donepezil)
│
├── parameter_profiles/         # Cryptographically verified parameter profiles
│   ├── ad4_standard_4.2/       # AD4_parameters.dat (SHA-256 verified)
│   └── ad4_1_bound/            # AD4.1_bound.dat (SHA-256 verified)
│
├── external/
│   └── autodock4zn/            # Isolated AutoDock4Zn specialized pipeline archive
│
├── bin/                        # Bundled Windows executables (Vina, AD4, AutoGrid4, Open Babel)
├── docs/                       # Comprehensive documentation (USER_GUIDE, DEVELOPER_GUIDE)
├── scripts/                    # Build, setup, and execution batch scripts
├── AutoDockSuitePro.spec       # Portable PyInstaller packaging specification
├── build_exe.py                # Standalone distribution packaging script
├── pyproject.toml              # Build system & packaging configuration
├── project_config.toml         # Portable configuration file
├── CITATION.cff                # Software citation metadata
├── LICENSE                     # GNU General Public License v2.0+
└── README.md
```

---

## Quick Start

### Prerequisites & Installation

ADSP requires Python 3.10+ (64-bit) on Windows.

1. Clone or download the repository:
   ```bash
   git clone https://github.com/your-org/ADSP.git
   cd ADSP
   ```

2. Install dependencies:
   ```bash
   scripts\install_dependencies.bat
   ```
   Or via pip:
   ```bash
   pip install -r requirements.txt
   ```

### Launching the Graphical Interface

Double-click `run_gui.bat` or run:
```bash
python main.py --gui
```

### CLI Docking Workflows

Initialize a project configuration:
```bash
python main.py --init
```

Validate environment, toolchains, and inputs:
```bash
python main.py --validate --config project_config.toml
```

Execute AutoDock Vina docking:
```bash
python main.py --config project_config.toml --engine VINA
```

Execute AutoDock4 docking:
```bash
python main.py --config project_config.toml --engine AUTODOCK4
```

---

## Scientific Parameter Profiles

ADSP bundles authoritative AutoDock4 parameter profiles verified by SHA-256 digests:

| Profile ID | Version | Parameter File | Reference / Source | SHA-256 Checksum |
|------------|---------|----------------|--------------------|------------------|
| `ad4_standard_4.2` | 4.2 | `AD4_parameters.dat` | AutoDock4 Force Field | `625DE5779B914382E21A135C776EFBC02B4221085BD0280118D103CCDD93EA7C` |
| `ad4_1_bound` | 4.1-bound | `AD4.1_bound.dat` | Huey et al. (2007) *J Comput Chem* | `6B98F7AB508F4882801938F8CED1C0BF38096496155A8005BAF941A201781CE8` |

### AutoDock4Zn Status & Isolation

The AutoDock4Zn zinc coordination pipeline (`external/autodock4zn/AutoDock4Zn-Pipeline-main.zip`) is isolated from standard AD4 runs. Standard AD4 parameter files include non-bonded Zn and Fe atom types; they do not apply specialized coordination pseudo-atoms unless explicitly configured through specialized legacy pipelines.

---

## Standalone Portable Packaging

ADSP can be compiled into a fully portable, zero-dependency Windows distribution:

```bash
python build_exe.py
```

This generates:
- Standalone folder: `dist\AutoDockSuitePro\`
- Portable ZIP: `dist\AutoDockSuitePro_Portable.zip`

The portable distribution contains all Python runtimes, Qt libraries, Open Babel binaries, Vina, AutoDock4, AutoGrid4, parameter profiles, and default configs. It runs out-of-the-box on clean Windows systems without Python installed and does not depend on developer paths.

---

## Testing & Quality Assurance

Run the comprehensive automated test suite (124 tests):

```bash
pytest
```

---

## Scientific Scope & Limitations

1. **Thermodynamics**: Vina calculates empirical binding affinities (kcal/mol), not thermodynamic free energies ($\Delta G$). Inhibition constants ($K_i$) are only derived when executing AutoDock4 with real DLG outputs.
2. **Charges**: Kollman partial charges assigned via ADSP fallback mode are structural approximations and should not be conflated with canonical quantum mechanical ESP charges.
3. **2NV6 / INH-NAD**: Crystallographic adduct extraction (`ZID.cif` / `2NV6 (2).pdb`, 52 heavy atoms / 82 explicit atoms) is deliberately kept distinct from independent SDF models (`2nv6_B_ZID.sdf`). False redocking validation is never claimed.

---

## Citation

If you use ADSP in academic research, please cite:

```bibtex
@software{adsp2026,
  author = {ADSP Contributors},
  title = {AutoDock Suite Platform (ADSP): Unified Molecular Docking Platform},
  version = {0.3.0},
  year = {2026},
  url = {https://github.com/your-org/ADSP}
}
```

Please also cite the underlying engine and toolchain publications:
- **AutoDock4**: Morris et al. (2009), *J Comput Chem* 30:2785-2791.
- **AutoDock Vina**: Trott & Olson (2010), *J Comput Chem* 31:455-461; Eberhardt et al. (2021), *J Chem Inf Model* 61:3891-3898.
- **Open Babel**: O'Boyle et al. (2011), *J Cheminform* 3:33.
- **RDKit**: RDKit: Open-source cheminformatics (https://www.rdkit.org).
- **Meeko**: Forli Lab, Scripps Research (https://github.com/forlilab/Meeko).

---

## License

ADSP is licensed under the **GNU General Public License v2.0 or later** ([`LICENSE`](LICENSE)).
Bundled third-party binaries and dependencies retain their respective licenses (see [`LICENSE`](LICENSE) and [`.licenses/`](.licenses/)).
