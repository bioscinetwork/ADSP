# AutoDock Suite Platform (ADSP)

ADSP is a research-oriented workflow scaffold for preparing receptors and ligands, running AutoDock4 or AutoDock Vina, parsing poses, validating against references, analysing interactions, generating complexes, and exporting reports.

## Scientific scope

ADSP preserves ordinary AutoDock4 and Vina workflows. Metal atoms can be retained and passed through standard AutoDock4 typing when the bundled parameter profile recognizes the atom type. This is structural/standard parameter support; it is not a specialized metal-coordination force field. AutoDock4Zn is detected as a separate, environment-dependent workflow and is not claimed to be executable unless its Python 2/MGLTools dependencies are available.

## Inputs and outputs

Curated benchmark inputs are kept at the repository root while the benchmark harness records SHA-256 checksums. Generated docking outputs belong in `results/`, `reports/`, and `logs/` and are ignored by Git by default. Scientific parameter files, executable assets, benchmark structures, small ligand fixtures, source code, and tests are intended to remain version controlled.

## Engines and preparation

- AutoDock Vina uses prepared PDBQT files and Vina configuration through `vina_workflow.py`.
- AutoDock4 uses GPF generation, AutoGrid4, DPF generation, AutoDock4, and DLG extraction through `autodock4_workflow.py` and `dlg_extract.py`.
- Receptor and ligand preparation is implemented in `prepare.py`, using available OpenBabel/RDKit/Meeko paths.
- CIF/mmCIF preparation prefers `gemmi` for coordinate and author chain/residue preservation, then uses OpenBabel for PDBQT generation.

## Benchmarks

The 2NV6 benchmark uses the crystallographic `ZID` INH–NAD adduct (`2NV6 (2).pdb` / `2NV6 (1).cif`) and the supplied independent `2nv6_B_ZID.sdf` as separate inputs. The benchmark records whether each scenario is available and never fabricates missing structures.

Run the regression suite with:

```text
python test_suite.py
python -m pytest -q
```

## Limitations

Charge mode named `kollman` is currently an ADSP fallback and must not be interpreted as guaranteed canonical AutoDockTools Kollman preparation. CIF conversion can still expose representation differences that require review. RMSD and interaction results are software outputs and require scientific interpretation.

## Citations

See `CITATION.cff`. Cite the underlying AutoDock4, AutoDock Vina, RDKit, Meeko, Open Babel, and any AutoDock4Zn reference implementation used for a particular run.
