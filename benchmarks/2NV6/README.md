# 2NV6 Benchmark: InhA / INH-NAD Adduct

- **Organism**: *Mycobacterium tuberculosis*
- **Ligand**: INH-NAD covalent adduct (ZID) bound to enoyl-acyl carrier protein reductase (InhA)
- **Scientific Role**:
  - Establishes clear separation between crystallographically extracted adduct (ZID from 2NV6 (2).pdb) and independently modeled SDF (2nv6_B_ZID.sdf).
  - Benchmarks coordinate validation and exposes explicit-hydrogen representation differences without claiming false redocking validation.
- **Files**:
  - 2NV6 (2).pdb / 2NV6 (1).cif: PDB and mmCIF entries
  - 2nv6_B_ZID.sdf: Independent SDF representation
  - ZID.cif: Chemical Component Dictionary entry for ZID
  - enchmark_2nv6.py: Verification and comparison script
