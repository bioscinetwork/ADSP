# ADSP Benchmark Suite

This directory contains curated crystallographic and structural benchmark fixtures used for verifying AutoDock4 and AutoDock Vina workflows in ADSP.

## Benchmark Index

| Benchmark | System | Engine Verification | Scientific Role |
|-----------|--------|---------------------|-----------------|
| 1CA2 | Human Carbonic Anhydrase II (Zn²⁺) | AutoDock4, Vina | Metalloprotein typing and grid centering validation |
| 1MBN | Sperm Whale Myoglobin (Fe²⁺/HEM) | AutoDock4, Vina | Heme-containing metalloprotein grid transparency and typing |
| 2NV6 | *Mycobacterium tuberculosis* InhA (INH-NAD) | Inspection / Verification | Antitubercular adduct benchmarking; isolates crystallographic extraction from independent SDF |

## Data Policy

- All benchmark structures are reproducible input fixtures.
- Generated docking outputs, temporary files, and grid maps are excluded from Git.
- Parameter checksums and coordinates are rigorously verified.
