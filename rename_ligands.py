#!/usr/bin/env python3
"""
Docking Automation Suite — Ligand Renaming Utility
====================================================
Reads compound names from SDF files (or fetches them from PubChem)
and renames ligand files (PDBQT or SDF) from PubChem CIDs to human-readable
compound names.

Problem:
    When SDF files are downloaded from PubChem (e.g. Structure3D_CID_2244.sdf)
    or imported into PyRx for PDBQT preparation (e.g. 2244_uff_e=-535.63.pdbqt),
    files are named with the numeric CID rather than meaningful chemical names
    like "Aspirin". Docking results and files then show CIDs instead of names.

Solution:
    This utility parses compound identifiers and properties directly from
    SDF files (or queries the PubChem API) to rename files:
      1. Rename SDF files in-place before docking (e.g. Structure3D_CID_2244.sdf -> aspirin.sdf)
      2. Rename PDBQT files prepared by PyRx before running docking
      3. Fully reversible with automatic rename_log.csv and --undo

Name Resolution Priority:
    1. Manual CSV mapping (--csv), if provided
    2. PubChem API common synonyms (--fetch-names)
    3. SDF property: PUBCHEM_IUPAC_TRADITIONAL_NAME
    4. SDF property: PUBCHEM_IUPAC_OPENEYE_NAME
    5. SDF property: PUBCHEM_IUPAC_NAME
    6. SDF molecule title (if it's not a pure number / CID)
    7. Original filename kept unchanged

Usage:
    # ── SDF Renaming (before docking or standalone) ──────────────────────
    # Rename SDF files downloaded from PubChem in-place
    python rename_ligands.py --sdf-dir ./my_sdfs/

    # Rename SDFs using common names from PubChem API (e.g. Aspirin instead of IUPAC)
    python rename_ligands.py --sdf-dir ./my_sdfs/ --fetch-names

    # Preview SDF renaming without changing any files
    python rename_ligands.py --sdf-dir ./my_sdfs/ --dry-run

    # ── PDBQT Renaming (PyRx outputs) ───────────────────────────────────
    # Extract names from original SDF and rename PyRx PDBQTs
    python rename_ligands.py --sdf compounds.sdf --pdbqt-dir ligands/

    # Fetch common names via PubChem API for PDBQTs
    python rename_ligands.py --sdf compounds.sdf --pdbqt-dir ligands/ --fetch-names

    # ── Common Workflows ────────────────────────────────────────────────
    # Use a manual CSV mapping file (overrides all other sources)
    python rename_ligands.py --csv mapping.csv --pdbqt-dir ligands/

    # Undo a previous rename using the generated log (works for both SDF and PDBQT)
    python rename_ligands.py --undo rename_log.csv --sdf-dir ./my_sdfs/
    python rename_ligands.py --undo rename_log.csv --pdbqt-dir ligands/
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.request
import urllib.error
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Configure UTF-8 output on Windows consoles
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════════════

BANNER = "Docking Automation Suite — Ligand Renaming Utility"

# SDF property names to try for compound names (in priority order)
NAME_PROPERTIES = [
    "PUBCHEM_IUPAC_TRADITIONAL_NAME",
    "PUBCHEM_IUPAC_OPENEYE_NAME",
    "PUBCHEM_IUPAC_NAME",
    "PUBCHEM_IUPAC_CAS_NAME",
    "PUBCHEM_IUPAC_SYSTEMATIC_NAME",
]

# SDF property for the CID
CID_PROPERTY = "PUBCHEM_COMPOUND_CID"

# PubChem REST API base URL
PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"

# Maximum CIDs per PubChem API request (to stay within URL length limits)
PUBCHEM_BATCH_SIZE = 100

# Maximum filename length (conservative for Windows paths)
MAX_FILENAME_LENGTH = 80

# Characters invalid in Windows filenames
INVALID_FILENAME_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')

# PyRx appends UFF minimization energy to filenames:
#   e.g. "2353_uff_e=-535.63.pdbqt" or "hrm_uff_e=-284.03.pdbqt"
# This regex captures the identifier before the suffix.
PYRX_UFF_SUFFIX = re.compile(r'^(.+?)_uff_e=[\-]?[\d.]+$')

# PubChem downloaded structure filename patterns:
# Examples:
#   Structure3D_CID_2244.sdf
#   Structure2D_CID_2244.sdf
#   Conformer3D_CID_2244.sdf
#   Conformer2D_CID_2244.sdf
#   Compound_2244.sdf
#   CID_2244.sdf, cid_2244.sdf, 2244.sdf
#   Structure3D_CID_2244_uff_e=-284.03.pdbqt (if converted by PyRx)
PUBCHEM_STEM_PATTERN = re.compile(
    r'^(?:Structure[23]D_CID_|Conformer[23]D_CID_|Compound_CID_|'
    r'Structure[23]D_|Conformer[23]D_|Compound_|CID_|cid_|lig_|ligand_)?'
    r'(\d+)'
    r'(?:_uff_e=[\-]?[\d.]+)?$',
    re.IGNORECASE
)


# ═══════════════════════════════════════════════════════════════════════════════
# Data Classes
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class MoleculeInfo:
    """Information extracted from a single molecule in an SDF file."""
    title: str = ""                           # First line of the MOL block (often the CID)
    cid: str = ""                             # PubChem CID (from property or title)
    properties: Dict[str, str] = field(default_factory=dict)
    source_file: str = ""                     # Which SDF file this came from

    @property
    def best_sdf_name(self) -> Optional[str]:
        """Return the best compound name available from SDF properties."""
        for prop in NAME_PROPERTIES:
            value = self.properties.get(prop, "").strip()
            if value:
                return value
        return None

    @property
    def identifier(self) -> str:
        """The identifier used for PDBQT filename matching (CID or title)."""
        return self.cid or self.title


@dataclass
class RenamePlan:
    """A planned file rename operation."""
    original_path: Path                       # Current file path
    new_path: Path                            # Desired file path
    original_name: str                        # Original filename (e.g., "Structure3D_CID_2244.sdf")
    new_name: str                             # New filename (e.g., "aspirin.sdf")
    name_source: str                          # Where the name came from
    cid: str = ""                             # PubChem CID (if known)
    status: str = "PENDING"                   # PENDING, RENAMED, SKIPPED, ERROR
    message: str = ""                         # Status message


# ═══════════════════════════════════════════════════════════════════════════════
# Filename & Identifier Extraction
# ═══════════════════════════════════════════════════════════════════════════════

def extract_compound_identifier(stem: str) -> str:
    """Extract compound identifier or CID from various naming conventions.

    Handles:
      - PubChem downloads:
          Structure3D_CID_2244  -> 2244
          Structure2D_CID_2244  -> 2244
          Conformer3D_CID_2244  -> 2244
          CID_2244              -> 2244
          cid_2244              -> 2244
          Compound_2244         -> 2244
          2244                  -> 2244
      - PyRx UFF energy suffix:
          2353_uff_e=-535.63    -> 2353
          Structure3D_CID_2244_uff_e=-284.03 -> 2244
          hrm_uff_e=-284.03     -> hrm
      - Common prefixes:
          lig_2244              -> 2244
          ligand_2244           -> 2244

    Args:
        stem: Filename stem (without extension).

    Returns:
        The extracted identifier or CID.
    """
    # 1. Match numeric PubChem / PyRx pattern
    match = PUBCHEM_STEM_PATTERN.match(stem)
    if match:
        return match.group(1)

    # 2. Match PyRx UFF suffix with non-numeric identifier (e.g. hrm_uff_e=-284.03)
    match_pyrx = PYRX_UFF_SUFFIX.match(stem)
    if match_pyrx:
        identifier = match_pyrx.group(1)
    else:
        identifier = stem

    # 3. Strip common prefixes if still present
    for prefix in [
        "Structure3D_CID_", "Structure2D_CID_",
        "Conformer3D_CID_", "Conformer2D_CID_",
        "Compound_CID_", "Compound_",
        "CID_", "cid_", "lig_", "ligand_"
    ]:
        if identifier.startswith(prefix):
            identifier = identifier[len(prefix):]
            break

    return identifier

# Alias for backwards compatibility
extract_pyrx_identifier = extract_compound_identifier


# ═══════════════════════════════════════════════════════════════════════════════
# SDF Parsing
# ═══════════════════════════════════════════════════════════════════════════════

def parse_sdf(sdf_path: Path) -> List[MoleculeInfo]:
    """Parse an SDF file and extract molecule information.

    Handles both single-molecule and multi-molecule (concatenated) SDF files.
    Each molecule record ends with '$$$$'.

    Args:
        sdf_path: Path to the SDF file.

    Returns:
        List of MoleculeInfo objects, one per molecule in the file.
    """
    molecules: List[MoleculeInfo] = []

    try:
        content = sdf_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        print(f"  [ERROR] Cannot read SDF file: {sdf_path} ({e})")
        return []

    # Split into individual molecule records
    records = content.split("$$$$")

    for record in records:
        record = record.strip()
        if not record:
            continue

        mol = MoleculeInfo(source_file=str(sdf_path))
        lines = record.splitlines()

        # First line is the molecule title/name
        if lines:
            mol.title = lines[0].strip()

        # Parse properties (lines after "M  END")
        in_properties = False
        current_prop = None

        for line in lines:
            if line.strip() == "M  END":
                in_properties = True
                continue

            if not in_properties:
                continue

            # Property header: "> <PROPERTY_NAME>"
            prop_match = re.match(r'^>\s*<(.+?)>\s*$', line)
            if prop_match:
                current_prop = prop_match.group(1).strip()
                continue

            # Property value (line after the header)
            if current_prop and line.strip():
                mol.properties[current_prop] = line.strip()
                current_prop = None

        # Determine CID
        cid_from_prop = mol.properties.get(CID_PROPERTY, "").strip()
        if cid_from_prop:
            mol.cid = cid_from_prop
        elif mol.title.isdigit():
            # Title is a pure number — likely a CID
            mol.cid = mol.title
        else:
            # Check if source filename has a CID (e.g. Structure3D_CID_2244.sdf)
            file_stem = Path(mol.source_file).stem
            stem_cid = extract_compound_identifier(file_stem)
            if stem_cid.isdigit():
                mol.cid = stem_cid

        molecules.append(mol)

    return molecules


def parse_multiple_sdfs(sdf_paths: List[Path]) -> List[MoleculeInfo]:
    """Parse multiple SDF files and return combined molecule list."""
    all_molecules: List[MoleculeInfo] = []
    for sdf_path in sdf_paths:
        print(f"  Parsing: {sdf_path.name}")
        molecules = parse_sdf(sdf_path)
        print(f"    -> Found {len(molecules)} molecule(s)")
        all_molecules.extend(molecules)
    return all_molecules


# ===============================================================================
# PubChem API Name Lookup
# ===============================================================================

def fetch_pubchem_names(cids: List[str]) -> Dict[str, str]:
    """Fetch common compound names from PubChem REST API.

    Queries the PubChem synonyms endpoint in batches. Returns the first
    (most common) synonym for each CID, which is typically the widely
    recognized name (e.g., "Aspirin" for CID 2244).

    Args:
        cids: List of PubChem CID strings.

    Returns:
        Dict mapping CID string -> common name.
    """
    if not cids:
        return {}

    # Filter to valid numeric CIDs only
    valid_cids = [c for c in cids if c.isdigit()]
    if not valid_cids:
        print("  [WARN] No valid numeric CIDs found for PubChem lookup.")
        return {}

    names: Dict[str, str] = {}
    total_batches = (len(valid_cids) + PUBCHEM_BATCH_SIZE - 1) // PUBCHEM_BATCH_SIZE

    print(f"\n  Fetching names from PubChem API ({len(valid_cids)} compounds, "
          f"{total_batches} batch(es))...")

    for batch_idx in range(0, len(valid_cids), PUBCHEM_BATCH_SIZE):
        batch = valid_cids[batch_idx:batch_idx + PUBCHEM_BATCH_SIZE]
        batch_num = batch_idx // PUBCHEM_BATCH_SIZE + 1

        cid_str = ",".join(batch)
        url = f"{PUBCHEM_API_BASE}/{cid_str}/synonyms/JSON"

        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent", "DockingAutomationSuite/1.0")

            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode("utf-8"))

            info_list = data.get("InformationList", {}).get("Information", [])
            for info in info_list:
                cid = str(info.get("CID", ""))
                synonyms = info.get("Synonym", [])
                if cid and synonyms:
                    # First synonym is typically the most common name
                    names[cid] = synonyms[0]

            found = sum(1 for c in batch if c in names)
            print(f"    Batch {batch_num}/{total_batches}: "
                  f"{found}/{len(batch)} names resolved")

        except urllib.error.HTTPError as e:
            print(f"    [WARN] Batch {batch_num} HTTP error: {e.code} {e.reason}")
            # Try individual lookups for this batch as fallback
            for cid in batch:
                name = _fetch_single_pubchem_name(cid)
                if name:
                    names[cid] = name

        except urllib.error.URLError as e:
            print(f"    [ERROR] Network error: {e.reason}")
            print("    Cannot reach PubChem API. Check your internet connection.")
            print("    Falling back to SDF property names only.")
            return names

        except Exception as e:
            print(f"    [WARN] Batch {batch_num} error: {e}")

        # Respect PubChem rate limits (max 5 requests/second)
        if batch_idx + PUBCHEM_BATCH_SIZE < len(valid_cids):
            time.sleep(0.25)

    resolved = len(names)
    print(f"  PubChem lookup complete: {resolved}/{len(valid_cids)} names resolved")

    return names


def _fetch_single_pubchem_name(cid: str) -> Optional[str]:
    """Fetch a single compound name from PubChem (fallback for batch failures)."""
    url = f"{PUBCHEM_API_BASE}/{cid}/synonyms/JSON"
    try:
        req = urllib.request.Request(url)
        req.add_header("User-Agent", "DockingAutomationSuite/1.0")
        with urllib.request.urlopen(req, timeout=15) as response:
            data = json.loads(response.read().decode("utf-8"))
        info_list = data.get("InformationList", {}).get("Information", [])
        if info_list and info_list[0].get("Synonym"):
            return info_list[0]["Synonym"][0]
    except Exception:
        pass
    return None


# ===============================================================================
# CSV Mapping
# ===============================================================================

def load_csv_mapping(csv_path: Path) -> Dict[str, str]:
    """Load a user-provided CID/filename-to-name mapping CSV.

    Expected CSV format (with or without header):
        filename,compound_name
        2244,Aspirin
        5090,Naproxen
        3672,Ibuprofen

    Or:
        cid,name
        2244,Aspirin

    The first column is the PDBQT filename stem (or CID).
    The second column is the desired compound name.

    Args:
        csv_path: Path to the CSV mapping file.

    Returns:
        Dict mapping filename stem / CID -> compound name.
    """
    mapping: Dict[str, str] = {}

    try:
        with open(csv_path, "r", encoding="utf-8-sig") as f:
            # Sniff for header
            sample = f.read(2048)
            f.seek(0)

            has_header = csv.Sniffer().has_header(sample)
            reader = csv.reader(f)

            if has_header:
                next(reader)  # Skip header

            for row_num, row in enumerate(reader, start=2 if has_header else 1):
                if len(row) < 2:
                    continue
                key = row[0].strip()
                name = row[1].strip()
                if key and name:
                    mapping[key] = name

    except Exception as e:
        print(f"  [ERROR] Cannot read CSV mapping: {csv_path} ({e})")
        return {}

    print(f"  Loaded {len(mapping)} name mappings from {csv_path.name}")
    return mapping


# ===============================================================================
# Name Sanitization
# ===============================================================================

def sanitize_filename(name: str, max_length: int = MAX_FILENAME_LENGTH) -> str:
    """Make a compound name safe for use as a filename.

    - Replaces invalid characters with underscores
    - Replaces spaces with underscores
    - Collapses consecutive underscores
    - Truncates to max_length
    - Strips leading/trailing underscores and dots

    Args:
        name: Raw compound name.
        max_length: Maximum filename length (stem, not including extension).

    Returns:
        Sanitized filename-safe string.
    """
    # Replace invalid filename characters
    safe = INVALID_FILENAME_CHARS.sub("_", name)

    # Replace spaces, tabs, and multiple whitespace with single underscore
    safe = re.sub(r'\s+', '_', safe)

    # Replace parentheses, brackets, commas with underscore for cleanliness
    # (these are technically valid on Windows but can cause issues with some tools)
    safe = re.sub(r'[(),\[\]{}]+', '_', safe)

    # Collapse consecutive underscores
    safe = re.sub(r'_+', '_', safe)

    # Strip leading/trailing underscores and dots
    safe = safe.strip('_.')

    # Truncate if too long
    if len(safe) > max_length:
        safe = safe[:max_length].rstrip('_')

    return safe or "unnamed_compound"


# ===============================================================================
# Rename Planning
# ===============================================================================

def resolve_names(
    molecules: List[MoleculeInfo],
    csv_mapping: Dict[str, str],
    pubchem_names: Dict[str, str],
) -> Dict[str, Tuple[str, str]]:
    """Resolve the best name for each molecule identifier.

    Applies the name resolution priority:
        1. CSV mapping (user override)
        2. PubChem API name
        3. SDF property names (in priority order)
        4. Non-numeric title
        5. Keep original

    Args:
        molecules: Parsed molecule info from SDF files.
        csv_mapping: User-provided name overrides.
        pubchem_names: Names fetched from PubChem API.

    Returns:
        Dict mapping identifier -> (resolved_name, source_label).
    """
    resolved: Dict[str, Tuple[str, str]] = {}

    for mol in molecules:
        identifier = mol.identifier
        if not identifier:
            continue

        # Already resolved (from an earlier SDF file or duplicate entry)
        if identifier in resolved:
            continue

        # Priority 1: CSV mapping
        if identifier in csv_mapping:
            resolved[identifier] = (csv_mapping[identifier], "CSV mapping")
            continue

        # Priority 2: PubChem API
        cid = mol.cid
        if cid and cid in pubchem_names:
            resolved[identifier] = (pubchem_names[cid], "PubChem API")
            continue

        # Priority 3: SDF properties
        sdf_name = mol.best_sdf_name
        if sdf_name:
            resolved[identifier] = (sdf_name, "SDF property")
            continue

        # Priority 4: Title is not a pure number (it's already a name)
        if mol.title and not mol.title.isdigit():
            resolved[identifier] = (mol.title, "SDF title")
            continue

        # Priority 5: No name found - skip (keep original)
        resolved[identifier] = (identifier, "unchanged")

    return resolved


def build_rename_plan(
    target_dir: Path,
    name_map: Dict[str, Tuple[str, str]],
    file_ext: str = ".pdbqt",
    max_length: int = MAX_FILENAME_LENGTH,
) -> List[RenamePlan]:
    """Build a list of rename operations by matching target files to resolved names.

    Scans the target directory for files whose stems match known identifiers
    (CIDs, molecule titles, or PubChem filenames). Handles:
      - PyRx naming convention: "2353_uff_e=-535.63.pdbqt"
      - PubChem downloads: "Structure3D_CID_2244.sdf", "Conformer3D_CID_2244.sdf", "CID_2244.sdf"
      - Raw CID stems: "2244.sdf", "2244.pdbqt"

    Args:
        target_dir: Directory containing files to rename.
        name_map: Mapping of identifier -> (resolved_name, source_label).
        file_ext: Extension of files to rename (".pdbqt" or ".sdf").
        max_length: Maximum filename length for sanitization.

    Returns:
        List of RenamePlan objects.
    """
    plan: List[RenamePlan] = []

    # Ensure extension starts with dot
    ext_pattern = f"*{file_ext}" if file_ext.startswith(".") else f"*.{file_ext}"
    target_files = sorted(target_dir.glob(ext_pattern))
    if not target_files:
        print(f"  [WARN] No {file_ext.upper()} files found in {target_dir}")
        return []

    # Track new names to detect conflicts
    used_names: Dict[str, int] = {}

    for file_path in target_files:
        stem = file_path.stem

        # Extract identifier from filename
        # e.g. "Structure3D_CID_2244" -> "2244", "2353_uff_e=-535.63" -> "2353"
        identifier = extract_compound_identifier(stem)

        # Try to match to a known identifier in the name map
        if identifier not in name_map:
            # Also try the raw stem
            if stem in name_map:
                identifier = stem
            else:
                continue  # Not a file we know about - skip

        raw_name, source = name_map[identifier]

        # Skip if name is unchanged (no rename needed)
        if source == "unchanged":
            # Still rename if the file has a PyRx suffix or PubChem prefix to clean up
            has_pyrx_suffix = PYRX_UFF_SUFFIX.match(stem) is not None
            has_pubchem_prefix = PUBCHEM_STEM_PATTERN.match(stem) is not None and stem != identifier
            if not has_pyrx_suffix and not has_pubchem_prefix:
                continue
            raw_name = identifier
            source = "prefix/suffix cleaned"

        # Sanitize the name for filesystem safety
        safe_name = sanitize_filename(raw_name, max_length=max_length)

        # Skip if sanitized name is the same as current stem
        if safe_name == stem:
            continue

        # Handle name conflicts (two CIDs resolving to the same name)
        if safe_name in used_names:
            used_names[safe_name] += 1
            # Append CID to disambiguate
            safe_name = f"{safe_name}_CID{identifier}"
        else:
            used_names[safe_name] = 1

        new_path = file_path.parent / f"{safe_name}{file_ext}"

        entry = RenamePlan(
            original_path=file_path,
            new_path=new_path,
            original_name=file_path.name,
            new_name=new_path.name,
            name_source=source,
            cid=identifier,
        )
        plan.append(entry)

    return plan


# ===============================================================================
# Execution
# ===============================================================================

def execute_renames(plan: List[RenamePlan], dry_run: bool = False) -> None:
    """Execute the rename plan.

    Args:
        plan: List of RenamePlan operations.
        dry_run: If True, only display what would happen without renaming.
    """
    if not plan:
        print("\n  No renames to perform.")
        return

    mode_str = "DRY-RUN PREVIEW" if dry_run else "EXECUTING RENAMES"

    print()
    print("=" * 70)
    print(f"  {mode_str}")
    print("=" * 70)

    # Display table header
    print(f"\n  {'#':<4} {'Original':<30} {'-->':<3} {'New Name':<30} {'Source'}")
    print(f"  {'-'*4} {'-'*30} {'-'*3} {'-'*30} {'-'*15}")

    renamed = 0
    skipped = 0
    errors = 0

    for i, entry in enumerate(plan, 1):
        orig_display = _truncate(entry.original_name, 28)
        new_display = _truncate(entry.new_name, 28)

        if dry_run:
            print(f"  {i:<4} {orig_display:<30} --> {new_display:<30} {entry.name_source}")
            entry.status = "DRY_RUN"
            continue

        # Check if target already exists
        if entry.new_path.exists():
            print(f"  {i:<4} {orig_display:<30} --> {new_display:<30} [SKIP: target exists]")
            entry.status = "SKIPPED"
            entry.message = "Target file already exists"
            skipped += 1
            continue

        # Check if source still exists
        if not entry.original_path.exists():
            print(f"  {i:<4} {orig_display:<30} --> {new_display:<30} [SKIP: source gone]")
            entry.status = "SKIPPED"
            entry.message = "Source file no longer exists"
            skipped += 1
            continue

        # Perform the rename
        try:
            entry.original_path.rename(entry.new_path)
            print(f"  {i:<4} {orig_display:<30} --> {new_display:<30} {entry.name_source}")
            entry.status = "RENAMED"
            renamed += 1
        except Exception as e:
            print(f"  {i:<4} {orig_display:<30} --> {new_display:<30} [ERROR: {e}]")
            entry.status = "ERROR"
            entry.message = str(e)
            errors += 1

    # Summary
    print()
    print("=" * 70)
    if dry_run:
        print(f"  DRY-RUN: {len(plan)} file(s) would be renamed.")
        print(f"  Run again without --dry-run to execute.")
    else:
        print(f"  Renamed  : {renamed}")
        if skipped:
            print(f"  Skipped  : {skipped}")
        if errors:
            print(f"  Errors   : {errors}")
    print("=" * 70)


def _truncate(text: str, max_len: int) -> str:
    """Truncate text with ellipsis if too long."""
    if len(text) <= max_len:
        return text
    return text[:max_len - 2] + ".."


# ===============================================================================
# Rename Log (for undo support)
# ===============================================================================

def save_rename_log(plan: List[RenamePlan], log_path: Path, dir_flag: str = "--pdbqt-dir") -> None:
    """Save a CSV log of all rename operations for undo capability.

    The log records original and new filenames, CIDs, name sources,
    and status for each operation.

    Args:
        plan: List of executed (or dry-run) RenamePlan operations.
        log_path: Path to save the CSV log.
        dir_flag: Flag name to display in the undo hint (--sdf-dir or --pdbqt-dir).
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with open(log_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)

        # Header with metadata
        writer.writerow([f"# Docking Automation Suite - Ligand Rename Log"])
        writer.writerow([f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"])
        writer.writerow([])

        # Column headers
        writer.writerow([
            "original_filename", "new_filename", "cid",
            "name_source", "status", "message",
        ])

        for entry in plan:
            writer.writerow([
                entry.original_path.name,
                entry.new_path.name,
                entry.cid,
                entry.name_source,
                entry.status,
                entry.message,
            ])

    print(f"\n  Rename log saved: {log_path}")
    print(f"  To undo: python rename_ligands.py --undo {log_path} {dir_flag} <dir>")


# ===============================================================================
# Undo
# ===============================================================================

def undo_renames(log_path: Path, target_dir: Path, dry_run: bool = False) -> None:
    """Undo a previous rename operation using the rename log.

    Reads the rename log CSV and reverses all RENAMED operations,
    restoring the original filenames.

    Args:
        log_path: Path to the rename_log.csv file.
        target_dir: Directory containing the renamed files (PDBQT or SDF).
        dry_run: If True, only preview without undoing.
    """
    if not log_path.exists():
        print(f"  [ERROR] Rename log not found: {log_path}")
        return

    entries: List[Tuple[str, str]] = []

    with open(log_path, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        for row in reader:
            # Skip comment lines and empty rows
            if not row or row[0].startswith("#"):
                continue
            # Skip header
            if row[0] == "original_filename":
                continue
            if len(row) >= 5 and row[4] == "RENAMED":
                entries.append((row[0], row[1]))  # original, new

    if not entries:
        print("  No RENAMED entries found in log. Nothing to undo.")
        return

    mode_str = "DRY-RUN UNDO PREVIEW" if dry_run else "UNDOING RENAMES"
    print()
    print("=" * 70)
    print(f"  {mode_str}")
    print("=" * 70)

    restored = 0
    skipped = 0

    for orig_name, new_name in entries:
        current_path = target_dir / new_name
        restore_path = target_dir / orig_name

        if not current_path.exists():
            print(f"  [SKIP] {new_name} not found (already undone or moved)")
            skipped += 1
            continue

        if restore_path.exists():
            print(f"  [SKIP] {orig_name} already exists (would overwrite)")
            skipped += 1
            continue

        if dry_run:
            print(f"  {new_name} -> {orig_name}")
        else:
            try:
                current_path.rename(restore_path)
                print(f"  [OK] {new_name} -> {orig_name}")
                restored += 1
            except Exception as e:
                print(f"  [ERROR] {new_name}: {e}")

    print()
    print("=" * 70)
    if dry_run:
        print(f"  DRY-RUN: {len(entries) - skipped} file(s) would be restored.")
    else:
        print(f"  Restored: {restored}, Skipped: {skipped}")
    print("=" * 70)


# ===============================================================================
# Summary Display
# ===============================================================================

def display_molecule_summary(
    molecules: List[MoleculeInfo],
    name_map: Dict[str, Tuple[str, str]],
) -> None:
    """Display a summary of parsed molecules and their resolved names.

    Shows each CID/identifier, the resolved name, and where the name
    came from. Useful for reviewing before committing to renames.

    Args:
        molecules: Parsed molecule information.
        name_map: Resolved names mapping.
    """
    print()
    print("=" * 70)
    print("  COMPOUND NAME RESOLUTION")
    print("=" * 70)
    print(f"\n  {'CID/ID':<15} {'Resolved Name':<40} {'Source'}")
    print(f"  {'-'*15} {'-'*40} {'-'*15}")

    seen = set()
    for mol in molecules:
        identifier = mol.identifier
        if identifier in seen:
            continue
        seen.add(identifier)

        name, source = name_map.get(identifier, (identifier, "not found"))
        safe_name = sanitize_filename(name) if source != "unchanged" else identifier

        id_display = _truncate(identifier, 13)
        name_display = _truncate(safe_name, 38)
        print(f"  {id_display:<15} {name_display:<40} {source}")

    print(f"\n  Total compounds: {len(seen)}")
    unchanged = sum(1 for _, (_, s) in name_map.items() if s == "unchanged")
    if unchanged:
        print(f"  No name found : {unchanged} (will keep original filename)")
    print("=" * 70)


# ═══════════════════════════════════════════════════════════════════════════════
# CLI Entry Point
# ═══════════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="rename_ligands.py",
        description=BANNER,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # ── Rename SDF files directly (e.g. downloaded from PubChem) ────────
  python rename_ligands.py --sdf-dir ./my_sdfs/
  python rename_ligands.py --sdf-dir ./my_sdfs/ --fetch-names
  python rename_ligands.py --sdf-dir ./my_sdfs/ --dry-run

  # ── Rename PDBQT files from PyRx ────────────────────────────────────
  python rename_ligands.py --sdf compounds.sdf --pdbqt-dir ligands/
  python rename_ligands.py --sdf compounds.sdf --pdbqt-dir ligands/ --fetch-names
  python rename_ligands.py --sdf compounds.sdf --pdbqt-dir ligands/ --dry-run

  # ── Use manual CSV mapping ──────────────────────────────────────────
  python rename_ligands.py --csv mapping.csv --sdf-dir ./my_sdfs/
  python rename_ligands.py --csv mapping.csv --pdbqt-dir ligands/

  # ── Undo a previous rename ──────────────────────────────────────────
  python rename_ligands.py --undo rename_log.csv --sdf-dir ./my_sdfs/
  python rename_ligands.py --undo rename_log.csv --pdbqt-dir ligands/
""",
    )

    parser.add_argument(
        "--sdf", nargs="+", type=Path,
        help="Path(s) to SDF file(s) containing compound information",
    )
    parser.add_argument(
        "--pdbqt-dir", type=Path,
        help="Directory containing PDBQT files to rename",
    )
    parser.add_argument(
        "--sdf-dir", type=Path,
        help="Directory containing SDF files to rename directly in-place",
    )
    parser.add_argument(
        "--fetch-names", action="store_true",
        help="Fetch common compound names from PubChem REST API "
             "(recommended for PubChem SDF files)",
    )
    parser.add_argument(
        "--csv", type=Path,
        help="Path to a CSV mapping file (columns: filename/CID, compound_name)",
    )
    parser.add_argument(
        "--undo", type=Path,
        help="Undo a previous rename using the specified rename_log.csv",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Preview rename operations without executing them",
    )
    parser.add_argument(
        "--max-length", type=int, default=MAX_FILENAME_LENGTH,
        help=f"Maximum filename length (default: {MAX_FILENAME_LENGTH})",
    )
    parser.add_argument(
        "--log", type=Path,
        help="Path to save the rename log CSV (default: <target-dir>/rename_log.csv)",
    )

    args = parser.parse_args()

    # Banner
    print()
    print("=" * 70)
    print(f"  {BANNER}")
    print("=" * 70)

    # Validate target directory
    if args.sdf_dir and args.pdbqt_dir:
        print("  [ERROR] Specify either --sdf-dir OR --pdbqt-dir, not both at once.")
        sys.exit(1)

    if not args.sdf_dir and not args.pdbqt_dir:
        print("  [ERROR] Provide a target directory:")
        print("    • Use --sdf-dir <dir> to rename SDF files directly")
        print("    • Use --pdbqt-dir <dir> to rename PDBQT files")
        sys.exit(1)

    target_dir: Path = args.sdf_dir if args.sdf_dir else args.pdbqt_dir
    file_ext: str = ".sdf" if args.sdf_dir else ".pdbqt"

    if not target_dir.exists():
        dir_type = "SDF" if args.sdf_dir else "PDBQT"
        print(f"  [ERROR] {dir_type} directory not found: {target_dir}")
        sys.exit(1)

    # ── Undo mode ─────────────────────────────────────────────────────────
    if args.undo:
        undo_renames(args.undo, target_dir, dry_run=args.dry_run)
        sys.exit(0)

    # If --sdf-dir is given and --sdf is not passed, auto-discover all .sdf in that dir
    if args.sdf_dir and not args.sdf:
        sdf_files_in_dir = sorted(args.sdf_dir.glob("*.sdf"))
        if not sdf_files_in_dir:
            print(f"  [ERROR] No .sdf files found in {args.sdf_dir}")
            sys.exit(1)
        args.sdf = sdf_files_in_dir
        print(f"\n  Auto-detected {len(args.sdf)} SDF file(s) in {args.sdf_dir}")

    # ── Normal mode: need either --sdf or --csv ──────────────────────────
    if not args.sdf and not args.csv:
        print("  [ERROR] Provide at least --sdf, --sdf-dir, or --csv.")
        print("  Run with --help for usage examples.")
        sys.exit(1)

    # Validate SDF files exist
    if args.sdf:
        for sdf_path in args.sdf:
            if not sdf_path.exists():
                print(f"  [ERROR] SDF file not found: {sdf_path}")
                sys.exit(1)

    # ── Step 1: Parse SDF files ───────────────────────────────────────────
    molecules: List[MoleculeInfo] = []
    if args.sdf:
        print("\n  Parsing SDF file(s)...")
        molecules = parse_multiple_sdfs(args.sdf)
        if not molecules:
            print("  [WARN] No molecules found in SDF file(s).")
            if not args.csv:
                sys.exit(1)

    # ── Step 2: Load CSV mapping ──────────────────────────────────────────
    csv_mapping: Dict[str, str] = {}
    if args.csv:
        if not args.csv.exists():
            print(f"  [ERROR] CSV mapping file not found: {args.csv}")
            sys.exit(1)
        csv_mapping = load_csv_mapping(args.csv)

    # ── Step 3: Fetch PubChem names (if requested) ────────────────────────
    pubchem_names: Dict[str, str] = {}
    if args.fetch_names:
        cids = list(set(
            mol.cid for mol in molecules
            if mol.cid and mol.cid not in csv_mapping
        ))
        if cids:
            pubchem_names = fetch_pubchem_names(cids)
        else:
            print("\n  [INFO] No CIDs to look up (all covered by CSV or no CIDs found).")

    # ── Step 4: Resolve names ─────────────────────────────────────────────
    name_map = resolve_names(molecules, csv_mapping, pubchem_names)

    # If only CSV was provided (no SDF), add CSV entries to name_map
    if args.csv and not args.sdf:
        for key, name in csv_mapping.items():
            if key not in name_map:
                name_map[key] = (name, "CSV mapping")

    # Display resolution summary
    display_molecule_summary(molecules, name_map)

    # ── Step 5: Build rename plan ─────────────────────────────────────────
    plan = build_rename_plan(
        target_dir=target_dir,
        name_map=name_map,
        file_ext=file_ext,
        max_length=args.max_length,
    )

    if not plan:
        target_type = "SDF" if args.sdf_dir else "PDBQT"
        print(f"\n  No {target_type} files matched or all names are already correct.")
        print("  Check that:")
        print(f"    • {target_type} filenames match CIDs or PubChem naming (e.g. Structure3D_CID_XXXX.sdf)")
        print(f"    • The target directory is correct ({target_dir})")
        sys.exit(0)

    # ── Step 6: Execute ───────────────────────────────────────────────────
    execute_renames(plan, dry_run=args.dry_run)

    # ── Step 7: Save rename log ───────────────────────────────────────────
    log_path = args.log or target_dir / "rename_log.csv"
    dir_flag = "--sdf-dir" if args.sdf_dir else "--pdbqt-dir"
    save_rename_log(plan, log_path, dir_flag=dir_flag)


if __name__ == "__main__":
    main()
