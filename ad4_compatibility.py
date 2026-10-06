"""AutoDock4 atom-type compatibility and parameter provenance primitives.

This module describes the ordinary AutoDock4 type vocabulary evidenced by the
bundled ``AD4_parameters.dat`` and the inspected AutoDockTools/MGLTools data.
It records parameter provenance and compatibility; it is not a claim that a
specialized metal model is active.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
import hashlib
import sys
import json
import tomllib
import zipfile
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, List, Mapping, Optional, Tuple


@dataclass(frozen=True)
class ParameterProfile:
    """Description of the parameter set against which atom types are checked."""

    profile_id: str
    name: str
    engine: str
    supported_atom_types: FrozenSet[str]
    version: str = "1"
    parameter_file: Optional[str] = None
    source: str = ""
    source_path: Optional[str] = None
    reference: str = ""
    doi: str = ""
    license: str = ""
    sha256: Optional[str] = None
    validation_status: str = "inventory_only"
    specialized_metal_workflow: bool = False
    specialized_workflow: Optional[str] = None
    compatible_preparation_workflow: str = "standard_ad4"
    runtime_parameter_resolution: bool = False
    import_timestamp: Optional[str] = None
    user_notes: str = ""
    notes: str = ""

    @property
    def profile_name(self) -> str:
        """Stable public alias used by profile manifests."""
        return self.name

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data["supported_atom_types"] = sorted(self.supported_atom_types)
        resolved = None
        if self.parameter_file:
            path = Path(self.parameter_file)
            if not path.is_absolute() and self.source_path:
                path = Path(self.source_path) / path
            resolved = str(path.resolve())
        data["parameter_filename"] = Path(self.parameter_file).name if self.parameter_file else None
        data["parameter_path"] = resolved
        return data


@dataclass(frozen=True)
class MetalRecord:
    element: str
    ad4_atom_types: Tuple[str, ...]
    standard_ad4_parameterized: bool
    specialized_workflow: Optional[str] = None
    aliases: Tuple[str, ...] = ()
    pdb_names: Tuple[str, ...] = ()
    parameter_profiles: Tuple[str, ...] = ("ad4_standard_4.2",)
    coordination_support: str = "not_specialized"
    validation_state: str = "validated_standard_type"
    notes: str = ""


class MetalRegistry:
    """Registry separates standard AD4 type availability from metal workflows."""

    def __init__(self, records: Iterable[MetalRecord]):
        self._records = {record.element.upper(): record for record in records}
        self._type_to_element = {
            atom_type.upper(): record.element
            for record in self._records.values()
            for atom_type in record.ad4_atom_types
        }

    def get(self, element_or_type: str) -> Optional[MetalRecord]:
        key = element_or_type.strip().upper()
        element = self._type_to_element.get(key, key)
        return self._records.get(element.upper())

    def element_for_type(self, atom_type: str) -> Optional[str]:
        return self._type_to_element.get(atom_type.strip().upper())

    def to_dict(self) -> Dict[str, Dict[str, Any]]:
        return {element: asdict(record) for element, record in sorted(self._records.items())}


# Atom types copied from the installed MGLTools 1.5.7 AD4_parameters.dat type
# declarations (including its upper/lower-case aliases and pseudo types).
STANDARD_AD4_ATOM_TYPES: FrozenSet[str] = frozenset({
    "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F",
    "Mg", "MG", "P", "SA", "S", "Cl", "CL", "Ca", "CA", "Mn", "MN",
    "Fe", "FE", "Zn", "ZN", "Br", "BR", "I", "Z", "G", "GA", "J", "Q",
})

_PROFILE_ROOT = Path(__file__).resolve().parent


def _resolve_bundled_file(filename: str, profile_subfolder: Optional[str] = None) -> Optional[Path]:
    candidates = []
    if profile_subfolder:
        candidates.append(_PROFILE_ROOT / "parameter_profiles" / profile_subfolder / filename)
        candidates.append(_PROFILE_ROOT / profile_subfolder / filename)
    candidates.append(_PROFILE_ROOT / "parameter_profiles" / filename)
    if filename.endswith(".zip"):
        candidates.append(_PROFILE_ROOT / "external" / "autodock4zn" / filename)
    candidates.append(_PROFILE_ROOT / filename)
    if hasattr(sys, "_MEIPASS"):
        meipass = Path(sys._MEIPASS)
        if profile_subfolder:
            candidates.append(meipass / "parameter_profiles" / profile_subfolder / filename)
            candidates.append(meipass / profile_subfolder / filename)
        candidates.append(meipass / "parameter_profiles" / filename)
        if filename.endswith(".zip"):
            candidates.append(meipass / "external" / "autodock4zn" / filename)
        candidates.append(meipass / filename)
    for c in candidates:
        if c.is_file():
            return c
    return None


def _bundled_checksum(filename: str, profile_subfolder: Optional[str] = None) -> Optional[str]:
    path = _resolve_bundled_file(filename, profile_subfolder)
    if not path or not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

_std_p = _resolve_bundled_file("AD4_parameters.dat", "ad4_standard_4.2")
_bound_p = _resolve_bundled_file("AD4.1_bound.dat", "ad4_1_bound")
_zn_p = _resolve_bundled_file("AutoDock4Zn-Pipeline-main.zip", "external/autodock4zn")

STANDARD_AD4_PROFILE = ParameterProfile(
    profile_id="ad4_standard_4.2",
    name="AutoDock4 standard parameters",
    engine="AUTODOCK4",
    supported_atom_types=STANDARD_AD4_ATOM_TYPES,
    version="4.2",
    parameter_file="AD4_parameters.dat",
    source_path=str(_std_p.parent if _std_p else _PROFILE_ROOT),
    source="Bundled AD4_parameters.dat; AutoDockTools/MGLTools-compatible AD4 parameter asset",
    reference="AutoDock4 force-field parameterization",
    license="AutoDock GPL header and AutoDock copyright notice",
    sha256=_bundled_checksum("AD4_parameters.dat", "ad4_standard_4.2"),
    validation_status="valid",
    specialized_metal_workflow=False,
    compatible_preparation_workflow="standard_ad4",
    runtime_parameter_resolution=False,
    notes="Standard non-specialized AutoDock4 parameters. Type recognition does not imply specialized metal coordination treatment.",
)

AD4_1_BOUND_PROFILE = ParameterProfile(
    profile_id="ad4_1_bound",
    name="AutoDock4.1 bound-state parameters",
    engine="AUTODOCK4",
    supported_atom_types=STANDARD_AD4_ATOM_TYPES,
    version="4.1-bound",
    parameter_file="AD4.1_bound.dat",
    source_path=str(_bound_p.parent if _bound_p else _PROFILE_ROOT),
    source="Bundled AD4.1_bound.dat; source header identifies Version 4.1 Bound",
    reference="Huey, Morris, Olson & Goodsell (2007), J Comput Chem 28:1145-1152",
    license="GNU General Public License v2 or later; AutoDock copyright notice",
    sha256=_bundled_checksum("AD4.1_bound.dat", "ad4_1_bound"),
    validation_status="valid",
    compatible_preparation_workflow="ad4_1_bound",
    notes="Distinct bound-state free-energy coefficients; not silently merged with standard AD4.",
)

AD4ZN_PROFILE = ParameterProfile(
    profile_id="ad4zn_specialized",
    name="AutoDock4Zn specialized parameters",
    engine="AUTODOCK4",
    supported_atom_types=STANDARD_AD4_ATOM_TYPES | frozenset({"TZ"}),
    version="reference",
    parameter_file="AutoDock4Zn-Pipeline-main.zip",
    source_path=str(_zn_p.parent if _zn_p else _PROFILE_ROOT),
    source="AutoDock4Zn-Pipeline-main.zip; AD4Zn.dat is an archive member",
    reference="AutoDock4Zn reference pipeline README and AD4Zn.dat",
    license="GPLv3; see archive LICENSE",
    sha256=_bundled_checksum("AutoDock4Zn-Pipeline-main.zip", "external/autodock4zn"),
    validation_status="reference_only",
    specialized_metal_workflow=True,
    specialized_workflow="AutoDock4Zn",
    compatible_preparation_workflow="ad4zn",
    notes="Requires isolated legacy preparation and dependency preflight; never selected by ordinary AD4 runs.",
)

METAL_REGISTRY = MetalRegistry([
    MetalRecord("Zn", ("Zn", "ZN"), True, specialized_workflow="AutoDock4Zn", aliases=("zinc",), pdb_names=("ZN",), parameter_profiles=("ad4_standard_4.2", "ad4_1_bound", "ad4zn_specialized"), notes="Ordinary non-bonded AD4 parameters; specialized reference workflow available, not selected."),
    MetalRecord("Fe", ("Fe", "FE"), True, aliases=("iron",), pdb_names=("FE",), parameter_profiles=("ad4_standard_4.2", "ad4_1_bound"), notes="Ordinary non-bonded AD4 parameters; no specialized heme chemistry."),
    MetalRecord("Mg", ("Mg", "MG"), True, aliases=("magnesium",), pdb_names=("MG",), parameter_profiles=("ad4_standard_4.2", "ad4_1_bound")),
    MetalRecord("Ca", ("Ca", "CA"), True, aliases=("calcium",), pdb_names=("CA",), parameter_profiles=("ad4_standard_4.2", "ad4_1_bound")),
    MetalRecord("Mn", ("Mn", "MN"), True, aliases=("manganese",), pdb_names=("MN",), parameter_profiles=("ad4_standard_4.2", "ad4_1_bound")),
    MetalRecord("Cu", (), False, aliases=("copper",), pdb_names=("CU",), parameter_profiles=(), validation_state="unsupported", notes="No standard AD4 parameter record in inspected assets."),
    MetalRecord("Co", (), False, aliases=("cobalt",), pdb_names=("CO",), parameter_profiles=(), validation_state="unsupported", notes="No standard AD4 parameter record in inspected assets."),
    MetalRecord("Ni", (), False, aliases=("nickel",), pdb_names=("NI",), parameter_profiles=(), validation_state="unsupported", notes="No standard AD4 parameter record in inspected assets."),
    MetalRecord("Mo", (), False, aliases=("molybdenum",), pdb_names=("MO",), parameter_profiles=(), validation_state="unsupported", notes="No standard AD4 parameter record in inspected assets."),
    MetalRecord("Hg", (), False, aliases=("mercury",), pdb_names=("HG",), parameter_profiles=(), validation_state="unsupported", notes="No standard AD4 parameter record in inspected assets."),
])

PARAMETER_PROFILES: Mapping[str, ParameterProfile] = {
    STANDARD_AD4_PROFILE.profile_id: STANDARD_AD4_PROFILE,
    AD4_1_BOUND_PROFILE.profile_id: AD4_1_BOUND_PROFILE,
    AD4ZN_PROFILE.profile_id: AD4ZN_PROFILE,
}


@dataclass(frozen=True)
class ProfileValidationResult:
    """Structured profile diagnostics suitable for configuration and metadata."""

    profile_id: str
    status: str
    errors: Tuple[str, ...] = ()
    warnings: Tuple[str, ...] = ()
    checked_parameter_file: Optional[str] = None
    observed_sha256: Optional[str] = None

    @property
    def valid(self) -> bool:
        return self.status == "VALID"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "profile_id": self.profile_id,
            "status": self.status,
            "valid": self.valid,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
            "checked_parameter_file": self.checked_parameter_file,
            "observed_sha256": self.observed_sha256,
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_parameter_atom_types(parameter_file: Path | str) -> List[str]:
    """Read ``atom_par`` declarations from a parameter asset."""
    path = Path(parameter_file)
    if not path.is_file():
        return []
    types: List[str] = []
    seen = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0] == "atom_par" and fields[1] not in seen:
            seen.add(fields[1])
            types.append(fields[1])
    return types


def _parameter_text_types(path: Path) -> List[str]:
    if path.suffix.lower() == ".zip" and zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as package:
            members = [name for name in package.namelist() if name.endswith("/AD4Zn.dat")]
            if members:
                return [line.split()[1] for line in package.read(members[0]).decode("utf-8", "replace").splitlines()
                        if len(line.split()) >= 2 and line.split()[0] == "atom_par"]
    return parse_parameter_atom_types(path)


def validate_profile(profile: ParameterProfile,
                     parameter_root: Path | str | None = None,
                     preparation_workflow: str = "standard_ad4") -> ProfileValidationResult:
    """Validate a profile without changing it or any bundled profile data."""
    errors: List[str] = []
    warnings: List[str] = []
    engine = profile.engine.upper()
    if engine not in {"AUTODOCK4", "VINA"}:
        errors.append(f"Unsupported engine '{profile.engine}'.")
    if not profile.profile_id.strip() or not profile.version.strip():
        errors.append("Profile ID and version are required.")
    if not profile.supported_atom_types:
        errors.append("Profile contains no atom types.")
    if len(profile.supported_atom_types) != len({t.strip() for t in profile.supported_atom_types}):
        errors.append("Profile atom types contain duplicates or blank values.")
    if profile.compatible_preparation_workflow != preparation_workflow:
        errors.append(
            f"Preparation workflow '{preparation_workflow}' is incompatible with "
            f"'{profile.compatible_preparation_workflow}'."
        )
    if profile.specialized_metal_workflow and not profile.specialized_workflow:
        errors.append("Specialized profiles must declare specialized_workflow.")

    parameter_path: Optional[Path] = None
    if profile.parameter_file:
        parameter_path = Path(profile.parameter_file)
        if not parameter_path.is_absolute():
            root = parameter_root or profile.source_path
            if root:
                parameter_path = Path(root) / parameter_path
        if not parameter_path.is_file():
            errors.append(f"Parameter file not found: {parameter_path}")
        else:
            observed = _sha256(parameter_path)
            if profile.sha256 and observed.lower() != profile.sha256.lower():
                errors.append(
                    f"Parameter file checksum mismatch: expected {profile.sha256}, observed {observed}."
                )
            declared_types = set(_parameter_text_types(parameter_path))
            if not declared_types:
                errors.append("Parameter file contains no atom_par declarations.")
            elif not set(profile.supported_atom_types).issubset(declared_types):
                missing = sorted(set(profile.supported_atom_types) - declared_types)
                errors.append("Profile atom types missing from parameter file: " + ", ".join(missing))
    elif not profile.runtime_parameter_resolution:
        errors.append("Profile requires parameter information but declares no parameter file or runtime resolution.")
    else:
        warnings.append("Parameter resolution is delegated to the docking runtime; ADSP does not bundle the file.")

    return ProfileValidationResult(
        profile_id=profile.profile_id,
        status="VALID" if not errors else "INVALID",
        errors=tuple(errors),
        warnings=tuple(warnings),
        checked_parameter_file=str(parameter_path) if parameter_path else None,
        observed_sha256=_sha256(parameter_path) if parameter_path and parameter_path.is_file() else None,
    )


def validate_profile_id(profile_id: str, **kwargs: Any) -> ProfileValidationResult:
    profile = PARAMETER_PROFILES.get(profile_id)
    if profile is None:
        return ProfileValidationResult(profile_id=profile_id, status="INVALID",
                                       errors=(f"Unknown parameter profile '{profile_id}'.",))
    return validate_profile(profile, **kwargs)


def register_user_profile(profile_directory: Path | str) -> ParameterProfile:
    """Read and register a user profile manifest without touching standard profiles."""
    directory = Path(profile_directory)
    manifest_path = directory / "profile.toml"
    if not manifest_path.is_file():
        raise ValueError(f"Profile manifest not found: {manifest_path}")
    data = tomllib.loads(manifest_path.read_text(encoding="utf-8"))
    raw = data.get("profile", data)
    atom_types = frozenset(str(value).strip() for value in raw.get("atom_types", []) if str(value).strip())
    parameter_file = raw.get("parameter_file")
    parameter_path = directory / str(parameter_file) if parameter_file else None
    checksum = str(raw.get("sha256", "")).strip() or None
    if parameter_path and parameter_path.is_file() and not checksum:
        checksum = _sha256(parameter_path)
    profile = ParameterProfile(
        profile_id=str(raw.get("profile_id", "")).strip(),
        name=str(raw.get("profile_name", raw.get("name", ""))).strip(),
        engine=str(raw.get("engine", "AUTODOCK4")).upper(),
        supported_atom_types=atom_types,
        version=str(raw.get("version", "1")),
        parameter_file=str(parameter_path) if parameter_path else None,
        source=str(raw.get("source", str(directory))),
        source_path=str(directory),
        reference=str(raw.get("reference", "")),
        doi=str(raw.get("doi", "")),
        license=str(raw.get("license", "")),
        sha256=checksum,
        validation_status="imported",
        specialized_metal_workflow=bool(raw.get("specialized", False)),
        specialized_workflow=raw.get("specialized_workflow"),
        compatible_preparation_workflow=str(raw.get("compatible_preparation_workflow", "standard_ad4")),
        runtime_parameter_resolution=bool(raw.get("runtime_parameter_resolution", False)),
        import_timestamp=datetime.now(timezone.utc).isoformat(),
        user_notes=str(raw.get("user_notes", "")),
        notes=str(raw.get("notes", "")),
    )
    if not profile.profile_id:
        raise ValueError("User profile requires profile_id.")
    result = validate_profile(profile, parameter_root=directory)
    profile = replace(profile, validation_status=result.status.lower())
    existing = PARAMETER_PROFILES.get(profile.profile_id)
    if existing is not None and existing != profile:
        raise ValueError(
            f"Profile ID '{profile.profile_id}' is already registered with different metadata; "
            "use a new versioned profile ID."
        )
    # This assignment only extends the in-memory registry; bundled standard data is untouched.
    PARAMETER_PROFILES[profile.profile_id] = profile  # type: ignore[index]
    return profile


@dataclass
class AtomTypeCheck:
    role: str
    path: str
    atom_types: List[str] = field(default_factory=list)
    unsupported_atom_types: List[str] = field(default_factory=list)
    metal_atoms: List[Dict[str, str]] = field(default_factory=list)
    compatibility_states: Dict[str, str] = field(default_factory=dict)
    sha256: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class AtomTypeCompatibilityError(ValueError):
    """Raised before AD4 file generation/execution for unsupported atom types."""

    def __init__(self, errors: List[str], provenance: Dict[str, Any]):
        self.errors = errors
        self.provenance = provenance
        super().__init__("AutoDock4 compatibility preflight failed: " + "; ".join(errors))


def extract_atom_types(pdbqt_path: Path | str) -> List[str]:
    """Return unique PDBQT AutoDock types, preserving file order."""
    path = Path(pdbqt_path)
    types: List[str] = []
    seen = set()
    if not path.is_file():
        return types
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 79:
            atom_type = line[77:79].strip()
            if atom_type and atom_type not in seen:
                seen.add(atom_type)
                types.append(atom_type)
    return types


def _metal_atoms(path: Path, known_types: Iterable[str]) -> List[Dict[str, str]]:
    known = set(known_types)
    atoms: List[Dict[str, str]] = []
    if not path.is_file():
        return atoms
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_type = line[77:79].strip() if len(line) >= 79 else ""
        element = line[76:78].strip() if len(line) >= 78 else ""
        record = METAL_REGISTRY.get(atom_type) or METAL_REGISTRY.get(element)
        if record:
            atoms.append({
                "element": record.element,
                "atom_type": atom_type,
                "residue_name": line[17:20].strip(),
                "atom_name": line[12:16].strip(),
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "standard_ad4_parameterized": str(record.standard_ad4_parameterized).lower(),
                "standard_ad4": "available" if record.standard_ad4_parameterized else "unavailable",
                "specialized_workflows": [record.specialized_workflow] if record.specialized_workflow else [],
                "specialized_workflow": record.specialized_workflow or "none",
                "specialized_selected": False,
                "validation_status": "standard_only" if record.standard_ad4_parameterized else "unsupported",
            })
    return atoms


def inspect_pdbqt(role: str, path: Path | str,
                   profile: ParameterProfile = STANDARD_AD4_PROFILE) -> AtomTypeCheck:
    path = Path(path)
    atom_types = extract_atom_types(path)
    unsupported = sorted({atom_type for atom_type in atom_types
                          if atom_type not in profile.supported_atom_types})
    supported_state = (
        "SUPPORTED STANDARD" if profile.profile_id == STANDARD_AD4_PROFILE.profile_id
        else ("SUPPORTED SPECIALIZED" if profile.specialized_metal_workflow
              else "CUSTOM / USER-PROVIDED")
    )
    states = {
        atom_type: (supported_state if atom_type in profile.supported_atom_types else "UNSUPPORTED")
        for atom_type in atom_types
    }
    return AtomTypeCheck(
        role=role,
        path=str(path),
        atom_types=atom_types,
        unsupported_atom_types=unsupported,
        metal_atoms=_metal_atoms(path, atom_types),
        compatibility_states=states,
        sha256=_sha256(path) if path.is_file() else None,
    )


def build_job_provenance(files: Iterable[Tuple[str, Path | str]],
                         profile: Optional[ParameterProfile] = STANDARD_AD4_PROFILE,
                         engine: str = "AUTODOCK4") -> Dict[str, Any]:
    checks = [inspect_pdbqt(role, path, profile or STANDARD_AD4_PROFILE) for role, path in files]
    if profile is None:
        for check in checks:
            check.unsupported_atom_types = []
    metals = [metal for check in checks for metal in check.metal_atoms]
    return {
        "engine": engine,
        "atom_typing": {
            "source": "PDBQT AutoDock atom-type column (originating preparer is not encoded)",
            "files": [check.to_dict() for check in checks],
        },
        "parameter_profile": profile.to_dict() if profile else None,
        "parameter_provenance": (
            "Standard AD4 profile type inventory used for compatibility preflight; "
            "the AutoDock4 executable retains its normal runtime parameter lookup."
            if profile else "Scoring and atom-type interpretation are managed by the Vina executable."
        ),
        "metals": {
            "detected": metals,
            "ordinary_atom_type_support": True if engine.upper() == "AUTODOCK4" else None,
            "specialized_metal_workflow_available": any(
                bool(METAL_REGISTRY.get(metal["element"]).specialized_workflow)
                for metal in metals if METAL_REGISTRY.get(metal["element"])
            ),
            "specialized_metal_workflow_selected": False,
            "interpretation": "Atom retention/type recognition is not specialized metal coordination treatment.",
        },
    }


def preflight_job(files: Iterable[Tuple[str, Path | str]],
                  profile_id: str = STANDARD_AD4_PROFILE.profile_id,
                  engine: str = "AUTODOCK4",
                  preparation_workflow: str = "standard_ad4",
                  required_map_types: Optional[Iterable[str]] = None) -> Dict[str, Any]:
    profile = PARAMETER_PROFILES.get(profile_id)
    if profile is None:
        raise AtomTypeCompatibilityError(
            [f"Unknown AutoDock4 parameter profile '{profile_id}'."],
            {"engine": engine, "parameter_profile_id": profile_id},
        )
    validation = validate_profile(profile, preparation_workflow=preparation_workflow)
    if not validation.valid:
        raise AtomTypeCompatibilityError(list(validation.errors), {
            "engine": engine, "parameter_profile": profile.to_dict(),
            "profile_validation": validation.to_dict(),
        })
    file_list = [(role, Path(path)) for role, path in files if path]
    provenance = build_job_provenance(file_list, profile, engine)
    provenance["profile_validation"] = validation.to_dict()
    provenance["compatibility_resolution"] = {
        "states": ["SUPPORTED STANDARD", "SUPPORTED SPECIALIZED", "CUSTOM / USER-PROVIDED",
                    "UNSUPPORTED", "UNKNOWN / NOT VALIDATED"],
        "required_map_types": list(required_map_types or []),
        "resolved_map_types": sorted({
            atom_type for check in provenance["atom_typing"]["files"]
            for atom_type in check["atom_types"]
        }),
    }
    errors: List[str] = []
    for check in provenance["atom_typing"]["files"]:
        if not check["atom_types"]:
            check["compatibility_states"] = {"<none>": "UNKNOWN / NOT VALIDATED"}
            errors.append(f"{check['role']} PDBQT has no readable AutoDock atom types: {check['path']}")
        if check["unsupported_atom_types"]:
            errors.append(
                f"{check['role']} PDBQT contains unsupported atom types for profile "
                f"'{profile.profile_id}': {', '.join(check['unsupported_atom_types'])}"
            )
    required = set(required_map_types or [])
    resolved = set(provenance["compatibility_resolution"]["resolved_map_types"])
    missing_maps = sorted(required - resolved)
    if missing_maps:
        errors.append("Required map atom types were not found: " + ", ".join(missing_maps))
    if errors:
        raise AtomTypeCompatibilityError(errors, provenance)
    return provenance
