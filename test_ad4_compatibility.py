"""Regression tests for AutoDock4 type compatibility and provenance."""

from __future__ import annotations

import tempfile
import unittest
import logging
import hashlib
from pathlib import Path
from unittest.mock import patch

from ad4_compatibility import (
    METAL_REGISTRY,
    STANDARD_AD4_PROFILE,
    AtomTypeCompatibilityError,
    extract_atom_types,
    preflight_job,
    PARAMETER_PROFILES,
    AD4_1_BOUND_PROFILE,
    ParameterProfile,
    register_user_profile,
    validate_profile,
    validate_profile_id,
    parse_parameter_atom_types,
    AD4ZN_PROFILE,
)
from autodock4_workflow import _extract_atom_types, generate_dpf, generate_gpf, run_ad4_job
from config import ProjectConfig
from models import DockingJob, Engine, ExecutionStatus, JobStatus
from prepare import parse_hetatm_records


def pdbqt_atom(serial: int, name: str, resname: str, atom_type: str,
               x: float = 0.0, y: float = 0.0, z: float = 0.0,
               charge: float = 0.0) -> str:
    return (
        f"ATOM  {serial:5d} {name:>4s} {resname:>3s} A{1:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}"
        f"{charge:10.3f} {atom_type:>2s}"
    )


class TestAD4Compatibility(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ad4_compatibility_")
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def _pdbqt(self, name: str, atoms: list[str]) -> Path:
        path = self.root / name
        path.write_text("\n".join(atoms + ["END"]) + "\n", encoding="utf-8")
        return path

    def test_type_extraction_is_shared_with_gpf_dpf(self):
        path = self._pdbqt("ligand.pdbqt", [
            pdbqt_atom(1, "C1", "LIG", "C"),
            pdbqt_atom(2, "O1", "LIG", "OA"),
            pdbqt_atom(3, "C2", "LIG", "C"),
        ])
        self.assertEqual(extract_atom_types(path), ["C", "OA"])
        self.assertEqual(_extract_atom_types(path), ["C", "OA"])

    def test_standard_metal_types_are_not_specialized_workflows(self):
        for metal in ("Zn", "Fe", "Mg", "Ca", "Mn"):
            record = METAL_REGISTRY.get(metal)
            self.assertIsNotNone(record)
            self.assertTrue(record.standard_ad4_parameterized)
            if metal == "Zn":
                self.assertEqual(record.specialized_workflow, "AutoDock4Zn")
            else:
                self.assertIsNone(record.specialized_workflow)
        zn = self._pdbqt("zn_receptor.pdbqt", [pdbqt_atom(1, "ZN", "ZN", "ZN")])
        provenance = preflight_job([("receptor", zn), ("ligand", self._pdbqt(
            "ligand2.pdbqt", [pdbqt_atom(1, "C1", "LIG", "C")]))])
        self.assertEqual(provenance["metals"]["detected"][0]["element"], "Zn")
        self.assertTrue(provenance["metals"]["specialized_metal_workflow_available"])
        self.assertFalse(provenance["metals"]["specialized_metal_workflow_selected"])

    def test_standard_profile_validates_with_runtime_provenance(self):
        result = validate_profile_id("ad4_standard_4.2")
        self.assertTrue(result.valid)
        self.assertEqual(result.observed_sha256, STANDARD_AD4_PROFILE.sha256)
        self.assertFalse(STANDARD_AD4_PROFILE.runtime_parameter_resolution)
        profile_data = STANDARD_AD4_PROFILE.to_dict()
        for key in ("profile_id", "version", "source", "source_path", "sha256", "specialized_workflow"):
            self.assertIn(key, profile_data)

    def test_bound_profile_is_distinct_and_valid(self):
        self.assertIn("ad4_1_bound", PARAMETER_PROFILES)
        self.assertNotEqual(AD4_1_BOUND_PROFILE.sha256, STANDARD_AD4_PROFILE.sha256)
        result = validate_profile(AD4_1_BOUND_PROFILE, preparation_workflow="ad4_1_bound")
        self.assertTrue(result.valid)

    def test_bundled_parameter_assets_are_parsed(self):
        root = Path(__file__).resolve().parent
        standard_types = set(parse_parameter_atom_types(root / "AD4_parameters.dat"))
        bound_types = set(parse_parameter_atom_types(root / "AD4.1_bound.dat"))
        self.assertIn("Zn", standard_types)
        self.assertIn("Fe", standard_types)
        self.assertEqual(standard_types, bound_types)
        self.assertTrue(AD4ZN_PROFILE.sha256)

    def test_invalid_profile_and_missing_parameter_file(self):
        self.assertFalse(validate_profile_id("does_not_exist").valid)
        profile = ParameterProfile(
            profile_id="missing", name="Missing", engine="AUTODOCK4",
            supported_atom_types=frozenset({"C"}), version="1",
            parameter_file="missing.dat", compatible_preparation_workflow="standard_ad4",
        )
        result = validate_profile(profile, parameter_root=self.root)
        self.assertFalse(result.valid)
        self.assertIn("not found", " ".join(result.errors))

    def test_checksum_validation(self):
        parameter_file = self.root / "parameters.dat"
        parameter_file.write_bytes(b"atom_par C 4.0 0.1 1.0 0.0 0 0 0 -1 -1 0\n")
        checksum = hashlib.sha256(parameter_file.read_bytes()).hexdigest()
        base = dict(profile_id="custom", name="Custom", engine="AUTODOCK4",
                    supported_atom_types=frozenset({"C"}), version="1",
                    parameter_file=str(parameter_file), compatible_preparation_workflow="standard_ad4")
        self.assertTrue(validate_profile(ParameterProfile(**base, sha256=checksum)).valid)
        self.assertFalse(validate_profile(ParameterProfile(**base, sha256="0" * 64)).valid)

    def test_user_profile_registration_is_isolated_from_standard(self):
        directory = self.root / "parameter_profiles" / "user" / "custom_cu_profile"
        directory.mkdir(parents=True)
        (directory / "parameters.dat").write_bytes(b"custom\n")
        (directory / "profile.toml").write_text(
            "[profile]\nprofile_id='custom_cu_profile'\nprofile_name='Custom Cu'\n"
            "engine='AUTODOCK4'\nversion='1'\natom_types=['C','Cu']\n"
            "parameter_file='parameters.dat'\ncompatible_preparation_workflow='standard_ad4'\n",
            encoding="utf-8")
        profile = register_user_profile(directory)
        self.assertEqual(profile.profile_id, "custom_cu_profile")
        self.assertIsNotNone(profile.import_timestamp)
        self.assertEqual(profile.source_path, str(directory))
        self.assertIn("custom_cu_profile", PARAMETER_PROFILES)
        self.assertEqual(PARAMETER_PROFILES["ad4_standard_4.2"], STANDARD_AD4_PROFILE)

    def test_receptor_ligand_and_flexible_types_are_resolved(self):
        receptor = self._pdbqt("rec.pdbqt", [pdbqt_atom(1, "ZN", "ZN", "ZN")])
        ligand = self._pdbqt("lig.pdbqt", [pdbqt_atom(1, "C1", "LIG", "C")])
        flexible = self._pdbqt("flex.pdbqt", [pdbqt_atom(1, "O1", "SER", "OA")])
        provenance = preflight_job(
            [("receptor", receptor), ("ligand", ligand), ("flex_receptor", flexible)],
            required_map_types=("ZN", "C", "OA"),
        )
        self.assertEqual(provenance["compatibility_resolution"]["required_map_types"], ["ZN", "C", "OA"])
        self.assertEqual(len(provenance["atom_typing"]["files"]), 3)

    def test_real_receptor_component_classification(self):
        root = Path(__file__).resolve().parent
        ca2 = parse_hetatm_records(root / "1CA2" / "1CA2_edited.pdb")
        mbn = parse_hetatm_records(root / "1MBN" / "1MBN_edited.pdb")
        zinc = next(item for item in ca2 if item["resname"] == "ZN")
        heme = next(item for item in mbn if item["resname"] == "HEM")
        self.assertEqual((zinc["chain"], zinc["resseq"], zinc["component_type"]), ("A", 262, "METAL"))
        self.assertEqual(heme["component_type"], "COFACTOR")
        self.assertIn("FE", heme["metal_elements"])

    def test_cif_preparation_preserves_metal_atom_types(self):
        from prepare import prepare_receptor_pdbqt
        root = Path(__file__).resolve().parent
        for filename, expected in (("1CA2.cif", "Zn"), ("1MBN.cif", "Fe")):
            output = self.root / (filename + ".pdbqt")
            prepare_receptor_pdbqt(root / filename, output, cleanup_water=False, add_hydrogens=False)
            atom_lines = [line for line in output.read_text(errors="replace").splitlines()
                          if line.startswith(("ATOM", "HETATM"))]
            self.assertIn(expected, {line[77:79].strip() for line in atom_lines})

    def test_cif_and_pdb_preparation_preserve_identity_and_coordinates(self):
        from prepare import prepare_receptor_pdbqt
        root = Path(__file__).resolve().parent

        def records(path):
            return {
                (line[12:16].strip(), line[17:20].strip(), line[21], line[22:26].strip(),
                 tuple(round(float(line[i:i + 8]), 3) for i in (30, 38, 46)))
                for line in path.read_text(errors="replace").splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }

        for stem in ("1CA2", "1MBN"):
            pdb_out = self.root / f"{stem}_pdb.pdbqt"
            cif_out = self.root / f"{stem}_cif.pdbqt"
            prepare_receptor_pdbqt(root / f"{stem}.pdb", pdb_out, cleanup_water=True, add_hydrogens=False)
            prepare_receptor_pdbqt(root / f"{stem}.cif", cif_out, cleanup_water=True, add_hydrogens=False)
            self.assertEqual(records(pdb_out), records(cif_out))
            pdb_out.unlink(missing_ok=True)
            cif_out.unlink(missing_ok=True)

    def test_2nv6_cif_and_pdb_preparation_preserve_zid(self):
        from prepare import prepare_receptor_pdbqt
        root = Path(__file__).resolve().parent
        pdb_out = self.root / "2nv6_pdb.pdbqt"
        cif_out = self.root / "2nv6_cif.pdbqt"
        prepare_receptor_pdbqt(root / "2NV6 (2).pdb", pdb_out, cleanup_water=True, add_hydrogens=False)
        prepare_receptor_pdbqt(root / "2NV6 (1).cif", cif_out, cleanup_water=True, add_hydrogens=False)

        def zid(path):
            return sorted((line[12:16].strip(), line[21], line[22:26].strip(), line[30:54])
                          for line in path.read_text(errors="replace").splitlines()
                          if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "ZID")

        self.assertEqual(zid(pdb_out), zid(cif_out))
        pdb_out.unlink(missing_ok=True)
        cif_out.unlink(missing_ok=True)

    def test_unknown_atom_type_fails_preflight(self):
        receptor = self._pdbqt("unknown_receptor.pdbqt", [pdbqt_atom(1, "X1", "UNK", "TZ")])
        ligand = self._pdbqt("ligand.pdbqt", [pdbqt_atom(1, "C1", "LIG", "C")])
        with self.assertRaises(AtomTypeCompatibilityError) as caught:
            preflight_job([("receptor", receptor), ("ligand", ligand)])
        self.assertIn("TZ", str(caught.exception))

    def test_ordinary_gpf_dpf_syntax_is_preserved(self):
        receptor = self._pdbqt("rec.pdbqt", [pdbqt_atom(1, "C", "ALA", "C")])
        ligand = self._pdbqt("lig.pdbqt", [pdbqt_atom(1, "C1", "LIG", "C")])
        cfg = ProjectConfig()
        cfg.ga_run = 1
        gpf = generate_gpf(receptor, ligand, self.root / "rec.gpf", cfg)
        dpf = generate_dpf(receptor, ligand, self.root / "lig.dpf", cfg)
        gpf_text = gpf.read_text(encoding="utf-8")
        dpf_text = dpf.read_text(encoding="utf-8")
        self.assertIn("receptor_types C", gpf_text)
        self.assertIn("ligand_types C", gpf_text)
        self.assertIn("map rec.C.map", gpf_text)
        self.assertIn("ligand_types C", dpf_text)
        self.assertIn("move lig.pdbqt", dpf_text)
        self.assertNotIn("nbp_r_eps", gpf_text)
        self.assertFalse(STANDARD_AD4_PROFILE.specialized_metal_workflow)

    def test_preflight_blocks_unsupported_type_before_autogrid(self):
        receptor = self._pdbqt("rec.pdbqt", [pdbqt_atom(1, "TZ", "ZN", "TZ")])
        ligand = self._pdbqt("lig.pdbqt", [pdbqt_atom(1, "C1", "LIG", "C")])
        config = ProjectConfig()
        config.log_directory = self.root / "logs"
        config.ad4_parameter_profile = STANDARD_AD4_PROFILE.profile_id
        job = DockingJob(
            job_id="unsupported_type",
            receptor_name="rec",
            ligand_name="lig",
            receptor_path=receptor,
            ligand_path=ligand,
            output_dir=self.root / "out",
            engine=Engine.AUTODOCK4,
        )
        with patch("autodock4_workflow._run_autogrid4") as autogrid:
            result = run_ad4_job(job, config)
        for handler in logging.getLogger("docking_automation.AD4.unsupported_type").handlers[:]:
            handler.close()
            logging.getLogger("docking_automation.AD4.unsupported_type").removeHandler(handler)
        autogrid.assert_not_called()
        self.assertEqual(result.status, JobStatus.INVALID_INPUT)
        self.assertEqual(result.execution_status, ExecutionStatus.FAILED)
        self.assertIn("compatibility", result.ad4_results)


if __name__ == "__main__":
    unittest.main()
