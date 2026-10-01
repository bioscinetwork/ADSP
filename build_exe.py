"""
AutoDock Suite Pro — Standalone Packaging Script
================================================
Compiles AutoDock Suite Pro into a standalone Windows executable distribution
using PyInstaller. Packages Python runtime, CustomTkinter assets, RDKit,
Meeko, OpenBabel, SciPy, and docking binaries (vina.exe, autodock4.exe,
autogrid4.exe, vina_split.exe).

Usage:
    python build_exe.py
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional


def ensure_pyinstaller() -> None:
    """Ensure pyinstaller is installed."""
    try:
        import PyInstaller
        print(f"[OK] PyInstaller is installed (v{PyInstaller.__version__})")
    except ImportError:
        print("[INFO] Installing PyInstaller via pip...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])


def find_binaries() -> Dict[str, Optional[Path]]:
    """Locate local, distribution, or system docking binaries to bundle."""
    bins: Dict[str, Optional[Path]] = {
        "vina.exe": None,
        "vina_split.exe": None,
        "autodock4.exe": None,
        "autogrid4.exe": None,
        "obabel.exe": None,
    }

    # Search candidates in priority order
    search_dirs = [
        Path("bin"),
        Path("dist/AutoDockSuitePro/_internal/bin"),
        Path("dist/AutoDockSuitePro/bin"),
    ]
    try:
        import openbabel
        ob_dir = Path(openbabel.__file__).parent / "bin"
        if ob_dir.is_dir():
            search_dirs.append(ob_dir)
    except Exception:
        pass

    if sys.prefix:
        search_dirs.append(Path(sys.prefix) / "Scripts")

    appdata = os.environ.get("APPDATA")
    if appdata:
        search_dirs.append(Path(appdata) / "Python" / f"Python{sys.version_info.major}{sys.version_info.minor}" / "site-packages" / "openbabel" / "bin")

    prog_files = os.environ.get("ProgramFiles", "C:/Program Files")
    prog_files_x86 = os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")
    for base in [prog_files, prog_files_x86]:
        search_dirs.append(Path(base) / "MGLTools-1.5.7")
        search_dirs.append(Path(base) / "The Scripps Research Institute" / "Vina")

    for name in bins.keys():
        for d in search_dirs:
            candidate = d / name
            if candidate.is_file():
                bins[name] = candidate.resolve()
                break
        if not bins[name]:
            which_p = shutil.which(name)
            if which_p:
                bins[name] = Path(which_p).resolve()

    return bins


def generate_spec_file(bins: Dict[str, Optional[Path]]) -> Path:
    """Generate PyInstaller spec file."""
    import customtkinter
    ctk_path = Path(customtkinter.__file__).parent.resolve()

    # Collect data files
    datas = [
        f"('{ctk_path.as_posix()}', 'customtkinter')",
        "('project_config.toml', '.')",
    ]
    assets_dir = Path("gui/assets").resolve()
    if assets_dir.is_dir():
        datas.append(f"('{assets_dir.as_posix()}', 'gui/assets')")

    # Bundle Meeko parameter data files
    try:
        import meeko
        meeko_data = Path(meeko.__file__).parent.resolve() / "data"
        if meeko_data.is_dir():
            datas.append(f"('{meeko_data.as_posix()}', 'meeko/data')")
            print(f"  + Bundling Meeko data: {meeko_data} -> meeko/data")
    except Exception as e:
        print(f"  - Warning: Could not locate Meeko data: {e}")

    # Bundle OpenBabel plugins (.obf) and data files
    try:
        import openbabel
        ob_pkg = Path(openbabel.__file__).parent.resolve()
        ob_bin = ob_pkg / "bin"
        if ob_bin.is_dir():
            datas.append(f"('{ob_bin.as_posix()}', 'openbabel/bin')")
            datas.append(f"('{ob_bin.as_posix()}', 'openbabel/lib/openbabel/3.2.1')")
            ob_data = ob_bin / "data"
            if ob_data.is_dir():
                datas.append(f"('{ob_data.as_posix()}', 'openbabel/bin/data')")
                datas.append(f"('{ob_data.as_posix()}', 'openbabel/share/openbabel/3.2.1')")
                datas.append(f"('{ob_data.as_posix()}', 'bin/data')")
            print(f"  + Bundling OpenBabel plugins & data from {ob_bin}")
    except Exception as e:
        print(f"  - Warning: Could not locate OpenBabel data: {e}")

    # Bundle local bin directory contents (plugins, DLLs, data)
    local_bin = Path("bin").resolve()
    if local_bin.is_dir():
        datas.append(f"('{local_bin.as_posix()}', 'bin')")
        datas.append(f"('{local_bin.as_posix()}', '_internal/bin')")
        if (local_bin / "data").is_dir():
            datas.append(f"('{(local_bin / 'data').as_posix()}', 'bin/data')")
            datas.append(f"('{(local_bin / 'data').as_posix()}', '_internal/bin/data')")

    # Add discovered binaries into bin/ folder of distribution
    binaries = []
    for name, path in bins.items():
        if path and path.is_file():
            binaries.append(f"('{path.as_posix()}', 'bin')")
            print(f"  + Bundling binary: {name} from {path}")
        else:
            print(f"  - Binary not found locally: {name} (will prompt user in GUI)")

    datas_str = ",\n        ".join(datas)
    binaries_str = ",\n        ".join(binaries) if binaries else ""

    spec_content = f"""# -*- mode: python ; coding: utf-8 -*-
# AutoDock Suite Pro — PyInstaller Spec

block_cipher = None

a = Analysis(
    ['main.py'],
    pathex=['.'],
    binaries=[
        {binaries_str}
    ],
    datas=[
        {datas_str}
    ],
    hiddenimports=[
        'customtkinter',
        'tkinterweb',
        'PySide6',
        'PySide6.QtCore',
        'PySide6.QtGui',
        'PySide6.QtWidgets',
        'PySide6.QtWebEngineWidgets',
        'PySide6.QtSvgWidgets',
        'PySide6.QtSvg',
        'openpyxl',
        'interactions',
        'depiction',
        'dlg_extract',
        'viewer_3d',
        'pymol_exporter',
        'gui',
        'gui.app',
        'gui.themes',
        'gui.widgets',
        'gui.settings',
        'gui.workspace_tab',
        'gui.preparation_tab',
        'gui.grid_tab',
        'gui.library_tab',
        'gui.runner_tab',
        'gui.results_tab',
        'gui.interaction_diagram',
        'gui_qt',
        'gui_qt.app',
        'gui_qt.styles',
        'gui_qt.workspace_tab',
        'gui_qt.preparation_tab',
        'gui_qt.grid_tab',
        'gui_qt.library_tab',
        'gui_qt.runner_tab',
        'gui_qt.results_tab',
        'gui_qt.comparison_tab',
        'gui_qt.settings_tab',
        'workflow_service',
        'complex_builder',
        'process_manager',
        'rename_ligands',
        'vina_workflow',
        'autodock4_workflow',
        'job_manager',
        'reporting',
        'validators',
        'vina_parser',
        'vina_splitter',
        'prepare',
        'executables',
        'config',
        'models',
        'logging_utils',
        'meeko',
        'rdkit',
        'rdkit.Chem',
        'rdkit.Chem.AllChem',
        'rdkit.Chem.Descriptors',
        'rdkit.Chem.Lipinski',
        'rdkit.Chem.Draw',
        'openbabel',
        'openbabel.pybel',
        'scipy',
        'scipy.spatial',
        'numpy',
        'PIL',
        'PIL.Image',
        'PIL.ImageTk',
    ],
    hookspath=[],
    hooksconfig={{}},
    runtime_hooks=[],
    excludes=['tkinter.test', 'matplotlib', 'torch', 'IPython'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='AutoDockSuitePro',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='gui/assets/app_icon.ico',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='AutoDockSuitePro',
)
"""
    spec_path = Path("AutoDockSuitePro.spec")
    spec_path.write_text(spec_content, encoding="utf-8")
    print(f"[OK] Generated {spec_path}")
    return spec_path


def copy_distribution_assets(dist_dir: Path, bins: Dict[str, Optional[Path]]) -> None:
    """Copies all configuration, binaries, and sample data into dist for 100% self-containment."""
    print("\n[INFO] Assembling self-contained portable distribution...")

    # 1. Ensure bin/ directories exist and copy all tools, DLLs, and data
    bin_dir = dist_dir / "bin"
    internal_bin_dir = dist_dir / "_internal" / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    internal_bin_dir.mkdir(parents=True, exist_ok=True)

    # Copy explicitly discovered named binaries
    for name, path in bins.items():
        if path and path.is_file():
            shutil.copy2(path, bin_dir / name)
            shutil.copy2(path, internal_bin_dir / name)
            print(f"  + Placed {name} in {bin_dir} and {internal_bin_dir}")

    # Copy full contents of local bin/ (including openbabel-3.dll, plugins, and data/)
    local_bin = Path("bin")
    if local_bin.is_dir():
        for item in local_bin.iterdir():
            dst_target = bin_dir / item.name
            dst_internal = internal_bin_dir / item.name
            if item.is_dir():
                if dst_target.exists():
                    shutil.rmtree(dst_target)
                shutil.copytree(item, dst_target)
                if dst_internal.exists():
                    shutil.rmtree(dst_internal)
                shutil.copytree(item, dst_internal)
                print(f"  + Bundled directory {item.name}/ into {bin_dir} and {internal_bin_dir}")
            elif item.is_file():
                shutil.copy2(item, dst_target)
                shutil.copy2(item, dst_internal)

    # Copy Meeko parameter data into _internal/meeko/data for standalone ligand preparation
    try:
        import meeko
        meeko_src = Path(meeko.__file__).parent / "data"
        if meeko_src.is_dir():
            dst_meeko = dist_dir / "_internal" / "meeko" / "data"
            dst_meeko.parent.mkdir(parents=True, exist_ok=True)
            if dst_meeko.exists():
                shutil.rmtree(dst_meeko)
            shutil.copytree(meeko_src, dst_meeko)
            # Also copy to dist/meeko/data for root-relative lookups
            dst_meeko_root = dist_dir / "meeko" / "data"
            dst_meeko_root.parent.mkdir(parents=True, exist_ok=True)
            if dst_meeko_root.exists():
                shutil.rmtree(dst_meeko_root)
            shutil.copytree(meeko_src, dst_meeko_root)
            print(f"  + Bundled Meeko parameter data into {dst_meeko}")
    except Exception as e:
        print(f"  - Warning copying Meeko data: {e}")

    # Copy OpenBabel plugins (.obf) and data files into openbabel module paths
    try:
        import openbabel
        ob_src = Path(openbabel.__file__).parent
        ob_bin_src = ob_src / "bin"
        if not ob_bin_src.is_dir():
            ob_bin_src = local_bin.resolve()

        plugin_targets = [
            dist_dir / "bin",
            dist_dir / "_internal" / "bin",
            dist_dir / "_internal" / "openbabel" / "bin",
            dist_dir / "_internal" / "openbabel" / "lib" / "openbabel" / "3.2.1",
        ]
        data_targets = [
            dist_dir / "bin" / "data",
            dist_dir / "_internal" / "bin" / "data",
            dist_dir / "_internal" / "openbabel" / "bin" / "data",
            dist_dir / "_internal" / "openbabel" / "share" / "openbabel" / "3.2.1",
        ]

        for pt in plugin_targets:
            pt.mkdir(parents=True, exist_ok=True)
            for obf in ob_bin_src.glob("*.obf"):
                shutil.copy2(obf, pt / obf.name)
            for dll in ob_bin_src.glob("*.dll"):
                shutil.copy2(dll, pt / dll.name)

        ob_data_src = ob_bin_src / "data" if (ob_bin_src / "data").is_dir() else (local_bin / "data")
        if ob_data_src.is_dir():
            for dt in data_targets:
                dt.mkdir(parents=True, exist_ok=True)
                for df in ob_data_src.iterdir():
                    if df.is_file():
                        shutil.copy2(df, dt / df.name)

        print("  + Bundled OpenBabel plugins and parameters into distribution")
    except Exception as e:
        print(f"  - Warning copying OpenBabel plugins: {e}")

    # 2. Copy and adapt project_config.toml to dist root with portable relative paths
    cfg_src = Path("project_config.toml")
    if cfg_src.is_file():
        cfg_text = cfg_src.read_text(encoding="utf-8")
        cfg_text = re.sub(r'vina\s*=\s*".*?[/\\]bin[/\\]vina\.exe"', 'vina       = "bin/vina.exe"', cfg_text)
        cfg_text = re.sub(r'vina_split\s*=\s*".*?[/\\]bin[/\\]vina_split\.exe"', 'vina_split = "bin/vina_split.exe"', cfg_text)
        cfg_text = re.sub(r'autogrid4\s*=\s*".*?[/\\]bin[/\\]autogrid4\.exe"', 'autogrid4  = "bin/autogrid4.exe"', cfg_text)
        cfg_text = re.sub(r'autodock4\s*=\s*".*?[/\\]bin[/\\]autodock4\.exe"', 'autodock4  = "bin/autodock4.exe"', cfg_text)
        cfg_text = re.sub(r'obabel\s*=\s*".*?[/\\]bin[/\\]obabel\.exe"', 'obabel     = "bin/obabel.exe"', cfg_text)
        cfg_text = re.sub(r'receptor_directory\s*=\s*".*?"', 'receptor_directory = "receptors"', cfg_text)
        cfg_text = re.sub(r'ligand_directory\s*=\s*".*?"', 'ligand_directory   = "ligands"', cfg_text)
        cfg_text = re.sub(r'result_directory\s*=\s*".*?"', 'result_directory = "results"', cfg_text)
        cfg_text = re.sub(r'log_directory\s*=\s*".*?"', 'log_directory    = "logs"', cfg_text)
        cfg_text = re.sub(r'report_directory\s*=\s*".*?"', 'report_directory = "reports"', cfg_text)
        if "[ui]" in cfg_text:
            cfg_text = re.sub(r'theme\s*=\s*".*?"', 'theme = "Noir"', cfg_text)
        else:
            cfg_text += '\n[ui]\ntheme = "Noir"\n'
        (dist_dir / "project_config.toml").write_text(cfg_text, encoding="utf-8")
        print("  + Generated portable project_config.toml in distribution root (theme=Noir)")

    # Copy dlg_extract.py for standalone AD4 result extraction
    dlg_py = Path("dlg_extract.py")
    if dlg_py.is_file():
        shutil.copy2(dlg_py, dist_dir / "dlg_extract.py")
        print("  + Copied dlg_extract.py to distribution root")

    # 3. Create / Copy working directories
    for folder_name in ("results", "reports", "logs"):
        (dist_dir / folder_name).mkdir(parents=True, exist_ok=True)

    # 4. Copy receptors and ligands
    for data_dir in ("receptors", "ligands"):
        src = Path(data_dir)
        dst = dist_dir / data_dir
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)
            print(f"  + Bundled {data_dir}/ into distribution")

    # 5. Create standalone executable launcher batch files for users without Python
    launch_gui_bat = dist_dir / "Launch_GUI.bat"
    launch_gui_bat.write_text(
        '@echo off\n'
        'title AutoDock Suite Pro\n'
        'start "" "%~dp0AutoDockSuitePro.exe"\n',
        encoding="utf-8"
    )
    print("  + Created Launch_GUI.bat in distribution")

    run_cli_bat = dist_dir / "Run_CLI.bat"
    run_cli_bat.write_text(
        '@echo off\n'
        'title AutoDock Suite Pro CLI\n'
        '"%~dp0AutoDockSuitePro.exe" %*\n'
        'pause\n',
        encoding="utf-8"
    )
    print("  + Created Run_CLI.bat in distribution")

    # 6. Copy documentation and configuration files
    for doc_file in ("USER_GUIDE.md", "DEVELOPER_GUIDE.md"):
        b_src = Path(doc_file)
        if b_src.is_file():
            shutil.copy2(b_src, dist_dir / doc_file)
            print(f"  + Copied {doc_file} to distribution")


def create_zip_package(dist_dir: Path) -> Path:
    """Create a single portable zip archive of the distribution folder."""
    zip_base = dist_dir.parent / "AutoDockSuitePro_Portable"
    print(f"\n[INFO] Creating portable ZIP archive: {zip_base}.zip ...")
    archive_path = Path(shutil.make_archive(str(zip_base), "zip", root_dir=str(dist_dir.parent), base_dir=dist_dir.name))
    size_mb = archive_path.stat().st_size / (1024 * 1024)
    print(f"[OK] Portable ZIP created successfully: {archive_path.resolve()} ({size_mb:.1f} MB)")
    return archive_path


def kill_running_instances() -> None:
    """Kill any running instances that could lock files in dist/."""
    if sys.platform == "win32":
        for proc in ("AutoDockSuitePro.exe", "autodock4.exe", "autogrid4.exe", "vina.exe", "vina_split.exe", "obabel.exe"):
            try:
                subprocess.run(["taskkill", "/F", "/IM", proc], capture_output=True, text=True)
            except Exception:
                pass


def build() -> None:
    """Execute PyInstaller build and distribution assembly."""
    print("=" * 60)
    print("  AutoDock Suite Pro — Standalone Executable Builder")
    print("=" * 60)

    kill_running_instances()
    ensure_pyinstaller()
    bins = find_binaries()
    spec_path = generate_spec_file(bins)

    print("\n[INFO] Running PyInstaller...")
    cmd = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", str(spec_path)]
    subprocess.check_call(cmd)

    dist_dir = Path("dist/AutoDockSuitePro")
    copy_distribution_assets(dist_dir, bins)
    zip_path = create_zip_package(dist_dir)

    print("\n" + "=" * 60)
    print("  BUILD SUCCESSFUL — FULL SELF-CONTAINED ENVIRONMENT READY!")
    print("=" * 60)
    print(f"  Distribution directory: {dist_dir.resolve()}")
    print(f"  Executable:             {dist_dir.resolve() / 'AutoDockSuitePro.exe'}")
    print(f"  Portable ZIP package:   {zip_path.resolve()}")
    print("\nTo run the portable application, double click AutoDockSuitePro.exe or run:")
    print(f'  "{dist_dir.resolve() / "AutoDockSuitePro.exe"}"')
    print(f"\nTo share with colleagues or students, simply send:")
    print(f'  "{zip_path.resolve()}"')
    print("=" * 60)


if __name__ == "__main__":
    build()
