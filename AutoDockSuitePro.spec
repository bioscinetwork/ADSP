# -*- mode: python ; coding: utf-8 -*-
# AutoDock Suite Pro - PyInstaller Spec (PORTABLE TEMPLATE)
#
# PURPOSE
# -------
# This file is a PORTABLE template.  All paths that are machine- or
# installation-specific (site-packages, Python interpreter location, etc.)
# are resolved dynamically at build time using SPECPATH (provided by
# PyInstaller) so this file can be committed and shared across machines.
#
# PREFERRED BUILD PATH
# --------------------
#     python build_exe.py
#
# DIRECT PYINSTALLER INVOCATION
# ------------------------------
#     pyinstaller --noconfirm --clean AutoDockSuitePro.spec
#
# NOTE ON AD4 PARAMETER FILES
# ----------------------------
# AD4_parameters.dat and AD4.1_bound.dat are resolved by ad4_compatibility.py
# relative to Path(__file__).resolve().parent.  PyInstaller collects them
# via the datas entry below, placing them alongside ad4_compatibility.py in
# _internal/ so the module-relative lookup continues to work in the bundle.

import sys
from pathlib import Path

# SPECPATH is provided by PyInstaller at build time.
_here = Path(SPECPATH).resolve()

_ad4_params  = _here / 'parameter_profiles' / 'ad4_standard_4.2' / 'AD4_parameters.dat'
_ad41_bound  = _here / 'parameter_profiles' / 'ad4_1_bound' / 'AD4.1_bound.dat'
_ad4zn_zip   = _here / 'external' / 'autodock4zn' / 'AutoDock4Zn-Pipeline-main.zip'
_local_bin   = _here / 'bin'
_assets_dir  = _here / 'gui' / 'assets'


def _site_pkg_dir(pkg_name):
    try:
        import importlib.util
        spec = importlib.util.find_spec(pkg_name)
        if spec and spec.origin:
            return Path(spec.origin).resolve().parent
    except Exception:
        pass
    return None


_ctk_path   = _site_pkg_dir('customtkinter') or (_here / 'customtkinter')
_meeko_path = _site_pkg_dir('meeko')
_ob_path    = _site_pkg_dir('openbabel')

_datas = [
    (str(_here / 'project_config.toml'), '.'),
]

if _assets_dir.is_dir():
    _datas.append((str(_assets_dir), 'gui/assets'))

if _ctk_path and Path(_ctk_path).is_dir():
    _datas.append((str(_ctk_path), 'customtkinter'))

if (_here / 'parameter_profiles').is_dir():
    _datas.append((str(_here / 'parameter_profiles'), 'parameter_profiles'))

if (_here / 'external').is_dir():
    _datas.append((str(_here / 'external'), 'external'))

if (_here / 'data').is_dir():
    _datas.append((str(_here / 'data'), 'data'))

# AD4 scientific parameter assets -- REQUIRED for packaged AD4 operation.
# These are placed at the root of _internal/ so ad4_compatibility.py's
# Path(__file__).resolve().parent lookup finds them after bundling.
_param_candidates = [
    (_ad4_params, 'AD4_parameters.dat'),
    (_ad41_bound, 'AD4.1_bound.dat'),
    (_ad4zn_zip, 'AutoDock4Zn-Pipeline-main.zip'),
    (_here / 'AD4_parameters.dat', 'AD4_parameters.dat'),
    (_here / 'AD4.1_bound.dat', 'AD4.1_bound.dat'),
    (_here / 'AutoDock4Zn-Pipeline-main.zip', 'AutoDock4Zn-Pipeline-main.zip'),
]
for _src, _dst_name in _param_candidates:
    if _src.is_file():
        _datas.append((str(_src), '.'))

if _meeko_path:
    _meeko_data = Path(_meeko_path) / 'data'
    if _meeko_data.is_dir():
        _datas.append((str(_meeko_data), 'meeko/data'))

if _ob_path:
    _ob_bin = Path(_ob_path) / 'bin'
    if _ob_bin.is_dir():
        _datas.append((str(_ob_bin), 'openbabel/bin'))
        _datas.append((str(_ob_bin), 'openbabel/lib/openbabel/3.2.1'))
        _ob_data = _ob_bin / 'data'
        if _ob_data.is_dir():
            _datas.append((str(_ob_data), 'openbabel/bin/data'))
            _datas.append((str(_ob_data), 'openbabel/share/openbabel/3.2.1'))
            _datas.append((str(_ob_data), 'bin/data'))

if _local_bin.is_dir():
    _datas.append((str(_local_bin), 'bin'))
    _datas.append((str(_local_bin), '_internal/bin'))
    _local_bin_data = _local_bin / 'data'
    if _local_bin_data.is_dir():
        _datas.append((str(_local_bin_data), 'bin/data'))
        _datas.append((str(_local_bin_data), '_internal/bin/data'))

_binaries = []
for _exe_name in ('vina.exe', 'vina_split.exe', 'autodock4.exe', 'autogrid4.exe', 'obabel.exe'):
    _candidate = _local_bin / _exe_name
    if _candidate.is_file():
        _binaries.append((str(_candidate), 'bin'))

block_cipher = None

a = Analysis(
    ['main.py'],
    pathex=[str(_here), str(_here / 'scripts')],
    binaries=_binaries,
    datas=_datas,
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
        'ad4_compatibility',
        'ad4zn_adapter',
        'rmsd_validation',
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
    hooksconfig={},
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
    icon=str(_assets_dir / 'app_icon.ico') if (_assets_dir / 'app_icon.ico').is_file() else None,
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
