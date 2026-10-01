@echo off
title AutoDock Suite Pro — Dependency Installer
color 0A
echo ======================================================================
echo   AutoDock Suite Pro — Automated Environment & Dependency Setup
echo ======================================================================
echo.

echo [1/3] Checking Python installation...
where python >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    color 0C
    echo [ERROR] Python was not found in PATH!
    echo Please install Python 3.10+ (64-bit) from https://www.python.org/
    pause
    exit /b 1
)

python --version
echo.

echo [2/3] Installing and upgrading required Python packages...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if %ERRORLEVEL% NEQ 0 (
    color 0C
    echo [ERROR] Failed to install dependencies via pip.
    pause
    exit /b 1
)
echo.

echo [3/3] Verifying environment & scientific toolchains...
python -c "import customtkinter, openpyxl, pandas, numpy, scipy; print('  [OK] GUI & Numerical stack (CustomTkinter, OpenPyXL, Pandas, NumPy, SciPy)')"
python -c "import meeko; print('  [OK] Scripps Meeko (v' + meeko.__version__ + ')')"
python -c "import rdkit; print('  [OK] RDKit (v' + rdkit.__version__ + ')')"
python -c "import openbabel; print('  [OK] OpenBabel (v' + openbabel.__version__ + ')')"

echo.
echo ======================================================================
echo   SUCCESS: All dependencies installed and ready!
echo ======================================================================
echo.
echo To launch the graphical interface, run: run_gui.bat
echo To run the docking pipeline from CLI, run: run_docking.bat
echo.
pause
