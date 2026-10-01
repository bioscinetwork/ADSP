@echo off
title Docking Automation Suite
echo ======================================================================
echo   Docking Automation Suite
echo   Cross-Engine Molecular Docking Automation Platform
echo ======================================================================
echo.

:: Check if python is available
where python >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Python is not found in system PATH.
    echo Please install Python 3.10+ or add it to PATH.
    echo.
    pause
    exit /b 1
)

:: Run the pipeline
python "%~dp0main.py" %*

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ======================================================================
    echo [NOTICE] Docking pipeline finished with exit code %ERRORLEVEL%.
    echo ======================================================================
    pause
) else (
    echo.
    echo ======================================================================
    echo [SUCCESS] Docking pipeline completed successfully!
    echo ======================================================================
    pause
)
