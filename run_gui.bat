@echo off
title AutoDock Suite Pro GUI
python main.py --gui
if %errorlevel% neq 0 (
    echo.
    echo An error occurred launching the GUI.
    pause
)
