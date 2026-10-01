"""
AutoDock Suite Pro — Unified GUI Package (gui/__init__.py)
==========================================================
Prioritizes the hardware-accelerated PySide6 (Qt6) interface,
with automatic fallback to CustomTkinter if PySide6 is unavailable.
"""

try:
    import PySide6
    from gui_qt.app import run_gui, MainWindow
    # Aliases
    launch_gui = run_gui
    AutoDockSuiteProApp = MainWindow
    AutoDockSuiteProGUI = MainWindow
    HAS_PYSIDE6 = True
except ImportError:
    from gui.app import run_gui, AutoDockSuiteProApp
    launch_gui = run_gui
    AutoDockSuiteProGUI = AutoDockSuiteProApp
    HAS_PYSIDE6 = False

from gui import themes

__all__ = ["run_gui", "launch_gui", "AutoDockSuiteProApp", "AutoDockSuiteProGUI", "themes", "HAS_PYSIDE6"]
