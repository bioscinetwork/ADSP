"""
AutoDock Suite Pro — Reusable Widget Helpers (gui/widgets.py)
=============================================================
Thin wrappers over customtkinter and ttk widgets that automatically apply
the current theme colors and standardized typography scale.
Keeps tab code clean, readable, and visually consistent.
"""

from __future__ import annotations
from typing import Any, Callable, Optional
import tkinter as tk
from tkinter import ttk
import customtkinter as ctk

from gui import themes


# ---------------------------------------------------------------------------
# Labelled Entry Row
# ---------------------------------------------------------------------------

def labelled_entry(
    parent: ctk.CTkFrame,
    label_text: str,
    row: int,
    default: str = "",
    width: int = 380,
) -> ctk.CTkEntry:
    """Creates a label + entry pair on a 2-column grid row with standardized typography."""
    t = themes.get()
    lbl = ctk.CTkLabel(
        parent,
        text=label_text,
        text_color=t.text_secondary,
        font=themes.font_body(),
        anchor="w",
    )
    lbl.grid(row=row, column=0, sticky="w", padx=(14, 8), pady=4)

    entry = ctk.CTkEntry(
        parent,
        width=width,
        height=32,
        corner_radius=6,
        fg_color=t.bg_tertiary,
        border_color=t.border,
        text_color=t.text_primary,
        font=themes.font_body(),
    )
    entry.insert(0, default)
    entry.grid(row=row, column=1, sticky="ew", padx=(0, 14), pady=4)
    return entry


# ---------------------------------------------------------------------------
# Section Header
# ---------------------------------------------------------------------------

def section_header(
    parent: Any,
    text: str,
    row: int,
    col: int = 0,
    colspan: int = 2,
) -> ctk.CTkLabel:
    """A bold section-header label inside a frame with standard typography."""
    t = themes.get()
    lbl = ctk.CTkLabel(
        parent,
        text=text.upper(),
        font=themes.font_section(),
        text_color=t.accent,
        anchor="w",
    )
    lbl.grid(
        row=row,
        column=col,
        columnspan=colspan,
        sticky="w",
        padx=14,
        pady=(14, 4),
    )
    return lbl


# ---------------------------------------------------------------------------
# Thin Horizontal Divider
# ---------------------------------------------------------------------------

def divider(parent: Any, row: int, col: int = 0, colspan: int = 2) -> tk.Frame:
    t = themes.get()
    sep = tk.Frame(parent, bg=t.border, height=1)
    sep.grid(
        row=row,
        column=col,
        columnspan=colspan,
        sticky="ew",
        padx=14,
        pady=6,
    )
    return sep


# ---------------------------------------------------------------------------
# Accent Button
# ---------------------------------------------------------------------------

def accent_button(
    parent: Any,
    text: str,
    command: Callable,
    width: int = 160,
    small: bool = False,
) -> ctk.CTkButton:
    """Accent-styled action button with proper touch/click targets."""
    t = themes.get()
    h = 28 if small else 34
    btn = ctk.CTkButton(
        parent,
        text=text,
        command=command,
        width=width,
        height=h,
        fg_color=t.btn_fg,
        hover_color=t.btn_hover,
        text_color=t.btn_text,
        font=themes.font_button(),
        corner_radius=6,
    )
    return btn


# ---------------------------------------------------------------------------
# Styled Treeview (ttk)
# ---------------------------------------------------------------------------

def styled_treeview(
    parent: Any,
    columns: tuple,
    height: int = 8,
) -> ttk.Treeview:
    """Creates a ttk.Treeview pre-styled with clean row height and legible typography."""
    t = themes.get()

    style = ttk.Style()
    style_name = "AutoDock.Treeview"
    try:
        style.theme_use("default")
    except Exception:
        pass

    style.configure(
        style_name,
        background=t.row_even,
        foreground=t.text_primary,
        fieldbackground=t.row_even,
        rowheight=28,  # Generous row height for high readability
        font=("Segoe UI", 11),
    )
    style.configure(
        f"{style_name}.Heading",
        background=t.bg_secondary,
        foreground=t.accent,
        font=("Segoe UI", 11, "bold"),
        relief="flat",
        padding=(6, 4),
    )
    style.map(
        style_name,
        background=[("selected", t.row_selected)],
        foreground=[("selected", "#ffffff" if t.ctk_mode == "dark" else t.text_primary)],
    )

    tree = ttk.Treeview(
        parent,
        columns=columns,
        show="headings",
        height=height,
        style=style_name,
    )
    tree.tag_configure("odd", background=t.row_odd)
    tree.tag_configure("even", background=t.row_even)
    return tree


def tree_insert_alternating(tree: ttk.Treeview, values: tuple) -> str:
    """Insert a row with alternating row background."""
    count = len(tree.get_children())
    tag = "even" if count % 2 == 0 else "odd"
    return tree.insert("", "end", values=values, tags=(tag,))


# ---------------------------------------------------------------------------
# Status Badge Label
# ---------------------------------------------------------------------------

def status_badge(
    parent: Any,
    text: str,
    ok: Optional[bool] = None,
) -> ctk.CTkLabel:
    """Clean status label with proper font weight.
    ok=True -> status_ok_text, ok=False -> status_fail_text, ok=None -> text_secondary.
    """
    t = themes.get()
    if ok is True:
        color = t.status_ok_text
    elif ok is False:
        color = t.status_fail_text
    else:
        color = t.text_secondary
    return ctk.CTkLabel(
        parent,
        text=text,
        text_color=color,
        font=themes.font_caption_bold(),
    )


# ---------------------------------------------------------------------------
# Console Text Redirector  (color-coded log output)
# ---------------------------------------------------------------------------

class TextRedirector:
    """Redirects stdout/stderr writes to a CTkTextbox widget (thread-safe).

    Lines containing specific markers are displayed in appropriate colors:
      [OK]   / ✔   → status_ok_text  (green)
      [WARN] / ⚠   → status_warn_text (amber)
      [FAIL] / [ERROR] / ✖ → status_fail_text (red)
    """

    def __init__(self, widget: ctk.CTkTextbox) -> None:
        self.widget = widget
        self._setup_tags()

    def _setup_tags(self) -> None:
        """Configure color tags on the underlying tk.Text widget."""
        try:
            t = themes.get()
            tw = self.widget._textbox  # access underlying tk.Text
            tw.tag_configure("ok",   foreground=t.status_ok_text)
            tw.tag_configure("warn", foreground=t.status_warn_text)
            tw.tag_configure("fail", foreground=t.status_fail_text)
            tw.tag_configure("accent", foreground=t.accent)
        except Exception:
            pass

    def _classify(self, text: str) -> str:
        """Return tag name for a line of text."""
        upper = text.upper()
        if any(m in upper for m in ("[OK]", "✔", "[SUCCESS]", "SUCCESS")):
            return "ok"
        if any(m in upper for m in ("[WARN]", "⚠", "[WARNING]", "WARNING")):
            return "warn"
        if any(m in upper for m in ("[FAIL]", "[ERROR]", "✖", "ERROR", "FAILED")):
            return "fail"
        if any(m in upper for m in ("PHASE", "DUAL-ENGINE", "ENGINE:", "AUTODOCK", "VINA")):
            return "accent"
        return ""

    def write(self, text: str) -> None:
        try:
            self.widget.after(0, self._do_write, text)
        except Exception:
            pass

    def _do_write(self, text: str) -> None:
        try:
            self.widget.configure(state="normal")
            tw = self.widget._textbox
            for line in text.splitlines(keepends=True):
                tag = self._classify(line)
                if tag:
                    tw.insert("end", line, (tag,))
                else:
                    tw.insert("end", line)
            tw.see("end")
            self.widget.configure(state="disabled")
        except Exception:
            try:
                self.widget.configure(state="normal")
                self.widget.insert("end", text)
                self.widget.see("end")
                self.widget.configure(state="disabled")
            except Exception:
                pass

    def flush(self) -> None:
        pass


# ---------------------------------------------------------------------------
# ToolTip — hover-activated popup for any widget
# ---------------------------------------------------------------------------

class ToolTip:
    """Lightweight tooltip that appears on mouse-over.

    Usage::

        ToolTip(some_button, "Click to start docking")
    """

    def __init__(self, widget: tk.Widget, text: str, delay_ms: int = 600) -> None:
        self.widget = widget
        self.text = text
        self.delay_ms = delay_ms
        self._tip_window: Optional[tk.Toplevel] = None
        self._after_id: Optional[str] = None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._cancel)
        widget.bind("<ButtonPress>", self._cancel)

    def _schedule(self, _event=None) -> None:
        self._cancel()
        self._after_id = self.widget.after(self.delay_ms, self._show)

    def _cancel(self, _event=None) -> None:
        if self._after_id:
            self.widget.after_cancel(self._after_id)
            self._after_id = None
        self._hide()

    def _show(self) -> None:
        if self._tip_window:
            return
        t = themes.get()
        try:
            x = self.widget.winfo_rootx() + 20
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
            self._tip_window = tw = tk.Toplevel(self.widget)
            tw.wm_overrideredirect(True)
            tw.wm_geometry(f"+{x}+{y}")
            tk.Label(
                tw, text=self.text,
                background=t.bg_tertiary,
                foreground=t.text_primary,
                relief="flat", borderwidth=0,
                font=("Segoe UI", 10),
                padx=8, pady=4,
            ).pack()
        except Exception:
            self._tip_window = None

    def _hide(self) -> None:
        if self._tip_window:
            try:
                self._tip_window.destroy()
            except Exception:
                pass
            self._tip_window = None
