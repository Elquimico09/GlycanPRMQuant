#!/usr/bin/env python3
"""
Unified Tkinter GUI for:

1.  classifyGlycan  – reads a consolidated AUC CSV, appends pos1…pos5, Class, Type
2.  plot_barplot    – plots mean ± SEM of glycan-class abundances

Both tools live in separate tabs (Classification  ▸  Barplot).

Dependencies
------------
pandas, numpy, matplotlib, scienceplots  (+ tkinter, bundled with CPython)

Adjust the two import lines below if your package paths differ.
"""

from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots                                       # noqa: F401

# --------------------------------------------------------------------------
#  IMPORT  ANALYTICAL ROUTINES
# --------------------------------------------------------------------------
from glycanPRMQuant.glycanClassification import classifyGlycan        # noqa: F401
from glycanPRMQuant.glycantypeBarplot import plot_barplot, plot_type_barplot             # noqa: F401
# --------------------------------------------------------------------------


# ----------------------------- classification helpers --------------------- #
def classify_glycan_file(in_path: Path, out_path: Path | None = None) -> Path:
    df = pd.read_csv(in_path)
    out_df = classifyGlycan(in_path) if isinstance(in_path, str) else classifyGlycan(str(in_path))
    if out_path is None:
        out_path = in_path
    out_df.to_csv(out_path, index=False)
    return out_path


# ----------------------------- GUI root ----------------------------------- #
class GlycanToolsGUI(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("Glycan Tools")
        self.geometry("580x290")
        self.resizable(False, False)

        nb = ttk.Notebook(self)
        nb.pack(expand=True, fill="both", padx=4, pady=4)

        self.classify_tab = ClassifyTab(nb)
        self.barplot_tab  = BarplotTab(nb)

        nb.add(self.classify_tab, text="Classification")
        nb.add(self.barplot_tab,  text="Barplot")


# ----------------------------- classification tab ------------------------- #
class ClassifyTab(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        self.csv_path   = tk.StringVar()
        self.status_msg = tk.StringVar()

        ttk.Label(self, text="Input consolidated CSV:").pack(anchor="w", padx=12, pady=(12, 0))

        entry = ttk.Entry(self, textvariable=self.csv_path, width=70)
        entry.pack(padx=12, fill="x")

        ttk.Button(self, text="Browse…", command=self.browse_csv).pack(pady=4)

        ttk.Button(self,
                   text="Run classification",
                   command=self.run_classification,
                   style="Accent.TButton").pack(pady=6)

        ttk.Separator(self, orient="horizontal").pack(fill="x", pady=4, padx=8)
        ttk.Label(self, textvariable=self.status_msg, foreground="green",
                  wraplength=520, justify="center").pack(pady=6)

    # ------------------------------------------------------------------ #
    def browse_csv(self):
        path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.csv_path.set(path)

    def run_classification(self):
        in_path = self.csv_path.get()
        if not in_path:
            messagebox.showwarning("Select file", "Please choose an input CSV first.")
            return

        out_path = filedialog.asksaveasfilename(
            title="Save classified CSV as…",
            initialfile=Path(in_path).stem + "_classified.csv",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )
        if not out_path:
            return

        try:
            written = classify_glycan_file(Path(in_path), Path(out_path))
            self.status_msg.set(f"Classification complete!\nSaved to:\n{written}")
        except Exception as exc:
            messagebox.showerror("Error", str(exc))


# ----------------------------- bar-plot tab -------------------------------- #
class BarplotTab(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        self.csv_path   = tk.StringVar()
        self.save_plot  = tk.BooleanVar(value=False)
        self.status_msg = tk.StringVar()
        self.by_var = tk.StringVar(value="Class")
        ttk.Radiobutton(self, text="Group by Class",  variable=self.by_var, value="Class").pack()
        ttk.Radiobutton(self, text="Group by Type",   variable=self.by_var, value="Type").pack()
        


        ttk.Label(self, text="Input consolidated CSV for barplot:").pack(anchor="w",
                                                                         padx=12, pady=(12, 0))
        entry = ttk.Entry(self, textvariable=self.csv_path, width=70)
        entry.pack(padx=12, fill="x")

        ttk.Button(self, text="Browse…", command=self.browse_csv).pack(pady=4)

        ttk.Checkbutton(self, variable=self.save_plot,
                        text="Save plot to file instead of displaying").pack(pady=(4, 2))

        ttk.Button(self,
                   text="Generate barplot",
                   command=self.run_plot,
                   style="Accent.TButton").pack(pady=6)

        ttk.Separator(self, orient="horizontal").pack(fill="x", pady=4, padx=8)
        ttk.Label(self, textvariable=self.status_msg, foreground="green",
                  wraplength=520, justify="center").pack(pady=6)

    # ------------------------------------------------------------------ #
    def browse_csv(self):
        path = filedialog.askopenfilename(
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if path:
            self.csv_path.set(path)

    def run_plot(self):
        csv = self.csv_path.get()
        if not csv:
            messagebox.showwarning("Select file", "Please choose an input CSV first.")
            return

        save_path = None
        if self.save_plot.get():
            save_path = filedialog.asksaveasfilename(
                title="Save barplot image as…",
                initialfile=Path(csv).stem + "_barplot.png",
                defaultextension=".png",
                filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("All files", "*.*")]
            )
            if not save_path:
                return

        try:
            if self.by_var.get() == "Class":
              plot_barplot(csv, save_path=save_path)
            else:
              plot_type_barplot(csv, save_path=save_path, figsize=(4.8, 4))
            msg = "Plot saved!" if save_path else "Plot displayed!"
            if save_path:
                msg += f"\n{save_path}"
            self.status_msg.set(msg)
        except Exception as exc:
            messagebox.showerror("Error", str(exc))


# ---------------------------- run it -------------------------------------- #
if __name__ == "__main__":
    GlycanToolsGUI().mainloop()
