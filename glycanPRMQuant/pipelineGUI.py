import multiprocessing
import threading
import os
import json
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import scrolledtext
from tkinter import ttk
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

class PipelineGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("glycanPRMQuant Parallel Pipeline")
        self.geometry("640x740")
        self.pipeline_proc = None
        self.log_queue = None
        self.progress_queue = None
        self.selected_files = []
        self._polling = False
        self._prefs_path = os.path.join(os.getcwd(), ".pipeline_gui_prefs.json")
        self._build_widgets()
        self._load_prefs()

    def _build_widgets(self):
        # Scrollable container for all controls
        container = tk.Frame(self)
        container.pack(fill="both", expand=True)
        canvas = tk.Canvas(container, borderwidth=0)
        vscroll = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)
        vscroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        self.inner = tk.Frame(canvas)
        self.inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.inner, anchor="nw")

        row = 0

        def add_label_entry(name, default, row, col):
            tk.Label(self.inner, text=name).grid(column=col, row=row, sticky="w", padx=5, pady=5)
            var = tk.StringVar(value=str(default))
            tk.Entry(self.inner, textvariable=var, width=12).grid(column=col + 1, row=row, padx=5, pady=5, sticky="w")
            return var

        # Input / Output
        tk.Label(self.inner, text="Input .mzML files:").grid(column=0, row=row, sticky="w", padx=5, pady=5)
        self.files_label = tk.StringVar(value="No files selected")
        tk.Label(self.inner, textvariable=self.files_label, wraplength=320, justify="left").grid(column=1, row=row, columnspan=2, padx=5, pady=5, sticky="w")
        tk.Button(self.inner, text="Browse", command=self._choose_files).grid(column=3, row=row, padx=5)
        row += 1

        tk.Label(self.inner, text="Output folder:").grid(column=0, row=row, sticky="w", padx=5, pady=5)
        self.out_dir = tk.StringVar()
        tk.Entry(self.inner, textvariable=self.out_dir, width=30).grid(column=1, row=row, columnspan=2, padx=5, pady=5, sticky="w")
        tk.Button(self.inner, text="Browse", command=self._choose_output).grid(column=3, row=row, padx=5)
        row += 1

        # Numeric parameters arranged in 3 columns
        params = [
            ("Workers", 2),
            ("MS1 ppm tol", 10),
            ("MS1 m/z min", 400),
            ("MS1 m/z max", 2000),
            ("MS2 intensity", 1e2),
            ("MS2 ppm tol", 10),
            ("MS2 m/z tol", 0.02),
            ("Smoothing window", 5),
            ("Mass offset", 0.0),
            ("m/z offset", 0.0),
            ("AUC rel. height", 0.7),
        ]
        self._param_vars = {}
        base_row = row
        for i, (label, default) in enumerate(params):
            r = base_row + i // 3
            c = (i % 3) * 2
            var = add_label_entry(label, default, r, c)
            self._param_vars[label] = var
        row = base_row + (len(params) + 2) // 3 + 1

        # Bind param vars to existing attributes
        self.n_workers = self._param_vars["Workers"]
        self.ppm_ms1_tol = self._param_vars["MS1 ppm tol"]
        self.mz_min = self._param_vars["MS1 m/z min"]
        self.mz_max = self._param_vars["MS1 m/z max"]
        self.intensity_threshold = self._param_vars["MS2 intensity"]
        self.ppm_ms2_tol = self._param_vars["MS2 ppm tol"]
        self.mz_tol = self._param_vars["MS2 m/z tol"]
        self.smoothing_window = self._param_vars["Smoothing window"]
        self.mass_offset = self._param_vars["Mass offset"]
        self.mz_offset = self._param_vars["m/z offset"]
        self.rel_height = self._param_vars["AUC rel. height"]

        # Toggles
        self.overwrite_var = tk.BooleanVar(value=False)
        self.dryrun_var = tk.BooleanVar(value=False)
        self.adduct_plot_var = tk.BooleanVar(value=True)
        self.total_plot_var = tk.BooleanVar(value=True)
        self.skyline_var = tk.BooleanVar(value=False)

        tk.Checkbutton(self.inner, text="Overwrite existing outputs", variable=self.overwrite_var).grid(column=0, row=row, columnspan=4, sticky="w", padx=5); row += 1
        tk.Checkbutton(self.inner, text="Dry run (plan only, no processing)", variable=self.dryrun_var).grid(column=0, row=row, columnspan=4, sticky="w", padx=5); row += 1
        tk.Checkbutton(self.inner, text="Generate precursor-adduct chromatograms", variable=self.adduct_plot_var).grid(column=0, row=row, columnspan=4, sticky="w", padx=5); row += 1
        tk.Checkbutton(self.inner, text="Generate total chromatograms", variable=self.total_plot_var).grid(column=0, row=row, columnspan=4, sticky="w", padx=5); row += 1
        tk.Checkbutton(self.inner, text="Skyline Transition list", variable=self.skyline_var).grid(column=0, row=row, columnspan=4, sticky="w", padx=5); row += 1

        # Run and Stop buttons
        self.run_btn = tk.Button(self.inner, text="Run Pipeline", command=self._on_run,
                                 bg="#4CAF50", fg="white", font=("Arial", 12, "bold"))
        self.run_btn.grid(column=0, row=row, columnspan=4, pady=10, sticky="ew")
        row += 1

        self.stop_btn = tk.Button(self.inner, text="Stop Run", command=self._on_stop,
                                  bg="#F44336", fg="white", font=("Arial", 12, "bold"),
                                  state="disabled")
        self.stop_btn.grid(column=0, row=row, columnspan=4, pady=10, sticky="ew")
        row += 1

        # Live log
        tk.Label(self.inner, text="Pipeline log:").grid(column=0, row=row, sticky="nw", padx=5)
        self.log_box = scrolledtext.ScrolledText(self.inner, height=12, width=80, state="disabled")
        self.log_box.grid(column=0, row=row + 1, columnspan=4, padx=5, pady=5, sticky="nsew")
        row += 2

        # Progress bar
        tk.Label(self.inner, text="Progress:").grid(column=0, row=row, sticky="w", padx=5)
        self.progress = ttk.Progressbar(self.inner, orient="horizontal", mode="determinate", length=450)
        self.progress.grid(column=1, row=row, columnspan=3, sticky="ew", padx=5, pady=5)
        row += 1

        # Status summary and open-output
        self.status_var = tk.StringVar(value="Idle")
        tk.Label(self.inner, textvariable=self.status_var, anchor="w").grid(column=0, row=row, columnspan=3, sticky="w", padx=5, pady=5)
        tk.Button(self.inner, text="Open output folder", command=self._open_output).grid(column=3, row=row, padx=5, pady=5, sticky="e")

    def _choose_files(self):
        files = filedialog.askopenfilenames(filetypes=[("mzML files","*.mzML *.mzml")])
        if files:
            self.selected_files = list(files)
            display = "\n".join(os.path.basename(f) for f in self.selected_files)
            self.files_label.set(display)
        else:
            self.selected_files = []
            self.files_label.set("No files selected")

    def _choose_output(self):
        d = filedialog.askdirectory()
        if d: self.out_dir.set(d)

    def _on_run(self):
        # gather parameters
        try:
            params = {
                "input_files":         self.selected_files,
                "input_dir":           None,
                "output_root":         self.out_dir.get(),
                "n_workers":           int(self.n_workers.get()),
                "ppm_ms1_tol":         float(self.ppm_ms1_tol.get()),
                "mz_min":              float(self.mz_min.get()),
                "mz_max":              float(self.mz_max.get()),
                "intensity_threshold": float(self.intensity_threshold.get()),
                "ppm_ms2_tol":         float(self.ppm_ms2_tol.get()),
                "mz_tol":              float(self.mz_tol.get()),
                "smoothing_window":    int(self.smoothing_window.get()),
                "mass_offset":         float(self.mass_offset.get()),
                "mz_offset":           float(self.mz_offset.get()),
                "rel_height":          float(self.rel_height.get()),
                "overwrite":           bool(self.overwrite_var.get()),
                "dry_run":             bool(self.dryrun_var.get()),
                "enable_adduct_plots": bool(self.adduct_plot_var.get()),
                "enable_total_plots":  bool(self.total_plot_var.get()),
                "skyline_transition":  bool(self.skyline_var.get()),
            }
        except ValueError as e:
            messagebox.showerror("Parameter error", f"Invalid number: {e}")
            return

        # Basic validation
        missing = [f for f in params["input_files"] if not os.path.isfile(f)]
        if not params["input_files"] or missing:
            messagebox.showerror("Input error", "No mzML files selected" if not params["input_files"] else f"Missing files: {missing[:3]}")
            return
        if not os.path.isdir(params["output_root"]):
            messagebox.showerror("Output error", "Output directory does not exist")
            return
        if params["mz_min"] >= params["mz_max"]:
            messagebox.showerror("Parameter error", "MS1 m/z min must be < m/z max")
            return
        if params["n_workers"] <= 0:
            params["n_workers"] = os.cpu_count() or 1
        if params["n_workers"] > 61:
            params["n_workers"] = 61
        self.status_var.set(f"Running… files={len(params['input_files'])}, workers={params['n_workers']}")

        self._save_prefs()

        # disable run, enable stop
        self.run_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
        self._clear_log()
        self.log_queue = multiprocessing.Queue()
        self.progress_queue = multiprocessing.Queue()
        self._polling = True
        self._init_counters(total=len(params["input_files"]))

        # launch pipeline in separate process with log queue
        self.pipeline_proc = multiprocessing.Process(
            target=run_parallel_pipeline,
            kwargs={**params, "log_queue": self.log_queue, "progress_queue": self.progress_queue}
        )
        self.pipeline_proc.start()

        # monitor process completion and log
        threading.Thread(target=self._monitor_proc, daemon=True).start()
        self.after(200, self._poll_log)

    def _on_stop(self):
        if self.pipeline_proc and self.pipeline_proc.is_alive():
            self.pipeline_proc.terminate()
            messagebox.showinfo("Stopped", "Pipeline process terminated.")
        self._finish_run()

    def _monitor_proc(self):
        if self.pipeline_proc:
            self.pipeline_proc.join()
            self.after(0, self._finish_run)

    def _finish_run(self):
        self.pipeline_proc = None
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self._polling = False
        self.status_var.set("Done")
        messagebox.showinfo("Done", "Pipeline run finished or stopped.")

    def _poll_log(self):
        if not self.log_queue:
            return
        try:
            while True:
                msg = self.log_queue.get_nowait()
                if msg is None:
                    return
                self._append_log(msg)
        except Exception:
            pass
        if self.pipeline_proc and self.pipeline_proc.is_alive():
            self.after(200, self._poll_log)

    def _poll_progress(self):
        if not self.progress_queue:
            return
        try:
            while True:
                item = self.progress_queue.get_nowait()
                if item is None:
                    return
                base, status, msg = item
                self._update_counters(status)
        except Exception:
            pass
        if self._polling:
            self.after(300, self._poll_progress)

    def _append_log(self, text):
        self.log_box.config(state="normal")
        self.log_box.insert(tk.END, text)
        self.log_box.see(tk.END)
        self.log_box.config(state="disabled")

    def _clear_log(self):
        self.log_box.config(state="normal")
        self.log_box.delete("1.0", tk.END)
        self.log_box.config(state="disabled")

    def _init_counters(self, total):
        self._counts = {"total": total, "done": 0, "skipped": 0, "error": 0, "dry-run": 0}
        self.progress["maximum"] = max(total, 1)
        self.progress["value"] = 0
        self._update_status_text()
        self._poll_progress()

    def _update_counters(self, status):
        if status in self._counts:
            self._counts[status] += 1
        # completed = done+skipped+error+dry-run
        completed = self._counts["done"] + self._counts["skipped"] + self._counts["error"] + self._counts["dry-run"]
        self.progress["value"] = completed
        self._update_status_text()

    def _update_status_text(self):
        c = getattr(self, "_counts", {"total":0,"done":0,"skipped":0,"error":0,"dry-run":0})
        self.status_var.set(
            f"Total {c['total']} | done {c['done']} | skipped {c['skipped']} | dry-run {c['dry-run']} | errors {c['error']}"
        )

    def _open_output(self):
        out_dir = self.out_dir.get()
        if out_dir and os.path.isdir(out_dir):
            try:
                os.startfile(out_dir)
            except Exception as e:
                messagebox.showerror("Open output", f"Cannot open folder: {e}")
        else:
            messagebox.showwarning("Open output", "No output folder set")

    def _load_prefs(self):
        try:
            if os.path.isfile(self._prefs_path):
                with open(self._prefs_path, "r") as fh:
                    prefs = json.load(fh)
                self.out_dir.set(prefs.get("out_dir",""))
                self.n_workers.set(str(prefs.get("n_workers",2)))
                self.ppm_ms1_tol.set(str(prefs.get("ppm_ms1_tol",10)))
                self.mz_min.set(str(prefs.get("mz_min",400)))
                self.mz_max.set(str(prefs.get("mz_max",2000)))
                self.intensity_threshold.set(str(prefs.get("intensity_threshold",1e2)))
                self.ppm_ms2_tol.set(str(prefs.get("ppm_ms2_tol",10)))
                self.mz_tol.set(str(prefs.get("mz_tol",0.02)))
                self.smoothing_window.set(str(prefs.get("smoothing_window",5)))
                self.mass_offset.set(str(prefs.get("mass_offset",0.0)))
                self.mz_offset.set(str(prefs.get("mz_offset",0.0)))
                self.rel_height.set(str(prefs.get("rel_height",0.7)))
                self.overwrite_var.set(prefs.get("overwrite", False))
                self.dryrun_var.set(prefs.get("dry_run", False))
                self.adduct_plot_var.set(prefs.get("enable_adduct_plots", True))
                self.total_plot_var.set(prefs.get("enable_total_plots", True))
                self.skyline_var.set(prefs.get("skyline_transition", False))
        except Exception:
            pass

    def _save_prefs(self):
        prefs = {
            "out_dir": self.out_dir.get(),
            "n_workers": self.n_workers.get(),
            "ppm_ms1_tol": self.ppm_ms1_tol.get(),
            "mz_min": self.mz_min.get(),
            "mz_max": self.mz_max.get(),
            "intensity_threshold": self.intensity_threshold.get(),
            "ppm_ms2_tol": self.ppm_ms2_tol.get(),
            "mz_tol": self.mz_tol.get(),
            "smoothing_window": self.smoothing_window.get(),
            "mass_offset": self.mass_offset.get(),
            "mz_offset": self.mz_offset.get(),
            "rel_height": self.rel_height.get(),
            "overwrite": self.overwrite_var.get(),
            "dry_run": self.dryrun_var.get(),
            "enable_adduct_plots": self.adduct_plot_var.get(),
            "enable_total_plots": self.total_plot_var.get(),
            "skyline_transition": self.skyline_var.get(),
        }
        try:
            with open(self._prefs_path, "w") as fh:
                json.dump(prefs, fh)
        except Exception:
            pass

if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = PipelineGUI()
    app.mainloop()
