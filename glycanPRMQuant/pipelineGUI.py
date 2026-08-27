import multiprocessing
import threading
import os
import json
import time
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import scrolledtext
from tkinter import ttk
from glycanPRMQuant.parallelProcess import run_parallel_pipeline
from glycanPRMQuant.spectra import validate_input_file_types


class PipelineGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("glycanPRMQuant Parallel Pipeline")
        self.geometry("980x1080")
        self.configure(bg="#0f172a")
        self.pipeline_proc = None
        self.queue_manager = None
        self.log_queue = None
        self.progress_queue = None
        self.selected_files = []
        self._polling = False
        self._process_done = False
        self._log_done = False
        self._progress_done = False
        self._run_started_at = None
        self._run_requested_files = 0
        self._summary_written = False
        self._prefs_path = os.path.join(os.getcwd(), ".pipeline_gui_prefs.json")
        self._init_style()
        self._build_widgets()
        self._load_prefs()

    def _init_style(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        self._colors = {
            "bg": "#0f172a",
            "panel": "#111827",
            "panel2": "#0b1220",
            "text": "#e5e7eb",
            "muted": "#9ca3af",
            "accent": "#14b8a6",
            "danger": "#ef4444",
            "outline": "#1f2937",
            "input": "#0b1020"
        }

        self.option_add("*Font", ("Segoe UI", 10))
        style.configure("TLabel", background=self._colors["panel"], foreground=self._colors["text"])
        style.configure("Muted.TLabel", background=self._colors["panel"], foreground=self._colors["muted"])
        style.configure("Title.TLabel", background=self._colors["bg"], foreground=self._colors["text"],
                        font=("Segoe UI Semibold", 16))
        style.configure("Section.TLabel", background=self._colors["panel"], foreground=self._colors["text"],
                        font=("Segoe UI Semibold", 11))
        style.configure("TFrame", background=self._colors["panel"])
        style.configure("Card.TFrame", background=self._colors["panel"], relief="flat")
        style.configure("TButton", padding=(10, 6))
        style.configure("Primary.TButton", background=self._colors["accent"], foreground="white")
        style.map("Primary.TButton",
                  background=[("active", "#0ea5a0"), ("disabled", "#1f2937")],
                  foreground=[("disabled", "#9ca3af")])
        style.configure("Danger.TButton", background=self._colors["danger"], foreground="white")
        style.map("Danger.TButton",
                  background=[("active", "#dc2626"), ("disabled", "#1f2937")],
                  foreground=[("disabled", "#9ca3af")])
        style.configure("TCheckbutton", background=self._colors["panel"], foreground=self._colors["text"])
        style.map("TCheckbutton",
                  background=[("active", self._colors["panel"])],
                  foreground=[("active", self._colors["text"])])
        style.configure("TEntry", fieldbackground=self._colors["input"], foreground=self._colors["text"])
        style.configure("TCombobox", fieldbackground=self._colors["input"], foreground=self._colors["text"])
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", self._colors["input"])],
            foreground=[("readonly", self._colors["text"])],
            selectbackground=[("readonly", self._colors["input"])],
            selectforeground=[("readonly", self._colors["text"])]
        )
        style.configure("TProgressbar", troughcolor=self._colors["panel2"], background=self._colors["accent"])

    def _build_widgets(self):
        # Scrollable container for all controls
        container = tk.Frame(self, bg=self._colors["bg"])
        container.pack(fill="both", expand=True)
        canvas = tk.Canvas(container, borderwidth=0, bg=self._colors["bg"], highlightthickness=0)
        vscroll = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)
        vscroll.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        self.inner = tk.Frame(canvas, bg=self._colors["bg"])
        self.inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=self.inner, anchor="nw")

        row = 0

        header = ttk.Frame(self.inner, style="TFrame")
        header.grid(column=0, row=row, columnspan=4, sticky="ew", padx=18, pady=(16, 6))
        ttk.Label(header, text="glycanPRMQuant Pipeline", style="Title.TLabel").grid(column=0, row=0, sticky="w")
        ttk.Label(header, text="Parallel processing for PRM glycan quantification", style="Muted.TLabel") \
            .grid(column=0, row=1, sticky="w")
        row += 1

        def section(title):
            nonlocal row
            frame = ttk.Frame(self.inner, style="Card.TFrame")
            frame.grid(column=0, row=row, columnspan=4, sticky="ew", padx=18, pady=8)
            ttk.Label(frame, text=title, style="Section.TLabel") \
                .grid(column=0, row=0, columnspan=4, sticky="w", padx=12, pady=(10, 6))
            row += 1
            return frame

        def add_label_entry(name, default, row, col, parent):
            ttk.Label(parent, text=name).grid(column=col, row=row, sticky="w", padx=8, pady=6)
            var = tk.StringVar(value=str(default))
            ttk.Entry(parent, textvariable=var, width=14).grid(column=col + 1, row=row, padx=8, pady=6, sticky="w")
            return var

        # Input / Output
        io_frame = section("Input & Output")
        tk.Label(io_frame, text="Input RAW or mzML files:", bg=self._colors["panel"], fg=self._colors["text"]) \
            .grid(column=0, row=1, sticky="w", padx=8, pady=6)
        self.files_label = tk.StringVar(value="No files selected")
        tk.Label(io_frame, textvariable=self.files_label, bg=self._colors["panel"], fg=self._colors["muted"],
                 wraplength=420, justify="left").grid(column=1, row=1, columnspan=2, padx=8, pady=6, sticky="w")
        ttk.Button(io_frame, text="Browse", command=self._choose_files).grid(column=3, row=1, padx=8, pady=6)

        tk.Label(io_frame, text="Output folder:", bg=self._colors["panel"], fg=self._colors["text"]) \
            .grid(column=0, row=2, sticky="w", padx=8, pady=6)
        self.out_dir = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.out_dir, width=40).grid(column=1, row=2, columnspan=2, padx=8, pady=6, sticky="w")
        ttk.Button(io_frame, text="Browse", command=self._choose_output).grid(column=3, row=2, padx=8, pady=6)

        tk.Label(io_frame, text="Precursor DB:", bg=self._colors["panel"], fg=self._colors["text"]) \
            .grid(column=0, row=3, sticky="w", padx=8, pady=6)
        self.precursor_db_path = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.precursor_db_path, width=40) \
            .grid(column=1, row=3, columnspan=2, padx=8, pady=6, sticky="w")
        ttk.Button(io_frame, text="Browse", command=self._choose_precursor_db).grid(column=3, row=3, padx=8, pady=6)

        tk.Label(io_frame, text="Structure DB:", bg=self._colors["panel"], fg=self._colors["text"]) \
            .grid(column=0, row=4, sticky="w", padx=8, pady=6)
        self.structure_db_path = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.structure_db_path, width=40) \
            .grid(column=1, row=4, columnspan=2, padx=8, pady=6, sticky="w")
        ttk.Button(io_frame, text="Browse", command=self._choose_structure_db).grid(column=3, row=4, padx=8, pady=6)
        row += 1

        # Numeric parameters arranged in 3 columns
        params_frame = section("Processing Parameters")
        params = [
            ("Workers", 2),
            ("MS1 ppm tol", 10),
            ("MS1 m/z min", 400),
            ("MS1 m/z max", 2000),
            ("MS2 intensity", 1e2),
            ("MS2 ppm tol", 10),
            ("MS2 m/z tol", 0.02),
            ("Max cleavages", 2),
            ("Smoothing window", 5),
            ("Mass offset", 0.0),
            ("m/z offset", 0.0),
            ("AUC rel. height", 0.7),
        ]
        self._param_vars = {}
        base_row = 1
        for i, (label, default) in enumerate(params):
            r = base_row + i // 3
            c = (i % 3) * 2
            var = add_label_entry(label, default, r, c, params_frame)
            self._param_vars[label] = var
        row += 1

        # Bind param vars to existing attributes
        self.n_workers = self._param_vars["Workers"]
        self.ppm_ms1_tol = self._param_vars["MS1 ppm tol"]
        self.mz_min = self._param_vars["MS1 m/z min"]
        self.mz_max = self._param_vars["MS1 m/z max"]
        self.intensity_threshold = self._param_vars["MS2 intensity"]
        self.ppm_ms2_tol = self._param_vars["MS2 ppm tol"]
        self.mz_tol = self._param_vars["MS2 m/z tol"]
        self.fragment_max_cleavages = self._param_vars["Max cleavages"]
        self.smoothing_window = self._param_vars["Smoothing window"]
        self.mass_offset = self._param_vars["Mass offset"]
        self.mz_offset = self._param_vars["m/z offset"]
        self.rel_height = self._param_vars["AUC rel. height"]

        # Smoothing method selector
        ttk.Label(params_frame, text="Smoothing method").grid(column=0, row=base_row + 4, sticky="w", padx=8, pady=6)
        self.smoothing_method = tk.StringVar(value="gaussian")
        ttk.Combobox(
            params_frame,
            textvariable=self.smoothing_method,
            values=["gaussian", "savgol"],
            state="readonly",
            width=12
        ).grid(column=1, row=base_row + 4, padx=8, pady=6, sticky="w")

        # AUC rel-height mode selector
        ttk.Label(params_frame, text="AUC rel height mode").grid(column=2, row=base_row + 4, sticky="w", padx=8, pady=6)
        self.rel_height_mode = tk.StringVar(value="prominence")
        ttk.Combobox(
            params_frame,
            textvariable=self.rel_height_mode,
            values=["prominence", "height"],
            state="readonly",
            width=12
        ).grid(column=3, row=base_row + 4, padx=8, pady=6, sticky="w")

        ttk.Label(params_frame, text="Fragment ion series").grid(column=4, row=base_row + 4, sticky="w", padx=8, pady=6)
        self.fragment_ion_series = tk.StringVar(value="ABCXYZ")
        ttk.Entry(params_frame, textvariable=self.fragment_ion_series, width=14) \
            .grid(column=5, row=base_row + 4, padx=8, pady=6, sticky="w")
        row += 1

        # Toggles
        toggles_frame = section("Outputs & Options")
        self.overwrite_var = tk.BooleanVar(value=False)
        self.dryrun_var = tk.BooleanVar(value=False)
        self.adduct_plot_var = tk.BooleanVar(value=True)
        self.total_plot_var = tk.BooleanVar(value=True)
        self.skyline_var = tk.BooleanVar(value=False)
        self.smoothing_var = tk.BooleanVar(value=True)
        self.isobaric_resolution_var = tk.BooleanVar(value=True)

        ttk.Checkbutton(toggles_frame, text="Overwrite existing outputs", variable=self.overwrite_var) \
            .grid(column=0, row=1, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Dry run (plan only, no processing)", variable=self.dryrun_var) \
            .grid(column=0, row=2, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Generate precursor-adduct chromatograms", variable=self.adduct_plot_var) \
            .grid(column=0, row=3, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Generate total chromatograms", variable=self.total_plot_var) \
            .grid(column=0, row=4, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Skyline Transition list", variable=self.skyline_var) \
            .grid(column=0, row=5, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Enable smoothing", variable=self.smoothing_var) \
            .grid(column=0, row=6, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(
            toggles_frame,
            text="Resolve isobaric precursor conflicts",
            variable=self.isobaric_resolution_var
        ) \
            .grid(column=0, row=7, columnspan=4, sticky="w", padx=8, pady=4)
        row += 1

        # Run and Stop buttons
        action_frame = section("Run")
        self.run_btn = ttk.Button(action_frame, text="Run Pipeline", command=self._on_run, style="Primary.TButton")
        self.run_btn.grid(column=0, row=1, columnspan=4, pady=8, padx=8, sticky="ew")
        self.stop_btn = ttk.Button(action_frame, text="Stop Run", command=self._on_stop,
                                   style="Danger.TButton", state="disabled")
        self.stop_btn.grid(column=0, row=2, columnspan=4, pady=(0, 8), padx=8, sticky="ew")
        row += 1

        # Live log
        log_frame = section("Pipeline Log")
        ttk.Label(log_frame, text="Live output", style="Muted.TLabel").grid(column=0, row=1, sticky="w", padx=8)
        self.log_box = scrolledtext.ScrolledText(
            log_frame,
            height=12,
            width=96,
            state="disabled",
            bg=self._colors["panel2"],
            fg=self._colors["text"],
            insertbackground=self._colors["text"]
        )
        self.log_box.grid(column=0, row=2, columnspan=4, padx=8, pady=6, sticky="nsew")
        row += 1

        # Progress bar
        progress_frame = section("Progress")
        ttk.Label(progress_frame, text="Status", style="Muted.TLabel").grid(column=0, row=1, sticky="w", padx=8)
        self.progress = ttk.Progressbar(progress_frame, orient="horizontal", mode="determinate", length=560)
        self.progress.grid(column=1, row=1, columnspan=3, sticky="ew", padx=8, pady=6)
        row += 1

        # Status summary and open-output
        self.status_var = tk.StringVar(value="Idle")
        status_frame = section("Output")
        ttk.Label(status_frame, textvariable=self.status_var, anchor="w") \
            .grid(column=0, row=1, columnspan=3, sticky="w", padx=8, pady=6)
        ttk.Button(status_frame, text="Open output folder", command=self._open_output) \
            .grid(column=3, row=1, padx=8, pady=6, sticky="e")

    def _choose_files(self):
        files = filedialog.askopenfilenames(
            filetypes=[
                ("Mass spectrometry files", "*.raw *.RAW *.mzML *.mzml"),
                ("Thermo RAW files", "*.raw *.RAW"),
                ("mzML files", "*.mzML *.mzml"),
            ]
        )
        if files:
            self.selected_files = list(files)
            count = len(self.selected_files)
            self.files_label.set(f"{count} file{'s' if count != 1 else ''} selected")
        else:
            self.selected_files = []
            self.files_label.set("No files selected")

    def _choose_output(self):
        d = filedialog.askdirectory()
        if d:
            self.out_dir.set(d)

    def _choose_precursor_db(self):
        path = filedialog.askopenfilename(
            filetypes=[("Database files", "*.csv *.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.precursor_db_path.set(path)

    def _choose_structure_db(self):
        path = filedialog.askopenfilename(
            filetypes=[("Database files", "*.csv *.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            self.structure_db_path.set(path)

    def _on_run(self):
        # gather parameters
        try:
            params = {
                "input_files": self.selected_files,
                "input_dir": None,
                "output_root": self.out_dir.get(),
                "precursor_db_path": self.precursor_db_path.get().strip() or None,
                "structure_db_path": self.structure_db_path.get().strip() or None,
                "n_workers": int(self.n_workers.get()),
                "ppm_ms1_tol": float(self.ppm_ms1_tol.get()),
                "mz_min": float(self.mz_min.get()),
                "mz_max": float(self.mz_max.get()),
                "intensity_threshold": float(self.intensity_threshold.get()),
                "ppm_ms2_tol": float(self.ppm_ms2_tol.get()),
                "mz_tol": float(self.mz_tol.get()),
                "fragment_ion_series": self.fragment_ion_series.get(),
                "fragment_max_cleavages": int(self.fragment_max_cleavages.get()),
                "smoothing_window": int(self.smoothing_window.get()),
                "smoothing_method": self.smoothing_method.get(),
                "mass_offset": float(self.mass_offset.get()),
                "mz_offset": float(self.mz_offset.get()),
                "rel_height": float(self.rel_height.get()),
                "rel_height_mode": self.rel_height_mode.get(),
                "overwrite": bool(self.overwrite_var.get()),
                "dry_run": bool(self.dryrun_var.get()),
                "enable_adduct_plots": bool(self.adduct_plot_var.get()),
                "enable_total_plots": bool(self.total_plot_var.get()),
                "skyline_transition": bool(self.skyline_var.get()),
                "enable_smoothing": bool(self.smoothing_var.get()),
                "resolve_isobaric_conflicts": bool(self.isobaric_resolution_var.get()),
            }
        except ValueError as e:
            messagebox.showerror("Parameter error", f"Invalid number: {e}")
            return

        # Basic validation
        missing = [f for f in params["input_files"] if not os.path.isfile(f)]
        if not params["input_files"] or missing:
            messagebox.showerror(
                "Input error",
                "No RAW or mzML files selected" if not params["input_files"] else f"Missing files: {missing[:3]}",
            )
            return
        try:
            validate_input_file_types(params["input_files"])
        except ValueError as exc:
            messagebox.showerror("Input error", str(exc))
            return
        if not os.path.isdir(params["output_root"]):
            messagebox.showerror("Output error", "Output directory does not exist")
            return
        for label, key in [("Precursor DB", "precursor_db_path"), ("Structure DB", "structure_db_path")]:
            if params[key] and not os.path.isfile(params[key]):
                messagebox.showerror("Input error", f"{label} does not exist: {params[key]}")
                return
        if params["mz_min"] >= params["mz_max"]:
            messagebox.showerror("Parameter error", "MS1 m/z min must be < m/z max")
            return
        if params["n_workers"] <= 0:
            params["n_workers"] = os.cpu_count() or 1
        if params["n_workers"] > 61:
            params["n_workers"] = 61
        if params["fragment_max_cleavages"] < 1:
            messagebox.showerror("Parameter error", "Max cleavages must be >= 1")
            return
        self.status_var.set(f"Running... files={len(params['input_files'])}, workers={params['n_workers']}")

        self._save_prefs()

        # disable run, enable stop
        self.run_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
        self._clear_log()
        self.queue_manager = multiprocessing.Manager()
        self.log_queue = self.queue_manager.Queue()
        self.progress_queue = self.queue_manager.Queue()
        self._polling = True
        self._process_done = False
        self._log_done = False
        self._progress_done = False
        self._run_started_at = time.monotonic()
        self._run_requested_files = len(params["input_files"])
        self._summary_written = False
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
            self.after(0, self._mark_process_done)

    def _mark_process_done(self):
        self._process_done = True
        self._poll_log()

    def _finish_run(self):
        if not self._polling and self.pipeline_proc is None:
            return
        self._drain_queues()
        self._append_completion_summary()
        self.pipeline_proc = None
        if self.queue_manager:
            self.queue_manager.shutdown()
            self.queue_manager = None
        self.log_queue = None
        self.progress_queue = None
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self._polling = False
        self._process_done = False
        self.status_var.set("Done")
        messagebox.showinfo("Done", "Pipeline run finished or stopped.")

    def _poll_log(self):
        if not self._polling and self.pipeline_proc is None:
            return
        self._drain_queues()

        if self._process_done:
            exitcode = self.pipeline_proc.exitcode if self.pipeline_proc else 0
            if self._log_done or exitcode not in (0, None):
                self._finish_run()
                return

        if self._polling:
            self.after(200, self._poll_log)

    def _drain_queues(self):
        if self.log_queue:
            try:
                while True:
                    msg = self.log_queue.get_nowait()
                    if msg is None:
                        self._log_done = True
                        break
                    self._append_log(msg)
            except Exception:
                pass

        if self.progress_queue:
            try:
                while True:
                    item = self.progress_queue.get_nowait()
                    if item is None:
                        self._progress_done = True
                        break
                    self._update_progress(item)
            except Exception:
                pass

    def _append_log(self, msg):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg)
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def _clear_log(self):
        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", "end")
        self.log_box.configure(state="disabled")

    def _format_duration(self, seconds: float) -> str:
        total_seconds = max(int(round(seconds)), 0)
        minutes, seconds = divmod(total_seconds, 60)
        return f"{minutes} minutes {seconds} seconds"

    def _append_completion_summary(self):
        if self._summary_written or self._run_started_at is None:
            return
        elapsed = time.monotonic() - self._run_started_at
        processed = (
            self._progress_done
            + self._progress_skipped
            + self._progress_error
            + self._progress_dryrun
        )
        total = self._run_requested_files or self._progress_total
        selected_note = "" if processed == total else f" ({total} selected)"
        self._append_log(
            f"Processed {processed} files in {self._format_duration(elapsed)}{selected_note}.\n"
        )
        self._summary_written = True

    def _init_counters(self, total: int):
        self._progress_total = max(int(total), 1)
        self._progress_done = 0
        self._progress_skipped = 0
        self._progress_error = 0
        self._progress_dryrun = 0
        self.progress["value"] = 0

    def _update_progress(self, item):
        base, status, msg = item
        if status == "done":
            self._progress_done += 1
        elif status == "skipped":
            self._progress_skipped += 1
        elif status == "error":
            self._progress_error += 1
        elif status == "dry-run":
            self._progress_dryrun += 1
        self.progress["value"] = 100.0 * (
            self._progress_done
            + self._progress_skipped
            + self._progress_error
            + self._progress_dryrun
        ) / self._progress_total
        self.status_var.set(
            f"Done={self._progress_done}, Skipped={self._progress_skipped}, "
            f"Dry-run={self._progress_dryrun}, Errors={self._progress_error}, "
            f"Total={self._progress_total}"
        )

    def _open_output(self):
        out_dir = self.out_dir.get()
        if not out_dir or not os.path.isdir(out_dir):
            messagebox.showerror("Output error", "Output directory does not exist")
            return
        try:
            os.startfile(out_dir)
        except Exception:
            messagebox.showerror("Open error", "Failed to open output folder.")

    def _save_prefs(self):
        prefs = {
            "output_root": self.out_dir.get(),
            "precursor_db_path": self.precursor_db_path.get(),
            "structure_db_path": self.structure_db_path.get(),
            "n_workers": self.n_workers.get(),
            "ppm_ms1_tol": self.ppm_ms1_tol.get(),
            "mz_min": self.mz_min.get(),
            "mz_max": self.mz_max.get(),
            "intensity_threshold": self.intensity_threshold.get(),
            "ppm_ms2_tol": self.ppm_ms2_tol.get(),
            "mz_tol": self.mz_tol.get(),
            "fragment_ion_series": self.fragment_ion_series.get(),
            "fragment_max_cleavages": self.fragment_max_cleavages.get(),
            "smoothing_window": self.smoothing_window.get(),
            "smoothing_method": self.smoothing_method.get(),
            "mass_offset": self.mass_offset.get(),
            "mz_offset": self.mz_offset.get(),
            "rel_height": self.rel_height.get(),
            "rel_height_mode": self.rel_height_mode.get(),
            "overwrite": self.overwrite_var.get(),
            "dry_run": self.dryrun_var.get(),
            "enable_adduct_plots": self.adduct_plot_var.get(),
            "enable_total_plots": self.total_plot_var.get(),
            "skyline_transition": self.skyline_var.get(),
            "enable_smoothing": self.smoothing_var.get(),
            "resolve_isobaric_conflicts": self.isobaric_resolution_var.get(),
        }
        try:
            with open(self._prefs_path, "w", encoding="utf-8") as f:
                json.dump(prefs, f, indent=2)
        except Exception:
            pass

    def _load_prefs(self):
        if not os.path.isfile(self._prefs_path):
            return
        try:
            with open(self._prefs_path, "r", encoding="utf-8") as f:
                prefs = json.load(f)
        except Exception:
            return
        if "output_root" in prefs:
            self.out_dir.set(prefs["output_root"])
        if "precursor_db_path" in prefs:
            self.precursor_db_path.set(prefs["precursor_db_path"])
        if "structure_db_path" in prefs:
            self.structure_db_path.set(prefs["structure_db_path"])
        for key, var in [
            ("n_workers", self.n_workers),
            ("ppm_ms1_tol", self.ppm_ms1_tol),
            ("mz_min", self.mz_min),
            ("mz_max", self.mz_max),
            ("intensity_threshold", self.intensity_threshold),
            ("ppm_ms2_tol", self.ppm_ms2_tol),
            ("mz_tol", self.mz_tol),
            ("fragment_ion_series", self.fragment_ion_series),
            ("fragment_max_cleavages", self.fragment_max_cleavages),
            ("smoothing_window", self.smoothing_window),
            ("smoothing_method", self.smoothing_method),
            ("mass_offset", self.mass_offset),
            ("mz_offset", self.mz_offset),
            ("rel_height", self.rel_height),
            ("rel_height_mode", self.rel_height_mode),
        ]:
            if key in prefs:
                var.set(prefs[key])
        if "overwrite" in prefs:
            self.overwrite_var.set(prefs["overwrite"])
        if "dry_run" in prefs:
            self.dryrun_var.set(prefs["dry_run"])
        if "enable_adduct_plots" in prefs:
            self.adduct_plot_var.set(prefs["enable_adduct_plots"])
        if "enable_total_plots" in prefs:
            self.total_plot_var.set(prefs["enable_total_plots"])
        if "skyline_transition" in prefs:
            self.skyline_var.set(prefs["skyline_transition"])
        if "enable_smoothing" in prefs:
            self.smoothing_var.set(prefs["enable_smoothing"])
        if "resolve_isobaric_conflicts" in prefs:
            self.isobaric_resolution_var.set(prefs["resolve_isobaric_conflicts"])


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = PipelineGUI()
    app.mainloop()
