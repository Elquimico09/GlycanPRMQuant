import multiprocessing
import threading
import os
import json
import time
from datetime import datetime
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import scrolledtext
from tkinter import ttk
from glycanPRMQuant.database_utils import validate_glycan_database
from glycanPRMQuant.parallelProcess import run_parallel_pipeline
from glycanPRMQuant.spectra import validate_input_file_types


def _ensure_output_directory(path: str) -> str:
    """Create and return a normalized GUI output directory."""
    requested = os.fspath(path).strip()
    if not requested:
        raise ValueError("Enter or select an output directory")

    output_directory = os.path.abspath(os.path.expanduser(requested))
    if os.path.exists(output_directory) and not os.path.isdir(output_directory):
        raise NotADirectoryError(
            f"Output path is not a directory: {output_directory}"
        )
    os.makedirs(output_directory, exist_ok=True)
    return output_directory


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
        self._log_file_handle = None
        self._log_file_path = None
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
        style.configure("Warning.TLabel", background=self._colors["panel"], foreground="#f59e0b")
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

        tk.Label(io_frame, text="Custom glycan database (optional):", bg=self._colors["panel"], fg=self._colors["text"]) \
            .grid(column=0, row=3, sticky="w", padx=8, pady=6)
        self.database_path = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.database_path, width=40) \
            .grid(column=1, row=3, columnspan=2, padx=8, pady=6, sticky="w")
        ttk.Button(io_frame, text="Browse", command=self._choose_database).grid(column=3, row=3, padx=8, pady=6)
        ttk.Label(
            io_frame,
            text=(
                "Warning: Custom databases require nonblank Condensed IUPAC, Composition, "
                "and Numerical Composition columns. Leave blank to use the bundled database."
            ),
            style="Warning.TLabel",
            wraplength=650,
        ).grid(column=1, row=4, columnspan=3, padx=8, pady=(0, 6), sticky="w")
        row += 1

        # Keep related acquisition, fragment, smoothing, and AUC controls
        # together instead of filling the grid strictly by input type.
        params_frame = section("Processing Parameters")
        self._param_vars = {}

        # General processing controls.
        self.n_workers = add_label_entry("Workers", 2, 1, 0, params_frame)
        self.ppm_ms1_tol = add_label_entry("MS1 ppm tol", 10, 1, 2, params_frame)
        self.intensity_threshold = add_label_entry(
            "MS2 intensity", 1e2, 1, 4, params_frame
        )

        # Fragment-generation and fragment-matching controls.
        self.fragment_mass_tol = add_label_entry(
            "Fragment tolerance value", 0.02, 2, 0, params_frame
        )
        ttk.Label(params_frame, text="Fragment tolerance unit").grid(
            column=2, row=2, sticky="w", padx=8, pady=6
        )
        self.fragment_mass_tol_unit = tk.StringVar(value="Da")
        ttk.Combobox(
            params_frame,
            textvariable=self.fragment_mass_tol_unit,
            values=["Da", "ppm"],
            state="readonly",
            width=12,
        ).grid(column=3, row=2, padx=8, pady=6, sticky="w")

        ttk.Label(params_frame, text="Fragment ion series").grid(
            column=4, row=2, sticky="w", padx=8, pady=6
        )
        self.fragment_ion_series = tk.StringVar(value="ABCXYZ")
        ttk.Entry(
            params_frame, textvariable=self.fragment_ion_series, width=14
        ).grid(column=5, row=2, padx=8, pady=6, sticky="w")
        self.fragment_max_cleavages = add_label_entry(
            "Max cleavages", 2, 3, 4, params_frame
        )

        # Chromatogram smoothing controls.
        self.smoothing_window = add_label_entry(
            "Smoothing window", 5, 4, 0, params_frame
        )
        ttk.Label(params_frame, text="Smoothing method").grid(
            column=2, row=4, sticky="w", padx=8, pady=6
        )
        self.smoothing_method = tk.StringVar(value="gaussian")
        ttk.Combobox(
            params_frame,
            textvariable=self.smoothing_method,
            values=["gaussian", "savgol"],
            state="readonly",
            width=12
        ).grid(column=3, row=4, padx=8, pady=6, sticky="w")

        # AUC boundary controls.
        self.rel_height = add_label_entry(
            "AUC rel. height", 0.7, 5, 0, params_frame
        )
        ttk.Label(params_frame, text="AUC rel height mode").grid(
            column=2, row=5, sticky="w", padx=8, pady=6
        )
        self.rel_height_mode = tk.StringVar(value="prominence")
        ttk.Combobox(
            params_frame,
            textvariable=self.rel_height_mode,
            values=["prominence", "height"],
            state="readonly",
            width=12
        ).grid(column=3, row=5, padx=8, pady=6, sticky="w")

        self._param_vars.update(
            {
                "Workers": self.n_workers,
                "MS1 ppm tol": self.ppm_ms1_tol,
                "MS2 intensity": self.intensity_threshold,
                "Fragment tolerance value": self.fragment_mass_tol,
                "Max cleavages": self.fragment_max_cleavages,
                "Smoothing window": self.smoothing_window,
                "AUC rel. height": self.rel_height,
            }
        )
        row += 1

        scoring_frame = section("Candidate Scoring")
        scoring_params = [
            ("Minimum fragments", 2),
            ("Minimum explained intensity", 0.01),
            ("Minimum candidate score", 35.0),
            ("Minimum evidence difference", 4.0),
            ("Mass-outlier minimum Δppm", 2.0),
            ("Maximum assignment q-value", 0.05),
        ]
        self._scoring_vars = {}
        for i, (label, default) in enumerate(scoring_params):
            r = 1 + i // 3
            c = (i % 3) * 2
            self._scoring_vars[label] = add_label_entry(label, default, r, c, scoring_frame)
        self.candidate_min_fragments = self._scoring_vars["Minimum fragments"]
        self.candidate_min_explained_intensity = self._scoring_vars["Minimum explained intensity"]
        self.candidate_min_score = self._scoring_vars["Minimum candidate score"]
        self.candidate_min_evidence_difference = self._scoring_vars["Minimum evidence difference"]
        self.candidate_mass_outlier_min_delta = self._scoring_vars["Mass-outlier minimum Δppm"]
        self.candidate_max_q_value = self._scoring_vars["Maximum assignment q-value"]
        row += 1

        consensus_frame = section("Cross-run Peak Consensus")
        self.consensus_rt_tolerance = add_label_entry(
            "Peak ΔRT tolerance (± min)", 0.3, 1, 0, consensus_frame
        )
        self.consensus_min_replicate_fraction = add_label_entry(
            "Minimum run coverage", 0.8, 1, 2, consensus_frame
        )
        ttk.Label(
            consensus_frame,
            text=(
                "Peaks are matched when their aligned apex is within ± the specified "
                "ΔRT of the group center (for example, 0.5 means ±0.5 min). "
                "Applied across all files in this batch; alternatives remain audited."
            ),
            style="Muted.TLabel",
            wraplength=650,
        ).grid(column=0, row=2, columnspan=4, sticky="w", padx=8, pady=(0, 6))
        row += 1

        # Optional mass corrections
        optional_frame = section("Optional")
        self.mz_offset = add_label_entry("m/z offset", 0.0, 1, 0, optional_frame)
        self.mass_offset = add_label_entry("Mass offset", 0.0, 1, 2, optional_frame)
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
        self.target_decoy_var = tk.BooleanVar(value=True)
        self.consensus_peak_var = tk.BooleanVar(value=True)

        ttk.Label(toggles_frame, text="Figure file type").grid(
            column=0, row=1, sticky="w", padx=8, pady=6
        )
        self.figure_filetype = tk.StringVar(value="pdf")
        ttk.Combobox(
            toggles_frame,
            textvariable=self.figure_filetype,
            values=["png", "pdf", "svg"],
            state="readonly",
            width=12,
        ).grid(column=1, row=1, sticky="w", padx=8, pady=6)

        ttk.Checkbutton(toggles_frame, text="Overwrite existing outputs", variable=self.overwrite_var) \
            .grid(column=0, row=2, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Dry run (plan only, no processing)", variable=self.dryrun_var) \
            .grid(column=0, row=3, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Generate precursor-adduct chromatograms", variable=self.adduct_plot_var) \
            .grid(column=0, row=4, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Generate total chromatograms", variable=self.total_plot_var) \
            .grid(column=0, row=5, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Skyline Transition list", variable=self.skyline_var) \
            .grid(column=0, row=6, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(toggles_frame, text="Enable smoothing", variable=self.smoothing_var) \
            .grid(column=0, row=7, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(
            toggles_frame,
            text="Resolve isobaric precursor conflicts",
            variable=self.isobaric_resolution_var
        ) \
            .grid(column=0, row=8, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(
            toggles_frame,
            text="Enable target-decoy validation",
            variable=self.target_decoy_var,
        ).grid(column=0, row=9, columnspan=4, sticky="w", padx=8, pady=4)
        ttk.Checkbutton(
            toggles_frame,
            text="Enable cross-run consensus peak selection",
            variable=self.consensus_peak_var,
        ).grid(column=0, row=10, columnspan=4, sticky="w", padx=8, pady=4)
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

    def _choose_database(self):
        path = filedialog.askopenfilename(
            filetypes=[("Database files", "*.csv *.xlsx *.xls"), ("All files", "*.*")]
        )
        if path:
            try:
                validate_glycan_database(path)
            except (FileNotFoundError, ValueError, OSError, ImportError) as exc:
                messagebox.showerror("Database error", str(exc))
                return
            self.database_path.set(path)

    def _on_run(self):
        # gather parameters
        try:
            params = {
                "input_files": self.selected_files,
                "input_dir": None,
                "output_root": self.out_dir.get(),
                "database_path": self.database_path.get().strip() or None,
                "n_workers": int(self.n_workers.get()),
                "ppm_ms1_tol": float(self.ppm_ms1_tol.get()),
                "intensity_threshold": float(self.intensity_threshold.get()),
                "fragment_mass_tol": float(self.fragment_mass_tol.get()),
                "fragment_mass_tol_unit": self.fragment_mass_tol_unit.get(),
                "fragment_ion_series": self.fragment_ion_series.get(),
                "fragment_max_cleavages": int(self.fragment_max_cleavages.get()),
                "smoothing_window": int(self.smoothing_window.get()),
                "smoothing_method": self.smoothing_method.get(),
                "mass_offset": float(self.mass_offset.get()),
                "mz_offset": float(self.mz_offset.get()),
                "rel_height": float(self.rel_height.get()),
                "rel_height_mode": self.rel_height_mode.get(),
                "figure_filetype": self.figure_filetype.get(),
                "overwrite": bool(self.overwrite_var.get()),
                "dry_run": bool(self.dryrun_var.get()),
                "enable_adduct_plots": bool(self.adduct_plot_var.get()),
                "enable_total_plots": bool(self.total_plot_var.get()),
                "skyline_transition": bool(self.skyline_var.get()),
                "enable_smoothing": bool(self.smoothing_var.get()),
                "resolve_isobaric_conflicts": bool(self.isobaric_resolution_var.get()),
                "candidate_min_fragments": int(self.candidate_min_fragments.get()),
                "candidate_min_explained_intensity": float(self.candidate_min_explained_intensity.get()),
                "candidate_min_score": float(self.candidate_min_score.get()),
                "candidate_min_evidence_difference": float(self.candidate_min_evidence_difference.get()),
                "candidate_mass_outlier_min_delta": float(self.candidate_mass_outlier_min_delta.get()),
                "candidate_max_q_value": float(self.candidate_max_q_value.get()),
                "enable_target_decoy": bool(self.target_decoy_var.get()),
                "enable_consensus_peak_selection": bool(self.consensus_peak_var.get()),
                "consensus_rt_tolerance": float(self.consensus_rt_tolerance.get()),
                "consensus_min_replicate_fraction": float(
                    self.consensus_min_replicate_fraction.get()
                ),
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
        try:
            params["output_root"] = _ensure_output_directory(
                params["output_root"]
            )
            self.out_dir.set(params["output_root"])
        except (OSError, ValueError) as exc:
            messagebox.showerror(
                "Output error", f"Could not create the output directory:\n{exc}"
            )
            return
        if params["database_path"]:
            try:
                validate_glycan_database(params["database_path"])
            except (FileNotFoundError, ValueError, OSError, ImportError) as exc:
                messagebox.showerror("Database error", str(exc))
                return
        if params["n_workers"] <= 0:
            params["n_workers"] = os.cpu_count() or 1
        if params["n_workers"] > 61:
            params["n_workers"] = 61
        if params["fragment_max_cleavages"] < 1:
            messagebox.showerror("Parameter error", "Max cleavages must be >= 1")
            return
        if params["fragment_mass_tol"] <= 0:
            messagebox.showerror(
                "Parameter error", "Fragment tolerance value must be positive"
            )
            return
        if params["figure_filetype"] not in {"png", "pdf", "svg"}:
            messagebox.showerror(
                "Parameter error", "Figure file type must be png, pdf, or svg"
            )
            return
        if params["candidate_min_fragments"] < 1:
            messagebox.showerror("Parameter error", "Minimum candidate fragments must be >= 1")
            return
        if not 0 <= params["candidate_min_explained_intensity"] <= 1:
            messagebox.showerror("Parameter error", "Minimum explained intensity must be between 0 and 1")
            return
        if (
            params["candidate_min_score"] < 0
            or params["candidate_min_evidence_difference"] < 0
            or params["candidate_mass_outlier_min_delta"] < 0
        ):
            messagebox.showerror("Parameter error", "Candidate score thresholds cannot be negative")
            return
        if not 0 < params["candidate_max_q_value"] <= 1:
            messagebox.showerror(
                "Parameter error", "Maximum assignment q-value must be in (0, 1]"
            )
            return
        if params["consensus_rt_tolerance"] <= 0:
            messagebox.showerror(
                "Parameter error", "Consensus RT tolerance must be positive"
            )
            return
        if not 0 < params["consensus_min_replicate_fraction"] <= 1:
            messagebox.showerror(
                "Parameter error", "Minimum run coverage must be in (0, 1]"
            )
            return
        database_path = params.pop("database_path")
        params["precursor_db_path"] = database_path
        params["structure_db_path"] = database_path

        self._save_prefs()
        self._clear_log()
        try:
            self._start_log_file(params["output_root"])
        except OSError as exc:
            messagebox.showerror(
                "Log error",
                f"Could not create the GUI pipeline log in the output folder:\n{exc}",
            )
            return
        self._append_run_header(params)
        self.status_var.set(f"Running... files={len(params['input_files'])}, workers={params['n_workers']}")

        # disable run, enable stop
        self.run_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
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
        if self._log_file_path:
            self._append_log(f"GUI pipeline log saved to: {self._log_file_path}\n")
        log_path = self._log_file_path
        self._close_log_file()
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
        if log_path:
            self.status_var.set(f"Done — log saved to {log_path}")
            messagebox.showinfo(
                "Done", f"Pipeline run finished or stopped.\n\nLog saved to:\n{log_path}"
            )
        else:
            self.status_var.set("Done — log file could not be saved")
            messagebox.showwarning(
                "Done", "Pipeline run finished or stopped, but the log file could not be saved."
            )

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
        msg = str(msg)
        if self._log_file_handle is not None:
            try:
                self._log_file_handle.write(msg)
                self._log_file_handle.flush()
            except (OSError, ValueError) as exc:
                failed_path = self._log_file_path
                self._close_log_file()
                self._log_file_path = None
                msg += (
                    f"\nERROR: Failed to write GUI pipeline log {failed_path}: {exc}\n"
                )
        self.log_box.configure(state="normal")
        self.log_box.insert("end", msg)
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def _start_log_file(self, output_root):
        self._close_log_file()
        timestamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
        stem = f"glycanPRMQuant_pipeline_log_{timestamp}"
        log_path = os.path.join(output_root, f"{stem}.txt")
        suffix = 2
        while os.path.exists(log_path):
            log_path = os.path.join(output_root, f"{stem}_{suffix}.txt")
            suffix += 1
        self._log_file_handle = open(log_path, "w", encoding="utf-8", buffering=1)
        self._log_file_path = log_path

    def _close_log_file(self):
        if self._log_file_handle is None:
            return
        try:
            self._log_file_handle.close()
        finally:
            self._log_file_handle = None

    def _append_run_header(self, params):
        started = datetime.now().astimezone().isoformat(timespec="seconds")
        input_files = params.get("input_files", [])
        lines = [
            "glycanPRMQuant GUI pipeline log\n",
            f"Started: {started}\n",
            f"Log file: {self._log_file_path}\n",
            "Input files:\n",
        ]
        lines.extend(f"  {path}\n" for path in input_files)
        lines.append("Parameters:\n")
        for name in sorted(params):
            if name in {"input_files", "input_dir"}:
                continue
            lines.append(f"  {name}: {params[name]}\n")
        lines.append("\n")
        self._append_log("".join(lines))

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
            "database_path": self.database_path.get(),
            "n_workers": self.n_workers.get(),
            "ppm_ms1_tol": self.ppm_ms1_tol.get(),
            "intensity_threshold": self.intensity_threshold.get(),
            "fragment_mass_tol": self.fragment_mass_tol.get(),
            "fragment_mass_tol_unit": self.fragment_mass_tol_unit.get(),
            "fragment_ion_series": self.fragment_ion_series.get(),
            "fragment_max_cleavages": self.fragment_max_cleavages.get(),
            "smoothing_window": self.smoothing_window.get(),
            "smoothing_method": self.smoothing_method.get(),
            "mass_offset": self.mass_offset.get(),
            "mz_offset": self.mz_offset.get(),
            "rel_height": self.rel_height.get(),
            "rel_height_mode": self.rel_height_mode.get(),
            "figure_filetype": self.figure_filetype.get(),
            "overwrite": self.overwrite_var.get(),
            "dry_run": self.dryrun_var.get(),
            "enable_adduct_plots": self.adduct_plot_var.get(),
            "enable_total_plots": self.total_plot_var.get(),
            "skyline_transition": self.skyline_var.get(),
            "enable_smoothing": self.smoothing_var.get(),
            "resolve_isobaric_conflicts": self.isobaric_resolution_var.get(),
            "candidate_min_fragments": self.candidate_min_fragments.get(),
            "candidate_min_explained_intensity": self.candidate_min_explained_intensity.get(),
            "candidate_min_score": self.candidate_min_score.get(),
            "candidate_min_evidence_difference": self.candidate_min_evidence_difference.get(),
            "candidate_mass_outlier_min_delta": self.candidate_mass_outlier_min_delta.get(),
            "candidate_max_q_value": self.candidate_max_q_value.get(),
            "enable_target_decoy": self.target_decoy_var.get(),
            "enable_consensus_peak_selection": self.consensus_peak_var.get(),
            "consensus_rt_tolerance": self.consensus_rt_tolerance.get(),
            "consensus_min_replicate_fraction": self.consensus_min_replicate_fraction.get(),
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
        if "database_path" in prefs:
            self.database_path.set(prefs["database_path"])
        for key, var in [
            ("n_workers", self.n_workers),
            ("ppm_ms1_tol", self.ppm_ms1_tol),
            ("intensity_threshold", self.intensity_threshold),
            ("fragment_mass_tol", self.fragment_mass_tol),
            ("fragment_mass_tol_unit", self.fragment_mass_tol_unit),
            ("fragment_ion_series", self.fragment_ion_series),
            ("fragment_max_cleavages", self.fragment_max_cleavages),
            ("smoothing_window", self.smoothing_window),
            ("smoothing_method", self.smoothing_method),
            ("mass_offset", self.mass_offset),
            ("mz_offset", self.mz_offset),
            ("rel_height", self.rel_height),
            ("rel_height_mode", self.rel_height_mode),
            ("figure_filetype", self.figure_filetype),
            ("candidate_min_fragments", self.candidate_min_fragments),
            ("candidate_min_explained_intensity", self.candidate_min_explained_intensity),
            ("candidate_min_score", self.candidate_min_score),
            ("candidate_min_evidence_difference", self.candidate_min_evidence_difference),
            ("candidate_mass_outlier_min_delta", self.candidate_mass_outlier_min_delta),
            ("candidate_max_q_value", self.candidate_max_q_value),
            ("consensus_rt_tolerance", self.consensus_rt_tolerance),
            (
                "consensus_min_replicate_fraction",
                self.consensus_min_replicate_fraction,
            ),
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
        if "enable_target_decoy" in prefs:
            self.target_decoy_var.set(prefs["enable_target_decoy"])
        if "enable_consensus_peak_selection" in prefs:
            self.consensus_peak_var.set(prefs["enable_consensus_peak_selection"])


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = PipelineGUI()
    app.mainloop()
