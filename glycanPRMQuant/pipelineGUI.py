import multiprocessing
import threading
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from glycanPRMQuant.parallelProcess import run_parallel_pipeline

class PipelineGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("glycanPRMQuant Parallel Pipeline")
        self.geometry("500x650")
        self.pipeline_proc = None
        self._build_widgets()

    def _build_widgets(self):
        row = 0
        def add_label_entry(name, default, row):
            tk.Label(self, text=name).grid(column=0, row=row, sticky="w", padx=5, pady=5)
            var = tk.StringVar(value=str(default))
            tk.Entry(self, textvariable=var, width=20).grid(column=1, row=row, padx=5, pady=5)
            return var

        # Input / Output
        tk.Label(self, text="Input .mzML folder:").grid(column=0,row=row,sticky="w",padx=5,pady=5)
        self.in_dir = tk.StringVar()
        tk.Entry(self, textvariable=self.in_dir, width=30).grid(column=1,row=row,padx=5,pady=5)
        tk.Button(self, text="Browse", command=self._choose_input).grid(column=2,row=row,padx=5)
        row += 1

        tk.Label(self, text="Output folder:").grid(column=0,row=row,sticky="w",padx=5,pady=5)
        self.out_dir = tk.StringVar()
        tk.Entry(self, textvariable=self.out_dir, width=30).grid(column=1,row=row,padx=5,pady=5)
        tk.Button(self, text="Browse", command=self._choose_output).grid(column=2,row=row,padx=5)
        row += 1

        # Numeric parameters
        self.n_workers           = add_label_entry("Workers",            2,    row); row+=1
        self.ppm_ms1_tol         = add_label_entry("MS1 ppm tol",       10,   row); row+=1
        self.mz_min              = add_label_entry("MS1 m/z min",       400,  row); row+=1
        self.mz_max              = add_label_entry("MS1 m/z max",       2000, row); row+=1
        self.intensity_threshold = add_label_entry("MS2 intensity",     1e2,  row); row+=1
        self.ppm_ms2_tol         = add_label_entry("MS2 ppm tol",       10,   row); row+=1
        self.mz_tol              = add_label_entry("MS2 m/z tol",       0.02, row); row+=1
        self.smoothing_window    = add_label_entry("Smoothing window",  5,   row); row+=1
        self.mass_offset         = add_label_entry("Mass offset",       0.0,  row); row+=1
        self.mz_offset           = add_label_entry("m/z offset",        0.0,  row); row+=1

        # Run and Stop buttons
        self.run_btn = tk.Button(self, text="Run Pipeline", command=self._on_run,
                                 bg="#4CAF50", fg="white", font=("Arial",12,"bold"))
        self.run_btn.grid(column=0,row=row,columnspan=3,pady=10,sticky="ew")
        row += 1

        self.stop_btn = tk.Button(self, text="Stop Run", command=self._on_stop,
                                  bg="#F44336", fg="white", font=("Arial",12,"bold"),
                                  state="disabled")
        self.stop_btn.grid(column=0,row=row,columnspan=3,pady=10,sticky="ew")

    def _choose_input(self):
        d = filedialog.askdirectory()
        if d: self.in_dir.set(d)

    def _choose_output(self):
        d = filedialog.askdirectory()
        if d: self.out_dir.set(d)

    def _on_run(self):
        # gather parameters
        try:
            params = {
                "input_dir":           self.in_dir.get(),
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
            }
        except ValueError as e:
            messagebox.showerror("Parameter error", f"Invalid number: {e}")
            return

        if not os.path.isdir(params["input_dir"]):
            messagebox.showerror("Input error", "Input directory does not exist")
            return
        if not os.path.isdir(params["output_root"]):
            messagebox.showerror("Output error", "Output directory does not exist")
            return

        # disable run, enable stop
        self.run_btn.config(state="disabled")
        self.stop_btn.config(state="normal")

        # launch pipeline in separate process
        self.pipeline_proc = multiprocessing.Process(target=run_parallel_pipeline, kwargs=params)
        self.pipeline_proc.start()

        # monitor process completion
        threading.Thread(target=self._monitor_proc, daemon=True).start()

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
        messagebox.showinfo("Done", "Pipeline run finished or stopped.")

if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = PipelineGUI()
    app.mainloop()
