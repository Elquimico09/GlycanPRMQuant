# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files


project_root = Path(SPEC).resolve().parent
datas = [
    (str(project_root / "glycanPRMQuant" / "database"), "glycanPRMQuant/database"),
]
datas += collect_data_files("scienceplots")
datas += collect_data_files("alpharaw", includes=["ext/thermo_fisher/**/*"])


a = Analysis(
    [str(project_root / "glycanPRMQuant" / "pipelineGUI.py")],
    pathex=[str(project_root)],
    binaries=[],
    datas=datas,
    hiddenimports=[
        "alpharaw.thermo",
        "alpharaw.raw_access.pythermorawfilereader",
        "clr",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "IPython",
        "PyQt5",
        "PyQt6",
        "PySide2",
        "PySide6",
        "conda",
        "jupyter",
        "notebook",
        "pyarrow",
        "sqlalchemy",
        "tensorflow",
        "torch",
        "transformers",
    ],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='GlycanPRMQuant',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='GlycanPRMQuant',
)
