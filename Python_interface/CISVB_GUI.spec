# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['CISVB_GUI.py'],
    pathex=['/home/sourav/Desktop/SYMM_VB_SOFTWARE/Python_interface'],
    binaries=[],
    datas=[('/home/sourav/Desktop/logistic_regression/myenv/lib/python3.12/site-packages/Pmw', 'Pmw')],
    hiddenimports=['Pmw'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='CisvbGen',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
