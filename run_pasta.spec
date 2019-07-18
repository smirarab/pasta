# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['run_pasta.py'],
             pathex=['/Users/smirarab/workspace/pasta'],
             binaries=[],
             datas=[('bin','bin'),('run_pasta_gui.py','.')],
             hiddenimports=['wx'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=True)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='run_pasta',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='run_pasta')
app = BUNDLE(coll,
             name='pasta.app',
             icon='pasta.icns',
             bundle_identifier=None)
