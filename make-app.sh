#!/bin/bash

rm -r build dist
python setup.py py2app --resources run_pasta.py,bin,run_seqtools.py --iconfile pasta.ico 2>&1 |tee py2app.log; 
chmod +x dist/PASTA.app/Contents/Resources/bin/*;
