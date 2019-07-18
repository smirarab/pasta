#!/bin/bash

set -u
set -e
set -x
set -o pipefail

rm -fr build dist

# python setup.py py2app --resources run_pasta.py,bin,run_seqtools.py --iconfile pasta.ico 2>&1 |tee py2app.log; 

pyinstaller --windowed --add-binary bin/:. run_pasta.spec 2>&1 |tee pyinstaller.log

#chmod +x dist/PASTA.app/Contents/Resources/bin/*;
