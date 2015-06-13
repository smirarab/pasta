#!/usr/bin/env python

#############################################################################
##  this file is part of pasta.
##  see "license.txt" for terms and conditions of usage.
#############################################################################


"""
Package setup and installation.
"""

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages

from datetime import datetime
import os
import platform
import sys
import pasta

script_name = 'run_pasta.py' 
gui_script_name = 'run_pasta_gui.py'

def compose_build_distribution_name(build_type):
    return "pasta%s-v%s-%s" % (build_type, pasta.PROGRAM_VERSION, datetime.now().strftime("%Y%b%d"))

param = {
    'name': pasta.PROGRAM_NAME,
    'version': pasta.PROGRAM_VERSION,
    'app': {gui_script_name},
    'description': pasta.PROGRAM_DESCRIPTION,
    'author': pasta.PROGRAM_AUTHOR,
    'author_email': ['pasta-users@googlegroups.com'],
    'url': pasta.PROGRAM_WEBSITE,
    'license': pasta.PROGRAM_LICENSE,
    'packages': find_packages(),
    'package_dir': {'pasta': 'pasta'},
    'test_suite': "pasta.test",
    'include_package_data': True,
    'install_requires': ['dendropy<=3.12'],
    'scripts' : [script_name,gui_script_name,'run_seqtools.py'],
    'zip_safe': True,
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': """A Python implementation of the Practical Alignment using SATe and Transitivity. 
    The package requires configuration to refer to third-party tools such as ClustalW2, MAFFT, MUCLE, OPAL, Prank, and RAxML, HMMMER,
    and the code is heavily based on SATe""",
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: GNU General Public License (GPL)",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    ],
    }

if sys.argv[1] == 'py2exe':
    PY2EXE_DIST_DIR = compose_build_distribution_name("win")
    if not platform.system() == 'Windows':
        raise ValueError('py2exe option only works on MS Windows.\n')
        from distutils.core import setup
    import glob
    import py2exe

    def find_data_files(source,target,patterns):
        if glob.has_magic(source) or glob.has_magic(target):
            raise ValueError("Magic not allowed in src, target")
        ret = {}
        for pattern in patterns:
            pattern = os.path.join(source,pattern)
            for filename in glob.glob(pattern):
                if os.path.isfile(filename):
                    targetpath = os.path.join(target,os.path.relpath(filename,source))
                    path = os.path.dirname(targetpath)
                    ret.setdefault(path,[]).append(filename)
        return sorted(ret.items())

    def extend_data_files(file_list, dirname, filenames):
        l = dirname.split(os.path.sep)
        target_list = l[l.index('data'):]
        target = os.path.join(*target_list)
        file_list.append((target,
                [os.path.join(dirname,
                        f) for f in filenames if not os.path.isdir(
                                os.path.join(dirname, f))]))

    def extend_output_files(file_list, dirname, filenames):
        l = dirname.split(os.path.sep)
        target_list = l[l.index('sample-output'):]
        target = os.path.join(*target_list)
        file_list.append((target,
                [os.path.join(dirname,
                        f) for f in filenames if not os.path.isdir(
                                os.path.join(dirname, f))]))

    bin_win_src = pasta.pasta_tools_dev_dir()
    bin_win_dest = pasta.pasta_tools_deploy_subpath()
    pasta_src_root = pasta.pasta_home_dir()
    data_dir = os.path.join(pasta_src_root, 'data')
    sample_output_dir = os.path.join(pasta_src_root, 'sample-output')
    my_files = []
    my_files.extend( find_data_files(
            bin_win_src,
            bin_win_dest,
            ['*'] ) )
    my_files.extend( find_data_files(
            os.path.join(bin_win_src, 'real_bin'),
            os.path.join(bin_win_dest, 'real_bin'),
            ['*'] ) )
    my_files.extend( find_data_files(
            os.path.join(pasta_src_root, 'doc'),
            'doc',
            ['*']))
    os.path.walk(data_dir, extend_data_files, my_files)
    os.path.walk(sample_output_dir, extend_output_files, my_files)

    PY2EXE_OPTIONS = {
        "unbuffered": True,
        "optimize": 2,
        "compressed": True,
        "bundle_files": 1,
        "includes": ['pasta'],
        "dll_excludes": ['w9xpopen.exe'],
        "dist_dir" : PY2EXE_DIST_DIR,
    }

    param.update({
        'console': [script_name, #os.path.join(pasta.PASTA_SCRIPT_RESOURCES, "mafft"),
                    os.path.join(pasta.PASTA_SCRIPT_RESOURCES, "hmmeralign")],
        'data_files': my_files,
        'zipfile': None,
        'options': {'py2exe': PY2EXE_OPTIONS},
        }
    )

### hack upon hack upon hack ...
if sys.argv[1] == 'py2exe':
    sys.stderr.write("\nMoving 'mafft.exe' into bundled binary directory ... \n")
    src_path = os.path.join(PY2EXE_DIST_DIR, "mafft.exe")
    dest_path = os.path.join(PY2EXE_DIST_DIR,
            bin_win_dest,
            "mafft.exe")
    if os.path.exists(dest_path):
        os.remove(dest_path)
    os.rename(src_path, dest_path)
    sys.stderr.write("OK\n")

# On Linux and OS X systems, sym-link all tool scripts
# to `bin` subdirectory, so PASTA can be run from the command-line
# I know this is ugly. Trust me, I hate it as much as you do.
if platform.system() != "Windows":

    DEST_DIR_ROOT = pasta.pasta_tools_deploy_dir(default_to_dev_dir=False)
    def create_symlink(src_path, subdir=None):
        if subdir:
            dest_dir = os.path.join(DEST_DIR_ROOT, subdir)
        else:
            dest_dir = DEST_DIR_ROOT
        dest_path = os.path.join(dest_dir, os.path.basename(src_path))
        sys.stderr.write("\nCreating link: '%s' => '%s'\n" % (src_path, dest_path))
        if os.path.exists(dest_path) and os.path.islink(dest_path):
            real_dest = os.path.abspath(os.path.realpath(dest_path))
            if real_dest != os.path.abspath(os.path.realpath(src_path)):
                msg = "ERROR: Symbolic link '%s' already exists, but points to different source: '%s'\n[Aborting]\nIf the old file was part of older PASTA versions, remove the old path manually and rerun." % (src_path, real_dest)
                sys.exit(msg)
            else:
                sys.stderr.write("Path already exists and is linked correctly.\n")
        elif os.path.exists(dest_path):
            msg = "ERROR: Path already exists: '%s'\n[Aborting]\n" % dest_path
            sys.exit(msg)
        else:
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            os.symlink(src_path, dest_path)

    # mafft
    #create_symlink(os.path.abspath(os.path.join(pasta.PASTA_SCRIPT_RESOURCES, "mafft")))
    create_symlink(os.path.abspath(os.path.join(pasta.PASTA_SCRIPT_RESOURCES, "hmmeralign")))

    # others
    tools_bin_srcdir = pasta.pasta_tools_dev_dir()
    tools_bin_subdirs = ['', 'mafftdir/bin','mafftdir/libexec']

    for subdir in tools_bin_subdirs:
        if subdir:
            tdir = os.path.join(tools_bin_srcdir, subdir)
        else:
            tdir = tools_bin_srcdir
        for fpath in os.listdir(tdir):
            src_path = os.path.join(tdir, fpath)
            if os.path.isfile(src_path) and not src_path.endswith('.txt'):
                create_symlink(src_path, subdir)
                


setup(**param)
