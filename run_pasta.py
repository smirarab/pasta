#! /usr/bin/env python

"""Main script of PASTA in command-line mode - this simply invokes the main
    function found in pasta/mainpasta.py
"""

# This file is part of PASTA which is forked from SATe

# PASTA like SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

if __name__ == "__main__":
    import os
    import sys
    from pasta.mainpasta import pasta_main
    from pasta import MESSENGER
    sys.setrecursionlimit(100000)
    _PASTA_DEBUG = os.environ.get('PASTA_DEBUG')
    _DEVELOPER = _PASTA_DEBUG and _PASTA_DEBUG != '0'

    if not _DEVELOPER:
        _PASTA_DEVELOPER = os.environ.get('PASTA_DEVELOPER')
        _DEVELOPER = _PASTA_DEVELOPER and _PASTA_DEVELOPER != '0'
    try:
        rc, temp_dir, temp_fs = pasta_main()
        if not rc:
            raise ValueError("Unknown PASTA execution error")
        if (temp_dir is not None) and (os.path.exists(temp_dir)):
            MESSENGER.send_info("Note that temporary files from the run have not been deleted, they can be found in:\n   '%s'\n" % temp_dir)
            if sys.platform.lower().startswith('darwin') and ("'" not in temp_dir):
                MESSENGER.send_info('''
If you cannot see this directory in the Finder application, you may want to use
the 'open' command executed from a Terminal.  You can do this by launching the
/Applications/Utilities/Terminal program and then typing

open '%s'

followed by a return at the prompt. If the argument to the open command is a
directory, then it should open a Finder window in the directory (even if that
directory is hidden by default).
''' % temp_dir)
    except Exception, x:
        if _DEVELOPER:
            raise
        message = "PASTA is exiting because of an error:\n%s " % str(x)
        try:
            from pasta import MESSENGER
            MESSENGER.send_error(message)
        except:
            sys.stderr.write(message)
        sys.exit(1)
