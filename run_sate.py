#! /usr/bin/env python

"""Main script of SATe in command-line mode - this simply invokes the main
    function found in satelib/mainsate.py
"""

# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
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
    from sate.mainsate import sate_main
    from sate import MESSENGER
    _SATE_DEBUG = os.environ.get('SATE_DEBUG')
    _DEVELOPER = _SATE_DEBUG and _SATE_DEBUG != '0'

    if not _DEVELOPER:
        _SATE_DEVELOPER = os.environ.get('SATE_DEVELOPER')
        _DEVELOPER = _SATE_DEVELOPER and _SATE_DEVELOPER != '0'
    try:
        rc, temp_dir, temp_fs = sate_main()
        if not rc:
            raise ValueError("Unknown SATe execution error")
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
        sys.exit("SATe is exiting because of an error:\n%s " % str(x))
