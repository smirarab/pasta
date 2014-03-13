#!/usr/bin/env python
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
# along with this program.      If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu, Mark Holder, and Jeet Sukumaran, University of Kansas

"""
Tests of the pasta/__init__.py
"""
import unittest

from pasta import get_logger
_LOG = get_logger(__name__)

class InitTest(unittest.TestCase):
    def testLogMessage(self):
        _LOG.debug("This is a debug message to test the debug message level")
        _LOG.info("This is an info message to test the info message level")
        _LOG.warn("This is a warn message to test the warn message level")
        _LOG.error("This is an error message to test the error message level")

if __name__ == "__main__":
    unittest.main()

