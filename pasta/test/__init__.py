#! /usr/bin/env python

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
SATe testing suite.
"""

import unittest
import re
import os

from pasta import get_logger
from pasta.configure import get_configuration
_LOG = get_logger("pasta.tests")

try:
    import pkg_resources
    TESTS_DIR = pkg_resources.resource_filename("pasta", "test")
    #SCRIPTS_DIR = pkg_resources.resource_filename("pasta", os.path.join(os.pardir, "scripts"))
    _LOG.info("using pkg_resources path mapping")
except:
    LOCAL_DIR = os.path.dirname(__file__)
    TESTS_DIR = os.path.join(LOCAL_DIR, os.path.pardir)
    PACKAGE_DIR = os.path.join(TESTS_DIR, os.path.pardir)
    #SCRIPTS_DIR = os.path.join(PACKAGE_DIR, os.path.pardir, "scripts")
    _LOG.info("using local filesystem path mapping")

TESTS_DATA_DIR = os.path.join(TESTS_DIR, "data")
TESTS_OUTPUT_DIR = os.path.join(TESTS_DIR, "output")
TESTS_COVERAGE_DIR = os.path.join(TESTS_DIR, "coverage")
TESTS_COVERAGE_REPORT_DIR = os.path.join(TESTS_COVERAGE_DIR, "report")
TESTS_COVERAGE_SOURCE_DIR = os.path.join(TESTS_COVERAGE_DIR, "source")

def get_test_file_names():
    """Get list of test file names."""
    path = os.path.dirname(__file__)
    files = os.listdir(path)
    test_file_pattern = re.compile("test.*\.py$", re.IGNORECASE)
    test_files = []
    for f in files:
        if test_file_pattern.search(f):
            test_files.append("pasta.test." + os.path.splitext(f)[0])
    return test_files

def get_test_suite(test_file_names=None):
    """
    Creates a unittest.TestSuite from all of the modules in
    `dendropy.test`. Right now, assumes (a) no subdirectories (though
    this can easily be accommodated) and (b) every test to be run is
    sitting in a module with a file name of 'test*.py', and, conversely,
    every file with a name of 'test*.py' has test(s) to be run.
    """
    if test_file_names is None:
        test_file_names = get_test_file_names()
    tests = unittest.defaultTestLoader.loadTestsFromNames(test_file_names)
    return unittest.TestSuite(tests)

def get_testing_configuration():
    """This function reads the users installation specific configuration files
    (so that we can get the path to the tools), but then strips all other
    settings out of the config object (so that users with different defaults
    settings will not experience failures of tests involving the configurable
    tools).
    """
    c = get_configuration()

    for sect in c._categories:
        g = getattr(c, sect)
        to_del = [opt for opt in list(g.options.keys()) if opt != 'path']
        for d in to_del:
            g.remove_option(d)
    return c

class TestLevel:
    FAST, NORMAL, SLOW, EXHAUSTIVE = 0, 10, 20, 30
    def name(i):
        if i <= TestLevel.FAST:
            return "FAST"
        if i <= TestLevel.NORMAL:
            return "NORMAL"
        if i <= TestLevel.SLOW:
            return "SLOW"
        return "EXHAUSTIVE"

    name = staticmethod(name)

    def name_to_int(l):
        try:
            return int(l)
        except:
            pass
        l = l.upper()
        if l == "FAST":
            return TestLevel.FAST
        if l == "NORMAL":
            return TestLevel.NORMAL
        if l == "SLOW":
            return TestLevel.SLOW
        if l == "EXHAUSTIVE":
            return TestLevel.EXHAUSTIVE
        raise ValueError("TestLevel %s unrecognized" % l)

    name_to_int = staticmethod(name_to_int)

def fast_testing_notification(logger, module_name, message=None, level=TestLevel.FAST):
    if message is None:
        message = "tests skipped"
    logger.warning('\nRunning in %s Testing Level. Skipping %s tests in %s: %s' % (TestLevel.name(get_current_testing_level()), TestLevel.name(level), module_name, message))

def get_current_testing_level():
    l = os.environ.get("SATELIB_TESTING_LEVEL")
    if l is None:
        if "SATELIB_FAST_TESTS" in os.environ:
            return TestLevel.FAST
        return TestLevel.NORMAL
    try:
        return TestLevel.name_to_int(l)
    except:
        _LOG.warn("the value %s for SATELIB_TESTING_LEVEL is not recognized.  Using NORMAL level" % l)
    return TestLevel.NORMAL

def is_test_enabled(level, logger=None, module_name="", message=None):
    tl = get_current_testing_level()
    if level > tl:
        if logger:
            fast_testing_notification(logger, module_name, message, level)
        return False
    return True

def data_source_path(filename=None):
    if filename is None:
        filename = ""
    elif isinstance(filename, list):
        filename = os.path.sep.join(filename)
    return os.path.join(TESTS_DATA_DIR, filename)

def data_target_path(filename=None):
    if filename is None:
        filename = ""
    elif isinstance(filename, list):
        filename = os.path.sep.join(filename)
    return os.path.join(TESTS_OUTPUT_DIR, filename)

def run():
    "Runs all of the unittests"
    runner = unittest.TextTestRunner()
    tests = get_test_suite()
    runner.run(tests)

if __name__ == "__main__":
    run()

