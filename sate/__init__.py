#!/usr/bin/env python

#############################################################################
##  this file is part of sate.
##  see "license.txt" for terms and conditions of usage.
#############################################################################

"""
Functions for configuring a runtime logger for the sate module.
"""

PROGRAM_NAME = "SATe"
PROGRAM_AUTHOR = ["Jiaye Yu", "Mark T. Holder", "Jeet Sukumaran"]
PROGRAM_LICENSE = "GNU General Public License, version 3"
PROGRAM_VERSION = "2.0.3"
PROGRAM_YEAR = "2009-2011"
PROGRAM_DESCRIPTION = "Simultaneous Alignment and Tree Estimation"
PROGRAM_WEBSITE = "http://phylo.bio.ku.edu/software/sate/sate.html"
PROGRAM_INSTITUTE = "Department of Ecology and Evolutionary Biology, Univesity of Kansas"
PROGRAM_LONG_DESCRIPTION = """
SATe performs iterative realignment and tree inference.

Minimally you must provide a sequence file (with the '--input' option) a starting tree is optional.

The command line allows you to alter the behavior of the algorithm (termination criteria, when the algorithm switches to "Blind" acceptance of new alignments, how the tree is decomposed to find subproblems to be used, and the external tools to use).

Options can also be passed in as configuration files.

With the format:
####################################################
[commandline]
option-name = value

[sate]
option-name = value
####################################################

If you tell sate to keep its temporary files (-k option), then the configuration for the run will be stored in a file "last_used.cfg" in the "temporary" directory.

If configuration files are read in the order they occur as arguments (with values in later files replacing previously read values). Options specified in the command line are read last. Thus these values "overwrite" any settings from the configuration files.
"""

__all__ = []

import os
import sys
import platform
import logging

_LOGGING_LEVEL_ENVAR = "SATE_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SATE_LOGGING_FORMAT"

# global debugging flag
if "SATE_DEBUG" in os.environ:
    if os.environ["SATE_DEBUG"]:
        if os.environ["SATE_DEBUG"].lower()[0] in ["1", "t", "y", "d"]:
            GLOBAL_DEBUG = True
        else:
            GLOBAL_DEBUG = False
    else:
        GLOBAL_DEBUG = False
else:
    GLOBAL_DEBUG = False

def get_logging_level():
    """Checks environment for SATE_LOGGING_LEVEL and returns a logging level
    integer.
    """
    import logging
    if _LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET
    return level

def get_logger(name="sate"):
    """
    Returns a logger with name set as given, and configured
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """
    logger_set = False
    logger = logging.getLogger(name)
    if not logger_set:
        level = get_logging_level()
        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (line %(lineno)d): %(levelname) 8s: %(message)s")
        simple_formatter = logging.Formatter("%(levelname) 8s: %(message)s")
        default_formatter = None
        logging_formatter = default_formatter
        if _LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter
        if logging_formatter is not None:
            logging_formatter.datefmt = '%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
    return logger

TIMING_LOG = logging.getLogger("TIMING_LOG")
TIMING_LOG.setLevel(logging.CRITICAL)

def set_timing_log_filepath(fp):
    global TIMING_LOG
    if not fp:
        TIMING_LOG.setLevel(logging.CRITICAL)
    else:
        TIMING_LOG.setLevel(logging.DEBUG)
        h = logging.FileHandler(fp)
        f = logging.Formatter("[%(asctime)s] : %(message)s")
        f.datefmt = '%D %H:%M:%S'
        h.setLevel(logging.DEBUG)
        h.setFormatter(f)
        TIMING_LOG.addHandler(h)

def log_exception(logger):
    '''Logs the exception trace to the logObj as an error'''
    import traceback, cStringIO
    s = cStringIO.StringIO()
    traceback.print_exc(None, s)
    logger.debug(s.getvalue())

class Messenger(object):
    """
    Wraps reporting of messages, progress, warnings and errors to users.
    Singleton (instantiated below).
    """

    def __init__(self):
        self.err_log_streams = [sys.stderr]
        self.run_log_streams = [sys.stdout]

    def _format_msg(self, msg, msg_type):
        return "SATe %s: %s\n" % (msg_type, msg)

    def _write_to_streams(self, streams, msg, flush=True):
        for s in streams:
            s.write(msg)
            if flush:
                s.flush()

    def _write_to_err_streams(self, msg, flush=True):
        self._write_to_streams(self.err_log_streams, msg, flush=flush)

    def _write_to_out_streams(self, msg, flush=True):
        self._write_to_streams(self.run_log_streams, msg, flush=flush)

    def send_error(self, msg):
        msg = self._format_msg(msg, "ERROR")
        self._write_to_err_streams(msg)

    def send_warning(self, msg):
        msg = self._format_msg(msg, "WARNING")
        self._write_to_err_streams(msg)

    def send_info(self, msg):
        msg = self._format_msg(msg, "INFO")
        self._write_to_out_streams(msg)

##############################################
## Instantiation Of Singleton
## more idiomatic way might be to implement
## this functionality as a module instead of
## a class
MESSENGER = Messenger()

##############################################
## Other globals
SATE_SYSTEM_PATHS_CFGFILE =  os.path.abspath(os.path.join(os.path.expanduser('~'), '.sate', 'sate_tool_paths.cfg'))
SATE_GUI_RESOURCES = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "resources", "images"))
SATE_SCRIPT_RESOURCES = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "resources", "scripts"))
DEFAULT_MAX_MB = os.environ.get("SATE_MAX_MB", 1024)

##############################################
## Path configuration

class ToolsDirNotFoundError(Exception):

    def __init__(self, paths_tried, env_var_name, *args, **kwargs):
        self.paths_tried = paths_tried
        self.env_var_name = env_var_name
        Exception.__init__(self, *args, **kwargs)

    def __str__(self):
        return "Could not find SATe tools bundle directory: '%s'; set $%s and try again" % (
                ", ".join("'" + s + "'" for s in self.paths_tried),
                self.env_var_name)

def sate_is_frozen():
    """Will return True if SATe is frozen.
    """
    import imp
    return (
        hasattr(sys, "frozen")          # new py2exe
        or hasattr(sys, "importers")    # old py2exe
        or imp.is_frozen("__main__")    # tools/freeze
    )

def sate_home_dir():
    """Attempts to return the directory that is the parent of the binary directories.
    """
    if sate_is_frozen():
        if platform.system() == "Darwin":
            retpath = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(sys.executable))), 'Resources')
        else:
            retpath = os.path.dirname(sys.executable)
        return retpath
        #retpath = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(sys.executable))), 'Resources') if platform.system() == "Darwin" else os.path.dirname(sys.executable)
        #return retpath
    else:
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def sate_tools_deploy_subpath():
    return "bin"

def sate_tools_deploy_dir(default_to_dev_dir=False):
    """
    If the environmental variable $SATE_TOOLS_RUNDIR is set, then this is returned.
    If not, and a '../bin' directory exists (relative to *this* file, i.e.
    __file__) *or* `default_to_dev_dir` is False, then that is returned.
    If not, then the return value of `sate_tools_dev_dir` is returned
    """
    if "SATE_TOOLS_RUNDIR" in os.environ:
        bin_path = os.path.expanduser(os.path.expandvars(os.environ["SATE_TOOLS_RUNDIR"]))
    else:
        home_path = sate_home_dir()
        bin_path = os.path.join(home_path, sate_tools_deploy_subpath())
    if os.path.exists(bin_path) or not default_to_dev_dir:
        return bin_path
    try:
        return sate_tools_dev_dir()
    except ToolsDirNotFoundError:
        #it is possible that correct paths to the tools will be
        #provided for in system configuration file
        #let client code deal with that ...
        return bin_path

def sate_tools_dev_dir(platform_name=None):
    """
    If the environmental variable $SATE_TOOLS_DEVDIR is set, then this is used.
    Otherwise, '../../sate-tools-linux', '../../sate-tools-mac', or
    '../../sate-tools-win' is used (path relative to this __file__).
    """
    if "SATE_TOOLS_DEVDIR" in os.environ and platform_name is None:
        bin_path = os.path.expanduser(os.path.expandvars(os.environ["SATE_TOOLS_DEVDIR"]))
    else:
        if platform_name is None:
            platform_name = platform.system()
        root_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        if platform_name == 'Windows':
            sub_path = "sate-tools-win"
        elif platform_name == 'Darwin':
            sub_path = "sate-tools-mac"
        elif platform_name == 'Linux':
            sub_path = "sate-tools-linux"
        else:
            raise OSError("SATe does not bundle tools for '%s' at this time!" % platform_name)
        bin_path = os.path.join(root_path, sub_path)
    if not os.path.exists(bin_path):
        raise ToolsDirNotFoundError(paths_tried=[bin_path],
                env_var_name="SATE_TOOLS_DEVDIR")
    return bin_path

