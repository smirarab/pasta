#!/usr/bin/env python

"""Accessor for runtime configuration settings.

In general one should be able to simply call get_configuration()

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

import ConfigParser
import os
import platform
import sys

from sate import get_logger, SATE_SYSTEM_PATHS_CFGFILE
from sate.settings import SateUserSettings
from sate.tools import get_external_tool_classes
import sate
from sate import filemgr

_LOG = get_logger(__name__)

def get_invoke_run_sate_command():
    """Used by GUI's this tries to return the python invocation for the run_sate
    program or script.
    """
    if sate.sate_is_frozen():
        if platform.system() == 'Windows':
            return ['run_sate.exe']
        elif platform.system() == 'Darwin':
            return [filemgr.quoted_file_path(sys.executable),
                    filemgr.quoted_file_path(os.path.join(sate.sate_home_dir(), 'run_sate.py'))]
        else:
            raise OSError('SATe is not frozen for %s' % platform.system())
    else:
        return [filemgr.quoted_file_path(sys.executable), filemgr.quoted_file_path('run_sate.py')]

_DEFAULT_TOOLS_PATH = None

def init_sate(sate_home=None):
    """
    Sets the _DEFAULT_TOOLS_PATH.
    """
    global _DEFAULT_TOOLS_PATH
    bin_dir = sate.sate_tools_deploy_dir(default_to_dev_dir=True)
    if _DEFAULT_TOOLS_PATH is None:
        _DEFAULT_TOOLS_PATH = {}
    for i in get_external_tool_classes():
        tool_name = i.section_name.split()[0]
        _DEFAULT_TOOLS_PATH[tool_name] = os.path.join(bin_dir, tool_name)
        if platform.system() == 'Windows' and tool_name != 'opal':
            _DEFAULT_TOOLS_PATH[tool_name] += '.exe'
    _DEFAULT_TOOLS_PATH['opal'] += '.jar'

def set_configuration_from_defaults(cfg):
    "Uses _DEFAULT_TOOLS_PATH to add paths to external tools to `cfg`"
    global _DEFAULT_TOOLS_PATH
    if _DEFAULT_TOOLS_PATH is None:
        init_sate()
    for name, path in _DEFAULT_TOOLS_PATH.items():
        x = getattr(cfg, name)
        x.path = path

def get_configuration(configfile=None):
    """Returns an instance of SateUserSettings that reflects the current
    defaults based on:
        1. paths inferred from SATE_SYSTEM_PATHS_CFGFILE
        2. any settings in `configfile` (these will overwrite the settings
            based on the defaults from sate_home).
    """
    cfg = SateUserSettings()
    set_configuration_from_defaults(cfg)

    if os.path.exists(SATE_SYSTEM_PATHS_CFGFILE):
        paths_parser = ConfigParser.RawConfigParser()
        try:
            paths_parser.read(SATE_SYSTEM_PATHS_CFGFILE)
        except:
            ### TODO: wrap up in messaging system
            sys.stderr.write("The specified configuration file %s cannot be read, the default settings are used instead.\n" % SATE_SYSTEM_PATHS_CFGFILE)
        else:
            for tool_class in get_external_tool_classes():
                try:
                    new_path = paths_parser.get(tool_class.section_name, 'path')
                    if new_path:
                        tool_name = tool_class.section_name.split()[0] # TODO: this (taking the first word) is a bad idea we can get clashes in _DEFAULT_TOOLS_PATH
                        getattr(cfg, tool_name).path = new_path
                except:
                    ### TODO: wrap up in messaging system
                    sys.stderr.write("Section '%s' not found: using default settings instead.\n" % tool_class.section_name)

    if configfile is not None:
        if os.path.isfile(configfile):
            cfg.read_config_filepath(configfile)
        else:
            ### TODO: wrap up in messaging system
            sys.stderr.write("The specified configuration file %s cannot be found, the default settings are used instead.\n" % configfile)
    return cfg

def get_input_source_directory(config):
    """
    Given a configuration object, returns the directory of the input file(s).
    """
    options = config.commandline
    if options.multilocus:
        # multilocus dataset: assume directory is given as input source
        return os.path.abspath(options.input)
    else:
        # single locus dataset: return directory nanme
        return os.path.dirname(os.path.abspath(options.input))


