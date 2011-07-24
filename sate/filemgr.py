#!/usr/bin/env python

#############################################################################
##  This file is part of SATe2.
##  See "LICENSE.txt" for terms and conditions of usage.
#############################################################################

"""
File and path management.
"""

import os
import re
import tempfile
import shutil
from threading import Lock
from sate import get_logger
_LOG = get_logger(__name__)

_ILLEGAL_FILENAME_PATTERN = re.compile(r'[^-_a-zA-Z0-9.]')
def get_safe_filename(filename):
    return "".join(_ILLEGAL_FILENAME_PATTERN.split(filename))

def quoted_file_path(path):
    if '"' not in path:
        return '"'+ path + '"'
    elif "'" not in path:
        return "'" + path + "'"
    else:
        path = path.replace('"', r'\"')
        return '"' + path + '"'

def open_with_intermediates(filepath, mode):
    """Opens a `filepath` in the `mode`, but will create intermediate directories if they are absent."""
    d = os.path.dirname(filepath)
    if d:
        if not os.path.exists(d):
            os.makedirs(d)
        elif not os.path.isdir(d):
            raise IOError('The file "%s" cannot be created because "%s" exists but is not a directory' % (filepath, d))
    return open(filepath, mode)

class TempFS(object):
    '''A wrapper for creating temporary directories that will safeguard against
    removing directories that were not created by SATe (there is no evidence
    that this has happened, but rmtree is a dangerous call and it would be a
    horrible bug to have).

    Note that this class does not protect against incorrect order of deletion.
        If you delete the top level tempdirectory then all of the subdirectories
        will be deleted.
    '''

    def __init__(self):
        self._directories_created = set()
        self._top_level_temp = None
        self._top_level_temp_real = None
        self._directories_created_lock = Lock()

    def _is_already_created(self, real_path):
        self._directories_created_lock.acquire()
        b = real_path in self._directories_created
        self._directories_created_lock.release()
        return b


    def create_subdir(self, dir):
        '''Creates a directory `dir`

        `dir` must be a subdirectory of self.top_level_temp and `dir`
        cannot exist.

        The canonical file path is returned.

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        '''
        rp = os.path.realpath(dir)
        if os.path.exists(rp):
            raise OSError("Path exists: '%s'" % rp)
        if not self._is_in_top_level_temp(rp):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % rp)
        self._directories_created_lock.acquire()
        try:
            already_in = rp in self._directories_created
            if already_in:
                raise OSError("Subdirectory is flagged as having already been created: '%s'" % rp)
            else:
                os.makedirs(rp)
                self._directories_created.add(rp)
                return rp
        finally:
            self._directories_created_lock.release()


    def create_top_level_temp(self, parent, prefix='temp'):
        '''Creates (and stores the path to) the top-level temporary
        directory under `parent`

        `parent` must already exist
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        assert(self._top_level_temp is None)
        assert(self._top_level_temp_real is None)
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            self._top_level_temp = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            self._top_level_temp_real = os.path.realpath(self._top_level_temp)
            self._directories_created.add(self._top_level_temp_real)
            return self._top_level_temp_real
        finally:
             self._directories_created_lock.release()


    def create_temp_subdir(self, parent, prefix='temp'):
        '''Creates (and stores the path to) a temporary directory under `parent`

        `parent` should be "below" the top_level_temp.
        `prefix` is passed to tempfile.mkdtemp

        Client code is responsible for deleting directory by a call to
        `remove_dir` using the same `TempFS` instance.

        The canonical file path is returned.

        '''
        r_parent = os.path.realpath(parent)
        if not os.path.exists(r_parent):
            raise OSError("Path does not exist: '%s'" % r_parent)
        if not os.path.isdir(r_parent):
            raise OSError("Path is not a directory: '%s'" % r_parent)
        if not self._is_in_top_level_temp(r_parent):
            raise OSError("Subdirectory is not under the top level temp dir: '%s'" % r_parent)
        self._directories_created_lock.acquire()
        try:
            d = tempfile.mkdtemp(prefix=prefix, dir=r_parent)
            rp = os.path.realpath(d)
            self._directories_created.add(rp)
            return rp
        finally:
            self._directories_created_lock.release()

    def remove_dir(self, real_path):
        '''Recursively removes `real_path` from the filesystem if it is
        listed as one of the directories created by this TempFS object (or raises
        a ValueError if it is not listed).
        '''
        self._directories_created_lock.acquire()
        try:
            if real_path in self._directories_created:
                self._directories_created.remove(real_path)
                if (real_path == self._top_level_temp) or (real_path == self._top_level_temp_real):
                    self._top_level_temp = None
                    self._top_level_temp_real = None
            else:
                raise ValueError("'%s' is not registered as a temporary directory that was created by this process!" % real_path)
        finally:
             self._directories_created_lock.release()
        _LOG.debug("Removing temp dir: '%s'" % real_path)
        # Because we raise an exception if real_path is not in _directories_created,
        #   we only get down here if real_path was in self._directories_created
        #   thus, this call should only delete a directory created by this
        #   TempFS instance.
        if os.path.exists(real_path):
            try:
                shutil.rmtree(real_path)
                return True
            except:
                return False
        return False

    def get_remaining_directories(self):
        '''Returns a copy of the set of directories that have been created but
        not deleted.
        '''
        self._directories_created_lock.acquire()
        c = set(self._directories_created)
        self._directories_created_lock.release()
        return c

    def get_top_level_temp(self):
        return self._top_level_temp_real
    top_level_temp = property(get_top_level_temp)

    def _is_in_top_level_temp(self, real_path):
        if self._top_level_temp_real is None:
            raise ValueError("_top_level_temp has not been set, yet!")
        in_common = os.path.commonprefix([self._top_level_temp_real, real_path])
        return in_common == self._top_level_temp_real

class SateProducts(object):
    """
    Handles paths to all (final) output produced by SATe.
    """

    meta_product_types = {
            "score": ".score",
            "tree": ".tre",
            "run_log" : ".out",
            "err_log" : ".err"
            }

    def __init__(self, sate_user_settings):
        """
        Configures self based on a fully-populated `SateUserSettings` config
        object.
        """
        self.sate_user_settings = sate_user_settings
        self.job_name = self.sate_user_settings.commandline.job
        if self.job_name is None:
            self.job_name = "satejob"
        self._job_file_name = get_safe_filename(self.job_name)
        self._disambiguator_idx = 0

        if not self.sate_user_settings.sate.output_directory:
            self._output_directory = self.get_input_source_directory()
        else:
            self._output_directory = self.sate_user_settings.sate.output_directory

        if os.path.exists(self._output_directory) and not os.path.isdir(self._output_directory):
            raise Exception("Requested output directory is not valid because a file of the same name already exists: %s" % self._output_directory)

        # ensure input sources populated
        assert self.sate_user_settings.input_seq_filepaths
        # alignment input/output names
        self._original_input_files = [os.path.abspath(f) for f in self.sate_user_settings.input_seq_filepaths]
        #self.output_alignment_suffixes = ["." + get_safe_filename(f) + ".aln" \
        #        for f in self.original_input_files]
        self._output_alignment_suffixes = []
        self._input_fpath_alignment_suffix_map = {}
        self._alignment_suffix_input_fpath_map = {}
        fasta_extension_pattern = re.compile(r"\.fast?a?$")
        for fidx, f in enumerate(self._original_input_files):
            safe_fn = get_safe_filename(os.path.basename(f))
            trunc_fn = fasta_extension_pattern.sub("", safe_fn)
            fn_stem = "marker%03d.%s" % (fidx+1, trunc_fn)
            suffix = "." + fn_stem + ".aln"
            sidx = 1
            while suffix in self._output_alignment_suffixes:
                suffix = "." + fn_stem + "-" + str(sidx) + ".aln"
                sidx += 1
            self._output_alignment_suffixes.append(suffix)
            self._input_fpath_alignment_suffix_map[f] = suffix
            self._alignment_suffix_input_fpath_map[suffix] = f

        # initialize/create attributes, setting to dummy values
        for stream_name in self.meta_product_types:
            self._set_stream(stream_name, None)
        self.alignment_streams = []
        self.input_fpath_alignment_stream_map = {}

        # dummy output prefix
        self.output_prefix = None

        # create working output streams
        self.setup()

    def _compose_stream_attr(self, stream_name):
        return stream_name + "_stream"

    def _get_stream(self, stream_name):
        return getattr(self, self._compose_stream_attr(stream_name), None)

    def _set_stream(self, stream_name, value):
        setattr(self, self._compose_stream_attr(stream_name), value)

    def setup(self):
        """
        Checks for file name clashes, disambiguates if neccessary, and creates
        output files.
        """
        self.create_output_prefix()
        self.create_product_paths()

    def create_product_paths(self):
        assert self.output_prefix
        for stream_name, product_extension in self.meta_product_types.items():
            output_path = self.output_prefix + product_extension
            stream = open_with_intermediates(output_path, "w")
            self._set_stream(stream_name, stream)
        for asi, sf in enumerate(self._output_alignment_suffixes):
            output_path = self.output_prefix + sf
            stream = open_with_intermediates(output_path, "w")
            self.alignment_streams.append(stream)
            self.input_fpath_alignment_stream_map[self._alignment_suffix_input_fpath_map[sf]] = stream

    def create_output_prefix(self):
        output_prefix_stem = os.path.join(self._output_directory, self._job_file_name)
        idx = 0
        disambiguator = ""
        while True:
            output_prefix = output_prefix_stem + disambiguator
            if not self.check_for_existing_files(output_prefix):
                break
            idx += 1
            disambiguator = "%d" % idx
        self.output_prefix = output_prefix

    def check_for_existing_files(self, output_prefix):
        for ext in self.meta_product_types.values():
            if os.path.exists(output_prefix + ext):
                return True
        for fn in self._output_alignment_suffixes:
            if os.path.exists(output_prefix + fn):
                return True
        return False

    def get_input_source_directory(self):
        """
        Given a configuration object, returns the directory of the input file(s).
        """
        options = self.sate_user_settings.commandline
        if options.multilocus:
            # multilocus dataset: assume directory is given as input source
            return os.path.abspath(options.input)
        else:
            # single locus dataset: return directory nanme
            return os.path.dirname(os.path.abspath(options.input))


