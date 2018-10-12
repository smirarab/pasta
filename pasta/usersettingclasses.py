"""Classes for validating user-controlled values."""

# This file is part of PASTA, which is forked from SATe

# PASTA, like SATe is free software: you can redistribute it and/or modify
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

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import os
import glob
from pasta import get_logger
from pasta.filemgr import open_with_intermediates

_LOG = get_logger(__name__)

def _underscores_to_dashes(s):
    return '-'.join(s.split('_'))

class UserSetting(object):
    def __init__(self, name, default, **kwargs):
        self.name = name
        self.default = default
        self.short_name = kwargs.get('short_name')
        self.subcategory = kwargs.get('subcategory')
        self.must_be_valid = kwargs.get('must_be_valid')
        if self.subcategory is None:
            self.subcategory = ''
        self.help = kwargs.get('help')
        self._value = self.default
    def is_valid(self):
        return self.value is not None

    def get_value(self):
        return self._value
    def set_value(self, v):
        self._value = v
    value = property(get_value, set_value)
    def as_config_parsable(self):
        v = self.value
        if v is None:
            return None
        return str(v)
    def _get_optparse_option_args_dict(self):
        kwargs = {'dest' : self.name,
                    'default' : None
                    }
        if self.help is not None:
            kwargs['help'] = self.help
        return kwargs

    def add_to_optparser(self, parser):
        kwargs = self._get_optparse_option_args_dict()
        long_name = '--%s' % _underscores_to_dashes(self.name)
        if self.short_name and len(self.short_name) == 1:
            parser.add_option('-%s' % self.short_name, long_name, **kwargs)
        else:
            parser.add_option(long_name, **kwargs)


class StringUserSetting(UserSetting):
    pass

class ChoiceUserSetting(UserSetting):
    def __init__(self, name, default, choices, multiple_choices=False, **kwargs):
        UserSetting.__init__(self, name, default, **kwargs)
        self.multiple_choices = multiple_choices
        self.choices = [i.lower() for i in choices]
        self.value = default
    def get_value(self):
        return self._value
    def set_value(self, v):
        if v is None:
            self._value = None
            return
        if self.multiple_choices:
            o = v.strip()
            if o.startswith('[') and o.endswith(']'):
                quoted_toks = [i.strip() for i in o[1:-1].split()]
                sq = "'"
                dq = '"'
                x = []
                if (i.startswith(sq) and i.endswith(sq)) or (i.startswith(dq) and i.endswith(dq)):
                    x.append(i[1:-1].lower())
            else:
                x = [i.lower() for i in o.split()]
            for el in x:
                if el not in self.choices:
                    raise ValueError("%s is not a valid setting for %s. Exepcting one of %s" % (el, self.name, str(self.choices)))
            self._value = x
        else:
            el = v.lower()
            if el not in self.choices:
                raise ValueError("%s is not a valid setting for %s" % (el, self.name))
            self._value = el
    value = property(get_value, set_value)
    def as_config_parsable(self):
        v = self.value
        if v is None:
            return None
        if self.multiple_choices:
            return ' '.join(v)
        return str(v)

class NumberUserSetting(UserSetting):
    def __init__(self, name, default, min=None, max=None, **kwargs):
        UserSetting.__init__(self, name, default, **kwargs)
        self.max = max
        self.min = min
    def get_value(self):
        return self._value
    def set_value(self, v):
        v = self.conv_type(v)
        if (self.max is not None) and v > self.max:
            raise ValueError("%s must be <= %s" % (self.name, str(self.max)))
        if (self.min is not None) and v < self.min:
            raise ValueError("%s must be >= %s" % (self.name, str(self.min)))
        self._value = v
    value = property(get_value, set_value)

class BoolUserSetting(UserSetting):
    def __init__(self, name, default, min=None, max=None, **kwargs):
        UserSetting.__init__(self, name, default, **kwargs)
        self.value = default
        self.default = self.value
    def get_value(self):
        return self._value
    def set_value(self, v):
        if isinstance(v, str):
            l = v.lower()
            self._value = l in ['1', 'true']
        elif isinstance(v, int):
            self._value = (v != 0)
        elif isinstance(v, bool):
            self._value = v
        else:
            raise TypeError("Expecting a boolean for %s" % self.name)
    value = property(get_value, set_value)

    def add_to_optparser(self, parser):
        kwargs = UserSetting._get_optparse_option_args_dict(self)
        if bool(self.default):
            kwargs['action'] = 'store_false'
            long_name = '--no-%s' % _underscores_to_dashes(self.name)
        else:
            kwargs['action'] = 'store_true'
            long_name = '--%s' % _underscores_to_dashes(self.name)
        if self.short_name and len(self.short_name) == 1:
            parser.add_option('-%s' % self.short_name, long_name, **kwargs)
        else:
            parser.add_option(long_name, **kwargs)

class FloatUserSetting(NumberUserSetting):
    def __init__(self, name, default, min=None, max=None, **kwargs):
        self.conv_type = float
        NumberUserSetting.__init__(self, name, default, min=min, max=max, **kwargs)
        self.value = default
    def _get_optparse_option_args_dict(self):
        kwargs = UserSetting._get_optparse_option_args_dict(self)
        kwargs['type'] = 'float'
        kwargs['metavar'] = '#.#'
        return kwargs

class IntUserSetting(NumberUserSetting):
    def __init__(self, name, default, min=None, max=None, **kwargs):
        self.conv_type = int
        NumberUserSetting.__init__(self, name, default, min=min, max=max, **kwargs)
        self.value = default
    def _get_optparse_option_args_dict(self):
        kwargs = UserSetting._get_optparse_option_args_dict(self)
        kwargs['type'] = 'int'
        kwargs['metavar'] = '#'
        return kwargs

class UserSettingGroup(object):
    def __init__(self, name):
        self.name = name
        self.options = {}

    def add_option(self, k, o):
        self.options[k] = o

    def remove_option(self, o):
        del self.options[o]

    def dict(self):
        d = {}
        for k, v in self.options.items():
            d[k] = v.value
        return d

    def keys(self):
        return list(self.dict().keys())

    def __getattr__(self, att_name):
        try:
            return self.__dict__['options'][att_name].value
        except:
            raise AttributeError("Attribute %s not found" % att_name)
    def __setattr__(self, att_name, value):
        if att_name in ['name', 'options']:
            self.__dict__[att_name] = value
        else:
            o = self.options.get(att_name)
            if o is None:
                raise AttributeError("Attribute %s not found" % att_name)
            o.value = value
    def set_config_parser_fields(self, p):
        if not p.has_section(self.name):
            p.add_section(self.name)
        for opt in self.all_options():
            v = opt.as_config_parsable()
            if v is not None:
                p.set(self.name, opt.name, v)

    def read_config_parser_fields(self, parsed):
        try:
            items = parsed.items(self.name)
        except configparser.NoSectionError:
            return
        for p in items:
            k, v = p
            opt = self.options.get(k)
            if opt is not None:
                opt.value = v
            else:
                _LOG.warn('Unknown option "%s" in section "%s" skipped!' % (k, self.name))

    def all_options(self):
        key_opt_list =[i for i in self.options.items()]
        key_opt_list.sort()
        return [i[1] for i in key_opt_list]

    def subcategories(self):
        s = set()
        for opt in self.all_options():
            s.add(opt.subcategory)
        l = list(s)
        l.sort()
        return l

    def add_to_optparser(self, parser):
        from optparse import OptionGroup
        s_list = self.subcategories()
        for s in s_list:
            n = self.name
            if n == 'sate': # how is this  for a hack
                n = 'SATe'
            if s:
                g = OptionGroup(parser, '%s %s options' % (n, s))
            else:
                g = OptionGroup(parser, '%s options' % n)
            parser.add_option_group(g)
            for opt in self.all_options():
                if opt.subcategory == s:
                    opt.add_to_optparser(g)

    def set_values_from_dict(self, d):
        """Walks through the keys of `d` and sets any matching UserSetting.name
        to the value in `d` ONLY IF the value in `d` is not None.
        """
        for o in self.all_options():
            v = d.get(o.name)
            if v is not None:
                o.value = v

def get_list_of_seq_filepaths_from_dir(dir_path):
    """
    Given a directory, finds all '*.fas' and '*.fasta' files and returns them 
    as a list of absolute paths.
    """
    if not dir_path:
        raise Exception("Sequence file directory must be specified for multilocus analysis.\n")
    if not os.path.exists(dir_path):
        raise Exception("The input sequence file directory '%s' does not exist.\n" % dir_path)
    if not os.path.isdir(dir_path):
        raise Exception("The input sequence files must be put into one directory for multilocus analysis, '%s' is not a directory.\n" % dir_path)
    dot_fas_glob = os.path.join(os.path.abspath(dir_path), '*.fas')
    dot_fasta_glob = os.path.join(os.path.abspath(dir_path), '*.fasta')
    seq_filename_list = glob.glob(dot_fas_glob) + glob.glob(dot_fasta_glob)
    if len(seq_filename_list) == 0:
        raise Exception("No files found with extension '.fas' or '.fasta' in directory: '%s'" % dir_path)
    seq_filename_list.sort()
    return seq_filename_list

class UserSettingsContainer(object):

    def __init__(self):
        self._categories = []
        self._config_parser = configparser.RawConfigParser()
        self.input_seq_filepaths = []

    def read_seq_filepaths_from_dir(self, dir_path):
        """
        Given a directory, finds all '*.fas' and '*.fasta' files and adds
        to `self.input_seq_filepaths`
        """
        seq_filename_list = get_list_of_seq_filepaths_from_dir(dir_path)
        self.input_seq_filepaths.extend(seq_filename_list)

    def read_seq_filepaths_from_delimited_string(self, dstr, delimiter=","):
        """
        Given a `delimiter`-delimited string of tokens, split strings into
        tokens and appends tokens to `self.input_seq_filepaths`. Assumes tokens
        do not include delimiter character or leading/trailing spaces.
        """
        fpaths = [s.strip() for s in dstr.split(delimiter)]
        self.input_seq_filepaths.extend(fpaths)

    def read_seq_filepaths(self, src, multilocus=False, delimiter=','):
        """
        If `multilocus` is True, then assumes `src` is a directory and scans
        directory for all files ending with `.fas` and `.fasta`. Otherwise,
        assumes that `src` is a filepath given directly.
        NOTE: this behavior is for compatibility with previous SATe behavior.
        TODO: NEW BEHAVIOR: `src` will be treated as a `delimiter`-delimited
        string.
        """
        if multilocus:
            #self.read_seq_filepaths_from_delimited_string(src)
            self.read_seq_filepaths_from_dir(src)
        else:
            self.input_seq_filepaths = [src]
        return self.input_seq_filepaths

    def read_config_filepath(self, filepath):
        self._config_parser.read(filepath)
        for gn in self._categories:
            g = getattr(self, gn)
            g.read_config_parser_fields(self._config_parser)

    def save_to_filepath(self, filepath):
        if filepath is None:
            filepath = os.path.expanduser(os.path.join( '~', '.pasta', 'pasta.cfg'))
        f = open_with_intermediates(filepath, 'w')
        for g in self.get_categories():
            g.set_config_parser_fields(self._config_parser)
        self._config_parser.write(f)
        f.close()

    def create_aligner(self, temp_fs, name=None):
        sate = self.sate
        if name is None:
            name = sate.aligner
        _LOG.debug("Creating aligner %s" % name)
        n = name.lower()
        try:
            g = getattr(self, n)
        except AttributeError:
            raise RuntimeError("Unrecognized aligner %s" % name)
        p = ''
        try:
            p = g.path
        except:
            pass
        if not p:
            raise RuntimeError("Cannot create a wrapper around the aligner %s because a path setting was not found in the configuration for this tool" % name)

        from pasta.tools import get_aligner_classes, CustomAligner
        d = g.dict()

        for c in get_aligner_classes():
            if c.section_name.startswith(n):
                _LOG.debug("Creating instance of %s" % str(c.__name__))
                return c(temp_fs=temp_fs, **d)
        return CustomAligner(name, temp_fs=temp_fs, **d)

    def create_treeshrink_wrapper(self, temp_fs):
        from pasta.tools import TreeShrink
        d = self.treeshrink.dict()
        return TreeShrink(temp_fs=temp_fs, **d)

    def create_merger(self, temp_fs, name=None):
        sate = self.sate
        if name is None:
            name = sate.merger
        _LOG.debug("Creating merger %s" % name)
        n = name.lower()
        try:
            g = getattr(self, n)
        except AttributeError:
            raise RuntimeError("Unrecognized merger %s" % name)
        p = ''
        try:
            p = g.path
        except:
            pass
        if not p:
            raise RuntimeError("Cannot create a wrapper around the merger %s because a path setting was not found in the configuration for this tool" % p)

        from pasta.tools import MergerClasses
        d = g.dict()
        for c in MergerClasses:
            if c.section_name.startswith(n):
                _LOG.debug("Creating instance of %s" % str(c.__name__))
                return c(temp_fs=temp_fs, **d)
        raise RuntimeError("The merger tool %s is not supported at this time" % name)

    def create_tree_estimator(self, temp_fs, name=None):
        sate = self.sate
        if name is None:
            name = sate.tree_estimator
        _LOG.debug("Creating tree estimator %s" % name)
        n = name.lower()

        try:
            g = getattr(self, n)
        except AttributeError:
            raise RuntimeError("Unrecognized tree_estimator %s" % name)
        p = ''
        try:
            p = g.path
        except:
            pass
        if not p:
            raise RuntimeError("Cannot create a wrapper around the tree_estimator %s because a path setting was not found in the configuration for this tool" % name)

        from pasta.tools import TreeEstimatorClasses
        d = g.dict()
        for c in TreeEstimatorClasses:
            if c.section_name.startswith(n):
                _LOG.debug("Creating instance of %s" % str(c.__name__))
                return c(temp_fs=temp_fs, **d)
        raise RuntimeError("The tree estimation tool %s is not supported at this time" % name)

    def get(self, category):
        "Returns a UserSettingGroup based on name or None"
        i = self.__dict__.get(category)
        if isinstance(i, UserSettingGroup):
            return i
        return None

    def get_categories(self):
        c = [i for i in self.__dict__.items() if isinstance(i[1], UserSettingGroup)]
        c.sort()
        return [i[1] for i in c]

    def set_values_from_dict(self, d):
        """Walks through the keys of `d` and sets any matching UserSetting.name
        to the value in `d` ONLY IF the value in `d` is not None.
        """
        for c in self.get_categories():
            c.set_values_from_dict(d)
    def dicts(self):
        d = {}
        for  c in self.get_categories():
            d[c.name] = c.dict()
        return d
