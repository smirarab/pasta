#!/usr/bin/env python

#############################################################################
##  this file is part of sate.
##  see "license.txt" for terms and conditions of usage.
#############################################################################

"""
Simple classes for reading and manipulating sequence data matrices
"""

import re
from copy import deepcopy
from sate import get_logger, log_exception, MESSENGER
from sate.filemgr import open_with_intermediates

_LOG = get_logger(__name__)
_INDEL = re.compile(r"[-]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")

DATASET_TAXA_ATTR = "taxon_sets"
DATASET_CHAR_ATTR = "char_matrices"

def is_sequence_legal(seq):
    """Check for illegal characters -- TODO, currently returns True"""
    return True

def read_fasta(src):
    """generator that returns (name, sequence) tuples from either a FASTA
    formatted file or file object.
    """
    file_obj = None
    if isinstance(src, str):
        try:
            file_obj = open(src, "rU")
        except IOError:
            print("The file `%s` does not exist, exiting gracefully" % src)
    elif isinstance(src, file):
            file_obj = src
    else:
        raise TypeError('FASTA reader cannot recognize the source of %s' % src)
    name = None
    seq_list = list()
    for line_number, i in enumerate(file_obj):
        if i.startswith('>'):
            if name:
                yield name, ''.join(seq_list)
                seq_list = list()
            name = i[1:].strip()
        else:
            seq = ''.join(i.strip().upper().split())
            if not is_sequence_legal(seq):
                raise Exception("Error: illegal characeters in sequence at line %d" % line_number)
            seq_list.append(seq)
    yield name, ''.join(seq_list)
    if isinstance(src, str):
        file_obj.close()

def read_nexus(src):
    "TODO: use dendropy to do this."
    raise NotImplementedError('Input of NEXUS file format is not supported yet.')

def read_phylip(src):
    "TODO: use dendropy to do this."
    raise NotImplementedError('Input of PHYLIP file format is not supported yet.')

def write_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name, seq in alignment.items():
        file_obj.write('>%s\n%s\n' % (name, seq) )
    if isinstance(dest, str):
        file_obj.close()

def write_phylip(alignment, dest):
    """Writes the `alignment` in relaxed PHYLIP format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    names = alignment.get_sequence_names()
    assert(names)
    ntax = len(names)
    seq = alignment[names[0]]
    nchar = len(seq)
    file_obj.write('%s\t%s\n' % (ntax, nchar) )
    for k in names:
        assert len(k.split()) == 1
        seq = alignment[k]
        assert len(seq) == nchar
        file_obj.write('%s\n%s\n' % (k, seq))
    if isinstance(dest, str):
        file_obj.close()

def write_nexus(alignment, file_obj):
    "TODO use dendropy"
    raise NotImplementedError('Output of NEXUS file format is not supported yet.')

class Alignment(dict, object):
    """A simple class that maps taxa names to sequences.
    TODO: switch to dendropy character_matrix
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None

    def get_datatype(self):
        return self._datatype

    def set_datatype(self, d):
        if d is None:
            self._datatype = None
        else:
            self._datatype = d.upper()

    datatype = property(get_datatype, set_datatype)

    def get_sequence_names(self):
        "returns a list of sequence names"
        return self.keys()

    def get_num_taxa(self):
        "returns the number sequences"
        return len(self.get_sequence_names())

    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        file_obj = open(filename, 'r')
        return self.read_file_object(file_obj, file_format=file_format)

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta
        elif ( file_format.upper() == 'NEXUS' ):
            read_func = read_nexus
        elif ( file_format.upper() == 'PHYLIP' ):
            read_func = read_phylip
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq

    def write_filepath(self, filename, file_format='FASTA'):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        file_obj = open_with_intermediates(filename,'w')
        self.write(file_obj, file_format=file_format)
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_fasta
        elif ( file_format.upper() == 'NEXUS' ):
            write_func = write_nexus
        elif ( file_format.upper() == 'PHYLIP' ):
            write_func = write_phylip
        else:
            write_func = write_fasta
        write_func(self, file_obj)

    def write_unaligned_fasta(self, filename):
        """Writes the sequence data without gaps as FASTA, but note that the
        lines may bet "ragged".
        """
        file_obj = open_with_intermediates(filename, 'w')
        for name, seq in self.items():
            file_obj.write('>%s\n%s\n' % (name, re.sub(_INDEL, '', seq)))
        file_obj.close()

    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            new_alignment[key] = self[key]
        return new_alignment

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            v = self.values()
            first_seq_len = len(v[0])
            return all([len(i) == first_seq_len for i in v])

    def partition_info(self, base=0):
        return (self.datatype, 1+base, self.sequence_length() + base)

    def sequence_length(self):
        if self.is_aligned():
            return len(self.values()[0])

    def max_sequence_length(self):
        return max(len(v) for v in self.values())


class SequenceDataset(object):
    """Class for creating a dendropy reader, validating the input, and
    keeping mapping of real taxa names to "safe" versions that will not
    cause problems for our alignment and tree tools.

    The general order of calls should be:

    ############################################################################
    # Initialization
    ############################################################################
    sd = SequenceDataset()
    sd.read(file_obj, file_format='FASTA', datatype=datatype)

    ############################################################################
    # Check matrix
    ############################################################################
    assert sd.sequences_are_valid(remap_missing, map_missing_to)

    ############################################################################
    # read trees before changing taxa names
    ############################################################################
    tree_list = sd.dataset.read_trees(tree_f, 'NEWICK', encode_splits=True)

    ############################################################################
    # Go to safe labels
    ############################################################################
    md = MultiLocusDataset([sd])
    md.relabel_for_sate()

    ############################################################################
    # use the dataset object
    ############################################################################
    job = SateJob(multilocus_dataset=md,
                    sate_team=sate_team,
                    name=options.jobname,
                    dataset=sd.dataset
                )
    job.tree = tree_list[0]

    job.run(tmp_dir_par=temporaries_dir)

    ############################################################################
    # restore the original names to change the dataset object held by the job
    ############################################################################
    sd.restore_taxon_names()

    ############################################################################
    # get the tree with the original labels
    ############################################################################
    tree_str = job.tree.as_newick_str()
    """

    def __init__(self):
        self.dataset = None
        self.alphabet = None
        self.safe_to_real_names = {}
        self.datatype = None
        self.filename = '<unknown>'

    def get_character_matrix(self):
        """Returns the first character matrix or raises IndexError if no
        characters have been read."""
        return self.dataset.char_matrices[0]
    character_matrix = property(get_character_matrix)

    def get_taxa_block(self):
        """Returns the list of taxa."""
        return getattr(self.dataset, DATASET_TAXA_ATTR)[0]
    taxa = property(get_taxa_block)

    def read(self, file_obj, file_format='FASTA', datatype=None, filename='<unknown>'):
        """If the datatype is fasta (or some other type that does not
        specify the type of data, then the datatype arg should be DNA, RNA
        or 'PROTEIN'
        """
        self.filename = filename
        fup = file_format.upper()
        amibig_formats = ['FASTA']
        if fup in amibig_formats:
            if not datatype:
                raise ValueError('datatype must be specified when the file_format is %s' % fup)
            dup = datatype.upper()
            datatype_list = ['DNA', 'RNA', 'PROTEIN']
            if not dup in datatype_list:
                raise ValueError('Expecting the datatype to be  DNA, RNA or PROTEIN')
            file_format = dup + file_format
        try:
            import dendropy
            self.dataset = dendropy.DataSet()
            self.dataset.read(file_obj, schema=file_format, row_type='str')
            n1 = len(self.dataset.taxon_sets[0].labels())
            n2 = len(set(self.dataset.taxon_sets[0].labels()))
            if n1 != n2:
                raise ValueError("There are redundant sequence names in your data set!")
        except:
            self.dataset = None
            raise
        try:
            tb = getattr(self.dataset, DATASET_TAXA_ATTR)[0]
            self.datatype = dup
        except:
            raise ValueError("No data was read from the file.")
        tb.lock()

    def sequences_are_valid(self, remap_missing=False, map_missing_to=None):
        """Check for ? in sequences"""
        try:
            taxa_block = getattr(self.dataset, DATASET_TAXA_ATTR)[0]
            char_block = getattr(self.dataset, DATASET_CHAR_ATTR)[0]
        except:
            raise ValueError("Data have not been read")
        try:
            sa_list = char_block.state_alphabets
            self.alphabet = sa_list[0]
            missing = self.alphabet.missing
        except:
            raise ValueError("Expecting a simple datatype with one state alphabet")
        if missing is None:
            raise ValueError("Expecting a DNA, RNA, or amino acid sequences")

        for taxon in taxa_block:
            char_vec = char_block[taxon]
            missing_inds = []
            for ind, s in enumerate(char_vec):
                if s == '?':
                    missing_inds.append(ind)
            if missing_inds:
                if remap_missing:
                    as_list = list(char_vec)
                    if map_missing_to:
                        for ind in missing_inds:
                            as_list[ind] = map_missing_to
                    else:
                        missing_inds.sort(reverse=True)
                        for ind in missing_inds:
                            as_list.pop(ind)
                    char_block[taxon] = ''.join(as_list)
                else:
                    return False
        return True

class MultiLocusDataset(list):
    def __init__(self, a=()):
        list.__init__(self, a)
        self.safe_to_real_names = {}
        self.filename_list = []
        self.taxa_label_to_taxon = {}
        self.dataset = None

    def new_with_shared_meta(self):
        m =  MultiLocusDataset()
        m.safe_to_real_names = self.safe_to_real_names
        m.filename_list = self.filename_list
        m.taxa_label_to_taxon = self.taxa_label_to_taxon
        m.dataset = self.dataset
        return m

    def read_files(self,
        seq_filename_list,
        datatype,
        missing=None,
        file_format='FASTA'):
        """
        Return a MultiLocusDataset object or raises an `Exception`

            - `seq_filename_list` should a be a list of paths to FASTA-formatted sequences
            - `datatype` should  be "DNA" or "PROTEIN"
            - `missing` should be "AMBIGUOUS" to "ABSENT" indicate whether these
            missing data symbols should be treated as "any residue" or "absent"

        """
        datatype = datatype.upper()
        if datatype not in ["DNA", "PROTEIN"]:
            raise Exception("Expecting the datatype to by 'DNA' or 'PROTEIN', but found: %s\n" % datatype)

        for seq_fn in seq_filename_list:
            MESSENGER.send_info("Reading input sequences from '%s'..." % seq_fn)
            sd = SequenceDataset()
            try:
                fileobj = open(seq_fn, 'rU')
                sd.read(fileobj,
                        file_format=file_format,
                        datatype=datatype,
                        filename=seq_fn)
                _LOG.debug("sd.datatype = %s" % sd.datatype)
            except Exception, x:
                raise Exception("Error reading file:\n%s\n" % str(x))

            try:
                if not sd.sequences_are_valid(remap_missing=False):
                    m = missing.upper() if missing is not None else "ABSENT"
                    if not m in ["AMBIGUOUS", "ABSENT"]:
                        if m:
                            msg = 'The value "%s" for --missing was not understood' % m
                        raise Exception('The missing data symbol ? was encountered.\nExpecting the "missing" command-line option to be either\n "Absent" to delete ? symbols, or\n "Ambiguous" to map them to "any residue".\n%s' % msg)
                    assert(sd.alphabet)
                    map_missing_to = (m == "AMBIGUOUS" and sd.alphabet.any_residue.symbol or None)
                    if not sd.sequences_are_valid(remap_missing=True, map_missing_to=map_missing_to):
                        raise Exception("Input sequences could not be prepared for SATe.  Please report this error\n")
            except Exception, x:
                raise Exception('Error in processing file "%s":\n%s\n' % (seq_fn, str(x)))
            self.append(sd)
        self.create_dendropy_dataset()

    def _register_safe_name(self, name, locus_index, filename):
        "Creates a unique entry in safe_to_real_names for name `n` (if needed)."
        ind = 0
        real_name = name
        safe_name_prefix = "".join(_DANGEROUS_NAME_CHARS.split(name))[:80].lower()
        safe_name = safe_name_prefix
        while True:
            if safe_name not in self.safe_to_real_names:
                self.safe_to_real_names[safe_name] = (real_name, set([locus_index]))
                return safe_name
            else:
                rn, loc_ind_set = self.safe_to_real_names[safe_name]
                if real_name == rn:
                    if locus_index in loc_ind_set:
                        raise ValueError("The taxon name '%s' was repeated in the file '%s'" % (real_name, filename))
                    loc_ind_set.add(locus_index)
                    return safe_name
            ind += 1
            safe_name = safe_name_prefix + str(ind)

    def create_dendropy_dataset(self):
        from dendropy import DataSet, TaxonSet, Taxon
        taxon_set = TaxonSet()
        self.taxa_label_to_taxon = {}
        for n, element in enumerate(self):
            if not isinstance(element, SequenceDataset):
                raise ValueError("Expecting all elements of MultiLocusDataset to be SequenceDataset objects when create_dataset is called!")
            taxa_block = element.taxa
            for taxon in taxa_block:
                if taxon.label not in self.taxa_label_to_taxon:
                    nt = Taxon(label=taxon.label)
                    self.taxa_label_to_taxon[taxon.label] = nt
                    taxon_set.add_taxon(nt)
        self.dataset = DataSet()
        self.dataset.attach_taxon_set(taxon_set)
        taxon_set.lock()

    def relabel_for_sate(self):
        self.safe_to_real_names = {}
        self.filename_list = []
        alignment_list = []
        for n, element in enumerate(self):
            if not isinstance(element, SequenceDataset):
                raise ValueError("Expecting all elements of MultiLocusDataset to be SequenceDataset objects when relabel_for_sate is called!")
            try:
                taxa_block = element.taxa
                char_block = element.character_matrix
            except:
                log_exception(_LOG)
                raise
            a = Alignment()
            a.datatype = element.datatype
            fn = element.filename
            self.filename_list.append(fn)
            for taxon in taxa_block:
                char_vec = char_block[taxon]
                safe_name = self._register_safe_name(taxon.label, n, fn)
                trees_taxon = self.taxa_label_to_taxon[taxon.label]
                trees_taxon.label = safe_name
                _LOG.debug("%s (%d) -> %s" % (taxon.label, id(taxon), safe_name))
                taxon.label = safe_name
                a[safe_name] = char_vec
            alignment_list.append(a)
        # replace the contents of the MultiLocusDataset with the newly
        #   created alignment dictionaries.
        del self[:]
        for a in alignment_list:
            self.append(a)

    def concatenate_alignments(self):
        _LOG.debug('Inside concatenate_alignments')
        combined_alignment = Alignment()
        partitions = []
        base = 0
        for a in self:
            assert(a.is_aligned())
            this_el_len = a.sequence_length()
            partitions.append( a.partition_info(base) )
            for k in a.keys():
                if combined_alignment.has_key(k):
                    combined_alignment[k] += a[k]
                else:
                    combined_alignment[k] = 'N'*base + a[k]
            for i in combined_alignment.keys():
                if not a.has_key(i):
                    combined_alignment[i] += 'N'*this_el_len
            base += this_el_len

        if len(set([a.datatype for a in self])) == 1:
            combined_alignment.datatype = self[0].datatype
        else:
            combined_alignment.datatype = "MIXED"
        return (combined_alignment, partitions)

    def restore_taxon_names(self):
        """Changes the labels in the contained alignment back to their original name"""
        for alignment in self:
            new_aln = {}
            for k, v in alignment.iteritems():
                real_name = self.safe_to_real_names[k][0]
                new_aln[real_name] = v
                self.taxa_label_to_taxon[real_name].label = real_name
            keys = alignment.keys()
            for k in keys:
                del alignment[k]
            for k, v in new_aln.iteritems():
                alignment[k] = v
        self.safe_to_real_names = {}
    def sub_alignment(self, taxon_names):
        m = self.new_with_shared_meta()
        for alignment in self:
            na = Alignment()
            na.datatype = alignment.datatype
            for k in taxon_names:
                if (alignment.has_key(k)):
                    na[k] = alignment[k]
            m.append(na)
        return m
    def get_num_taxa(self):
        t = set()
        for el in self:
            t.update(set(el.keys()))
        return len(t)
    def get_num_loci(self):
        return len(self)
