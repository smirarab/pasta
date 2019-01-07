#! /usr/bin/env python

import unittest
import os
import sys
import tempfile
import random
from pasta import filemgr
from pasta import get_logger
from pasta.configure import get_configuration

_LOG = get_logger(__name__)

class SateProductsTest(unittest.TestCase):

    def setUp(self):
        self.top_dir = tempfile.mkdtemp()
        self.output_prefixes = []
        self.product_results = []

    def tearDown(self):
        pass

    def create_input_files(self,
            job_subdir,
            input_subdir=None):
        src_paths = []
        seq_dir = os.path.join(self.top_dir, job_subdir)
        if input_subdir is not None:
            seq_dir = os.path.join(seq_dir, input_subdir)
        for i in range(5):
            fp = os.path.join(seq_dir, "data%d.fasta" % (i+1))
            src_paths.append(fp)
            f = filemgr.open_with_intermediates(fp, "w")
            f.close()
        return src_paths

    def generate_random_result(self):
        return "%d%d%d\n" % (random.randint(0, sys.maxsize),
                random.randint(0, sys.maxsize),
                random.randint(0, sys.maxsize))

    def create_and_verify(self,
            job_name="satejob",
            input_subdir=None,
            output_subdir=None,
            expected_index=None):

        ## create directories and files

        # job subdirectory
        job_subdir = "test-%s" % job_name

        # basic set of input sequences
        input_seq_filepaths = self.create_input_files(job_subdir=job_subdir,
                input_subdir=input_subdir)

        # check if we can handle redundant input files without overwriting output
        input_seq_filepaths.extend(list(input_seq_filepaths))

        # output directory
        if output_subdir is not None:
            output_dir = os.path.join(self.top_dir, job_subdir, output_subdir)
            expected_output_dir = output_dir
        else:
            output_dir = None
            expected_output_dir = os.path.dirname(input_seq_filepaths[0])

        ## create the product manager
        user_config = get_configuration()
        user_config.input_seq_filepaths = input_seq_filepaths
        user_config.commandline.input = input_seq_filepaths[0]
        sp = filemgr.PastaProducts(sate_user_settings=user_config)

        ## job prefix: must be unique
        output_prefix = sp.output_prefix
        self.assertTrue(output_prefix not in self.output_prefixes)
        self.output_prefixes.append(output_prefix)

        ## meta products (score, tree, and log files)
        self.assertTrue(hasattr(sp, "score_stream"))
        self.assertTrue(hasattr(sp, "tree_stream"))
        self.assertTrue(hasattr(sp, "run_log_stream"))
        self.assertTrue(hasattr(sp, "err_log_stream"))
        for stream_name, product_extension in list(filemgr.PastaProducts.meta_product_types.items()):
            expected_fn = output_prefix + product_extension
            self.assertTrue(os.path.exists(expected_fn))
            stream_attr_name = stream_name + "_stream"
            self.assertTrue(hasattr(sp, stream_attr_name))
            stream = getattr(sp, stream_attr_name)
            self.assertEqual(
                    os.path.abspath(stream.name),
                    os.path.abspath(expected_fn))
            random_result = self.generate_random_result()
            self.product_results.append((expected_fn, random_result,))
            stream.write(random_result)
            stream.flush()
            stream.close()

        ## final alignment files
        self.assertEqual(len(sp.alignment_streams), len(input_seq_filepaths))
        align_fnames = []
        for alignment_stream in sp.alignment_streams:
            fn = os.path.abspath(alignment_stream.name)
            self.assertTrue(os.path.exists(fn))
            align_fnames.append(fn)
            random_result = self.generate_random_result()
            self.product_results.append((os.path.abspath(alignment_stream.name), random_result,))
            alignment_stream.write(random_result)
            alignment_stream.flush()
            alignment_stream.close()
        self.assertEqual(len(set(align_fnames)), len(align_fnames))

        ## return sp, for futher tests if needed
        return sp

    def testDefault(self):
        sate_products = []
        for output_subdir in [None, "out"]:
            for expected_index in [None, 1, 2, 3]:
                sp = self.create_and_verify(job_name="satejob",
                        input_subdir=None,
                        output_subdir=output_subdir,
                        expected_index=expected_index)
                sate_products.append(sp)

        ## meta-level testing: ensure uniqueness across replicates
        self.assertEqual(len(set(self.output_prefixes)), len(self.output_prefixes))
        align_fnames = []
        for sp in sate_products:
            align_fnames.extend([os.path.abspath(a.name) for a in sp.alignment_streams])
        self.assertEqual(len(set(align_fnames)), len(align_fnames))

        ## meta-level testing: test results read == test results written
        for fpath, expected_result in self.product_results:
            f = open(fpath, "rU")
            actual_result = f.read()
            self.assertEqual(actual_result, expected_result)

if __name__ == "__main__":
    unittest.main()
