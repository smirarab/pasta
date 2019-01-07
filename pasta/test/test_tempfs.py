#! /usr/bin/env python
import unittest
import os
from pasta.filemgr import TempFS
from pasta import get_logger

_LOG = get_logger(__name__)

class TempFSTest(unittest.TestCase):
    def setUp(self):
        self.ts = TempFS()
    def tearDown(self):
        dir_list = self.ts.get_remaining_directories()
        for dir in dir_list:
            self.ts.remove_dir(dir)

    def testCreate(self):
        self.assertTrue(True)

    def testNoTopLevel(self):
        self.assertRaises(ValueError, self.ts.create_subdir, 'bogus')
        cur_parent = os.path.abspath(os.curdir)
        self.assertRaises(ValueError, self.ts.create_temp_subdir, prefix='bogus', parent=cur_parent)
        self.assertEqual(self.ts.get_remaining_directories(), set())

    def testBadParForTop(self):
        fn = 'THANKS.txt'
        if os.path.exists(fn) and (not os.path.isdir(fn)):
            self.assertRaises(OSError, self.ts.create_top_level_temp, prefix='bogus', parent=fn)
        else:
            _LOG.warn("test of create_top_level_temp with file as parent skipped because '%s' does not exist" % fn)
        bogus_par = 'bogus_par'
        if os.path.exists(bogus_par):
            _LOG.warn("test of create_top_level_temp with non existent parent skipped because '%s' exists" % bogus_par)
        else:
            self.assertRaises(OSError, self.ts.create_top_level_temp, prefix='bogus', parent=bogus_par)

    def testTop(self):

        d = self.ts.create_top_level_temp(prefix='bogus 4 testing', parent=os.curdir)
        self.assertEqual(os.path.realpath(d), d)
        self.assertEqual(self.ts.top_level_temp, d)
        self.assertEqual(os.path.abspath(d), d)
        self.assertTrue(os.path.exists(d))
        self.assertTrue(os.path.isdir(d))
        self.assertTrue(d in self.ts.get_remaining_directories())

        # there can be only on top
        self.assertRaises(AssertionError, self.ts.create_top_level_temp, prefix='bogus_for_testing', parent=os.curdir)

        # subdirectories cannot be created outside of  the top...
        self.assertRaises(OSError, self.ts.create_subdir, 'bogus')
        # subdirectories cannot be created outside of  the top...
        self.assertRaises(OSError, self.ts.create_temp_subdir, prefix='bogus', parent=os.curdir)

        # but  be created inside
        ssd = os.path.join(d, 'bogussd')
        sd = self.ts.create_subdir(ssd)
        self.assertEqual(sd, ssd)
        self.assertTrue(sd in self.ts.get_remaining_directories())
        self.assertTrue(d in self.ts.get_remaining_directories())
        self.assertTrue(os.path.exists(sd))
        self.assertTrue(os.path.isdir(sd))

        nssd = os.path.join(ssd, 'nested')
        nsd = self.ts.create_subdir(nssd)
        self.assertEqual(nsd, nssd)
        self.assertTrue(nsd in self.ts.get_remaining_directories())
        self.assertTrue(d in self.ts.get_remaining_directories())
        self.assertTrue(os.path.exists(nsd))
        self.assertTrue(os.path.isdir(nsd))

        tsd = self.ts.create_temp_subdir(prefix='bogus', parent=d)
        self.assertEqual(os.path.realpath(tsd), tsd)
        self.assertEqual(os.path.abspath(tsd), tsd)
        self.assertTrue(os.path.exists(tsd))
        self.assertTrue(os.path.isdir(tsd))
        self.assertTrue(tsd in self.ts.get_remaining_directories())

        self.assertTrue(sd in self.ts.get_remaining_directories())
        self.assertTrue(d in self.ts.get_remaining_directories())
        self.assertEqual(len(self.ts.get_remaining_directories()), 4)

        # create tempdir in nested
        tnsd = self.ts.create_temp_subdir(prefix='tempinnested', parent=nsd)
        self.assertEqual(os.path.realpath(tnsd), tnsd)
        self.assertEqual(os.path.abspath(tnsd), tnsd)
        self.assertTrue(os.path.exists(tnsd))
        self.assertTrue(os.path.isdir(tnsd))
        self.assertTrue(tnsd in self.ts.get_remaining_directories())

        # subdirectories within create_temp_subdir should work...
        innermost = os.path.join(tnsd, 'innermost')
        innermostsd = self.ts.create_subdir(innermost)
        self.assertEqual(innermostsd, innermost)
        self.assertTrue(innermostsd in self.ts.get_remaining_directories())
        self.assertTrue(d in self.ts.get_remaining_directories())
        self.assertTrue(os.path.exists(innermostsd))
        self.assertTrue(os.path.isdir(innermostsd))


        self.assertRaises(ValueError, self.ts.remove_dir, 'THANKS.txt')

        self.assertEqual(self.ts.remove_dir(sd), True)
        self.assertFalse(os.path.exists(sd))
        self.assertFalse(os.path.exists(innermostsd))
        # removing sd will remove nsd (because it is inside sd), so
        #   trying to created directories in that location should fail
        self.assertRaises(OSError, self.ts.create_temp_subdir, prefix='tempinnested', parent=nsd)
        self.assertTrue(not os.path.exists(sd))
        self.assertTrue(os.path.exists(tsd))
        self.assertEqual(self.ts.remove_dir(d), True)
        self.assertFalse(os.path.exists(tsd))
        self.assertRaises(ValueError, self.ts.remove_dir, tsd)
        self.assertRaises(ValueError, self.ts.remove_dir, d)
        self.assertRaises(ValueError, self.ts.create_subdir, 'bogus')
        self.assertRaises(OSError, self.ts.create_temp_subdir, prefix='bogus', parent=d)



if __name__ == "__main__":
    unittest.main()
