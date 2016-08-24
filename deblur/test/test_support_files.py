# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main

from deblur.support_files import pos_db, neg_db

import os.path


class supportFilesTests(TestCase):
    """Test the supporting data files
    """
    def test_reference_(self):
        """Test if the positive and negative filtering
        reference fasta files exist
        """
        # the positive filtering fasta file
        self.assertTrue(os.path.isfile(pos_db))
        # the negative filtering fasta file
        self.assertTrue(os.path.isfile(neg_db))

if __name__ == '__main__':
    main()
