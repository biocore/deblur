# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from StringIO import StringIO

from skbio.parse.sequences import parse_fasta
import numpy as np


class workflowTests(TestCase):
	def SetUp(self):
		pass

	def test_demultiplex_seqs(self):
		pass

	def test_split_seqs_on_sample_ids(self):
		pass

	def test_trim_seqs(self):
		pass

	def test_dereplicate_seqs(self):
		pass

	def test_remove_singletons_seqs(self):
		pass

	def test_remove_artifacts_seqs(self):
		pass

	def test_multiple_sequence_alignment(self):
		pass

	def test_remove_chimeras_denovo_from_seqs(self):
		pass

	def test_generate_biom_table(self):
		pass

	def test_assign_taxonomy(self):
		pass

	def test_launch_workflow(self):
		pass


if __name__ == '__main__':
    main()
