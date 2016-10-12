# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkdtemp
from os.path import join, dirname, abspath
import logging
from os import makedirs, remove

from deblur.parallel_deblur import parallel_deblur
from deblur.workflow import (build_index_sortmerna,
                             start_log, sequence_generator)


class parallelDeblurTests(TestCase):
    """Test parallel deblur methods functionality.
    """

    def setUp(self):
        """ Create working directory and two FASTA input
            files corresponding to two samples (s1 and s2).
            Each input file contains 120 sequences, of which
            100 are 16S amplicons (ART simulator), 10 are
            chimeric sequences (Grinder) and 10 are PhiX
            artifacts (ART). The 100 amplicon sequences
            intend to evenly represent a community of 10
            species.
        """
        # test output can be written to this directory
        self.working_dir = mkdtemp()
        deblur_working_dir = join(self.working_dir, "deblur_working_dir")
        makedirs(deblur_working_dir)

        # the data directory for the workflow test files
        self.test_data_dir = join(dirname(abspath(__file__)), 'data')
        self.seqs_s1_fp = join(self.test_data_dir, 'seqs_s1.fasta')
        self.seqs_s2_fp = join(self.test_data_dir, 'seqs_s2.fasta')
        self.seqs_s3_fp = join(self.test_data_dir, 'seqs_s3.fasta')
        self.orig_s1_fp = join(self.test_data_dir, 'simset.s1.fasta')
        self.orig_s2_fp = join(self.test_data_dir, 'simset.s2.fasta')
        self.orig_s3_fp = join(self.test_data_dir, 'simset.s3.fasta')

        self.files_to_remove = []

        logfilename = join(self.working_dir, "log.txt")
        start_log(level=logging.DEBUG, filename=logfilename)

    def tearDown(self):
        for f in self.files_to_remove:
            remove(f)
        rmtree(self.working_dir)

    def compare_result(self, simfilename, origfilename, trim_length):
        """Compare the results of deblurring to the original mixture

        Parameters
        ----------
        simfilename : str
            name of the simulated reads fasta file
        origfilename : str
            name of the fasta file with the ground truth sequences
        trim_length : int
            the length used for trimming the sequences in deblurring
        """

        # get the trimmed ground truth sequences
        orig_seqs = [item[1] for item in sequence_generator(origfilename)]
        orig_seqs = [item[:trim_length].upper() for item in orig_seqs]

        # get the deblurred fasta file sequences
        out_seqs = [item[1] for item in sequence_generator(simfilename)]
        out_seqs = [item.upper() for item in out_seqs]

        out_seqs.sort()
        orig_seqs.sort()

        # test we see all ground truth sequences and no other
        self.assertEqual(out_seqs, orig_seqs)

    def test_parallel_deblur(self):
        """Test parallel deblur using 3 simulated sequence files.
        seqs1 - 100 reads using art, original sequences are >0.5 identical.
        seqs2 - 200 reads using grinder, original sequences are >0.9 identical,
        0.1 chimeras, 35 phix reads
        seqs3 - simple - 15 reads from seqs1 (10 reads for 1001203,
        5 reads for 694276) for manual test validation
        """
        # index the 70% rep. set database
        ref_fp = join(self.test_data_dir, '70_otus.fasta')
        ref_db_fp = build_index_sortmerna(
            ref_fp=(ref_fp,),
            working_dir=self.working_dir)

        trim_length = 100
        params = ['deblur', 'workflow', '--seqs-fp', 'ignorethis',
                  '--output-dir', self.working_dir, '--pos-ref-fp', ref_fp,
                  '-d', '1,0.06,0.02,0.02,0.01,0.005,0.005,0.005,0.001'
                  ',0.001,0.001,0.0005',
                  '-t', str(trim_length)]
        parallel_deblur([self.seqs_s1_fp, self.seqs_s2_fp, self.seqs_s3_fp],
                        params, ref_db_fp, None, jobs_to_start=2)

        deblur_working_dir = join(self.working_dir, "deblur_working_dir")

        deb1res = join(deblur_working_dir,
                       'seqs_s1.fasta.trim.derep.no_artifacts'
                       '.msa.deblur.no_chimeras')
        self.compare_result(deb1res, self.orig_s1_fp, trim_length=trim_length)

        deb2res = join(deblur_working_dir,
                       'seqs_s2.fasta.trim.derep.no_artifacts'
                       '.msa.deblur.no_chimeras')
        self.compare_result(deb2res, self.orig_s2_fp, trim_length=trim_length)

        deb3res = join(deblur_working_dir,
                       'seqs_s3.fasta.trim.derep.no_artifacts'
                       '.msa.deblur.no_chimeras')
        self.compare_result(deb3res, self.orig_s3_fp, trim_length=trim_length)

if __name__ == '__main__':
    main()
