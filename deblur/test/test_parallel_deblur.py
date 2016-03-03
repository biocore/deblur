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
from os.path import join
from glob import glob

from bfillings.sortmerna_v2 import build_database_sortmerna

from deblur.parallel_deblur import parallel_deblur


class parallelDeblurTests(TestCase):
    """Test parallel deblur methods functionality.
    """
    def test_parallel_deblur(self):
        """Test parallel_deblur() functionality.
        """
        working_dir = mkdtemp()
        seqs = {"s1": [
                ("s1_seq1",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s1_seq2",
                 "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
                ("s1_seq3",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s1_seq4",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s1_seq5",
                 "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT")],
                "s2": [
                ("s2_seq6",
                 "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("s2_seq7",
                 "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("s2_seq8",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s2_seq9",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s2_seq10",
                 "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("s2_phix1",
                 "TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCC"),
                ("s2_phix2",
                 "CTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGC"),
                ("s2_phix3",
                 "GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAA")]}
        outputs = {}
        split_dir = mkdtemp()
        print "split_dir = ", split_dir
        for sample in seqs:
            if sample not in outputs:
                outputs[sample] = open(join(split_dir, sample + ".fa"), 'w')
            for seq in seqs[sample]:
                outputs[sample].write(">%s\n%s\n" % seq)
        inputs = glob('%s/*' % split_dir)
        ref = [("ref1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTA"
                        "GTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
               ("ref2", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                        "TCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
               ("ref3", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                        "TCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
               ("ref4", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                        "TCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
               ("ref5", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                        "TCGGCTTTGTAAATCCCTGGGTAAATAGGGT"),
               ("ref6", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                        "TCGGCTTTGTAAATCCCTGGGTAAATCGGGT")]
        ref_fp = join(working_dir, "ref.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        # build index
        ref_db_fp, files_to_remove = \
            build_database_sortmerna(
                fasta_path=ref_fp,
                max_pos=10000,
                output_dir=working_dir)
        params = {}
        params['output-dir'] = working_dir
        params['ref-db-fp'] = tuple([ref_db_fp])
        params['ref-fp'] = tuple([ref_fp])
        params['file-type'] = 'fasta'
        params['read-error'] = 0.05
        params['mean-error'] = None
        params['error-dist'] = None
        params['indel-prob'] = 0.01
        params['indel-max'] = 3
        params['trim-length'] = 40
        params['min-size'] = 2
        params['negate'] = False
        params['delim'] = '_'
        all_result_paths = parallel_deblur(inputs=inputs,
                                           params=params,
                                           jobs_to_start=1)
        rmtree(working_dir)


if __name__ == '__main__':
    main()
