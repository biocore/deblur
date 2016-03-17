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
from os import listdir
from types import GeneratorType
from os.path import join, isfile, basename, abspath, dirname, splitext

from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta
from skbio.alignment import SequenceCollection
from skbio.sequence import DNA
from bfillings.sortmerna_v2 import build_database_sortmerna
from biom.table import Table
from biom import load_table

from deblur.workflow import (dereplicate_seqs,
                             remove_chimeras_denovo_from_seqs,
                             remove_artifacts_seqs,
                             parse_deblur_output,
                             generate_biom_data,
                             generate_biom_table,
                             trim_seqs,
                             multiple_sequence_alignment,
                             launch_workflow,
                             split_sequence_file_on_sample_ids_to_files,
                             build_index_sortmerna)


class workflowTests(TestCase):
    """ Test deblur pipeline and individual methods functionality """

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

        # the data directory for the workflow test files
        self.test_data_dir = join(dirname(abspath(__file__)), 'data')
        self.seqs_s1_fp = join(self.test_data_dir, 'seqs_s1.fa')
        self.seqs_s2_fp = join(self.test_data_dir, 'seqs_s2.fa')
        self.orig_s1_fp = join(self.test_data_dir, 'simset.s1.fasta')
        self.orig_s2_fp = join(self.test_data_dir, 'simset.s2.fasta')

        self.files_to_remove = []

    def tearDown(self):
        remove_files(self.files_to_remove)
        rmtree(self.working_dir)

    def test_trim_seqs(self):
        seqs = [("seq1", "tagggcaagactccatggtatga"),
                ("seq2", "cggaggcgagatgcgtggta"),
                ("seq3", "tactagcaagattcctggtaaagga"),
                ("seq4", "aggatgcgagatgcgtg"),
                ("seq5", "gagtgcgagatgcgtggtgagg"),
                ("seq6", "ggatgcgagatgcgtggtgatt"),
                ("seq7", "agggcgagattcctagtgga--")]
        obs = trim_seqs(seqs, 20)

        self.assertTrue(isinstance(obs, GeneratorType))

        exp = [("seq1", "tagggcaagactccatggta"),
               ("seq2", "cggaggcgagatgcgtggta"),
               ("seq3", "tactagcaagattcctggta"),
               ("seq5", "gagtgcgagatgcgtggtga"),
               ("seq6", "ggatgcgagatgcgtggtga"),
               ("seq7", "agggcgagattcctagtgga")]
        self.assertEqual(list(obs), exp)

    def test_dereplicate_seqs_remove_singletons(self):
        """ Test dereplicate_seqs() method functionality with
            removing singletons
        """
        seqs = [("seq1", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq2", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq3", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq4", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCT"),
                ("seq5", "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCG"),
                ("seq6", "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT"),
                ("seq7", "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT")]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)

        output_fp = join(self.working_dir, "seqs_derep.fasta")
        log_fp = join(self.working_dir, "seqs_derep.log")

        dereplicate_seqs(seqs_fp=seqs_fp,
                         output_fp=output_fp)
        self.assertTrue(isfile(output_fp))
        self.assertTrue(isfile(log_fp))

        exp = [("seq1;size=3;",
                "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
               ("seq6;size=2;",
                "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT")]

        with open(output_fp, 'U') as out_f:
            act = [item for item in parse_fasta(out_f)]

        self.assertEqual(act, exp)

    def test_dereplicate_seqs(self):
        """ Test dereplicate_seqs() method functionality,
            keep singletons
        """
        seqs = [("seq1", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq2", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq3", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                ("seq4", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCT"),
                ("seq5", "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCG"),
                ("seq6", "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT"),
                ("seq7", "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT")]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)

        output_fp = join(self.working_dir, "seqs_derep.fasta")
        log_fp = join(self.working_dir, "seqs_derep.log")

        dereplicate_seqs(seqs_fp=seqs_fp,
                         output_fp=output_fp,
                         min_size=1)
        self.assertTrue(isfile(output_fp))
        self.assertTrue(isfile(log_fp))

        exp = [("seq1;size=3;",
                "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
               ("seq6;size=2;",
                "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT"),
               ("seq4;size=1;",
                "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCT"),
               ("seq5;size=1;",
                "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCG")]

        with open(output_fp, 'U') as out_f:
            act = [item for item in parse_fasta(out_f)]

        self.assertEqual(act, exp)

    def test_remove_artifacts_seqs(self):
        """ Test remove_artifacts_seqs() function for removing
            sequences not matching to a reference database
            using SortMeRNA. This test forces a new index
            construction for the reference sequences.
        """
        seqs = [("seq1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq2", "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
                ("seq3", "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("seq4", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq5", "CTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATAGGGTC"),
                ("seq6", "TTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAAT"),
                ("phix1", "TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCC"),
                ("phix2", "CTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGC"),
                ("phix3", "GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAA")]
        exp_seqs = ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6"]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
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
        ref_fp = join(self.working_dir, "ref2.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        ref_db_fp, files_to_remove = build_index_sortmerna(
            ref_fp=(ref_fp,),
            working_dir=self.working_dir)
        output_fp = remove_artifacts_seqs(seqs_fp=seqs_fp,
                                          ref_fp=(ref_fp,),
                                          working_dir=self.working_dir,
                                          ref_db_fp=ref_db_fp,
                                          negate=False,
                                          threads=1)
        obs_seqs = []
        with open(output_fp, 'U') as output_f:
            for label, seq in parse_fasta(output_f):
                obs_seqs.append(label)
        self.assertEqual(obs_seqs, exp_seqs)

    def test_remove_artifacts_seqs_index_prebuilt(self):
        """ Test remove_artifacts_seqs() function for removing
            sequences not matching to a reference database
            using SortMeRNA. This test passes a built index.
        """
        seqs = [("seq1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq2", "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
                ("seq3", "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("seq4", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq5", "CTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATAGGGTC"),
                ("seq6", "TTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAAT"),
                ("phix1", "TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCC"),
                ("phix2", "CTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGC"),
                ("phix3", "GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAA")]
        exp_seqs = ["seq1", "seq2", "seq3", "seq4", "seq5", "seq6"]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
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
        ref_fp = join(self.working_dir, "ref3.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        # build index
        sortmerna_db, files_to_remove = \
            build_database_sortmerna(
                fasta_path=ref_fp,
                max_pos=10000,
                output_dir=self.working_dir)
        self.files_to_remove.extend(files_to_remove)
        output_fp = join(self.working_dir, "seqs_filtered.fasta")
        output_fp = remove_artifacts_seqs(seqs_fp=seqs_fp,
                                          ref_fp=(ref_fp,),
                                          working_dir=self.working_dir,
                                          ref_db_fp=(sortmerna_db,),
                                          negate=False,
                                          threads=1)

        obs_seqs = []
        with open(output_fp, 'U') as output_f:
            for label, seq in parse_fasta(output_f):
                obs_seqs.append(label)
        self.assertEqual(obs_seqs, exp_seqs)

    def test_remove_artifacts_seqs_negate(self):
        """ Test remove_artifacts_seqs() function for removing
            sequences matching to a reference database
            using SortMeRNA (negate option).
        """
        seqs = [("seq1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq2", "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
                ("seq3", "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("seq4", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq5", "CTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATAGGGTC"),
                ("seq6", "TTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAAT"),
                ("phix1", "TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCC"),
                ("phix2", "CTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGC"),
                ("phix3", "GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAA")]
        exp_seqs = ["phix1", "phix2", "phix3"]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
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
        ref_fp = join(self.working_dir, "ref4.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        ref_db_fp, files_to_remove = build_index_sortmerna(
            ref_fp=(ref_fp,),
            working_dir=self.working_dir)
        output_fp = join(self.working_dir, "seqs_filtered.fasta")
        output_fp = remove_artifacts_seqs(seqs_fp=seqs_fp,
                                          ref_fp=(ref_fp,),
                                          working_dir=self.working_dir,
                                          ref_db_fp=ref_db_fp,
                                          negate=True,
                                          threads=1)
        obs_seqs = []
        with open(output_fp, 'U') as output_f:
            for label, seq in parse_fasta(output_f):
                obs_seqs.append(label)
        self.assertEqual(obs_seqs, exp_seqs)

    def test_remove_artifacts_seqs_mismatch_ref_index(self):
        """ Test remove_artifacts_seqs() function for removing
            sequences not matching to a reference database
            using SortMeRNA. A ValueError() should be raised
            when a user passes a reference sequence and an index
            database that do not match.
        """
        seqs = [("seq1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq2", "CCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"),
                ("seq3", "TCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCC"),
                ("seq4", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCC"),
                ("seq5", "CTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATAGGGTC"),
                ("seq6", "TTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAAT"),
                ("phix1", "TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCC"),
                ("phix2", "CTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGC"),
                ("phix3", "GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAA")]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
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
        ref_fp = join(self.working_dir, "ref5.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        ref_bis = [("ref7", "attaaatcagttatcgtttatttgatagttcctttactacatgga"
                            "tatc"),
                   ("ref8", "accttacgagaaatcaaagtctttgggttctggggggagtatggt"
                            "cgcaaggctgaaacttaaaggaattgacggaaggg"),
                   ("ref9", "aattgcgataacgaacgagaccttaacctactaaatagtgctgct"
                            "agcatttgc"),
                   ("ref10", "gacgggtgacggagaattagggttcgattccggagagggagcct"
                             "gagaaacggctaccacatccaag")]
        ref_bis_fp = join(self.working_dir, "ref_bis.fasta")
        with open(ref_bis_fp, 'w') as ref_bis_f:
            for seq in ref_bis:
                ref_bis_f.write(">%s\n%s\n" % seq)
        # build index
        sortmerna_db, files_to_remove = \
            build_database_sortmerna(
                fasta_path=ref_bis_fp,
                max_pos=10000,
                output_dir=self.working_dir)
        self.files_to_remove.extend(files_to_remove)
        self.assertRaises(ValueError,
                          remove_artifacts_seqs,
                          seqs_fp=seqs_fp,
                          ref_fp=(ref_fp,),
                          working_dir=self.working_dir,
                          ref_db_fp=(sortmerna_db,),
                          negate=False,
                          threads=1)

    def test_remove_chimeras_denovo_from_seqs(self):
        """ Test remove_chimeras_denovo_from_seqs() method functionality.
            Remove chimeric sequences from a FASTA file using the UCHIME
            algorithm, implemented in VSEARCH.
        """
        seqs = [("s1_104;size=2;", "GTGCCAGCCGCCGCGGTAATACCCGCAGCTCAAGTGGTG"
                                   "GTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTT"
                                   "GTAAATCCCTGGGTAAATCGGGAAGCTTAACTTTCCGAC"
                                   "TTCCGAGGAGACTGTCAAACTTGGGACCGGGAG"),
                ("s1_106;size=2;", "GTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTG"
                                   "TGGATGTTTATTGGGCCTAAAGCGTCCGTAGCCGGCTGC"
                                   "GCAAGTCTGTCGGGAAATCCGCACGCCTAACGTGCGGGC"
                                   "GTCCGGCGGAAACTGCGTGGCTTGGGACCGGAA"),
                ("s1_1;size=9;", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAA"
                                 "ACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT"
                                 "CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGA"
                                 "CCGGGCGAGGTTAGAGGTACTCTCGGG"),
                ("s1_20;size=9;", "TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAA"
                                  "AACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGG"
                                  "GAAGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTG"
                                  "GGACCGGGAGAGGCTAGAGGTACTTCTGGG"),
                ("s1_40;size=8;", "TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAA"
                                  "AGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCA"
                                  "CCGAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAG"
                                  "GGGACGAGAGAGGCAGACGGTATTTCCGGG"),
                ("s1_60;size=8;", "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAA"
                                  "AGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCG"
                                  "CACGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTG"
                                  "GGACCGGAAGACTCGAGGGGTACGTCAGGG")]
        seqs_non_chimera = ["s1_1;size=9;", "s1_20;size=9;",
                            "s1_40;size=8;", "s1_60;size=8;"]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
        output_fp = remove_chimeras_denovo_from_seqs(
            seqs_fp=seqs_fp,
            working_dir=self.working_dir)
        seqs_obs = []
        with open(output_fp, 'U') as output_f:
            for label, seq in parse_fasta(output_f):
                label = label.split()[0]
                seqs_obs.append(label)
        self.assertEqual(seqs_non_chimera, seqs_obs)

    def test_parse_deblur_output(self):
        """ Test parsing of deblur output into a dictionary
        """
        derep_clusters = {'s1_80': ['s1_80', 's1_81', 's1_82'],
                          's1_0': ['s1_0', 's1_1'],
                          's1_10': ['s1_10', 's1_12', 's1_13'],
                          's1_118': ['s1_118', 's1_119']}
        deblur_clusters_exp = {
            'AGTCGTACGTGCATGCA': ['s1_80', 's1_81', 's1_82'],
            'TGTGTAGCTGTGCTGAT': ['s1_0', 's1_1'],
            'CGGGTGCATGTCGTGAC': ['s1_10', 's1_12', 's1_13'],
            'TGGCTAGTGCTAGCTCC': ['s1_118', 's1_119']}
        seqs = [("s1_80;size=3;", "AGTCGTACGTGCATGCA"),
                ("s1_0;size=2;", "TGTGTAGCTGTGCTGAT"),
                ("s1_10;size=3;", "CGGGTGCATGTCGTGAC"),
                ("s1_118;size=2;", "TGGCTAGTGCTAGCTCC")]
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
        deblur_clusters = parse_deblur_output(seqs_fp, derep_clusters)
        self.assertDictEqual(deblur_clusters, deblur_clusters_exp)

    def test_generate_biom_data(self):
        """ Test generating BIOM table data for deblur output
        """
        data_exp = {(3, 0): 3, (2, 0): 2, (1, 0): 2, (0, 0): 3}
        deblur_clusters = {
            'AGTCGTACGTGCATGCA': ['s1_80', 's1_81', 's1_82'],
            'TGTGTAGCTGTGCTGAT': ['s1_0', 's1_1'],
            'CGGGTGCATGTCGTGAC': ['s1_10', 's1_12', 's1_13'],
            'TGGCTAGTGCTAGCTCC': ['s1_118', 's1_119']}
        otu_ids_exp = ['AGTCGTACGTGCATGCA', 'TGTGTAGCTGTGCTGAT',
                       'CGGGTGCATGTCGTGAC', 'TGGCTAGTGCTAGCTCC']
        sample_ids_exp = ['s1']
        data, otu_ids, sample_ids = generate_biom_data(deblur_clusters)
        self.assertDictEqual(data, data_exp)
        self.assertItemsEqual(otu_ids, otu_ids_exp)
        self.assertItemsEqual(sample_ids, sample_ids_exp)

    def test_generate_biom_table(self):
        """ Test generating BIOM table
        """
        seqs = [("s1_80;size=3;", "AGTCGTACGTGCATGCA"),
                ("s1_0;size=3;", "TGTGTAGCTGTGCTGAT"),
                ("s1_10;size=3;", "CGGGTGCATGTCGTGAC")]
        uc_output = """S\t0\t100\t*\t*\t*\t*\t*\ts1_80\t*
H\t0\t100\t100.0\t*\t0\t0\t*\ts1_81\ts1_80
H\t0\t100\t100.0\t*\t0\t0\t*\ts1_82\ts1_80
S\t1\t100\t*\t*\t*\t*\t*\ts1_0\t*
H\t1\t100\t100.0\t*\t0\t0\t*\ts1_1\ts1_0
H\t1\t100\t100.0\t*\t0\t0\t*\ts1_60\ts1_0
S\t2\t100\t*\t*\t*\t*\t*\ts1_10\t*
H\t2\t100\t100.0\t*\t0\t0\t*\ts1_12\ts1_10
H\t2\t100\t100.0\t*\t0\t0\t*\ts1_13\ts1_10
"""
        data = {(2, 0): 3, (1, 0): 3, (0, 0): 3}
        otu_ids = ['CGGGTGCATGTCGTGAC', 'TGTGTAGCTGTGCTGAT',
                   'AGTCGTACGTGCATGCA']
        sample_ids = ['s1']
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for seq in seqs:
                seqs_f.write(">%s\n%s\n" % seq)
        # temporary file for .uc output
        uc_output_fp = join(self.working_dir, "derep.uc")
        with open(uc_output_fp, 'w') as uc_output_f:
            uc_output_f.write(uc_output)
        table_exp = Table(data, otu_ids, sample_ids, sample_metadata=None)
        clusters, table = generate_biom_table(seqs_fp,
                                              uc_output_fp)
        self.assertEqual(table, table_exp)

    def test_multiple_sequence_alignment(self):
        """Test multiple sequence alignment.
        """
        seqs = [DNA('caccggcggcccggtggtggccattattattgggtctaaag', id='seq_1'),
                DNA('caccggcggcccgagtggtggccattattattgggtcaagg', id='seq_2'),
                DNA('caccggcggcccgagtgatggccattattattgggtctaaag', id='seq_3'),
                DNA('aaccggcggcccaagtggtggccattattattgggtctaaag', id='seq_4'),
                DNA('caccgggcccgagtggtggccattattattgggtctaaag', id='seq_5')]
        seqs_col = SequenceCollection(seqs)
        seqs_fp = join(self.working_dir, "seqs.fna")
        with open(seqs_fp, 'w') as o:
            o.write(seqs_col.to_fasta())
        alignment = multiple_sequence_alignment(seqs_fp)
        align_exp = [
            DNA(
                'caccggcggcccg-gtggtggccattattattgggtctaaag', id='seq_1'),
            DNA(
                'caccggcggcccgagtggtggccattattattgggtcaagg-', id='seq_2'),
            DNA(
                'caccggcggcccgagtgatggccattattattgggtctaaag', id='seq_3'),
            DNA(
                'aaccggcggcccaagtggtggccattattattgggtctaaag', id='seq_4'),
            DNA(
                'caccg--ggcccgagtggtggccattattattgggtctaaag', id='seq_5')]
        self.assertItemsEqual(alignment, align_exp)

    def test_build_index_sortmerna(self):
        """Test functionality of build_index_sortmerna()
        """
        ref1 = [("ref1", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTA"
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
        ref2 = [("ref1", "GTCGTAGCTAGCTGCCCACGATCGTAGCTAGCTAGCTACGTAGCTCATCAC"
                 "TCGCCGACCCACGTCCCACTGATGCTGTGGG"),
                ("ref2", "GCGGCGCCCAAAAATGTCGTGTAAAATTTTCTCGTACCCACTTGCTACCCA"
                 "TGGCCGCCATGCTGCTAACGCAATATATATA"),
                ("ref3", "TGTGAAAGCGCGCGAGAGAGTCGTATATATGGGCGCGGCGCGATGCTGCCC"
                 "GTCGATGCTGATCCCCCACGTACGTAGCCCC"),
                ("ref4", "GTGTGCTCGCGTAGCTAGCTTATATATCGGCGCGTAGTGCTAGCCCCAAAA"
                 "GTGTCCCCCCCCTCCTTTTTTATATATGCAA"),
                ("ref5", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                 "TCGGCTTTGTAAATCCCTGGGTAAATAGGGT"),
                ("ref6", "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAG"
                 "TCGGCTTTGTAAATCCCTGGGTAAATCGGGT")]
        ref1_fp = join(self.working_dir, "ref1.fasta")
        with open(ref1_fp, 'w') as ref_f:
            for seq in ref1:
                ref_f.write(">%s\n%s\n" % seq)
        ref2_fp = join(self.working_dir, "ref2.fasta")
        with open(ref2_fp, 'w') as ref_f:
            for seq in ref2:
                ref_f.write(">%s\n%s\n" % seq)
        ref_fps = tuple([ref1_fp, ref2_fp])
        ref_db_fp, files_to_remove = build_index_sortmerna(
            ref_fp=ref_fps,
            working_dir=self.working_dir)
        files_to_remove_exp = ['ref1.stats', 'ref1.pos_0.dat',
                               'ref1.kmer_0.dat', 'ref1.bursttrie_0.dat',
                               'ref2.stats', 'ref2.pos_0.dat',
                               'ref2.kmer_0.dat', 'ref2.bursttrie_0.dat']
        files_to_remove_act = [basename(f) for f in files_to_remove]
        self.assertListEqual(files_to_remove_exp, files_to_remove_act)

    def run_workflow_try(self, simfilename, origfilename, ref_fp, ref_db_fp):
        """Test launching the complete workflow using simulated sequences
        and compare to original ground truth.

        input:
        simfilename : string
            name of the simulated reads fasta file
        origfilename : string
            name of the fasta file with the ground truth sequences
        """
        seqs_fp = simfilename
        output_fp = self.working_dir
        read_error = 0.05
        mean_error = 0.005
        error_dist = None
        indel_prob = 0.01
        indel_max = 3
        trim_length = 100
        min_size = 2
        negate = False
        threads = 1
        delim = '_'
        biom_fp = launch_workflow(seqs_fp, output_fp, read_error, mean_error,
                                  error_dist, indel_prob, indel_max,
                                  trim_length, min_size, (ref_fp,),
                                  ref_db_fp, negate, threads, delim)

        # get the trimmed ground truth sequences
        with open(origfilename, 'U') as f:
            orig_seqs = [item[1] for item in parse_fasta(f)]
        orig_seqs = [item[:trim_length].upper() for item in orig_seqs]

        table_obs = load_table(biom_fp)
        outseqs = table_obs.ids(axis='observation')

        # test we see all ground truth sequences and no other
        self.assertItemsEqual(outseqs, orig_seqs)

    def test_launch_workflow(self):
        """Test launching complete workflow using 2 simulated sequence files.
        seqs1 - 100 reads using art, original sequences are >0.5 identical.
        seqs2 - 200 reads using grinder, original sequences are >0.9 identical,
        0.1 chimeras, 35 phix reads
        """
        # index the 70% rep. set database
        ref_fp = join(self.test_data_dir, '70_otus.fasta')
        ref_db_fp, files_to_remove = build_index_sortmerna(
            ref_fp=(ref_fp,),
            working_dir=self.working_dir)

        self.run_workflow_try(self.seqs_s1_fp,
                              self.orig_s1_fp, ref_fp, ref_db_fp)
        self.run_workflow_try(self.seqs_s2_fp,
                              self.orig_s2_fp, ref_fp, ref_db_fp)

    def get_seqs_act_split_sequence_on_sample_ids(self, output_dir):
        """Parse output of split_sequence_file_on_sample_ids_to_files()

        Parameters
        ----------
        output_dir: string
            output directory path storing FASTA files

        Returns
        -------
        seqs_act: dict
            dictionary with keys being sample IDs and values list of
            sequences belonging to sample ID
        """
        seqs_act = {}
        for fn in listdir(output_dir):
            input_fp = join(output_dir, fn)
            sample_file = splitext(fn)[0]
            with open(input_fp, 'U') as input_f:
                for label, seq in parse_fasta(input_f):
                    sample = label.split('_')[0]
                    self.assertEqual(sample_file, sample)
                    if sample not in seqs_act:
                        seqs_act[sample] = [(label, seq)]
                    else:
                        seqs_act[sample].append((label, seq))
        return seqs_act

    def test_split_sequence_file_on_sample_ids_to_files(self):
        """Test functionality of split_sequence_file_on_sample_ids_to_files()
        """
        seqs_fasta = {"s1": [
                      ("s1_seq1",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                      ("s1_seq2",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG")],
                      "s2": [
                      ("s2_seq3",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG"),
                      ("s2_seq4",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCT")],
                      "s3": [
                      ("s3_seq5",
                       "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCG"),
                      ("s3_seq6",
                       "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT")],
                      "s4": [
                      ("s4_seq7",
                       "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT")]}
        # Test FASTA file split on sample IDs to multiple FASTA files
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for sample in seqs_fasta:
                for seq in seqs_fasta[sample]:
                    seqs_f.write(">%s\n%s\n" % seq)
        output_dir = mkdtemp()
        with open(seqs_fp, 'U') as seqs_f:
            split_sequence_file_on_sample_ids_to_files(seqs=seqs_f,
                                                       outdir=output_dir)
        seqs_act = self.get_seqs_act_split_sequence_on_sample_ids(
            output_dir=output_dir)
        self.assertDictEqual(seqs_fasta, seqs_act)


if __name__ == '__main__':
    main()
