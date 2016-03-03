# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from os import close, listdir
from types import GeneratorType
from os.path import join, isfile, splitext

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
                             split_sequence_file_on_sample_ids_to_files)


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

        # create temp FASTA file for s1 sequences
        f, self.seqs_s1_fp = mkstemp(prefix='seqs_s1_',
                                     suffix='.fasta')
        close(f)
        with open(self.seqs_s1_fp, 'w') as t:
            t.write(seqs_s1)

        # create temp FASTA file for s2 sequences
        f, self.seqs_s2_fp = mkstemp(prefix='seqs_s2_',
                                     suffix='.fasta')
        close(f)
        with open(self.seqs_s2_fp, 'w') as t:
            t.write(seqs_s2)

        self.files_to_remove = [self.seqs_s1_fp,
                                self.seqs_s2_fp]

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
        ref_fp = join(self.working_dir, "ref.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        output_fp = remove_artifacts_seqs(seqs_fp=seqs_fp,
                                          ref_fp=(ref_fp,),
                                          working_dir=self.working_dir,
                                          ref_db_fp=None,
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
        ref_fp = join(self.working_dir, "ref.fasta")
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
        ref_fp = join(self.working_dir, "ref.fasta")
        with open(ref_fp, 'w') as ref_f:
            for seq in ref:
                ref_f.write(">%s\n%s\n" % seq)
        output_fp = join(self.working_dir, "seqs_filtered.fasta")
        output_fp = remove_artifacts_seqs(seqs_fp=seqs_fp,
                                          ref_fp=(ref_fp,),
                                          working_dir=self.working_dir,
                                          ref_db_fp=None,
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
        ref_fp = join(self.working_dir, "ref.fasta")
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

    def test_launch_workflow(self):
        """Test launching complete workflow using simulated sequences 1.
        """
        seqs_fp = self.seqs_s1_fp
        output_dir = self.working_dir
        read_error = 0.05
        mean_error = 0.05
        error_dist = None
        indel_prob = 0.01
        indel_max = 3
        trim_length = 100
        min_size = 2
        ref_fp = join(self.working_dir, "db.fasta")
        with open(ref_fp, 'w') as ref_f:
            ref_f.write(database_16S)
        ref_db_fp = None
        negate = False
        threads = 1
        delim = '_'
        output_fp = launch_workflow(seqs_fp, output_dir, read_error,
                                    mean_error, error_dist, indel_prob,
                                    indel_max, trim_length,
                                    min_size, tuple([ref_fp]), ref_db_fp,
                                    negate, threads, delim)
        data = {(0, 0): 9, (1, 0): 9, (2, 0): 10, (3, 0): 8, (4, 0): 2,
                (5, 0): 2, (6, 0): 8, (7, 0): 2, (8, 0): 8, (9, 0): 9,
                (10, 0): 7, (11, 0): 9}
        otu_ids = ["TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGG"
                   "TTTGATAAATCCTTGGGTAAATCGGGAAGCTTAACTTTCCGATTC",
                   "CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGG"
                   "TTTGATCAGTCTTCCGGGAAATCTGACAGCTCAACTGTTAGGTCT",
                   "CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGG"
                   "ACTGGTAAGTCCCTTGTGAAATCGGTCGGCTCAACCGTTCGGTGC",
                   "TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGG"
                   "CCTTGTAAGTTCTCCCTTAAATGGTATGGCTCAACCATACCGTGG",
                   "TTCCGGTCCCAAGCCGAACAGTTTCCGCCGGACGCCCATCCGTTGAGCGGATGGA"
                   "TTTCCCGACAGACTTGTTCGGCCAGCTACGGACGCTTTAGGCCCA",
                   "GTGCCAGCCGCCGCGGTAATACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGC"
                   "CTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGG",
                   "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGG"
                   "CTTAGTGAGTCTTTCGTTAAATCCGGCGACTTAATCGTCGTTTGC",
                   "CGCCCGGTCCCAAGCTTTGCAGTCTCCCCAGAAATCGGATCGTTAAGCGACCCGA"
                   "TTTACCCAGGGATTTACAAAGCCGACTACGGACGTTTTAGGCTCA",
                   "TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGG"
                   "CCGAACAAGTCTGTCGGGAAATCCATCCGCTCAACGGATGGGTCC",
                   "TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGG"
                   "CTCAACAAGTCTCCCGTTAAATCCAGCGATCTAATCGTTGGATGC",
                   "TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGG"
                   "CTAGGTTAGTCCCCTGTTAAATCCACCGAATTAATCGTTGGATGC",
                   "TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGG"
                   "CTTTGTAAATCCCTGGGTAAATCGGGTCGCTTAACGATCCGATTC"]
        sample_ids = ["s1"]
        table_exp = Table(data, otu_ids, sample_ids, sample_metadata=None)
        table_obs = load_table(output_fp)
        self.assertEqual(table_obs, table_exp)

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
        seqs_fastq = {"s1": [
                      ("s1_seq1",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG",
                       "ZZ[P[]^^PHGOLZ]_^^^N\^^^^TR\^\]^^^^^^^Z^^^[^BBBB"),
                      ("s1_seq2",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG",
                       "^__eceee^WYecgghhhhHbfhhU^afhgfgfhhW`_eghgghfheg")],
                      "s2": [
                      ("s2_seq3",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCG",
                       "___eeVccceeeefhdcghhfhhhhgfhhegfhihiicgh]ffiihii"),
                      ("s2_seq4",
                       "TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCT",
                       "^^^cYacc\JQbZddhhhhbeJbeddWbcdhhhhhcYcchhhX`_`cc")],
                      "s3": [
                      ("s3_seq5",
                       "TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCG",
                       "___VacaceSbeehfhhhhYbghhhhfhhhhhhgh^ecadfhaghhhh"),
                      ("s3_seq6",
                       "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT",
                       "__beeeeeecegghihhfhhhiiiifgfhfgiifihfgfh`fdhif]c")],
                      "s4": [
                      ("s4_seq7",
                       "CTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGGGAGAAT",
                       "aaaeceeeeggggiiiiiihfhiiidghhegfhfihihhihhiiiige")]}
        # Test FASTA file split on sample IDs to multiple FASTA files
        seqs_fp = join(self.working_dir, "seqs.fasta")
        with open(seqs_fp, 'w') as seqs_f:
            for sample in seqs_fasta:
                for seq in seqs_fasta[sample]:
                    seqs_f.write(">%s\n%s\n" % seq)
        output_dir = mkdtemp()
        with open(seqs_fp, 'U') as seqs_f:
            split_sequence_file_on_sample_ids_to_files(seqs=seqs_f,
                                                       filetype='fasta',
                                                       outdir=output_dir)
        seqs_act = self.get_seqs_act_split_sequence_on_sample_ids(
            output_dir=output_dir)
        self.assertDictEqual(seqs_fasta, seqs_act)
        # Test FASTQ file split on multiple sample IDs to multiple FASTA files
        seqs_fp = join(self.working_dir, "seqs.fastq")
        with open(seqs_fp, 'w') as seqs_f:
            for sample in seqs_fastq:
                for seq in seqs_fastq[sample]:
                    seqs_f.write("@%s\n%s\n+\n%s\n" % seq)
        output_dir = mkdtemp()
        with open(seqs_fp, 'U') as seqs_f:
            split_sequence_file_on_sample_ids_to_files(seqs=seqs_f,
                                                       filetype='fastq',
                                                       outdir=output_dir)
        seqs_act = self.get_seqs_act_split_sequence_on_sample_ids(
            output_dir=output_dir)
        self.assertDictEqual(seqs_fasta, seqs_act)


# 2 samples, each representing 10 species, 8/10 species
#   common between both samples
# Greengenes IDs species s1 = [813079, 630311, 565570, 565152,
#                              518927, 330533, 241952, 198216,
#                              191626, 187525]
# Greengenes IDs species s2 = [4483264, 4483458, 565570, 565152,
#                              518927, 330533, 241952, 198216,
#                              191626, 187525]
# 200 simulated reads capturing 12 species
# ((8 common species * 20 reads/species) + \
#  (4 unique species * 10 reads/species))
# art_illumina -amp -sam -i s1_10_species_V4.fna -l 150 \
#   -f 10 -o s1_art_even
# art_illumina -amp -sam -i s2_10_species_V4.fna -l 150 \
#   -f 10 -o s2_art_even
# 20 (10/sample) chimeric sequences belonging to any of the
#   12 species: using Grinder
# 20 (10/sample) phiX/mRNA reads as artifacts
# art_illumina -sam -i PhiX.fasta -l 150 -o phix_sim -f 20
# split_libraries_fastq.py -i s1_art_even.fq,s2_art_even.fq \
#   -o split_libraries_fastq --sample_ids s1,s2 --barcode_type='not-barcoded' \
#   --phred_offset 33
# split_sequence_file_on_sample_ids.py -i split_libraries_fastq/seqs.fna \
#  -o split_sequence_file_on_sample_ids

seqs_s1 = """>s1_0 813079-10
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAGAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_1 813079-9
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_2 813079-8
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_3 813079-7
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_4 813079-6
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATAGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_5 813079-5
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_6 813079-4
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCGAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_7 813079-3
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_8 813079-2
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_9 813079-1
TACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAAATCCCTGGGTAAATCGGGT
CGCTTAACGATCCGATTCTGGGGAGACTGCAAAGCTTGGGACCGGGCGAGGTTAGAGGTACTCTCGGG
>s1_10 630311-10
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_11 630311-9
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTCTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_12 630311-8
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_13 630311-7
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_14 630311-6
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_15 630311-5
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_16 630311-4
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_17 630311-3
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_18 630311-2
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_19 630311-1
CACCGGCGGCCCGAGTGGTGGCCATTATTATTGGGTCTAAAGGGTCCGTAGCCGGTTTGATCAGTCTTCCGGGAAATCTGAC
AGCTCAACTGTTAGGTCTGGGGGATACTGTCAGACTTGGAACCGGGAGAGGTAAGAGGTACTACAGGG
>s1_20 565570-10
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_21 565570-9
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_22 565570-8
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_23 565570-7
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_24 565570-6
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_25 565570-5
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGAAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_26 565570-4
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_27 565570-3
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_28 565570-2
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGATTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_29 565570-1
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s1_30 565152-10
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_31 565152-9
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAAGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_32 565152-8
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_33 565152-7
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTAGTCCGGGG
>s1_34 565152-6
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAGAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_35 565152-5
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGTGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_36 565152-4
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_37 565152-3
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_38 565152-2
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_39 565152-1
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s1_40 518927-10
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_41 518927-9
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTACATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_42 518927-8
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_43 518927-7
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_44 518927-6
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCAGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_45 518927-5
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_46 518927-4
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_47 518927-3
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_48 518927-2
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTGGGTTAGTCCCCTGTTAAATCCAAC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_49 518927-1
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s1_50 330533-10
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_51 330533-9
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_52 330533-8
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_53 330533-7
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_54 330533-6
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTAGGTCTGGG
>s1_55 330533-5
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGTG
>s1_56 330533-4
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCATAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_57 330533-3
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGCAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_58 330533-2
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_59 330533-1
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s1_60 241952-10
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCAGGG
>s1_61 241952-9
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_62 241952-8
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_63 241952-7
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTATCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGATTGGGACCAGAAGACTCGAGGGGTACGTCTGGG
>s1_64 241952-6
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_65 241952-5
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_66 241952-4
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGCCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_67 241952-3
TACCGGCAGCTCAAGTGATGACCGCTATTATTCGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_68 241952-2
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_69 241952-1
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s1_70 198216-10
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGTAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_71 198216-9
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_72 198216-8
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_73 198216-7
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_74 198216-6
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACGAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_75 198216-5
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_76 198216-4
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_77 198216-3
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATAATGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_78 198216-2
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_79 198216-1
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s1_80 191626-10
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATTGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_81 191626-9
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_82 191626-8
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_83 191626-7
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_84 191626-6
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_85 191626-5
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_86 191626-4
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_87 191626-3
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_88 191626-2
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_89 191626-1
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s1_90 187525-10
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTTGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_91 187525-9
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_92 187525-8
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_93 187525-7
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCAAACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_94 187525-6
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_95 187525-5
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_96 187525-4
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_97 187525-3
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTTAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_98 187525-2
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_99 187525-1
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s1_100 chimera1-1 ref=518927,330533 amp=479..770,449..741 pos=compl(1..150)
TTCCGGTCCCAAGCCGAACAGTTTCCGCCGGACGCCCATCCGTTGAGCGGATGGATTTCCCGACAGACTTGTTCGGCCAGCT
ACGGACGCTTTAGGCCCAATAAACATCCACACCACTCGGAGAGCTGGTATTACCGCGGCGGCTGACAC
>s1_101 chimera1-2 ref=518927,330533,241952 pos=compl(1..150)
TTCCGGTCCCAAGCCGAACAGTTTCCGCCGGACGCCCATCCGTTGAGCGGATGGATTTCCCGACAGACTTGTTCGGCCAGCT
ACGGACGCTTTAGGCCCAATAAACATCCACACCACTCGGAGAGCTGGTATTACCGCGGCGGCTGACAC
>s1_102 chimera2-1 ref=330533,813079 amp=449..741,448..739 pos=compl(1..150)
CGCCCGGTCCCAAGCTTTGCAGTCTCCCCAGAAATCGGATCGTTAAGCGACCCGATTTACCCAGGGATTTACAAAGCCGACT
ACGGACGTTTTAGGCTCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s1_103 chimera2-2 ref=330533,813079 amp=449..741,448..739 pos=compl(1..150)
CGCCCGGTCCCAAGCTTTGCAGTCTCCCCAGAAATCGGATCGTTAAGCGACCCGATTTACCCAGGGATTTACAAAGCCGACT
ACGGACGTTTTAGGCTCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s1_104 chimera3-1 ref=813079,565570 amp=448..739,429..720 pos=1..150
GTGCCAGCCGCCGCGGTAATACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAA
ATCCCTGGGTAAATCGGGAAGCTTAACTTTCCGACTTCCGAGGAGACTGTCAAACTTGGGACCGGGAG
>s1_105 chimera3-2 ref=813079,565570 amp=448..739,429..720 pos=1..150
GTGCCAGCCGCCGCGGTAATACCCGCAGCTCAAGTGGTGGTCGCTATTATTGAGCCTAAAACGTCCGTAGTCGGCTTTGTAA
ATCCCTGGGTAAATCGGGAAGCTTAACTTTCCGACTTCCGAGGAGACTGTCAAACTTGGGACCGGGAG
>s1_106 chimera4-1 ref=518927,241952 amp=479..770,452..744 pos=1..150
GTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAA
GTCTGTCGGGAAATCCGCACGCCTAACGTGCGGGCGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAA
>s1_107 chimera4-2 ref=518927,241952 amp=479..770,452..744 pos=1..150
GTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAA
GTCTGTCGGGAAATCCGCACGCCTAACGTGCGGGCGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAA
>s1_108 chimera5-1 ref=518927,330533 amp=479..770,449..741 pos=1..150
GTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCCGAACAA
GTCTGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAA
>s1_109 chimera5-2 ref=518927,330533 amp=479..770,449..741 pos=1..150
GTGTCAGCCGCCGCGGTAATACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCCGAACAA
GTCTGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAA
>s1_110 gi|304442482|gb|HM753704.1|-10-1
GTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATAT
TGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTA
>s1_111 gi|304442482|gb|HM753704.1|-10-2
GTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATAT
TGCTTTTGATGCCGACCCTAAATTTTTTGCCTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTA
>s1_112 gi|304442482|gb|HM753704.1|-9-1
GCGTGTAGCGAACTGCGATGGGCATACTGTAACCATAAGGCCACGTATTTTGCAAGCTATTTAACTGGCGGCGATTGCGTAC
CCGACGACCAAAATTAGGGTCAACGCTACCTGTAGGAAGTGTCCGCATAAAATGCACCGCATGGAAAT
>s1_113 gi|304442482|gb|HM753704.1|-9-2
GCGTGTAGCGAACTGCGATGGGCATACTGTAACCATAAGGCCACGTATTTTGCAAGCTATTTAACTGGCGGCGATTGCGTAC
CCGACGACCAAAATTAGGGTCAACGCTACCTGTAGGAAGTGTCCGCATAAAATGCACCGCATGGAAAT
>s1_114 gi|304442482|gb|HM753704.1|-8-1
GGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCC
TGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTT
>s1_115 gi|304442482|gb|HM753704.1|-8-2
GGATACATCTGTCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGCC
TGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTT
>s1_116 gi|304442482|gb|HM753704.1|-7-1
GTGCTGGTGCTAATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGA
CAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGA
>s1_117 gi|304442482|gb|HM753704.1|-7-2
GTGCTGGTGCTAATGCTTCCTCTGCTGGTATGGTTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGA
CAATCAGAAAGAGATTGCCGAGATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGA
>s1_118 gi|304442482|gb|HM753704.1|-6-1
GTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTT
ATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCG
>s1_119 gi|304442482|gb|HM753704.1|-6-2
GTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCATCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTT
ATTACCATTTCAACTACTCCGGTTATCGCTGGCGACTCCTTCGAGATGGACGCCGTTGGCGCTCTCCG
"""

seqs_s2 = """>s2_120 4483264-10
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_121 4483264-9
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_122 4483264-8
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACCTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_123 4483264-7
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_124 4483264-6
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_125 4483264-5
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGACAGGCGGAATTCCTAG
>s2_126 4483264-4
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTCGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_127 4483264-3
TACCTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATA
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGATGCAGGCGGAATTCCTAG
>s2_128 4483264-2
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_129 4483264-1
TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGTTGTCAAGTCAGATGTGAAATCTATG
GGCTTAACCCATAACGTGCATTTGAAACTGACAGTCTTGAGTGATGGAGAGGCAGGCGGAATTCCTAG
>s2_130 4483458-10
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGATGAGGAAGGTGGAATTTCCGG
>s2_131 4483458-9
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_132 4483458-8
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_133 4483458-7
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_134 4483458-6
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_135 4483458-5
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_136 4483458-4
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_137 4483458-3
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAGTTTCCGG
>s2_138 4483458-2
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_139 4483458-1
TACAGGGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGCAGGCGGATCAATAAGTCGGGGGTTAAATCCATG
TGCTTAACACATGCACGGCTTCCGATACTGTTGATCTAGAGTCTCGAAGAGGAAGGTGGAATTTCCGG
>s2_140 565570-10
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_141 565570-9
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_142 565570-8
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCCTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_143 565570-7
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCTGGAGAGGCTAGAGGTACTTCTGGG
>s2_144 565570-6
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_145 565570-5
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_146 565570-4
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_147 565570-3
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_148 565570-2
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_149 565570-1
TACCTGCAGCCCAAGTGGTGGTCGATTTTATTGAGTCTAAAACGTTCGTAGCCGGTTTGATAAATCCTTGGGTAAATCGGGA
AGCTTAACTTTCCGATTCCGAGGAGACTGTCAAACTTGGGACCGGGAGAGGCTAGAGGTACTTCTGGG
>s2_150 565152-10
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAAAGGTCGACGGTATTCCGGGG
>s2_151 565152-9
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_152 565152-8
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_153 565152-7
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_154 565152-6
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_155 565152-5
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTCAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_156 565152-4
GACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_157 565152-3
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_158 565152-2
TACAAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGTTCGACGGTATTCCGGGG
>s2_159 565152-1
TACCAGCCCCTTAAGTGGTAGGGACGATTATTTGGCCTAAAGCGTCCGTAGCCGGCTTAGTGAGTCTTTCGTTAAATCCGGC
GACTTAATCGTCGTTTGCGAGAGATACTGCTAGGCTTGAGGACGGGAGAGGTCGACGGTATTCCGGGG
>s2_160 518927-10
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_161 518927-9
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_162 518927-8
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_163 518927-7
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_164 518927-6
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTTTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_165 518927-5
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_166 518927-4
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_167 518927-3
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_168 518927-2
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_169 518927-1
TACCAGCTCTCCGAGTGGTGTGGATGTTTATTGGGCCTAAAGCATCCGTAGCTGGCTAGGTTAGTCCCCTGTTAAATCCACC
GAATTAATCGTTGGATGCGGGGGATACTGCTTGGCTAGGGGACGAGAGAGGCAGACGGTATTTCCGGG
>s2_170 330533-10
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_171 330533-9
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_172 330533-8
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_173 330533-7
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_174 330533-6
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_175 330533-5
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGAAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_176 330533-4
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_177 330533-3
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAATCTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_178 330533-2
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_179 330533-1
TACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCCGAACAAGTCTGTCGGGAAATCCATC
CGCTCAACGGATGGGTCCGGCGGAAACTGTTCGGCTTGGGACCGGAAGACCTGAGGGGTACGTCTGGG
>s2_180 241952-10
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_181 241952-9
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_182 241952-8
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCATGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGCGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_183 241952-7
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_184 241952-6
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGACGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_185 241952-5
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGGGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGCGTACGTCTGGG
>s2_186 241952-4
TACCGGCAGCTCCAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAGCGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAAACTCGAGGGGTACGTCTGGG
>s2_187 241952-3
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_188 241952-2
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_189 241952-1
TACCGGCAGCTCAAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCCGGCTGCGCAAGTCTGTCGGGAAATCCGCA
CGCCTAACGTGCGGGTCCGGCGGAAACTGCGTGGCTTGGGACCGGAAGACTCGAGGGGTACGTCTGGG
>s2_190 198216-10
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_191 198216-9
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_192 198216-8
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_193 198216-7
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_194 198216-6
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_195 198216-5
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGTGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_196 198216-4
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGTCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_197 198216-3
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_198 198216-2
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_199 198216-1
TACCAGCTCTCCGAGTGGTGGGGACGATTATTGGGCTTAAAGCGTTCGTAGCCGGCTCAACAAGTCTCCCGTTAAATCCAGC
GATCTAATCGTTGGATGCGGAAGATACTGTTGGGCTAGGAGGCGGGAGAAGCCGACGGTATTCTCGGG
>s2_200 191626-10
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_201 191626-9
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_202 191626-8
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_203 191626-7
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_204 191626-6
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_205 191626-5
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGCACCTCGGGG
>s2_206 191626-4
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_207 191626-3
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_208 191626-2
CACCGGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGG
>s2_209 191626-1
CACCAGCAGCCCAAGTGGTGGCCGGGTTTATTGGGCCTAAAGCGCTCGTAGCCGGACTGGTAAGTCCCTTGTGAAATCGGTC
GGCTCAACCGTTCGGTGCAAGGGATACTATCGGTCTTGAGGCCGGGAAGGGTCGGAGGTACCTCGGGT
>s2_210 187525-10
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGCTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGTTGGTACTTGAGGG
>s2_211 187525-9
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_212 187525-8
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_213 187525-7
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGACGCGGGTGGTACTTGAGGG
>s2_214 187525-6
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGGGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_215 187525-5
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_216 187525-4
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_217 187525-3
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_218 187525-2
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTAGTTGAGGG
>s2_219 187525-1
TACCAGCCCCGCGAGTGGTCAGGACGATTATTAAGCCTAAGGCGTCCGTAGCCGGCCTTGTAAGTTCTCCCTTAAATGGTAT
GGCTCAACCATACCGTGGGGAGAATACTGCAAGGCTAGGGGGCGGGAGAGGCGGGTGGTACTTGAGGG
>s2_220 chimera1-1 ref=191626,241952 amp=354..646,452..744 pos=compl(1..150)
TTCCGGTCCCAAGCCACGCAGTTTCCGCCGGACGCCCGCACGTTAGGCGTGCGGATTTCCCGACAGACTTGCGCAGCCGGCT
ACGGACGCTTTAGGCCCAATAAACCCGGCCACCACTTGGGCTGCTGGTGTTACCGCGGCGGCTGGCAC
>s2_221 chimera1-2 ref=191626,241952 amp=354..646,452..744 pos=compl(1..150)
TTCCGGTCCCAAGCCACGCAGTTTCCGCCGGACGCCCGCACGTTAGGCGTGCGGATTTCCCGACAGACTTGCGCAGCCGGCT
ACGGACGCTTTAGGCCCAATAAACCCGGCCACCACTTGGGCTGCTGGTGTTACCGCGGCGGCTGGCAC
>s2_222 chimera2-1 ref=241952,191626 amp=452..744,354..646 pos=compl(1..150)
TCCCGGCCTCAAGACCGATAGTATCCCTTGCAAGCCGAACGGTTGAGCCGACCGATTTCACAAGGGACTTACCAGTCCGGCT
ACGAGCGCTTTAGGCCCAATAATAGCGGTCATCACTTGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_223 chimera2-2 ref=241952,191626 amp=452..744,354..646 pos=compl(1..150)
TCCCGGCCTCAAGACCGATAGTATCCCTTGCAAGCCGAACGGTTGAGCCGACCGATTTCACAAGGGACTTACCAGTCCGGCT
ACGAGCGCTTTAGGCCCAATAATAGCGGTCATCACTTGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_224 chimera3-1 ref=330533,518927 amp=449..741,479..770 pos=1..150
GTGCCAGCCGCCGCGGTAATACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCTAGGTTA
GTCCCCTGTTAAATCCACCGAATTAATCGTTGGATTGCGGGGGATACTGCTTGGCTAGGGGACGAGAG
>s2_225 chimera3-2 ref=330533,518927 amp=449..741,479..770 pos=1..150
GTGCCAGCCGCCGCGGTAATACCGGCAGCTCGAGTGATGACCGCTATTATTGGGCCTAAAGCGTCCGTAGCTGGCTAGGTTA
GTCCCCTGTTAAATCCACCGAATTAATCGTTGGATTGCGGGGGATACTGCTTGGCTAGGGGACGAGAG
>s2_226 chimera4-1 ref=330533,191626 amp=449..741,354..646 pos=compl(1..150)
TCCCGGCCTCAAGACCGATAGTATCCCTTGCAAGCCGAACGGTTGAGCCGACCGATTTCACAAGGGACTTACCAGTCCGGCT
ACGAGCGCTTTAGGCCCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_227 chimera4-2 ref=330533,191626 amp=449..741,354..646 pos=compl(1..150)
TCCCGGCCTCAAGACCGATAGTATCCCTTGCAAGCCGAACGGTTGAGCCGACCGATTTCACAAGGGACTTACCAGTCCGGCT
ACGAGCGCTTTAGGCCCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_228 chimera5-1 ref=330533,518927 amp=449..741,479..770 pos=compl(1..150)
CTCTCGTCCCCTAGCCAAGCAGTATCCCCCGCAATCCAACGATTAATTCGGTGGATTTAACAGGGGACTAACCTAGCCAGCT
ACGGATGCTTTAGGCCCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_229 chimera5-2 ref=330533,518927 amp=449..741,479..770 pos=compl(1..150)
CTCTCGTCCCCTAGCCAAGCAGTATCCCCCGCAATCCAACGATTAATTCGGTGGATTTAACAGGGGACTAACCTAGCCAGCT
ACGGATGCTTTAGGCCCAATAATAGCGGTCATCACTCGAGCTGCCGGTATTACCGCGGCGGCTGGCAC
>s2_230 gi|304442482|gb|HM753704.1|-35-1
TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTC
GTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGTCCCTTACTTAAGGATAAATT
>s2_231 gi|304442482|gb|HM753704.1|-35-2
TCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACTAAAGGCAAGCGTAAAGGCGCTC
GTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGTCCCTTACTTAAGGATAAATT
>s2_232 gi|304442482|gb|HM753704.1|-34-1
GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAACGTCGGCTACAGTAACTTTTCCCAGCCTCAATCTC
ATCTCTCTTTTTGCGTTCTGCTTCAATATCTGGTTGAACGGCGTCGCGTCGTAACCCAGCTTGGTAAG
>s2_233 gi|304442482|gb|HM753704.1|-34-2
GCGCATAAATTTGAGCAGATTTGTCGTCACAGGTTGCGCCGCCAAAACGTCGGCTACAGTAACTTTTCCCAGCCTCAATCTC
ATCTCTCTTTTTGCGTTCTGCTTCAATATCTGGTTGAACGGCGTCGCGTCGTAACCCAGCTTGGTAAG
>s2_234 gi|304442482|gb|HM753704.1|-33-1
CTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGC
TAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACA
>s2_235 gi|304442482|gb|HM753704.1|-33-2
CTTATGGTTACAGTATGCCCATCGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGC
TAAAGGTGAGCCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACA
>s2_236 gi|304442482|gb|HM753704.1|-32-1
GCCGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAAT
CAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCT
>s2_237 gi|304442482|gb|HM753704.1|-32-2
GCCGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGTTTCGCTGAAT
CAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTGCTATTGCT
>s2_238 gi|304442482|gb|HM753704.1|-31-1
ACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATC
GTCAACGTTATATCTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGA
>s2_239 gi|304442482|gb|HM753704.1|-31-2
ACTTTTTTTCTGATAAGCTGGTTCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATC
GTCAACGTTATATCTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGA
"""

database_16S = """>Unc00pqt count=1; cluster_weight=1; cluster=Unc00pqt; \
cluster_score=1.000000; cluster_center=True;
gagtttgatcctggctcagaacgaacgctggcggcgtggataagacatgcaagtcgaacgagaactttcgctgt\
agcaatacagtggaagttttagtggcgaacgggtgcgtaacacgtggacaatctgccttaaagtgggggatagc\
tcggcgaaagccgaattaataccgcatgtgatcagccaaaacattttggcgaaattaaagacggcgcaagctgt\
cgctttttgaggagtccgcggcctatcagctagttggtgaggtaacggctcaccaaggctttgacgggtagctg\
gtctgagaggacgaccagccacactggaactgagacacggtccagacacctacgggtggcagcagtcgagaatt\
tttctcaatgggggaaaccctgaaggagcgacgccgcgtggaggatgaaggtcttcgggttgtaaactcctgtc\
atttgggaacaagtcgtacgtagtaactatacgtgcgttgatagtaccagaagaggaagagacggctaactctg\
tgccagcagccgcggtaacacgtagggaccaagcgttgtccggatttattgggcgtaaagcgctcgtaggcggt\
tcggtaagtcgggtgtgaaacccctgggctcaacccggggacgccacccgatactgctgtgactcgagttcggt\
aggggagcggggaattcccggtgtagcggtgaaatgcgtagatatcgggaggaacaccgttggcgaaggcgact\
ttctgggctgctactgacgctgaagcgcgaaagccaggggagcgaacaggattagataccctggtagtcctggc\
cgtaaacgttggacactaggtgttggaggcatcgacccctccagtgccgcagccaacgcattaagtgtcccgcc\
tggggagtacgaccgcaaggttgaaactcaaaggaattgacgggggcccgcacaagcggtggagcatgtggttt\
aattcgacgcaacgcgaagaaccctaccggggtttgacatgcacgggacaggtgtggaaacacaccctcccttc\
ggggtccgtgcacaggtgctgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgag\
cgcaacccttgtcgctagttgccagcgagtaatgtcgggcactctagcgagactgccggtggtaagccggagga\
aggtggggatgacgtcaagtcctcatggcccttacatcccgggctacacacgtgctacaatggccggtacaaag\
ggttccgatatcgcgagatggaggcaatcccaaaaagccggtctcagttcggattgcagtctgcaactcgactg\
catgaaggtggaatcgctagtaatcgcaaatcagcaggttgcggtgaatacgttcccgggccttgtacacac
>Unc40897 count=1; cluster_weight=3; cluster=Unc40897; cluster_score=1.000\
000; cluster_center=True;
gatgaacgctagcgagacgcttaacacatgcaagtcgaatgctgatttatcagcatggcagacgggtgagtaac\
acgtgcttaatttgccccggagtgggggacaacatcggaaacggtgataatcccgcatacgttccataaaatga\
aaatatatggaagaaagatttatcgctctgggaggagggcgcgcctgactagttagttggtgaggtaatggctc\
accaagacgatgatcagtagctggtctgagaggacgaacagccacaatgggactgagacaaggcccatacccct\
acggggggcagcagtggggaatcttgcgcaatggacgaaagtctgacgcagcgatctcgcgtgaaggatgaagc\
tttgcggagtgtaaacttcttttctgggagaaaacgactgatggtatctcaggaataagggacggctaactacg\
tgccagcagccgcggtaatacgtaggtcccaagcgttatccggatttactgggcgtaaagcgtccgtagccggc\
tttgtaagtccgctttcaaaacccgaagctcaacttcggtccgggagtagatactgcaaagctagagaaagata\
gaggtatgtagaatttccggtggaggggtgaaatccgttgatatcggaaggaataccaaaagcgaaggcagcat\
gctatatcttttctgacggtcagggacgaaagtatggggagcgaacgggattagataccccggtagtccatacc\
gtaaatgatgcccgctaggtgtgcatctcgctttgcgagatatgtgccgtaagctaacgcgttaagcgggtcgc\
ctgagtagtatagtcgcaaggctgaaactcaaagggataggcgggggaacacacaagcagaggattgtctcgat\
taattggataataagccaagaatcttaccccggtttgacatccttcaaattctgcagaaacgtagaagtgctcg\
gttcgccgagaatgaagtgacaggtgctgcatggccgtcgtcagctcgtgccgcaaggtgtctggttaagtcca\
ggaacgagcgcaaccctcatcgttagttagaatgtctaacgagactgcctcggtaacgaggaggaaggagagga\
tgacgtcaggtcattatggcccttatgccgggggcgtcagagacaatacaatgggtggtacagcaggtagcaat\
gcagcgatgcggagccaatccctaaaaccatcctcagttcggattgtagtctggaactcgactacatgaagttg\
gatttgctagtaatggtgggtcagccacaccaccgtgaatatgtccctgttcctt
>Unc03ol6 count=1; cluster_weight=177; cluster=Unc03ol6; cluster_score=1.0\
00000; cluster_center=True;
gatgaacgctagcgacaggcctaacacatgcaagtcgaggggtaacattggtgcttgcaccagatgacgaccgg\
cgcacgggtgagtaacgcgtatgcaacctgccttatactgggggatagcccgaggaaactcggattaataccgc\
ataatgatttgatttcgcatgagattaaatttaaaactccggtggtataagatgggcatgcgtaacattagcta\
gttggtgaggtaacggctcaccaaggcgatgatgtttagggggtctgagaggatagtcccccacactggtactg\
agacacggaccagactcctacgggaggcagcagtgaggaatattggtcaatggacgcaagtctgaaccagccat\
gtcgcgtgcaggaagactgccctatgggttgtaaactgcttttatataggaataatggtaattacgtgtagtta\
tttgaatgtactatatgaataaggatcggctaactccgtgccagcagccgcggtaatacggaggattcgagcgt\
tatccggatttattgggtttaaagggtccgtaggcgggctgataagtcagtggtgaaatctcacagctcaactg\
tgaaactgccattgaaactgtcggtcttgaatttggttgaggtaggcggaatacgatatgtagcggtgaaatgc\
atagatatatcgtagaacaccgattgcgaaggcagcttactaaaccaacattgacgctgatggacgaaagcgtg\
ggtagcgaacaggattagataccctggtagtccacgccgtaaacgatgattactcgttgtgcacgatacaatgt\
gcgtgactgagcgaaagcattaagtaatccacctggggagtacgttggcaacaatgaaactcaaaggaattgac\
gggggcccgcacaagcggaggaacatgtggtttaattcgatgatacgcgaggaaccttacctgggcttaaatgc\
actttgacagtttgggaaaccgaatttctcttcggagcaaagtacaaggtgctgcatggttgtcgtcagctcgt\
gccgtgaggtgtcgggttaagtcccataacgagcgcaacccctatcgttagttgctaacaggtaaagctgagga\
ctctagcgagactgccaccgtaaggtgtgaggaaggtggggatgacgtcaaatcagcacggcccttacgtccag\
ggctacacacgtgttacaatggccggtacaaagggcagctacctggtgacaggatgctaatctctaaagccggt\
ctcagttcggatcggagtctgcaactcgactccgtgaagttggattcgctagtaatcgcatatcagcaacgatg\
cggtgaatacgttcccgggcct
>Unc94190 count=1; cluster_weight=3; cluster=Unc94190; cluster_score=1.000\
000; cluster_center=True;
agagtttgatcctggctcaggatgaacgctggcggcgtggctaaggcatgcaagtcgaacgaagacattggtgt\
agcaatacattgatgtacttagtggcaaacgggagcgtaacacattggtaacgtacctcgatctcacgaataac\
tcagggaaacttgagctaatacgcgatatggatggggggccgcaaggcattcccattgaaagatttatcggatc\
gagagcgacctatgttccatcagctagttggcggggtaaaggcccaccaaggctatgacggataccaggcgcga\
gagcgtgacctggctcactgggactgagacactgcccagacacctacgggtggctgcagtcgcgaatattccac\
aatgcacgaaagtgtgatggagcgacgccgcgtgtcggaagaaggccttcgggtcgtaaacgacttttagcagg\
gacgaaccatgacggtacctgcagaataagaggttgctaactctgtgccagcagcagcggtaatacagagacct\
caagcgttatccggatttattgggcgtaaagcgttcgtaggcggctttgtaagtcatttgttaaatcctcaggc\
ttaaccctgggggctgcgggtgatactgcaaggcttgagtgtgggagaggcaagcggaatgcccggtgtagcgg\
taaaatgcgttaatatcgtggtagaacacccaaggcgaaggcagcttgctggaacacaactgacgctcagtgaa\
cgaaagcgtggggagcgaaagggattagatacccctgtagtccacgctgtaaactatggctactagatttcggg\
agtttcgaccctctcggagtcgacgaaacaagctaacgcgttaagtagcccgcctgggaagtacggtcgcaaga\
ctaaaactcaaaggaatagacggggacccgcacaagcggtggatcatgcggtttaattcgatgcaacgcgagaa\
acctcacctgggcttgaaatatagctgtgggcgattggaaacaaacgccgtcctcgagatcgctatgtaggtgc\
tgcatggctgtcgtcagctcgtgccgtgaggtgtacccttaagtggggtaacgagcgcaacccctgttgtgcgt\
tatactttcgcacaagactgccctcgcaagagggaggaaggaggggacgacgtcaagtcagcatggcccttacg\
tccagggcaacacgcatgatacaatggtcggtataatgggacgccaaagcgcaagctggagcaaatccctcaaa\
accgatcccagttcggatcggggtctgcaattcgaccccgtgaagccggaatcgctagtaaaggcggatcatca\
cgccgccttgaatacgttctcgggtcttgtacacaccgcccgtc
>Unc02sd5 count=1; cluster_weight=1; cluster=Unc02sd5; cluster_score=1.000\
000; cluster_center=True;
gattcgtggcaaacgggtgagtaacacgtgcttaatctgccctggagtggggaatacctgtcgaaaggcaacta\
attccccatacgttccggtaaatgagaaacaaccggaagaaagatttatcgctccaggaggagggcgcggccta\
ttagctagttggtgaggtaacggctcaccaaggcgatgatgggtagctggtctgagaggatgatcagccacaat\
gggactgagatacggcccatacccctacggggggcagcagtggggaatattgcgcaatggacgaaagtctgacg\
cagcgatctcgcgtggaggatgaagcattaaggtgtgtaaactccttttctggaggaagacgactgacggtact\
ccaggaataagggacggctaactacgtgccagcagccgcggtaatacgtaggtctcaagcgttatccggattta\
ctgggcgtaaagcgtccgtagccggctttgtaagtctgcattcaaataccggagctcaacttcggagagggtgt\
agatactgcattgctagagaaaaacaggggttagtagaatttccggtggaggggtgaaatccgttgatatcgga\
aggaataccaaaagcgaaggcagctaactatgttttttctgacggtaagggacgaaagcttgggtagcaaacag\
gattagataccctggtagtccacgccctaaacgatgctcactcgatatgcgatccgtaggattgcgtgtccaag\
tgaaagcgttaagtgagccacctggggagtacgtcggcaacgatgaaactcaaaggaattgacgggggtccgca\
caagcggtggagcatgtggtttaattcgatgatacgcgaggaaccttacctgggctagaatgcgagtgaccggc\
cctgaaaggggcttttccttcgggacacaaagcaaggtgctgcatggctgtcgtcagctcgtgccgtgaggtgt\
tgggttaagtcccgcaacgagcgcaacccttgtccttagttgccagcacttcgggtggggactctaaggagact\
gccggcgcaagccgcgaggaaggtggggatgacgtcaagtcatcatggcctttatgcccagggctacacacgtg\
ctacaatggccagtacagagggtagcgaagccgcgaggtgaagccaatcccagaaagctggtcccagttcggat\
tggagtctggaactcgactccatgaaggtggaatcgctagtaatcgcgcatcagccatggtgcggtgaatacgt\
tcccggaccttgtacacaccgcccgtcaaaccatgggagccgggggtgcctgaaggtggtggcc
>LynMarte count=1; cluster_weight=2; cluster=LynMarte; cluster_score=1.000\
000; cluster_center=True;
cttaagaatggggacaacagagggaaactgctgctaatacccgatgtgccggaaggtgaaagatttattgcctg\
aagatgagctcgcgtccgattagctagttggtagggtaaaagcctaccaaggcgacgatcggtagctggtctga\
gaggatgaccagccacactgggactgagacacggcccagactcctacgggaggcagcagtggggaattttccgc\
aatgggcgaaagcctgacggagcaagaccgcgtgagggaggaaggctcttgggttgtaaacctcttttctctgg\
gaagaagctctgacggtaccagaggaatcagcatcggctaactccgtgccagcagccgcggtaatacggaggat\
gcaagcgttatccggaattattgggcgtaaagcgtccgtaggtggctgttcaagtctgccgttaaaaccagtgg\
cttaaccactgaatagcggtggaaactgaatagctagagtgtggtaggggtagagggaattcccagtgtagcgg\
tgaaatgcgtagagattgggaagaacaccggtggcgaaagcgctctgctggaccacaactgacactcacaggac\
gaaagctaggggagcgaaagggattagatacccctgtagtcctagccgtaaacgatggatactaggtgttgtga\
gtatcgaccctcacagtgccggagccaacgcgttaagtatcccgcctggggagtacgcacgcaagtgtgaaact\
caaaggaattgacgggggcccgcacaagcggtggagtatgtggtttaattcgatgcaacgcgaagaaccttacc\
agggcttgacatgtcgcgaatcctcttgaaagggaggagtgccttcgggagcgcgaacacaggtggtgcatggc\
tgtcgtcagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaaccctcgtttttagttgccagca\
ttcagttgggcactctaaagagactgccggtgacaaaccggaggaaggtggggatgacgtcaagtcagcatgcc\
ccttacgtcctgggctacacacgtactacaatgcgacggacaaagagcagctaaacagcgatgtcttgctaatc\
tcgtaaaccgtggctcagttcagattgtaggctgcaactcgcctacatgaaggcggaatcgctagtaatcccag\
gtcagcatactggggtgaattcgttcccgggccttgtacacaccgcccgtcacaccatgggagttggccacgcc\
cgaagtcgttatcctaacctgcaaaggaaggagatgccgaaggtggggctgatgactggggtgaagtcgtaaca\
aggtagccgtaccggaaggtgtggctggatcacctccttttcagggagaccattccccaattgatgaccgaaaa\
atagtaattaggtcacttttgaggtcatacccaggtcgaacgagattggaacaattggctttcaaactatgatt\
aggttgggttaaatgggctattagctcaggtggttagagcgcacccctgataagggtgaggtccctggttcaag\
tccaggatggccc
>Unc03k80 count=1; cluster_weight=1; cluster=Unc03k80; cluster_score=1.000\
000; cluster_center=True;
agagtttgatcctggctcaggatgaacattggcggcgtggataaggcatgcaagtcgaacgggaacataatgtt\
ttgcgtaagtaaaattgtatgttccagtggcgaacgggttagtaacgcgtagatacatcccatagagtttggca\
tagcccatcgaaaggtgggataattccaaatagtccctgagagcaatcgaagggtaaaggagcaatccgctcta\
tgattggtctgcgtcccatcagcttgttggtgaggtaatggctcaccaaggcaatgacgggtagggggtgtgag\
agcacgacccccaacagggaaactgagacacggttcccactcccacgggaggcagcagtcgagaatcttcgaca\
atgggcgaaagcctgatcgagcgacgccgcgtgcgggacgaaggtcttcggattgtaaaccgcttttctagagg\
aagaaaccaatgtttacattggttgacggtactctaggaataagaggtgactaaactcgtgccagcagtagcgg\
taatacgagtgcctcaaatgttatccggaattattgggcgtaaagcgtccgtaggcggctctgcaagtttcttg\
ttaaaaactttggcttaaccaaagaattgcgaggaaaactgcagagcttgagggtgttagaggttggtggaact\
cacggtgtaggggtgaaatccgttgataccgtggggaacaccgaaagcgaaagcagccaaccagggcattcctg\
acgctgagggacgaaagcgtgggtagcgatacggattagatacccgtgtagtccacgccctaaacgatgctgac\
tagctatttggagtgtcgacccttcaagtggcaaagctaacgcgttaagtcagccgcctgggaagtacggccgc\
aaggctaaaactcaaaggaatagacgggggttcgcacaagcggtggatcacgtggcttaattcgacgacaagcg\
aagaaccttaccagggcttgacatgccaaggtgtagtcttactgaaaggcaaggcgaccgttaattcggagctt\
ggcacaggtgctgcatggttgtcgtcagcatgagtcttgagatgctcccttaaatggggtatcatgcgcaaccc\
ttattgcttgttttacgtatcaagcgagactgcctgtgttttcacaggaggaaggaggggatgacgccaaatcc\
gcatggccctcacgccctgggctgcacacgtgatacaatggccggtacaatgggtttgctaagcggtaacgcgg\
agctaatcccaccaaaaccggtttcagttcggattgaggtctgcaacccgacctcatgaagctggaatcgctag\
taatcgcggatcagccacgccgcggtgaatatgtccccgaaccttgtactcaccgcccgtcacaccagggaagt\
cggcagtagtttaagtcccaacttagattgggcccaaactaaggtcgataaccaaggtgaagtcgtaacaaggt\
aacc
>UncLept4 count=1; cluster_weight=1; cluster=UncLept4; cluster_score=1.000\
000; cluster_center=True;
gagtttgatcatggctcaggatgaacgctggcggcgtgcttaacacatgcaagtcgaacgaaccttcgggttag\
tggcggacgggtgagtaacgcgtgagaatctgccctcaggagggggataacagttggaaacgactgctaatacc\
ccatatgccgagaggtgaaacgtattttgcctggggatgagctcgcgtctgattagcttgtaggtagggtaacg\
gcctacctaggctacgatccgtagctggtctgagaggacgatcagccacactgggactgagacacggcccagac\
tcctacgggaggcagcagtggggaattttccgcaatgggcgcaagcctgacggagcaacgccgcgtgagggatg\
aaggcctttgggctgtaaacctcttttctcaaggaagaagacctgacggtacttgaggaataagccacggctaa\
ttccgtgccagcagccgcggtaatacgggagtggcaagcgttatccggaattattgggcgtaaagcgtccgcag\
gcggccgtttaagtctgttgttaaatcgtggagctcaactccatcacggcaatggaaactgttcggcttgagta\
tggtaggggcagagggaatttccggtgtagcggtgaaatgcgtagatatcgggaagaacaccagtggcaaaagc\
gctctgctgggccattactgacgctgaggcacgaaagcgtggggagcaaacaggattagataccctggtagtcc\
acgccgtaaacgatggatgctcggtgttggccctcattgtggggtcagcgcccaagctaacgcgttaagcatcc\
cacctggggagtacgctcgcaagagtgaaactcaaaggaattgacgggggcccgcacaagcggtggagcatgtg\
gcttaattcgatgcaacgcgaagaacctcacccaggctcgaacgcggtcggacaggtcttgaaagggatcctcc\
ttcgggtcgatcgcgagttggtgcatggctgtcgtcagttggtgtcgtgagatgttgggttaagtcccgcaacg\
agcgcaacccctgtcactagttaccagcggatcatgccggggactctagtgaaaccgcctgcgcaagcagtgag\
gaaggtggggacgacgtcaagtcatcatggcccttacgcctggggctgcacacgtgctacaatgggtgatacaa\
cgggcagccacctcgcgagagggagccaatccataaaatcactccaagttcggattggagtgtgcaactcgact\
ccatgaagccggaatcggtagtaatcgcgtatcagcaacgacgcggtgaatacgttcccgggccttgtacacac\
cgcccgta
>UncFi681 count=1; cluster_weight=2; cluster=UncFi681; cluster_score=1.000\
000; cluster_center=True;
agagtttgatcctggctcaggacgaacgctggcggcgtgcttaacacatgcaagttgtacgctccccatcggcc\
aaaaccttccggttgcggctgatgatcggagagtagcggacgggtgagtaacgcgtaaagaatctgtcttaaag\
tttgggataacctgacgaaagttgggctaatcctgaataagctaacagataggcatctatcagttagaaaagat\
ggtgcaagccatcgctttatgaggactttgcgtcggattagcttgctggtagggtaatggcctaccagggcgac\
gatccgtagttggtctgagaggacgatcaaccacactgggactgagacacggcccagactcctacgggaggctg\
cagtggggaatcttccgcaatgggtgaaaacctgacggagcaacgtcgcgtgagtgaagaaggccttcgggttg\
taaagctctgtccctggggacgaactggttgtacaggaaatggtacagctttgacggtacccggggaggaagcc\
ctggctaactacgtgccagccgccgcggtaatacgtagagggcaagcgctgtccggaatcattgagcgtaaagg\
gtacgcacgcggatatgcaagtcaaatgtgaaaggtactggcttaaccagtacagtgcatttgaaactgtatat\
cttgagtacagcacaggagaacggaattcctggagtagcggtgaaatgcgtagatatcagaagaacaccagtgg\
cgaaagcggttctctgggctgttactgacgctgaggtacgaaagctggggtagcgaacgggattagataccccg\
gtagtcccagctgtaaacactggacactaggtgttgggggttcgattccttcagtgccgcagttaacgcattaa\
gtgtcccgcctggggattacggtcgcaagactgaaactcaaaggaattgacgggggcccgcacaagcggtggag\
catgtggtttaattcgatgcaacgcgaagaaccttaccgaggtttgacatcccgtgaccacctgtgagagcagg\
atttggctttttagccacacggagacaggtggtgcatggctgtcgtcagctcgtgtcgtgagatgttgggttaa\
gtcccgcaacgagcgcaacccttattcctagttgccagcacatagtggtggggactctagggatactgccgatg\
aaagtcggaggaaggtggggatgacgtcaagtcctcatgccctttatgcctcgggctacacacgtgctacaatg\
gctgatacagagggcagcgatatcgcgagataaagctaatccttcaaaatcagtcatagttcggattgcaggct\
gcaattcgcctgcatgaagctggaatcgctagtaatcgcaggtcagcatactgcggtgaatacgttcccgggcc\
ttgtacacaccgcccgtcacaccacccgagttggatgcaccagaagtcacctgagggtgccaaaggtgtgcccg\
gtgaggggggtgaagtcgtaacaaggtagccg
>UncR6785 count=1; cluster_weight=1; cluster=UncR6785; cluster_score=1.000\
000; cluster_center=True;
ttcggcggtaaaagatgagtatgcgtcctattagctagatggtaaggtaacggcttaccatggctacgataggt\
aggggtcctgagagggagatcccccacactggtactgagacacggaccagactcctacgggaggcagcagtgag\
gaatattggacaatgggcgaaagcctgatccagcaatgccgcgtgagtgatgaaggtcttcggattgtaaagct\
ctttcgacggggaagataatgacggtacccgtagaagaagccccggctaacttcgtgtcagcagccgcggtaat\
acgaagggggctagcgttgttcggatttactgggcgtaaagcgtgcgtaggcggcttgtcaagtcagaagtgaa\
agcccggggctcaactccggaactgctttttgagactggccaggcttgagtttcgggaggaggggtaagtggga\
atttcccagtgtagaaggtgaaaattcgtagatattgggaaagaaccaccgggtggcgaaaggcggctgcctgg\
ccccgaaactgacgctgagggcacgaaagcgtgggggagcaaacaggatcagaataccctggtagtccacgccc\
gtaaacgatgagtgctagacgtttggggcaatttaatggtccgccagtgtccgaaagcttaacgcgttcaagcc\
actcccggccttgggggagttaccgggccgcaaaggttaaaaacctcaaagggaaattgtactgggggccccgc\
actaagcgggtggagcatggtggtttaattcgaagccaaccgccgcagaaaccttacctagcccttgacatgga\
cgttcgcggatcttcagagatgaagttttcggttcggccggaaccgtccacacaggtgctgtcatggctgtcgt\
cagctcgtgtcgtgagatgttgggttaagtcccgcaacgagcgcaaccctcatccttagttgccatcagttcgg\
ctgggcactctagggaaactgccggtgataagccggaggaaggtggggatgacgtcaagtcctcatggccctta\
tgggctgggctacacacgtgctacaatggcgactacagtgggcagtgacatcgcgaggtgaagctaatatccaa\
aagtcgtctcagttcggattgcactctgcaactcgagtgcatgaagttggaatcgatagtaatcgcggaacagc\
atgccgcggtgaatacgttcccgggccttgtacacaccgcccgtcacaccatgggagttggttttacccgaagc\
cggtgagctgacccgaaaggggggcagccgtccacggtaaggtcagcgactggggtgaagtggtaacaaggtaa\
ccgta
>UncCr944 count=1; cluster_weight=3; cluster=UncCr944; cluster_score=1.000\
000; cluster_center=True;
ccggttgatcctgccgggaccccactgctatcgggataggactaagacatgctagtcatgcgcttcctagccaa\
tttgggagcgcggcacacagctcagtaacacgtggctaacctgcccttgggacaagaacacccccgggaaactg\
gggctaattctcgataggtgaagaactctggaatgagtcttcgcttaaaagacgctgcgctatgcttgcaggcg\
ccgcccaaggatggggccgcgaccgatcaggttgttggtgaggtaatggctcaccaagccttttaccggtgcgg\
gccgtgagagcgggagcccggagatgggcactgagacaagggcccaggccctacggggcgcagcagtcgcgaaa\
actttgcaataagcgaaagcttgacagggctatcccgagtgccatccgctgaggaaggcttttacccagtctag\
aacgctgggggaataaggagagggcaagtctggtgtcagccgccgcggtaataccagctcttcgagtggtgtgg\
atgtttattgggcctaaagcatccgtagctggctaggttagtcccctgttaaatccaccgaattaatcgttgga\
ttgcgggggatactgcttggctaggggacgagagaggcagacggtatttccggggtaggggtgaaatcctataa\
tcccgggaggaccaccagtggcgaaggctgtctgctagaacgcgcccgacggtgagggatgaaagctgggggag\
cgaaccggattagatacccgggtagtcccagctgtaaacgatgcaggctaggtgtttggacggccacgtgccgt\
tctagtgccgcagggaaactgttaagcctgccgcctggggagtacgatcgcaagattgaaacttaaaggaattg\
gcgggggagcac
>UncA1400 count=1; cluster_weight=5; cluster=UncA1400; cluster_score=1.000\
000; cluster_center=True;
ctggttgatcctgccggaggccactgctatcggtttccgactaagccatgcgagtcgggtgtcgcaagacgccg\
gcgcacggctcggtaacacgcggattatctatcctctggtgggggataacctcgggaaactgaggctaataccc\
cataagcatttgaagtttgaaaactcagatgctgaaggcgcgagcgccagagggtgagtctgcggcctatcagg\
tagtaggtggtgtaacggaccacctagcctaagacgggtacgggccttgagagagggagcccggagatggattc\
tgagacacgaatccaggccctacggggcgcagcaggcgcgaaaccttcacactgcgcgtaagcgcgatgagggg\
agaccgagtgccttctctttacgggaaggcttttcacaagcctaaaaagcttgtggaataagggctgggcaaga\
cgggtgccagccgccgcggtaatacccgcagctcaagtggtggtcgctattattgagcctaaaacgtccgtagt\
cggctttgtaaatccctgggtaaatcgggtcgcttaacgatccgatttctggggagactgcaaagcttgggacc\
gggcgaggttagaggtactctcggggtaggggtgaaatcctgtaatcccgaggggacgacctgtggcgaaagcg\
tctaacttgaacggctccgacgatgagggacgaaggctaggggagcaaaccggattagatacccgggtagtcct\
agctgtaaacactgcccgcttggtgtgacctgtcctccgggggcaggtcgtgccggagcgaaggtgttaagcgg\
gctgcttggggagtacggccgcaaggttgaaacttaaaggaattggcgggggagcaccgcaacgggaggatgcg\
tgcggtttaattggattcaacgccggaaaactcaccaggggagaccagtggatgtgagtcaagctgacgacttt\
actcgaacaagagctggagaggtggtgcatggccatcgtcagctcgtaccgtagggcattcacttaagtgtgat\
aacgagcgagacccccattcctagttgccaacctccttgaaaaagggggcgcactctaggaagaccgcttccgc\
taaggaagaggaaggtggggtcgacggtaggtcagtatgccccgaatcccctgggctacacgcgcgctacaaag\
gacaggacaatgagttccgaccccgaaaggggaaggtaatctcgaaacctgttcgtagtccggattgagggctg\
taactcgccctcatgaaggtggattccgtagtaatcgcggatcaacagtccgcggtgaatatgcccctgctcct\
tgcacacaccgcccgtc
>UncA5959 count=1; cluster_weight=55; cluster=UncA5959; cluster_score=1.00\
0000; cluster_center=True;
actgctcagtaacacgtggagaacgtgcccttaagtggaggataatctcgggaaattgaggataatactccata\
gatcatgggatctggaatgacccatggttgaaagttccggcgcttaaggatcgctctgcggcttatcaggttgt\
agtgggtgtaacgtaccccctagccagtgacgagtatgggccttgagagagggagcccagagttggattctgag\
acacgaatccaggccctacggggcgcagcagtcgcgaaaacttcacactgggggcaaccccgatgagggaattc\
ctagtgctaggacatttgttctagcttttctctagcgtagataactagaggaataagggctgggtaagacgggt\
gccagccgccgcggtaatacctgcagcccaagtggtggtcgattttattgagtctaaagcgttcgtagccggtt\
tgataaatccttgggtaaatcgggaagcttaactttccgacttccgaggagactgtcaaacttgggaccgggag\
aggctagaggtacttctggggtaggggtaaaatcctgtaatcctagaaggaccaccggtggcgaaggcgtctag\
ctagaacggatccgacggtgagggacgaagccctgggtcgcaaacgggattagataccccggtagtccagggtg\
taaacgctgccgacttggtgttggaggcccttcgggggcattcagtgccggagagaagttgttaagtcggccac\
ttggggagtacgtccgcaaggatgaaacttaaaggaattggcgggggagcaccgcaacgggaggagcgtgcggt\
ttaattggattcaacaccggaaaactcaccaagggcgactgttacatgaaagccaggctaatgaccttgcttga\
ttttcagagaggtggtgcatggccgtcgtcagttcgtaccgtaaggcgttctcttaagtgagataacgaacgag\
accctcactaataattgctacttcgatctccggatcggaggcacattattgggaccgctggcgctaagtcagag\
gaaggagaggtcaacggtaggtcagtatgccccgaatctcttgggctacacgcgcgctacaaagggcgggacaa\
tgggctccgacgccgagaggcgaaggtaatctcgaaacccgtccgtagttcggattgagggttgtaactcaccc\
tcatgaagctggattccgtagtaatcgcgaatcaacaactcgcggtgaatatgcccctgctccttgcacacacc\
gcccgtcaaaccatccgagttgggtttcagtgaggtcacctctaattagggtgttcgaactgagatttagcaag\
gaaggttaagtcgtaacaaggtagccgt
"""

if __name__ == '__main__':
    main()
