# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from deblur.deblurring import get_sequences_stats, seq_to_array


class DeblurringTests(TestCase):
    def setUp(self):
        self.seqs = [("151_9240;size=170;",  "---tagggcaagactccatggt-----"),
                     ("151_6640;size=1068;", "---cggaggcgagatgcgtggt-----"),
                     ("151_64278;size=200;", "---tactagcaagattcctggt-----"),
                     ("151_4447;size=1812;", "---aggatgcgagatgcgtggt-----"),
                     ("151_14716;size=390;", "---gagtgcgagatgcgtggtg-----"),
                     ("151_3288;size=1337;", "---ggatgcgagatgcgtggtg-----"),
                     ("151_41690;size=157;", "---agg-gcgagattcctagtgg----"),
                     ("151_5155;size=998;",  "---gaggatgcgagatgcgtgg-----"),
                     ("151_527;size=964;",   "---acggaggatgatgcgcggt-----"),
                     ("151_5777;size=305;",  "---ggagtgcaagattccaggt-----")]

    def tearDown(self):
        pass

    def test_seq_to_array(self):
        exp = np.array([4, 4, 4, 0, 2, 2, 4, 2, 1, 2, 0, 2, 0, 3, 3, 1,
                        1, 3, 0, 2, 3, 2, 2, 4, 4, 4, 4])
        # all lowercase
        obs = seq_to_array("---agg-gcgagattcctagtgg----")
        npt.assert_equal(obs, exp)

        # all uppercase
        obs = seq_to_array("---AGG-GCGAGATTCCTAGTGG----")
        npt.assert_equal(obs, exp)

        # mixed case
        obs = seq_to_array("---AgG-gCgaGATtccTAGtgG----")
        npt.assert_equal(obs, exp)

        # non-ACTG chars
        obs = seq_to_array("---FOO-gcgagattcctagfoo----")
        exp = np.array([4, 4, 4, 5, 5, 5, 4, 2, 1, 2, 0, 2, 0, 3, 3, 1,
                        1, 3, 0, 2, 5, 5, 5, 4, 4, 4, 4])
        npt.assert_equal(obs, exp)

    def test_get_sequences_stats_error_real_length(self):
        seqs = [("151_9240;size=170;", "----tagggcaagactccatg-----"),
                ("151_6640;size=1068;", "--acggaggcga-atgcgtggt----"),
                ("151_64278;size=200;", "---tactagcaagattcctgg-----")]

        with self.assertRaises(ValueError):
            get_sequences_stats(seqs)

    def test_get_sequences_stats_error_length(self):
        seqs = [("151_9240;size=170;", "----tagggcaagactccatg----"),
                ("151_6640;size=1068;", "--acggaggcgagatgcgtggt----"),
                ("151_64278;size=200;", "---tactagcaagattcctgg-----")]

        with self.assertRaises(ValueError):
            get_sequences_stats(seqs)

    def test_get_sequences_stats(self):
        (obs_seq_freqs, obs_seqs, obs_seq_labels,
            obs_seqs_np) = get_sequences_stats(self.seqs)

        exp_seq_freqs = {
            '---aggatgcgagatgcgtggt-----': 1812,
            '---ggatgcgagatgcgtggtg-----': 1337,
            '---cggaggcgagatgcgtggt-----': 1068,
            '---gaggatgcgagatgcgtgg-----': 998,
            '---acggaggatgatgcgcggt-----': 964,
            '---gagtgcgagatgcgtggtg-----': 390,
            '---ggagtgcaagattccaggt-----': 305,
            '---tactagcaagattcctggt-----': 200,
            '---tagggcaagactccatggt-----': 170,
            '---agg-gcgagattcctagtgg----': 157,
        }

        exp_seqs = ('---aggatgcgagatgcgtggt-----',
                    '---ggatgcgagatgcgtggtg-----',
                    '---cggaggcgagatgcgtggt-----',
                    '---gaggatgcgagatgcgtgg-----',
                    '---acggaggatgatgcgcggt-----',
                    '---gagtgcgagatgcgtggtg-----',
                    '---ggagtgcaagattccaggt-----',
                    '---tactagcaagattcctggt-----',
                    '---tagggcaagactccatggt-----',
                    '---agg-gcgagattcctagtgg----')

        exp_seq_names = ('151_4447;size=1812;',
                         '151_3288;size=1337;',
                         '151_6640;size=1068;',
                         '151_5155;size=998;',
                         '151_527;size=964;',
                         '151_14716;size=390;',
                         '151_5777;size=305;',
                         '151_64278;size=200;',
                         '151_9240;size=170;',
                         '151_41690;size=157;')

        exp_seqs_np = [
            np.array([4, 4, 4, 0, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0, 3,
                      2, 1, 2, 3, 2, 2, 3, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 1, 2, 2, 0, 2, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 2, 0, 2, 2, 0, 3, 2, 1, 2, 0, 2,
                      0, 3, 2, 1, 2, 3, 2, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 0, 1, 2, 2, 0, 2, 2, 0, 3, 2, 0,
                      3, 2, 1, 2, 1, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 2, 0, 2, 3, 2, 1, 2, 0, 2, 0, 3,
                      2, 1, 2, 3, 2, 2, 3, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 2, 2, 0, 2, 3, 2, 1, 0, 0, 2, 0,
                      3, 3, 1, 1, 0, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 3, 0, 1, 3, 0, 2, 1, 0, 0, 2, 0,
                      3, 3, 1, 1, 3, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 3, 0, 2, 2, 2, 1, 0, 0, 2, 0, 1,
                      3, 1, 1, 0, 3, 2, 2, 3, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 0, 2, 2, 4, 2, 1, 2, 0, 2, 0, 3,
                      3, 1, 1, 3, 0, 2, 3, 2, 2, 4, 4, 4, 4])]

        self.assertEqual(obs_seq_freqs, exp_seq_freqs)
        self.assertEqual(obs_seqs, exp_seqs)
        self.assertEqual(obs_seq_labels, exp_seq_names)
        npt.assert_equal(obs_seqs_np, exp_seqs_np)

if __name__ == '__main__':
    main()
