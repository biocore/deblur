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

from deblur.sequence import Sequence


class SequenceTests(TestCase):
    def setUp(self):
        self.label = "151_4447;size=1812;"
        self.sequence = "---aggatgcgagatgcgtggt-----"
        self.exp_seq = "---AGGATGCGAGATGCGTGGT-----"
        self.exp_np_seq = np.array([4, 4, 4, 0, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0,
                                    3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4, 4])

    def test_init(self):
        obs = Sequence(self.label, self.sequence)

        self.assertEqual(obs.label, self.label)
        self.assertEqual(obs.sequence, self.exp_seq)
        self.assertEqual(obs.length, 27)
        self.assertEqual(obs.unaligned_length, 19)
        self.assertEqual(obs.frequency, 1812)
        npt.assert_equal(obs.np_sequence, self.exp_np_seq)

    def test_init_uppercase(self):
        sequence = "---AGGATGCGAGATGCGTGGT-----"

        obs = Sequence(self.label, sequence)

        self.assertEqual(obs.label, self.label)
        self.assertEqual(obs.sequence, self.exp_seq)
        self.assertEqual(obs.length, 27)
        self.assertEqual(obs.unaligned_length, 19)
        self.assertEqual(obs.frequency, 1812)
        npt.assert_equal(obs.np_sequence, self.exp_np_seq)

    def test_init_mixed_case(self):
        sequence = "---AggATgcGAgatGCgtgGT-----"

        obs = Sequence(self.label, sequence)

        self.assertEqual(obs.label, self.label)
        self.assertEqual(obs.sequence, self.exp_seq)
        self.assertEqual(obs.length, 27)
        self.assertEqual(obs.unaligned_length, 19)
        self.assertEqual(obs.frequency, 1812)
        npt.assert_equal(obs.np_sequence, self.exp_np_seq)

    def test_init_non_actg_chars(self):
        sequence = "---FoOatgcgagatgcgtfOo-----"
        exp_seq = "---FOOATGCGAGATGCGTFOO-----"
        exp_np_seq = np.array([4, 4, 4, 5, 5, 5, 0, 3, 2, 1, 2, 0, 2, 0,
                               3, 2, 1, 2, 3, 5, 5, 5, 4, 4, 4, 4, 4])

        obs = Sequence(self.label, sequence)

        self.assertEqual(obs.label, self.label)
        self.assertEqual(obs.sequence, exp_seq)
        self.assertEqual(obs.length, 27)
        self.assertEqual(obs.unaligned_length, 19)
        self.assertEqual(obs.frequency, 1812)
        npt.assert_equal(obs.np_sequence, exp_np_seq)

    def test_eq_true(self):
        s1 = Sequence(self.label, self.sequence)
        s2 = Sequence(self.label, self.sequence)
        self.assertTrue(s1 == s2)

        s2 = Sequence("152_4447;size=1812;", self.sequence)
        self.assertTrue(s1 == s2)

    def test_eq_false(self):
        s1 = Sequence(self.label, self.sequence)
        s2 = Sequence(self.label, "---aggatgcgagatgcgtgtt-----")
        self.assertFalse(s1 == s2)

        self.assertFalse(s1 == 2)
        self.assertFalse(s1 == "---aggatgcgagatgcgtgtt-----")

    def test_ne_true(self):
        s1 = Sequence(self.label, self.sequence)
        s2 = Sequence(self.label, "---aggatgcgagatgcgtgtt-----")
        self.assertTrue(s1 != s2)

        self.assertTrue(s1 != 2)
        self.assertTrue(s1 != "---aggatgcgagatgcgtgtt-----")

    def test_ne_false(self):
        s1 = Sequence(self.label, self.sequence)
        s2 = Sequence(self.label, self.sequence)
        self.assertFalse(s1 != s2)

        s2 = Sequence("152_4447;size=1812;", self.sequence)
        self.assertFalse(s1 != s2)

    def test_to_fasta(self):
        seq = Sequence(self.label, self.sequence)
        obs = seq.to_fasta()
        exp = ">151_4447;size=1812;\n---AGGATGCGAGATGCGTGGT-----\n"
        self.assertEqual(obs, exp)

        seq.frequency = 1811.1
        obs = seq.to_fasta()
        exp = ">151_4447;size=1811;\n---AGGATGCGAGATGCGTGGT-----\n"
        self.assertEqual(obs, exp)

        seq.frequency = 1811.5
        obs = seq.to_fasta()
        exp = ">151_4447;size=1812;\n---AGGATGCGAGATGCGTGGT-----\n"
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
