# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkstemp
from os import close, remove

import numpy as np
import numpy.testing as npt

from deblur.parse import parse_tuni


class ParseTests(TestCase):
    """Tests the parser functions"""
    def setUp(self):
        fd, self.tuni_fp = mkstemp(suffix='.fasta.tuni')
        close(fd)

        with open(self.tuni_fp, 'w') as f:
            f.write(TUNI_CONTENTS)

        self._clean_up_files = []

    def tearDown(self):
        for f in self._clean_up_files:
            remove(f)

    def test_parse_tuni(self):
        """.tuni files are correctly parsed"""
        (obs_seq_freqs, obs_seqs, obs_seq_names, obs_seqs_np,
         obs_sort_freq) = parse_tuni(self.tuni_fp)

        exp_seq_freqs = {
            '---aggatgcgagatgcgtggt----': 1812,
            '----ggatgcgagatgcgtggt----': 1337,
            '--acggaggcgagatgcgtggt----': 1068,
            '--gaggatgcgagatgcgtggt----': 998,
            '---acggaggatgatgcgcggt----': 964,
            '----gagtgcgagatgcgtggt----': 390,
            '---ggagtgcaagattccagg-----': 305,
            '---tactagcaagattcctgg-----': 200,
            '----tagggcaagactccatg-----': 170,
            '----agg-gcgagattcctag-----': 157,
        }

        exp_seqs = ['---aggatgcgagatgcgtggt----',
                    '----ggatgcgagatgcgtggt----',
                    '--acggaggcgagatgcgtggt----',
                    '--gaggatgcgagatgcgtggt----',
                    '---acggaggatgatgcgcggt----',
                    '----gagtgcgagatgcgtggt----',
                    '---ggagtgcaagattccagg-----',
                    '---tactagcaagattcctgg-----',
                    '----tagggcaagactccatg-----',
                    '----agg-gcgagattcctag-----']

        exp_seq_names = ['151_4447;size=1812;',
                         '151_3288;size=1337;',
                         '151_6640;size=1068;',
                         '151_5155;size=998;',
                         '151_527;size=964;',
                         '151_14716;size=390;',
                         '151_5777;size=305;',
                         '151_64278;size=200;',
                         '151_9240;size=170;',
                         '151_41690;size=157;']

        exp_seqs_np = [
            np.array([4, 4, 4, 0, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 4, 4, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 0, 1, 2, 2, 0, 2, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 2, 0, 2, 2, 0, 3, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 4, 0, 1, 2, 2, 0, 2, 2, 0, 3, 2, 0,
                      3, 2, 1, 2, 1, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 4, 4, 2, 0, 2, 3, 2, 1, 2, 0, 2, 0,
                      3, 2, 1, 2, 3, 2, 2, 3, 4, 4, 4, 4]),
            np.array([4, 4, 4, 2, 2, 0, 2, 3, 2, 1, 0, 0, 2, 0,
                      3, 3, 1, 1, 0, 2, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 3, 0, 1, 3, 0, 2, 1, 0, 0, 2, 0,
                      3, 3, 1, 1, 3, 2, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 4, 3, 0, 2, 2, 2, 1, 0, 0, 2, 0,
                      1, 3, 1, 1, 0, 3, 2, 4, 4, 4, 4, 4]),
            np.array([4, 4, 4, 4, 0, 2, 2, 4, 2, 1, 2, 0, 2, 0,
                      3, 3, 1, 1, 3, 0, 2, 4, 4, 4, 4, 4])]

        exp_sort_freq = [1812, 1337, 1068, 998, 964, 390, 305, 200, 170, 157]

        self.assertEqual(obs_seq_freqs, exp_seq_freqs)
        self.assertEqual(obs_seqs, exp_seqs)
        self.assertEqual(obs_seq_names, exp_seq_names)
        npt.assert_equal(obs_seqs_np, exp_seqs_np)
        self.assertEqual(obs_sort_freq, exp_sort_freq)

TUNI_CONTENTS = """>151_4447;size=1812;
---aggatgcgag
atgcgtggt----
>151_3288;size=1337;
----ggatgcgag
atgcgtggt----
>151_6640;size=1068;
--acggaggcgag
atgcgtggt----
>151_5155;size=998;
--gaggatgcgag
atgcgtggt----
>151_527;size=964;
---acggaggatg
atgcgcggt----
>151_14716;size=390;
----gagtgcgag
atgcgtggt----
>151_5777;size=305;
---ggagtgcaag
attccagg-----
>151_64278;size=200;
---tactagcaag
attcctgg-----
>151_9240;size=170;
----tagggcaag
actccatg-----
>151_41690;size=157;
----agg-gcgag
attcctag-----
"""

if __name__ == '__main__':
    main()
