# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from collections import defaultdict
import re

import numpy as np

trans_dict = defaultdict(lambda: 5, A=0, C=1, G=2, T=3)
trans_dict['-'] = 4


class Sequence(object):

    def __init__(self, label, sequence):
        self.label = label
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.unaligned_length = self.length - self.sequence.count('-')
        self.frequency = float(re.search('(?<=size=)\w+', self.label).group(0))
        self.np_sequence = np.array(
            [trans_dict[b] for b in self.sequence], dtype=np.int8)
