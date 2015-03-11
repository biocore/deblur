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
    """Sequence object to represent the aligned reads

    Attributes
    ----------
    label : str
        The sequence label
    sequence : str
        The sequence string
    length : int
        The sequence length (aligned)
    unaligned_length : int
        The unaligned sequence length
    frequency : float
        The number of times the sequence have been seen in the dataset
    np_sequence : numpy array of int8
        An int8 numpy representation of the sequence string

    Methods
    -------
    to_fasta
    """

    def __init__(self, label, sequence):
        self.label = label
        self.sequence = sequence.upper()
        self.length = len(self.sequence)
        self.unaligned_length = self.length - self.sequence.count('-')
        self.frequency = float(re.search('(?<=size=)\w+', self.label).group(0))
        self.np_sequence = np.array(
            [trans_dict[b] for b in self.sequence], dtype=np.int8)

    def __eq__(self, other):
        return (type(self) == type(other) and
                self.sequence == other.sequence and
                self.frequency == other.frequency)

    def __ne__(self, other):
        return not self.__eq__(other)

    def to_fasta(self):
        """Returns a string with the sequence in fasta format

        Returns
        -------
        str
            The FASTA representation of the sequence
        """
        prefix, suffix = re.split('(?<=size=)\w+', self.label, maxsplit=1)
        new_count = int(round(self.frequency))
        new_label = "%s%d%s" % (prefix, new_count, suffix)
        return ">%s\n%s\n" % (new_label, self.sequence)
