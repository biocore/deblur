# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import splitext

import numpy as np
from skbio import parse_fasta


def SeqToArray(seq):
    """ convert a string sequence to a numpy array"""
    seqa = np.zeros(len(seq), dtype=np.int8)
    for ind, base in enumerate(seq):
        if base == 'A':
            seqa[ind] = 0
        elif base == 'a':
            seqa[ind] = 0
        elif base == 'C':
            seqa[ind] = 1
        elif base == 'c':
            seqa[ind] = 1
        elif base == 'G':
            seqa[ind] = 2
        elif base == 'g':
            seqa[ind] = 2
        elif base == 'T':
            seqa[ind] = 3
        elif base == 't':
            seqa[ind] = 3
        elif base == '-':
            seqa[ind] = 4
        else:
            seqa[ind] = 5
    return(seqa)


def parse_tuni(filepath):
    """"""
    # Check that the input file is a .tuni file
    fname, ext = splitext(filepath)
    if ext != '.tuni':
        raise ValueError("The input file %s is not a '.tuni' file." % filepath)

    seq_freqs = {}
    seqs = []
    seq_names = []
    seqs_np = []
    sort_freq = []

    with open(filepath, 'U') as f:
        for label, seq in parse_fasta(f):
            # Get the number of reads from the header string
            # TODO: try to find a regex, rather than using magic numbers
            num_seqs = float(label[label.find(';size=')+6:-1])
            # Convert sequences to numpy array
            np_seq = SeqToArray(seq)
            # Hash the number of reads
            seq_freqs[seq] = num_seqs
            # Store the list of sequences (needs to be sorted)
            # We use a hash for the frequencies and a numpy array for the
            # hamming comparisons
            seqs.append(seq)
            seq_names.append(label)
            seqs_np.append(np_seq)
            # Store the ordered list of frequencies for sorting
            sort_freq.append(num_seqs)

    # Check if the file actually has sequences
    if len(seqs) == 0:
        raise ValueError("No sequences found in file: %s" % filepath)

    return seq_freqs, seqs, seq_names, seqs_np, sort_freq
