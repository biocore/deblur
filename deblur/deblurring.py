# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from collections import defaultdict
from operator import itemgetter
import re

import numpy as np


def seq_to_array(seq):
    """Convert a string sequence to a numpy array

    Parameters
    ----------
    seq : str
        The sequence string

    Returns
    -------
    np.array
        The numpy representation of `seq`
    """
    trans_dict = defaultdict(lambda: 5, A=0, C=1, G=2, T=3)
    trans_dict['-'] = 4
    return np.array([trans_dict[base] for base in seq.upper()], dtype=np.int8)


def get_sequences_stats(input_seqs):
    """Gets the sequence information

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format

    Returns
    -------
    (seq_freqs, seqs, seq_labels, np_seqs)
        seq_freqs: dictionary of {sequence: num_seqs}
        seqs: list of sequences, sorted by num_seqs
        seq_labels: list of sequence labels, sorted by num_seqs
        np_seqs: list of numpy-sequences, sorted by num_seqs

    Raises
    ------
    ValueError
        If no sequences where found in `input_seqs`
        If all the sequences do not have the same length either aligned or
        unaligned.
    """
    seq_freqs = {}
    seq_info = []
    seq_lengths = []
    real_lengths = []
    for label, seq in input_seqs:
        # Store the sequence length
        seq_lengths.append(len(seq))
        real_lengths.append(len(seq) - seq.count('-'))
        # Get the number of reads from the label
        num_seqs = float(re.search('(?<=size=)\w+', label).group(0))
        # Convert the sequences to numpy array
        np_seq = seq_to_array(seq)
        # Hash the number of reads
        seq_freqs[seq] = num_seqs
        # Store the list of sequences
        seq_info.append((num_seqs, seq, label, np_seq))

    # Check if we actually had sequences
    if len(seq_info) == 0:
        raise ValueError("No sequences found!")

    # The code assumes that all sequences are of equal length, raising a
    # ValueError in case that is not true
    seq_lengths = set(seq_lengths)
    real_lengths = set(real_lengths)
    if len(seq_lengths) != 1 or len(real_lengths) != 1:
        raise ValueError(
            "Not all sequence have the same length. Aligned lengths: %s, "
            "sequence lengths: %s"
            % (", ".join(map(str, seq_lengths)),
               ", ".join(map(str, real_lengths))))

    # Get 3 lists: the sequence, the sequence label and the numpy sequence
    # all lists sorted by the num_seqs value
    _, seqs, seq_labels, np_seqs = zip(
        *sorted(seq_info, key=itemgetter(0), reverse=True))

    return seq_freqs, seqs, seq_labels, np_seqs


def deblur(input_seqs, read_error=0.05, mean_error=None, error_dist=None,
           indel_prob=0.01, indel_max=3):
    """ Deblur the reads

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format
    read_error : float, optional
        the maximal read error expected (fraction - typically 0.01)
    mean_error :
        the mean read error used for peak spread normalization, typically 0.01
    error_dist :
        the error distribution array, or 0 if use default
    indel_prob :
        the probability for an indel (currently constant for number of indels
        until max is reached)
    indel_max :
        the maximal number of indels expected by errors (error cutoff)

    Results
    -------
    seq_freqs - the deblurred number of reads for each sequence
            (0 if not present)

    Notes
    -----
    mean_error is used only for normalizing the peak height before deblurring,
    whereas read_error is used for calculating the expected number of errors
    for each position.
    The error distribution array 'error_dist' should be of length >10, where
    Xi = max frequency of error hamming. If it is 0, we use the default
    distribution
    """
    # Step 1: get the sequence information
    seq_freqs, seqs, seq_labels, np_seqs = get_sequences_stats(input_seqs)

    num_real = len(seqs[0]) - seqs[0].count('-')
    mod_factor = pow((1 - mean_error), num_real)

    # if error_list not supplied, use the default (22 mock mixture setup)
    if not error_dist:
        error_dist = [1.0 / mod_factor, pow(read_error, 1) / mod_factor, 0.01,
                      0.01, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005,
                      0.001, 0.001, 0.001, 0.001]
    else:
        error_dist[:] = [val / mod_factor for val in error_dist]

    max_h_dist = len(error_dist)-1

    for idx, curr_seq in enumerate(seqs):
        curr_seq_freq = seq_freqs[curr_seq]

        # no need to remove neighbors if freq. is <=0
        if curr_seq_freq <= 0:
            continue

        # correct for the fact that many reads are expected to be mutated
        num_err = [val * curr_seq_freq for val in error_dist]

        # if it's low level, just continue
        if num_err[1] < 0.1:
            continue

        # compare to all other sequences and calculate hamming dist
        curr_np_seq = np_seqs[idx]

        curr_seq_len = len(seqs[idx].rstrip('-'))

        for idx_tmp, np_seq_tmp in enumerate(np_seqs):
            # Ignore current sequence
            if idx_tmp == idx:
                continue

            # calculate the hamming distance
            h_dist = np.count_nonzero(np.not_equal(np_seq_tmp, curr_np_seq))

            # if far away, don't need to correct
            if h_dist > max_h_dist:
                continue

            # close, so lets calculate exact distance
            # num_sub -> num substitutions
            num_sub = 0
            num_indel = 0

            # We stop checking in the shortest "real" sequence. We need to do
            # this in order to avoid double counting the insertion/deletions
            l = len(seqs[idx_tmp].rstrip('-'))
            stop_length = min(curr_seq_len, l)

            for curr_pos in range(stop_length):
                if not curr_np_seq[curr_pos] == np_seq_tmp[curr_pos]:
                    # 4 is '-'
                    if np_seq_tmp[curr_pos] == 4:
                        num_indel += 1
                    else:
                        if curr_np_seq[curr_pos] == 4:
                            num_indel += 1
                        else:
                            num_sub += 1

            correction_value = num_err[num_sub]

            if num_indel > indel_max:
                correction_value = 0
            elif num_indel > 0:
                # remove errors due to (PCR?) indels (saw in 22 mock mixture)
                correction_value = correction_value * indel_prob

            # met all the criteria - so correct the frequency of the neighbor
            seq_freqs[seqs[idx_tmp]] -= correction_value

    return seq_freqs
