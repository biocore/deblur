# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from operator import attrgetter

import numpy as np
import logging

from deblur.sequence import Sequence


def get_default_error_profile():
    """Return the default error profile for deblurring
    based on illumina run data
    """
    error_dist = [1, 0.06, 0.02, 0.02, 0.01,
                  0.005, 0.005, 0.005, 0.001, 0.001,
                  0.001, 0.0005]
    return error_dist


def get_sequences(input_seqs):
    """Returns a list of Sequences

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format

    Returns
    -------
    list of Sequence

    Raises
    ------
    ValueError
        If no sequences where found in `input_seqs`
        If all the sequences do not have the same length either aligned or
        unaligned.
    """
    seqs = [Sequence(id, seq) for id, seq in input_seqs]

    if len(seqs) == 0:
        logger = logging.getLogger(__name__)
        logger.warn('No sequences found in fasta file!')
        return None

    # Check that all the sequence lengths (aligned and unaligned are the same)
    aligned_lengths = set(s.length for s in seqs)
    unaligned_lengths = set(s.unaligned_length for s in seqs)

    if len(aligned_lengths) != 1 or len(unaligned_lengths) != 1:
        raise ValueError(
            "Not all sequence have the same length. Aligned lengths: %s, "
            "sequence lengths: %s"
            % (", ".join(map(str, aligned_lengths)),
               ", ".join(map(str, unaligned_lengths))))

    seqs = sorted(seqs, key=attrgetter('frequency'), reverse=True)
    return seqs


def deblur(input_seqs, mean_error=0.005,
           error_dist=None,
           indel_prob=0.01, indel_max=3):
    """Deblur the reads

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format. The label
        should include the sequence count in the 'size=X' format.
    mean_error : float, optional
        The mean illumina error, used for original sequence estimate.
        Default: 0.005
    error_dist : list of float, optional
        A list of error probabilities. The length of the list determines the
        amount of hamming distances taken into account. Default: None, use
        the default error profile (from get_default_error_profile() )
    indel_prob : float, optional
        Indel probability (same for N indels). Default: 0.01
    indel_max : int, optional
        The maximal number of indels expected by errors. Default: 3

    Results
    -------
    list of Sequence
        The deblurred sequences

    Notes
    -----
    mean_error is used only for normalizing the peak height before deblurring.
    The array 'error_dist' represents the error distribution, where
    Xi = max frequency of error hamming. The length of this array - 1 limits
    the hamming distance taken into account, i.e. if the length if `error_dist`
    is 10, sequences up to 10 - 1 = 9 hamming distance will be taken into
    account
    """
    logger = logging.getLogger(__name__)

    if error_dist is None:
        error_dist = get_default_error_profile()
    logger.debug('Using error profile %s' % error_dist)

    # Get the sequences
    seqs = get_sequences(input_seqs)
    if seqs is None:
        logger.warn('no sequences deblurred')
        return None
    logger.info('deblurring %d sequences' % len(seqs))

    # fix the original frequencies of each read error using the
    # mean error profile
    mod_factor = pow((1 - mean_error), seqs[0].unaligned_length)
    error_dist = np.array(error_dist) / mod_factor

    max_h_dist = len(error_dist)-1

    for seq_i in seqs:
        # no need to remove neighbors if freq. is <=0
        if seq_i.frequency <= 0:
            continue

        # Correct for the fact that many reads are expected to be mutated
        num_err = error_dist * seq_i.frequency

        # if it's low level, just continue
        if num_err[1] < 0.1:
            continue

        # Compare to all other sequences and calculate hamming dist
        seq_i_len = len(seq_i.sequence.rstrip('-'))
        for seq_j in seqs:
            # Ignore current sequence
            if seq_i == seq_j:
                continue

            # Calculate the hamming distance
            h_dist = np.count_nonzero(np.not_equal(seq_i.np_sequence,
                                                   seq_j.np_sequence))

            # If far away, don't need to correct
            if h_dist > max_h_dist:
                continue

            # Close, so lets calculate exact distance

            # We stop checking in the shortest sequence after removing trailing
            # indels. We need to do this in order to avoid double counting
            # the insertions/deletions
            l = min(seq_i_len, len(seq_j.sequence.rstrip('-')))
            sub_seq_i = seq_i.np_sequence[:l]
            sub_seq_j = seq_i.np_sequence[:l]

            mask = (sub_seq_i != sub_seq_j)
            # find all indels
            mut_is_indel = np.logical_or(sub_seq_i[mask] == 4,
                                         sub_seq_j[mask] == 4)
            num_indels = mut_is_indel.sum()
            num_substitutions = h_dist - num_indels

            correction_value = num_err[num_substitutions]

            if num_indels > indel_max:
                correction_value = 0
            elif num_indels > 0:
                # remove errors due to (PCR?) indels (saw in 22 mock mixture)
                correction_value = correction_value * indel_prob

            # met all the criteria - so correct the frequency of the neighbor
            seq_j.frequency -= correction_value

    result = [s for s in seqs if round(s.frequency) > 0]
    logger.info('%d unique sequences left following deblurring' % len(result))
    return result
