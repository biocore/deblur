# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

"""
Deblur sequences
Created on Thu Sep 19 17:43:09 2013
@author: amnon
"""

from os import listdir
from os.path import isfile, join, splitext

import numpy as np
from skbio import parse_fasta
# from skbio.alignment import global_pairwise_align
# for the pairwise alignment
from cogent.align.align import global_pairwise
from cogent.align.align import make_dna_scoring_dict
from cogent import DNA

__version__ = "0.0.1-dev"


def RemoveError(seqs, seqsnp, sfreq, readerror, meanerror, ofracerr,
                indelprob, indelmax, pyroseq):
    """ Deblur the reads

    Parameters
    ----------
    seqs :
        the list of sequences
    seqsnp :
        a list of numpy arrays of the sequences (for faster comparison)
    sfreq :
        dictionary (based on the sequence) of the number of reads for each
        sequence
    readerror :
        the maximal read error expected (fraction - typically 0.01)
    meanerror :
        the mean read error used for peak spread normalization, typically 0.01
    ofracerr :
        the error distribution array, or 0 if use default
    indelprob :
        the probability for an indel (currently constant for number of indels
        until max is reached)
    indelmax :
        the maximal number of indels expected by errors (error cutoff)
    pyroseq :
        if set, use pairwise alignment for pyrosequencing data

    Results
    -------
    sfreq - the deblurred number of reads for each sequence
            (0 if not present)
    debugdata - a list of strings

    Notes
    -----
    meanerror is used only for normalizing the peak height before deblurring,
    whereas readerror is used for calculating the expected number of errors for
    each position error distribution array X should be of length >10, where
    Xi = max frequency of error hamming if it is 0, we use the default
    distribution
    """
    # take the list values so it won't change
    fracerr = list(ofracerr)

    # we assume all sequences are of equal length
    commonlen = len(seqs[0])
    for cseq in seqs:
        if not(commonlen == len(cseq)):
            print("Not all sequences are same length!!!!")
            print(commonlen)
            print(len(cseq))
            print(cseq)
    print ("processing", len(seqs), "sequences")

    numreal = 0
    for cchar in seqs[0]:
        if not (cchar == '-'):
            numreal += 1
    modfactor = pow((1-meanerror), numreal)

    # create the error profile from the read error
    # exponential independent
    #   fracerr=[]
    #   for a in range(10):
    #       fracerr.append(pow(readerror,a)/modfactor)

    # empirical
    #    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,2*pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor]

    # used for the 22 mock mixture
    #    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,0.01,0.01,0.01,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.0005,0.0001,0.0001]

    # used for the 44 mock mixture
    #   e1=pow(readerror,1)/modfactor
    #   fracerr=[1.0/modfactor,e1,e1/4,e1/5,e1/6,e1/8,e1/10,e1/15,e1/20,e1/30,e1/40,e1/50,e1/50,e1/50,e1/50,e1/50,e1/50,e1/100,e1/500,e1/500]

    # if fracerr not supplied, use the default (22 mock mixture setup)
    print("original fracer parameter:", fracerr)
    if fracerr == 0:
        fracerr = [1.0/modfactor,
                   pow(readerror, 1)/modfactor, 0.01, 0.01, 0.01, 0.005, 0.005,
                   0.005, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.001]
        print("modified fracerr because it was 0")
    else:
        for idx, val in enumerate(fracerr):
            fracerr[idx] = fracerr[idx]/modfactor

    maxhdist = len(fracerr)-1

    print "fracerr", fracerr
    print "readerror", readerror
    print "modfactor", modfactor

    print("indel prob:", indelprob)
    print("indel max:", indelmax)
    print("readerror:", readerror)
    print("meanerror:", meanerror)
    print("mod factor:", modfactor)
    print("fracerr:", fracerr)

    # for pairwise alignment:
    # Match 10, Transition -8, Transversion -8
    DNAm = make_dna_scoring_dict(10, -8, -8)

    for idx, cseq in enumerate(seqs):
        csfreq = sfreq[cseq]
        # no need to remove neighbors if freq. is <=0
        if csfreq <= 0:
            continue
        # correct for the fact that many reads are expected to be mutated
        numerr = []
        for a in range(len(fracerr)):
            numerr.append(fracerr[a]*csfreq)

        # if it's low level, just continue
        if numerr[1] < 0.1:
            continue

        # compare to all other sequences and calculate hamming dist
        cseqnp = seqsnp[idx]
        oseqlen = len(seqs[idx].rstrip('-'))
        for idxtmp, seqnptmp in enumerate(seqsnp):
            # don't compare to ourselves (dist=0)
            if idxtmp == idx:
                continue
            # calculate the hamming distance
            hdist = np.count_nonzero(np.not_equal(seqnptmp, cseqnp))
            # if far away, don't need to correct
            if hdist > maxhdist:
                continue
            # close, so lets calculate exact distance

            numsub = 0
            numindel = 0
            # experimental try 2
            # s1=seqs[idx].replace('-','')
            # s2=seqs[idxtmp].replace('-','')
            # cseq1,cseq2=nw_align(s1,s2)

            # experimental: pairwise align the sequences
            if pyroseq:
                s0 = DNA.makeSequence(seqs[idx])
                s0 = s0.degap()
                s1 = DNA.makeSequence(seqs[idxtmp])
                s1 = s1.degap()
                print s0._seq
                print s1._seq
                align = global_pairwise(s0, s1, DNAm, 10, 9)
                a0 = align.getGappedSeq('seq_0')
                a1 = align.getGappedSeq('seq_1')
                cseq1 = a0._seq
                cseq2 = a1._seq

                len1 = len(cseq1.rstrip('-'))
                len2 = len(cseq2.rstrip('-'))
                oseqlen = len(cseq1)
                for cpos in range(oseqlen):
                    if not (cseq1[cpos] == cseq2[cpos]):
                        if cseq1[cpos] == '-':
                            if cpos < len1:
                                numindel += 1
                        else:
                            if cseq2[cpos] == '-':
                                if cpos < len2:
                                    numindel += 1
                            else:
                                numsub += 1

            # not pyrosequencing so use the faster global alignment
            else:
                for cpos in range(oseqlen):
                    if not (cseqnp[cpos] == seqnptmp[cpos]):
                        # 4 is '-'
                        if seqnptmp[cpos] == 4:
                            numindel += 1
                        else:
                            if cseqnp[cpos] == 4:
                                numindel += 1
                            else:
                                numsub += 1

            nerr = numerr[numsub]

            # remove errors due to (PCR?) indels (saw in 22 mock mixture)
            if numindel > 0:
                nerr = nerr*indelprob
            if numindel > indelmax:
                nerr = 0

            # if the effect is small - don't do anything
            if nerr < 0.1:
                continue
            # met all the criteria - so correct the frequency of the neighbor
            sfreq[seqs[idxtmp]] -= nerr

    return(sfreq)


def deblur_file(filepath, read_error, mean_error, error_dist, indel_prob,
                indel_max, pyroseq):
    """
    Parameters
    ----------
    filepath : str
        Path to the input filename
    read_error :
    mean_error :
    error_dist :
    indel_prob :
    indel_max :
    pyroseq : bool

    Returns
    -------
    """
    # Parse the input file
    seq_freqs, seqs, seq_names, seqs_np, sort_freq = parse_tuni(filepath)

    # Sort the sequences according to frequencies (descending)
    sort_order = sorted(range(len(sort_freq)), key=sort_freq.__getitem__,
                        reverse=True)
    seqs = [seqs[i] for i in sort_order]
    seq_names = [seq_names[i] for i in sort_order]
    seqs_np = [seqs_np[i] for i in sort_order]

    # After loading the file - remove the read errors (MAIN FUNCTION)
    c_freq = RemoveError(seqs, seqs_np, seq_freqs, read_error, mean_error,
                         error_dist, indel_prob, indel_max, pyroseq)

    # And finally save the new fasta as a '.clean' file
    with open(join(filepath, '.clean'), 'w') as out_f:
        c_id = 1
        for c_seq in seqs:
            freq_rounded = round(c_freq[c_seq])
            if freq_rounded > 0:
                out_f.write('>aa%d;size=%d;\n' % (c_id, int(freq_rounded)))
                for c in c_seq.upper():
                    # We remove the indels
                    if not c == '-':
                        out_f.write(c)
                out_f.write('\n')
