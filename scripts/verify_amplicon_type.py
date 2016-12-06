#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click
import numpy as np
import pandas as pd

def count_starting_kmers(input_fasta_fp, num_seqs, seed):
    """Generate value_counts dataframe of starting k-mers for random subsample
    of a fasta file"""
    kmer_length = 4
    if seed:
        np.random.seed(seed)
    starting_kmers = []
    with open(input_fasta_fp) as handle:
        lines = pd.Series(handle.readlines())
        num_lines = len(lines)
        if num_lines/2 < num_seqs:
            rand_line_nos = np.random.choice(np.arange(1,num_lines,2), 
                                             size=num_seqs, replace=True)
        else:
            rand_line_nos = np.random.choice(np.arange(1,num_lines,2), 
                                             size=num_seqs, replace=False)
        rand_lines = lines[rand_line_nos]
    for sequence in rand_lines:
        starting_kmers.append(sequence[:kmer_length])
    starting_kmer_value_counts = pd.Series(starting_kmers).value_counts()
    return(starting_kmer_value_counts)

@click.command()
@click.option('--input_fasta_fp', '-f', required=True,
              type=click.Path(resolve_path=True, readable=True, exists=True,
              file_okay=True), 
              help="Input fasta file from Deblur (.fa, .fna, .fasta)")
@click.option('--num_seqs', '-n', required=False, type=int, default=10000,
              help="Number of sequences to randomly subsample [default: 10000]")
@click.option('--cutoff', '-c', required=False, type=float, default=0.5,
              help="Minimum fraction of sequences required to match "
                   "a diagnostic k-mer [default: 0.5]")
@click.option('--seed', '-s', required=False, type=int, 
              help="Random number seed [default: None]")

def verify_amplicon_type(input_fasta_fp, num_seqs, cutoff, seed):
    """Determine the most likely amplicon type of a fasta file based on the
    first few nucleotides.

    The most frequent starting k-mer in a random subsample of sequences must 
    match, above a given cutoff fraction of sequences, one of the following  
    diagnostic k-mers:

      k-mer\tAmplicon\tForward primer                             
      TACG\t16S rRNA\t515f                                       
      NNNN\t18S rRNA\tEuk1391f                                    
      NNNN\tITS rRNA\tITS1f                                       
    """
    starting_kmer_value_counts = count_starting_kmers(input_fasta_fp, num_seqs, 
      seed)
    top_kmer = starting_kmer_value_counts.index[0]
    top_kmer_count = starting_kmer_value_counts[0]
    top_kmer_frac = top_kmer_count/num_seqs
    if (top_kmer == 'TACG') & (top_kmer_frac > cutoff):
        print('Amplicon type: 16S/515f (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    elif (top_kmer == 'NNNN') & (top_kmer_frac > cutoff):
        print('Amplicon type: 18S/Euk1391f (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    elif (top_kmer == 'NNNN') & (top_kmer_frac > cutoff):
        print('Amplicon type: ITS/ITS1f (%s%% of sequences start with %s)' %
              (round(top_kmer_frac*100, 1), top_kmer))
    else:
        print('Could not determine amplicon type'),
        print('(most frequent starting k-mer was %s with %s%%)' %
              (top_kmer, round(top_kmer_frac*100, 1)))

if __name__ == '__main__':
    verify_amplicon_type()
