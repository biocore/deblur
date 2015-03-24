# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


def trim_seqs(input_seqs, trim_len):
    """Step 1: trim FASTA sequences to specified length

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format

    Returns
    -------
    Generator of (str, str)
        The trimmed sequences in (label, sequence) format
    """
    for label, seq in input_seqs:
        if len(seq) >= trim_len:
            yield label, seq[:trim_len]


def dereplicate_seqs(seqs_fp):
    """Step 2a: dereplicate FASTA sequences using VSEARCH"""
    pass


def remove_singletons_seqs(seqs_fp):
    """Step 2b: remove singletons from dereplicated FASTA
       file using VSEARCH"""
    pass


def remove_artifacts_seqs(seqs_fp):
    """Step 3: remove artifacts from FASTA file"""
    pass


def multiple_sequence_alignment(seqs_fp):
    """Step 4: perform multiple sequence alignment on FASTA
       file using MAFFT"""
    pass


def remove_chimeras_denovo_from_seqs(seqs_fp):
    """Step 5: remove chimeras de novo using UCHIME"""
    pass


def generate_biom_table(seqs_fp):
    """Step 6: generate BIOM table and representative
       FASTA set"""
    pass


def launch_workflow(seqs_fp, mapping_fp):
    """Launch full deblur workflow"""
    pass
