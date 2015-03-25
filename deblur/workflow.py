# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from os.path import splitext
from bfillings.vsearch import vsearch_dereplicate_exact_seqs


def trim_seqs(seqs_fp, trim_len):
    """Step 1: trim FASTA sequences to specified length"""
    pass


def dereplicate_seqs(seqs_fp,
                     output_fp,
                     minuniquesize=2):
    """Step 2a: dereplicate FASTA sequences and remove
       singletons using VSEARCH

    Parameters
    ----------
    seqs_fp : string
        filepath to FASTA sequence file
    output_fp : string
        filepath to dereplicated FASTA file

    Returns
    -------
    None
    """
    log_name = "%s.log" % splitext(output_fp)[0]

    vsearch_dereplicate_exact_seqs(
        fasta_filepath=seqs_fp,
        output_filepath=output_fp,
        minuniquesize=minuniquesize,
        log_name=log_name)


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
