# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------


from os.path import splitext, join, dirname, exists
from bfillings.vsearch import (vsearch_dereplicate_exact_seqs,
                               vsearch_chimera_filter_de_novo)
from os import rename, makedirs


def trim_seqs(seqs_fp, trim_len):
    """Step 1: trim FASTA sequences to specified length"""
    pass


def dereplicate_seqs(seqs_fp,
                     output_fp,
                     min_size=2):
    """Step 2a: dereplicate FASTA sequences and remove
       singletons using VSEARCH

    Parameters
    ----------
    seqs_fp : string
        filepath to FASTA sequence file
    output_fp : string
        filepath to dereplicated FASTA file
    min_size : integer
        discard sequences with an abundance value smaller
        than integer
    """
    log_name = "%s.log" % splitext(output_fp)[0]

    vsearch_dereplicate_exact_seqs(
        fasta_filepath=seqs_fp,
        output_filepath=output_fp,
        minuniquesize=min_size,
        log_name=log_name)


def remove_artifacts_seqs(seqs_fp):
    """Step 3: remove artifacts from FASTA file"""
    pass


def multiple_sequence_alignment(seqs_fp):
    """Step 4: perform multiple sequence alignment on FASTA
       file using MAFFT"""
    pass


def remove_chimeras_denovo_from_seqs(seqs_fp, output_fp):
    """Step 5: remove chimeras de novo using UCHIME
       (VSEARCH implementation)

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    output_fp: string
        file path to store chimera-free results
    """
    working_dir = join(dirname(output_fp), "working_dir")
    if not exists(working_dir):
        makedirs(working_dir)

    output_chimera_filepath, output_non_chimera_filepath,\
        output_alns_filepath, output_tabular_filepath, log_filepath =\
        vsearch_chimera_filter_de_novo(
            fasta_filepath=seqs_fp,
            working_dir=working_dir,
            output_chimeras=False,
            output_nonchimeras=True,
            output_alns=False,
            output_tabular=False,
            log_name="vsearch_uchime_de_novo_chimera_filtering.log")

    rename(output_non_chimera_filepath, output_fp)


def generate_biom_table(seqs_fp):
    """Step 6: generate BIOM table and representative
       FASTA set"""
    pass


def launch_workflow(seqs_fp, mapping_fp):
    """Launch full deblur workflow"""
    pass
