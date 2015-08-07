# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os.path import splitext, dirname, join, exists, basename
from os import makedirs, stat

from bfillings.vsearch import vsearch_dereplicate_exact_seqs
from bfillings.sortmerna_v2 import (build_database_sortmerna,
                                    sortmerna_map)
from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta


def trim_seqs(input_seqs, trim_len):
    """Trim FASTA sequences to specified length

    Parameters
    ----------
    input_seqs : iterable of (str, str)
        The list of input sequences in (label, sequence) format
    trim_len : int
        Sequence trimming length

    Returns
    -------
    Generator of (str, str)
        The trimmed sequences in (label, sequence) format
    """
    for label, seq in input_seqs:
        if len(seq) >= trim_len:
            yield label, seq[:trim_len]


def dereplicate_seqs(seqs_fp,
                     output_fp,
                     min_size=2):
    """Dereplicate FASTA sequences and remove singletons using VSEARCH

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


def remove_artifacts_seqs(seqs_fp,
                          ref_fp,
                          output_fp,
                          ref_db_fp=None,
                          negate=False,
                          threads=1):
    """Remove artifacts from FASTA file using SortMeRNA

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    ref_fp: tuple
        file path(s) to FASTA database file
    output_fp: string
        file path to store output results
    ref_db_fp: string or tuple, optional
        file path(s) to indexed FASTA database
    negate: boolean, optional
        if True, discard all input sequences aligning
        to reference database
    threads: integer, optional
        number of threads to use for SortMeRNA
    """
    working_dir = join(dirname(output_fp), "working_dir")
    if not exists(working_dir):
        makedirs(working_dir)

    aligned_seq_ids = set()
    files_to_remove = []

    for i, db in enumerate(ref_fp):
        # create working directory for each
        # reference database
        db_dir_base = splitext(basename(db))[0]
        db_dir = join(working_dir, db_dir_base)
        if not exists(db_dir):
            makedirs(db_dir)

        if ref_db_fp:
            sortmerna_db = ref_db_fp[i]
        else:
            # build index
            sortmerna_db, files_to_remove = \
                build_database_sortmerna(
                    fasta_path=db,
                    max_pos=10000,
                    output_dir=db_dir)

        # run SortMeRNA
        app_result = sortmerna_map(
            seq_path=seqs_fp,
            output_dir=db_dir,
            refseqs_fp=db,
            sortmerna_db=sortmerna_db,
            threads=threads,
            best=1)

        # Print SortMeRNA errors
        stderr_fp = app_result['StdErr'].name
        if stat(stderr_fp).st_size != 0:
            with open(stderr_fp, 'U') as stderr_f:
                for line in stderr_f:
                    print line
            raise ValueError("Could not run SortMeRNA.")

        for line in app_result['BlastAlignments']:
            line = line.strip().split('\t')
            if line[1] == '*':
                continue
            else:
                aligned_seq_ids.add(line[0])

        # remove indexed database files
        remove_files(files_to_remove, error_on_missing=False)

    if negate:
        def op(x): return x not in aligned_seq_ids
    else:
        def op(x): return x in aligned_seq_ids

    # if negate = False, only output sequences
    # matching to at least one of the databases
    with open(seqs_fp, 'U') as seqs_f:
        with open(output_fp, 'w') as out_f:
            for label, seq in parse_fasta(seqs_f):
                label = label.split()[0]
                if op(label):
                        out_f.write(">%s\n%s\n" % (label, seq))


def multiple_sequence_alignment(seqs_fp):
    """Perform multiple sequence alignment on FASTA
       file using MAFFT"""
    pass


def remove_chimeras_denovo_from_seqs(seqs_fp):
    """Remove chimeras de novo using UCHIME"""
    pass


def generate_biom_table(seqs_fp):
    """Generate BIOM table and representative FASTA set"""
    pass


def launch_workflow(seqs_fp, mapping_fp):
    """Launch full deblur workflow"""
    pass
