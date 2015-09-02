# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os.path import splitext, dirname, join, exists, basename
from collections import defaultdict
from datetime import datetime
from os import makedirs, stat, rename

from bfillings.vsearch import (vsearch_dereplicate_exact_seqs,
                               vsearch_chimera_filter_de_novo,
                               parse_uc_to_clusters)
from bfillings.sortmerna_v2 import (build_database_sortmerna,
                                    sortmerna_map)
from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta
from skbio import Alignment

from biom.table import Table


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
                     min_size=2,
                     dereplication_map=True):
    """Dereplicate FASTA sequences and remove singletons using VSEARCH

    Parameters
    ----------
    seqs_fp : string
        filepath to FASTA sequence file
    output_fp : string
        file path to dereplicated sequences (FASTA format)
    min_size : integer
        discard sequences with an abundance value smaller
        than integer
    dereplication_map: boolean
        output the dereplication map in .uc format
    """
    log_name = "%s.log" % splitext(output_fp)[0]

    vsearch_dereplicate_exact_seqs(
        fasta_filepath=seqs_fp,
        output_filepath=output_fp,
        output_uc=dereplication_map,
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


def remove_chimeras_denovo_from_seqs(seqs_fp, output_fp):
    """Remove chimeras de novo using UCHIME (VSEARCH implementation)

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


def parse_deblur_output(seqs_fp, derep_clusters):
    """ Parse deblur output file into an OTU map

    For each deblurred sequence in seqs_fp, use the sequence label to
    obtain all dereplicated sequence labels belonging to it
    (from derep_clusters) to create entries in a new dictionary where the keys
    are actual sequences (not the labels). Note not all sequences
    in derep_clusters will be in seqs_fp since they could have been removed in
    the artifact filtering step. 

    Parameters
    ----------
    seqs_fp: string
        file path to deblurred sequences
    derep_clusters: dictionary
        dictionary of dereplicated sequences map

    Returns
    -------
    clusters: dictionary
        dictionary of clusters including dereplicated sequence labels
    """
    clusters = {}
    # Replace representative sequence name with actual sequence in cluster
    msa_fa = Alignment.read(seqs_fp, format='fasta')
    for label, seq in Alignment.iteritems(msa_fa):
        cluster_id = label.split(';')[0]
        seq2 = str(seq.degap())
        if seq2 not in clusters:
            clusters[seq2] = []
        if cluster_id not in derep_clusters:
            raise ValueError(
                'Seed ID %s does not exist in .uc file' % cluster_id)
        else:
            clusters[seq2].extend(derep_clusters[cluster_id])
    return clusters


def generate_biom_data(clusters, delim='_'):
    """ Parse OTU map dictionary into a sparse dictionary
    {(cluster_idx,sample_idx):count}

    Parameters
    ----------
    clusters: dictionary
        OTU map as dictionary
    delim: string, optional
        delimiter for splitting sample and sequence IDs in sequence label

    Returns
    -------
    data: dictionary
        sprase dictionary {(otu_idx,sample_idx):count}
    cluster_ids: list
        list of cluster IDs
    sample_ids: list
        list of sample IDs
    """
    sample_ids = []
    sample_id_idx = {}
    data = defaultdict(int)
    sample_count = 0
    for cluster_idx, otu in enumerate(clusters):
        for seq_id in clusters[otu]:
            seq_id_sample = seq_id.split(delim)[0]
            try:
                sample_idx = sample_id_idx[seq_id_sample]
            except KeyError:
                sample_idx = sample_count
                sample_id_idx[seq_id_sample] = sample_idx
                sample_ids.append(seq_id_sample)
                sample_count += 1
            data[(cluster_idx, sample_idx)] += 1
    cluster_ids = clusters.keys()

    return data, cluster_ids, sample_ids


def generate_biom_table(seqs_fp,
                        uc_fp,
                        output_biom_fp,
                        delim='_'):
    """Generate BIOM table and representative FASTA set

    Parameters
    ----------
    seqs_fp: string
        file path to deblurred sequences
    uc_fp: string
        file path to dereplicated sequences map (.uc format) 
    output_biom_fp: string
        file path to output BIOM table
    delim: string, optional
        delimiter for splitting sample and sequence IDs in sequence label
    """
    # parse clusters in dereplicated sequences map (.uc format)
    derep_clusters = parse_uc_to_clusters(uc_fp)
    # parse clusters in deblur file, set observation ID to be the sequence
    deblur_clusters = parse_deblur_output(seqs_fp, derep_clusters)
    # create sparse dictionary of observation and sample ID counts
    data, otu_ids, sample_ids = generate_biom_data(deblur_clusters, delim)
    # build BIOM table
    return Table(data, otu_ids, sample_ids,
                   observation_metadata=None, 
                   sample_metadata=None, table_id=None,
                   generated_by="deblur",
                   create_date=datetime.now().isoformat())


def launch_workflow(seqs_fp, mapping_fp):
    """Launch full deblur workflow"""
    pass
