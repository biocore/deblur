# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, dirname, join, exists, basename
from collections import defaultdict
from datetime import datetime
from os import mkdir, stat, rename

from bfillings.vsearch import (vsearch_dereplicate_exact_seqs,
                               vsearch_chimera_filter_de_novo)
from bfillings.uclust import clusters_from_uc_file
from bfillings.sortmerna_v2 import (build_database_sortmerna,
                                    sortmerna_map)
from bfillings.mafft_v7 import align_unaligned_seqs
from skbio.util import remove_files
from skbio.parse.sequences import parse_fasta, parse_fastq
from skbio import Alignment
from biom.table import Table
from biom import load_table
from biom.util import biom_open, HAVE_H5PY

from deblur.deblurring import deblur


def trim_seqs(input_seqs, trim_len):
    """Trim FASTA sequences to specified length.

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
                     uc_output=True):
    """Dereplicate FASTA sequences and remove singletons using VSEARCH.

    Parameters
    ----------
    seqs_fp : string
        filepath to FASTA sequence file
    output_fp : string
        file path to dereplicated sequences (FASTA format)
    min_size : integer, optional
        discard sequences with an abundance value smaller
        than integer
    uc_output: boolean, optional
        output the dereplication map in .uc format
    """
    log_name = "%s.log" % splitext(output_fp)[0]

    vsearch_dereplicate_exact_seqs(
        fasta_filepath=seqs_fp,
        output_filepath=output_fp,
        output_uc=uc_output,
        minuniquesize=min_size,
        log_name=log_name)


def remove_artifacts_seqs(seqs_fp,
                          ref_fp,
                          working_dir,
                          ref_db_fp=None,
                          negate=False,
                          threads=1,
                          verbose=False):
    """Remove artifacts from FASTA file using SortMeRNA.

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    ref_fp: tuple
        file path(s) to FASTA database file
    working_dir: string
        working directory path
    ref_db_fp: tuple, optional
        file path(s) to indexed FASTA database
    negate: boolean, optional
        if True, discard all input sequences aligning
        to reference database
    threads: integer, optional
        number of threads to use for SortMeRNA
    verbose: boolean, optional
        If true, output SortMeRNA errors
    """
    output_fp = join(working_dir,
                     "%s.no_artifacts" % basename(seqs_fp))
    aligned_seq_ids = set()
    files_to_remove = []

    for i, db in enumerate(ref_fp):
        # create working directory for each
        # reference database
        db_dir_base = splitext(basename(db))[0]
        db_dir = join(working_dir, db_dir_base)
        if not exists(db_dir):
            mkdir(db_dir)

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
            if verbose:
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
    return output_fp


def multiple_sequence_alignment(seqs_fp, threads=1):
    """Perform multiple sequence alignment on FASTA file using MAFFT.

    Parameters
    ----------
    seqs_fp: string
        filepath to FASTA file for multiple sequence alignment
    threads: integer, optional
        number of threads to use

    Returns
    -------
    Alignment object
        The aligned sequences.

    See Also
    --------
    skbio.Alignment
    """
    return align_unaligned_seqs(seqs_fp=seqs_fp, params={'--thread': threads})


def remove_chimeras_denovo_from_seqs(seqs_fp, working_dir):
    """Remove chimeras de novo using UCHIME (VSEARCH implementation).

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    output_fp: string
        file path to store chimera-free results
    """
    output_fp = join(
        working_dir, "%s.no_chimeras" % basename(seqs_fp))
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
    return output_fp


def parse_deblur_output(seqs_fp, derep_clusters):
    """ Parse deblur output file into an OTU map.

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

    Notes
    -----
    For each deblurred sequence in seqs_fp, use the sequence label to
    obtain all dereplicated sequence labels belonging to it
    (from derep_clusters) to create entries in a new dictionary where the keys
    are actual sequences (not the labels). Note not all sequences
    in derep_clusters will be in seqs_fp since they could have been removed in
    the artifact filtering step.
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
    """ Parse OTU map dictionary into a sparse dictionary.

    Parameters
    ----------
    clusters: dictionary
        OTU map as dictionary
    delim: string, optional
        delimiter for splitting sample and sequence IDs in sequence label
        default: '_'

    Returns
    -------
    data: dictionary
        sprase dictionary {(otu_idx,sample_idx):count}
    cluster_ids: list
        list of cluster IDs
    sample_ids: list
        list of sample IDs

    Notes
    -----
    Sparse dictionary format is {(cluster_idx,sample_idx):count}.
    This function is based on QIIME's parse_otu_map() function found at
    https://github.com/biocore/qiime/blob/master/qiime/parse.py.

    QIIME is a GPL project, but we obtained permission from the authors of this
    function to port it to deblur (and keep it under deblur's BSD license).
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


def split_sequence_file_on_sample_ids_to_files(seqs,
                                               filetype,
                                               outdir):
    """Split FASTA/Q file on sample IDs.

    Parameters
    ----------
    seqs: file handler
        file handler to demultiplexed FASTA/Q file
    filetype: string
        file type, FASTA or FASTQ
    outdir: string
        dirpath to output split FASTA files
    """
    if filetype == 'fasta':
        parser = parse_fasta
        ext = '.fna'
    else:
        parser = parse_fastq
        ext = '.fq'
    outputs = {}
    for bits in parser(seqs):
        sample = bits[0].split('_', 1)[0]
        if sample not in outputs:
            outputs[sample] = open(join(outdir, sample + ext), 'w')
        outputs[sample].write(">%s\n%s\n" % (bits[0], bits[1]))


def generate_biom_table(seqs_fp,
                        uc_fp,
                        delim='_'):
    """Generate BIOM table and representative FASTA set.

    Parameters
    ----------
    seqs_fp: string
        file path to deblurred sequences
    uc_fp: string
        file path to dereplicated sequences map (.uc format)
    delim: string, optional
        delimiter for splitting sample and sequence IDs in sequence label
        default: '_'

    Returns
    -------
    deblur_clusters: dictionary
        dictionary of clusters including dereplicated sequence labels
    Table: biom.table
        an instance of a BIOM table
    """
    # parse clusters in dereplicated sequences map (.uc format)
    with open(uc_fp, 'U') as uc_f:
        derep_clusters, failures, seeds = clusters_from_uc_file(uc_f)
    # parse clusters in deblur file, set observation ID to be the sequence
    deblur_clusters = parse_deblur_output(seqs_fp, derep_clusters)
    # create sparse dictionary of observation and sample ID counts
    data, otu_ids, sample_ids = generate_biom_data(deblur_clusters, delim)
    # build BIOM table
    return deblur_clusters, Table(data, otu_ids, sample_ids,
                                  observation_metadata=None,
                                  sample_metadata=None, table_id=None,
                                  generated_by="deblur",
                                  create_date=datetime.now().isoformat())


def write_biom_table(table, biom_fp):
    """Write BIOM table to file.

    Parameters
    ----------
    table: biom.table
        an instance of a BIOM table
    biom_fp: string
        filepath to output BIOM table
    """
    with biom_open(biom_fp, 'w') as f:
        if HAVE_H5PY:
            table.to_hdf5(h5grp=f, generated_by="deblur")
        else:
            table.to_json(direct_io=f, generated_by="deblur")


def merge_otu_tables(output_fp, all_tables):
    """Merge multiple BIOM tables into one.

    Code taken from QIIME's merge_otu_tables.py

    Parameters
    ----------
    output_fp: string
        filepath to final BIOM table
    all_tables: list
        list of filepaths for BIOM tables to merge
    """
    master = load_table(all_tables[0])
    for input_fp in all_tables[1:]:
        master = master.merge(load_table(input_fp))
    write_biom_table(master, output_fp)


def launch_workflow(seqs_fp, output_dir, read_error, mean_error, error_dist,
                    indel_prob, indel_max, trim_length, min_size, ref_fp,
                    ref_db_fp, negate, threads=1, delim='_'):
    """Launch full deblur workflow.

    Parameters
    ----------
    seqs_fp: string
        post split library sequences for debluring
    output_dir: string
        dirpath to output file
    read_error: float
        read error rate
    mean_error: float
        mean error for original sequence estimate
    error_dist: list
        list of error probabilities for each hamming distance
    indel_prob: float
        insertion/deletion (indel) probability
    indel_max: integer
        maximal indel number
    trim_length: integer
        sequence trim length
    min_size: integer
        upper limit on sequence abundance (discard sequences below limit)
    ref_fp: tuple
        filepath(s) to FASTA reference database for artifact removal
    ref_db_fp: tuple
        filepath(s) to SortMeRNA indexed database for artifact removal
    negate: boolean
        discard all sequences aligning to the ref_fp database
    threads: integer, optional
        number of threads to use for SortMeRNA
    delim: string, optional
        delimiter in FASTA labels to separate sample ID from sequence ID

    Return
    ------
    biom_fp: string
        filepath to BIOM table
    """
    # Create working directory
    working_dir = join(output_dir, "deblur_working_dir")
    if not exists(working_dir):
        mkdir(working_dir)
    # Step 1: Trim sequences to specified length
    output_trim_fp = join(working_dir, "%s.trim" % basename(seqs_fp))
    with open(seqs_fp, 'U') as in_f, open(output_trim_fp, 'w') as out_f:
        for label, seq in trim_seqs(
                input_seqs=parse_fasta(in_f), trim_len=trim_length):
            out_f.write(">%s\n%s\n" % (label, seq))
    # Step 2: Dereplicate sequences
    output_derep_fp = join(working_dir,
                           "%s.derep" % basename(output_trim_fp))
    dereplicate_seqs(seqs_fp=output_trim_fp,
                     output_fp=output_derep_fp,
                     min_size=min_size,
                     uc_output=True)
    # Step 3: Remove artifacts
    output_artif_fp = remove_artifacts_seqs(seqs_fp=output_derep_fp,
                          ref_fp=ref_fp,
                          working_dir=working_dir,
                          ref_db_fp=ref_db_fp,
                          negate=negate,
                          threads=threads)
    # Step 4: Multiple sequence alignment
    output_msa_fp = join(working_dir,
                         "%s.msa" % basename(output_artif_fp))
    with open(output_msa_fp, 'w') as f:
        alignment = multiple_sequence_alignment(seqs_fp=output_artif_fp,
                                                threads=threads)
        f.write(alignment.to_fasta())
    # Step 5: Launch deblur
    output_deblur_fp = join(working_dir,
                            "%s.deblur" % basename(output_msa_fp))
    with open(output_deblur_fp, 'w') as f:
        seqs = deblur(parse_fasta(output_msa_fp), read_error, mean_error,
                      error_dist, indel_prob, indel_max)
        for s in seqs:
            # remove '-' from aligned sequences
            s.sequence = s.sequence.replace('-', '')
            f.write(s.to_fasta())
    # Step 6: Chimera removal
    output_no_chimeras_fp = remove_chimeras_denovo_from_seqs(
        output_deblur_fp, working_dir)
    # Step 7: Generate BIOM table
    deblur_clrs, table = generate_biom_table(seqs_fp=output_no_chimeras_fp,
                                             uc_fp="%s.uc" % output_trim_fp,
                                             delim=delim)
    # Step 8: Write BIOM table to file
    if table.is_empty():
        raise ValueError(
            "Attempting to write an empty BIOM table.")
    biom_fp = join(working_dir, "%s.biom" % basename(seqs_fp))
    write_biom_table(table, biom_fp)

    return biom_fp
