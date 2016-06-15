# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, join, basename, isfile, split
from os import listdir
from collections import defaultdict
from datetime import datetime
from os import stat
import logging
import re
import scipy
import numpy as np
import subprocess
import time

#from bfillings.vsearch import vsearch_dereplicate_exact_seqs
#from bfillings.sortmerna_v2 import (build_database_sortmerna,
#                                    sortmerna_map)
# from bfillings.mafft_v7 import align_unaligned_seqs
from skbio.parse.sequences import parse_fasta
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
                     use_log=False):
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
    use_log: boolean, optional
        save the vsearch logfile as well (to output_fp.log)
        default=False
    """
    logger = logging.getLogger(__name__)
    logger.debug('dereplicate seqs file %s' % seqs_fp)

    log_name = "%s.log" % output_fp

    params = ['vsearch','--derep_fulllength', seqs_fp, '--output', output_fp, '--sizeout']
    params.extend(['--fasta_width', '0', '--threads','1', '--minuniquesize', str(min_size)])
    if use_log:
        params.extend(['--log', log_name])
    res=subprocess.call(params)
    if not res==0:
        logger.error('Problem running vsearch dereplication on file %s' % seqs_fp)
        logger.debug('parameters used:\n%s' % params)
        return


def build_index_sortmerna(ref_fp, working_dir):
    """Build a SortMeRNA index for all reference databases.

    Parameters
    ----------
    ref_fp: tuple
        filepaths to FASTA reference databases
    working_dir: string
        working directory path where to store the indexed database

    Returns
    -------
    all_db: tuple
        filepaths to SortMeRNA indexed reference databases
    """
    logger = logging.getLogger(__name__)
    logger.info('build_index_sortmerna files %s to'
                ' dir %s' % (ref_fp, working_dir))

    all_db = []
    for db in ref_fp:
        fasta_dir, fasta_filename = split(db)
        index_basename = splitext(fasta_filename)[0]
        db_output = join(working_dir, index_basename)
        logger.debug('processing file %s into location %s' % (db, db_output))
        params = ['indexdb_rna', '--ref', '%s,%s' % (db,db_output), '--tmpdir', working_dir]
        res=subprocess.call(params)
        if not res==0:
            logger.error('Problem running indexdb_rna on file %s to dir %s. database not indexed' % (db, db_output))
            continue
        logger.debug('file %s indexed' % db)
        all_db.append(db_output)
    return all_db


def remove_artifacts_seqs(seqs_fp,
                          ref_fp,
                          working_dir,
                          ref_db_fp,
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
    ref_db_fp: tuple
        file path(s) to indexed FASTA database
    negate: boolean, optional
        if True, discard all input sequences aligning
        to reference database
    threads: integer, optional
        number of threads to use for SortMeRNA
    verbose: boolean, optional
        If true, output SortMeRNA errors
    """
    logger = logging.getLogger(__name__)
    logger.info('remove_artifacts_seqs file %s' % seqs_fp)

    if stat(seqs_fp).st_size == 0:
        logger.warn('file %s has size 0, continuing' % seqs_fp)
        return

    output_fp = join(working_dir,
                     "%s.no_artifacts" % basename(seqs_fp))
    blast_output = join(working_dir,
                        '%s.sortmerna' % basename(seqs_fp))
    aligned_seq_ids = set()
    for i, db in enumerate(ref_fp):
        logger.debug('running on ref_fp %s working dir %s refdb_fp %s seqs %s'
                     % (db, working_dir, ref_db_fp[i], seqs_fp))
        # run SortMeRNA
        params = ['sortmerna', '--reads', seqs_fp, '--ref', '%s,%s' % (db, ref_db_fp[i])]
        params.extend(['--aligned', blast_output])
        params.extend(['--blast','3','--best','1','--print_all_reads'])
        params.extend(['-v'])

        with open(blast_output + '.log', "w") as f:
            res=subprocess.call(params, stdout=f)
        if not res==0:
            logger.critical('sortmerna error on file %s' % seqs_fp)
            return output_fp

        bfl = open(blast_output+'.blast','r')
        for line in bfl:
            line = line.strip().split('\t')
            if line[1] == '*':
                continue
            else:
                aligned_seq_ids.add(line[0])
        bfl.close()

    if negate:
        def op(x): return x not in aligned_seq_ids
    else:
        def op(x): return x in aligned_seq_ids

    # if negate = False, only output sequences
    # matching to at least one of the databases
    totalseqs = 0
    okseqs = 0
    badseqs = 0
    with open(seqs_fp, 'U') as seqs_f:
        with open(output_fp, 'w') as out_f:
            for label, seq in parse_fasta(seqs_f):
                totalseqs += 1
                label = label.split()[0]
                if op(label):
                        out_f.write(">%s\n%s\n" % (label, seq))
                        okseqs += 1
                else:
                        badseqs += 1
    logger.info('total sequences %d, passing sequences %d, '
                'failing sequences %d' % (totalseqs, okseqs, badseqs))
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
        The aligned sequences. False if file does not exist

    See Also
    --------
    skbio.Alignment
    """
    logger = logging.getLogger(__name__)
    logger.debug('multiple_sequence_alignment seqs file %s' % seqs_fp)

    if stat(seqs_fp).st_size == 0:
        logger.info('msa failed. file %s has no reads' % seqs_fp)
        return False
    try:
        #        aligned_seqs = align_unaligned_seqs(seqs_fp=seqs_fp,
        #                                            params={'--thread': threads})
        #        return aligned_seqs
        msa_fp = seqs_fp + '.msa'
        params = ['mafft', '--quiet', '--parttree', '--auto', seqs_fp]
        with open(msa_fp, "w") as f:
            subprocess.call(params, stdout=f)
        return msa_fp
    except:
        # alignment can fail if only 1 sequence present
        logger.info('msa failed for file %s (maybe only 1 read?)' % seqs_fp)
        return False


def remove_chimeras_denovo_from_seqs(seqs_fp, working_dir):
    """Remove chimeras de novo using UCHIME (VSEARCH implementation).

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    output_fp: string
        file path to store chimera-free results

    Returns
    -------
    output_fp
        the chimera removed fasta file name
    """
    logger = logging.getLogger(__name__)
    logger.debug('remove_chimeras_denovo_from_seqs seqs file %s'
                 'to working dir %s' % (seqs_fp, working_dir))

    output_fp = join(
        working_dir, "%s.no_chimeras" % basename(seqs_fp))

    # we use the parameters dn=0.000001, xn=1000, minh=10000000
    # so 1 mismatch in the A/B region will cancel it being labeled as chimera
    # and ~3 unique reads in each region will make it a chimera if
    # no mismatches
    params = ['vsearch', '--uchime_denovo', seqs_fp]
    params.extend(['--nonchimeras', output_fp])
    params.extend(['-dn', '0.000001', '-xn', '1000'])
    params.extend(['-minh', '10000000', '--mindiffs', '5'])
    params.extend(['--fasta_width', '0'])
    with open(output_fp+'.log', "w") as f:
        subprocess.call(params, stdout=f, stderr=f)
    return output_fp


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
    logger = logging.getLogger(__name__)
    logger.debug('generate_biom_data')

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
                                               outdir):
    """Split FASTA file on sample IDs.

    Parameters
    ----------
    seqs: file handler
        file handler to demultiplexed FASTA file
    outdir: string
        dirpath to output split FASTA files
    """
    logger = logging.getLogger(__name__)
    logger.info('split_sequence_file_on_sample_ids_to_files'
                ' for file %s into dir %s' % (seqs, outdir))

    outputs = {}
    for bits in parse_fasta(seqs):
        sample = bits[0].split('_', 1)[0]
        if sample not in outputs:
            outputs[sample] = open(join(outdir, sample + '.fasta'), 'w')
        outputs[sample].write(">%s\n%s\n" % (bits[0], bits[1]))
    for sample in outputs:
        outputs[sample].close()
    logger.info('split to %d files' % len(outputs))


def write_biom_table(table, biom_fp):
    """Write BIOM table to file.

    Parameters
    ----------
    table: biom.table
        an instance of a BIOM table
    biom_fp: string
        filepath to output BIOM table
    """
    logger = logging.getLogger(__name__)
    logger.debug('write_biom_table to file %s' % biom_fp)
    with biom_open(biom_fp, 'w') as f:
        if HAVE_H5PY:
            logger.debug('saving with h5py')
            table.to_hdf5(h5grp=f, generated_by="deblur")
            logger.debug('wrote to hdf5 file %s' % biom_fp)
        else:
            logger.debug('no h5py. saving to json')
            table.to_json(direct_io=f, generated_by="deblur")
            logger.debug('wrote to json file %s' % biom_fp)


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
    logger = logging.getLogger(__name__)
    logger.debug('merge_otu_tables for file %d tables'
                 ' into file %s' % (len(all_tables), output_fp))

    master = load_table(all_tables[0])
    for input_fp in all_tables[1:]:
        master = master.merge(load_table(input_fp))
    write_biom_table(master, output_fp)


def get_files_for_table(input_dir,
                        file_end='.fasta.trim.derep.no_artifacts'
                        '.msa.deblur.no_chimeras'):
    """Get a list of files to add to the output table

    Parameters:
    -----------
    input_dir : string
        name of the directory containing the deblurred fasta files
    file_end : string
        the ending of all the fasta files to be added to the table
        (default '.trim.derep.no_artifacts.msa.deblur.no_chimeras')

    Returns
    -------
    names : list of tuples of (string,string)
        list of tuples of:
            name of fasta files to be added to the biom table
            sampleid (file names without the file_end and path)
    """
    logger = logging.getLogger(__name__)
    logger.debug('get_files_for_table input dir %s, '
                 'file-ending %s' % (input_dir, file_end))

    names = []
    for cfile in listdir(input_dir):
        cname = join(input_dir, cfile)
        if not isfile(cname):
            continue
        if len(cfile) < len(file_end):
            continue
        if cfile[-len(file_end):] == file_end:
            sampleid = cfile[:-len(file_end)]
            names.append((cname, sampleid))
    logger.debug('found %d files' % len(names))
    return names


def create_otu_table(output_fp, deblurred_list):
    """Create a biom table out of all files in a directory

    Parameters
    ----------
    output_fp : string
        filepath to final BIOM table
    deblurred_list : list of string
        list of file names (including path) of all deblurred
        fasta files to add to the table
    """
    logger = logging.getLogger(__name__)
    logger.debug('create_otu_table for %d samples, '
                 'into output table %s' % (len(deblurred_list), output_fp))

    seqdict = {}
    seqlist = []
    sampdict = {}
    samplist = []
    obs = scipy.sparse.dok_matrix((1E9, 1E6), dtype=np.int)
    for (cfilename, csampleid) in deblurred_list:
        if csampleid in sampdict:
            logger.error('sample %s already in table!' % csampleid)
            continue
        sampdict[csampleid] = len(sampdict)-1
        samplist.append(csampleid)
        csampidx = len(sampdict)-1
        for chead, cseq in parse_fasta(open(cfilename, 'U')):
            if cseq not in seqdict:
                seqdict[cseq] = len(seqlist)
                seqlist.append(cseq)
            cseqidx = seqdict[cseq]
            cfreq = float(re.search('(?<=size=)\w+', chead).group(0))
            obs[cseqidx, csampidx] = cfreq
    logger.debug('loaded %d samples, %d unique sequences'
                 % (len(samplist), len(seqlist)))
    obs.resize((len(seqlist), len(samplist)))
    table = Table(obs, seqlist, samplist,
                  observation_metadata=None,
                  sample_metadata=None, table_id=None,
                  generated_by="deblur",
                  create_date=datetime.now().isoformat())
    logger.debug('converted to biom table')
    write_biom_table(table, output_fp)
    logger.debug('saved to biom file %s' % output_fp)


def launch_workflow(seqs_fp, working_dir, read_error, mean_error, error_dist,
                    indel_prob, indel_max, trim_length, min_size, ref_fp,
                    ref_db_fp, negate, threads=1, delim='_'):
    """Launch full deblur workflow.

    Parameters
    ----------
    seqs_fp: string
        post split library sequences for debluring
    working_dir: string
        working directory path
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
    output_no_chimers_fp : string
        filepath to fasta file with no chimeras of False if error encountered
    """
    logger = logging.getLogger(__name__)
    logger.info('--------------------------------------------------------')
    logger.info('launch_workflow for file %s' % seqs_fp)

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
                     min_size=min_size)
    # Step 3: Remove artifacts
    output_artif_fp = remove_artifacts_seqs(seqs_fp=output_derep_fp,
                                            ref_fp=ref_fp,
                                            working_dir=working_dir,
                                            ref_db_fp=ref_db_fp,
                                            negate=negate,
                                            threads=threads)
    if not output_artif_fp:
        logger.debug('remove artifacts failed, aborting')
        return
    # Step 4: Multiple sequence alignment
    output_msa_fp = join(working_dir,
                         "%s.msa" % basename(output_artif_fp))
    # with open(output_msa_fp, 'w') as f:
    #     alignment = multiple_sequence_alignment(seqs_fp=output_artif_fp,
    #                                             threads=threads)
    #     if not alignment:
    #         logger.debug('msa failed. aborting')
    #         return False
    #     f.write(alignment.to_fasta())
    alignment = multiple_sequence_alignment(seqs_fp=output_artif_fp,
                                            threads=threads)
    if not alignment:
        logger.debug('msa failed. aborting')
        return False

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
    return output_no_chimeras_fp


def start_log(level=logging.DEBUG, filename=''):
    """
    start the logger for the run

    Parameters
    ----------
    levele : logging.DEBUG, logging.INFO etc. for the log level.
    filename : str
      name of the filename to save the log to or
      empty (default) to use deblur.log.TIMESTAMP
    """
    if not filename:
        tstr = time.ctime()
        tstr = tstr.replace(' ', '.')
        tstr = tstr.replace(':', '.')
        filename = 'deblur.log.%s' % tstr
    logging.basicConfig(filename=filename, level=level,
                        format='%(levelname)s:%(asctime)s:%(message)s')
    logger = logging.getLogger(__name__)
    logger.debug('deblurring started')
