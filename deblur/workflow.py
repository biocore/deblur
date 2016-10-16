# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import splitext, join, basename, isfile, split
from datetime import datetime
from os import stat
from glob import glob
import logging
import re
import scipy
import numpy as np
import subprocess
import time
import warnings
import io

import skbio
from biom.table import Table
from biom.util import biom_open
from biom import load_table

from deblur.deblurring import deblur


sniff_fasta = skbio.io.io_registry.get_sniffer('fasta')
sniff_fastq = skbio.io.io_registry.get_sniffer('fastq')


def sequence_generator(input_fp):
    """Yield (id, sequence) from an input file

    Parameters
    ----------
    input_fp : filepath
        A filepath, which can be any valid fasta or fastq file within the
        limitations of scikit-bio's IO registry.

    Notes
    -----
    The use of this method is a stopgap to replicate the existing `parse_fasta`
    functionality while at the same time allowing for fastq support.

    Raises
    ------
    skbio.io.FormatIdentificationWarning
        If the format of the input file cannot be determined.

    Returns
    -------
    (str, str)
        The ID and sequence.

    """
    kw = {}
    if sniff_fasta(input_fp)[0]:
        format = 'fasta'
    elif sniff_fastq(input_fp)[0]:
        format = 'fastq'

        # WARNING: the variant is currently forced to illumina 1.8 as the
        # quality scores are _not_ used in downstream processing. However, if
        # in the future, quality scores are to be interrogated, it is critical
        # that this variant parameter be exposed to the user at the command
        # line. The list of allowable paramters can be found here:
        # http://scikit-bio.org/docs/latest/generated/skbio.io.format.fastq.html#format-parameters
        kw['variant'] = 'illumina1.8'
    else:
        # usually happens when the fasta file is empty
        # so need to return no sequences (and warn)
        logger = logging.getLogger(__name__)
        msg = "input file %s does not appear to be FASTA or FASTQ" % input_fp
        logger.warn(msg)
        warnings.warn(msg, UserWarning)
        return

    # some of the test code is using file paths, some is using StringIO.
    if isinstance(input_fp, io.TextIOBase):
        input_fp.seek(0)

    for record in skbio.read(input_fp, format=format, **kw):
        yield (record.metadata['id'], str(record))


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
    # counters for the number of trimmed and total sequences
    logger = logging.getLogger(__name__)

    okseqs = 0
    totseqs = 0

    for label, seq in input_seqs:
        totseqs += 1
        if len(seq) >= trim_len:
            okseqs += 1
            yield label, seq[:trim_len]

    if okseqs < 0.01*totseqs:
        logger = logging.getLogger(__name__)
        errmsg = 'Vast majority of sequences (%d / %d) are shorter ' \
                 'than the trim length (%d). ' \
                 'Are you using the correct -t trim length?' \
                 % (totseqs-okseqs, totseqs, trim_len)
        logger.warn(errmsg)
        warnings.warn(errmsg, UserWarning)
    else:
        logger.debug('trimmed to length %d (%d / %d remaining)'
                     % (trim_len, okseqs, totseqs))


def dereplicate_seqs(seqs_fp,
                     output_fp,
                     min_size=2,
                     use_log=False,
                     threads=1):
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
    threads : int, optional
        number of threads to use (0 for all available)
    """
    logger = logging.getLogger(__name__)
    logger.info('dereplicate seqs file %s' % seqs_fp)

    log_name = "%s.log" % output_fp

    params = ['vsearch', '--derep_fulllength', seqs_fp,
              '--output', output_fp, '--sizeout',
              '--fasta_width', '0', '--minuniquesize', str(min_size),
              '--quiet', '--threads', str(threads)]
    if use_log:
        params.extend(['--log', log_name])
    sout, serr, res = _system_call(params)
    if not res == 0:
        logger.error('Problem running vsearch dereplication on file %s' %
                     seqs_fp)
        logger.debug('parameters used:\n%s' % params)
        logger.debug('stdout: %s' % sout)
        logger.debug('stderr: %s' % serr)
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
        params = ['indexdb_rna', '--ref', '%s,%s' %
                  (db, db_output), '--tmpdir', working_dir]
        sout, serr, res = _system_call(params)
        if not res == 0:
            logger.error('Problem running indexdb_rna on file %s to dir %s. '
                         'database not indexed' % (db, db_output))
            logger.debug('stdout: %s' % sout)
            logger.debug('stderr: %s' % serr)
            logger.critical('execution halted')
            raise RuntimeError('Cannot index database file %s' % db)
        logger.debug('file %s indexed' % db)
        all_db.append(db_output)
    return all_db


def filter_minreads_samples_from_table(table, minreads=1, inplace=True):
    """Filter samples from biom table that have less than
    minreads reads total

    Paraneters
    ----------
    table : biom.Table
        the biom table to filter
    minreads : int (optional)
        the minimal number of reads in a sample in order to keep it
    inplace : bool (optional)
        if True, filter the biom table in place, if false create a new copy

    Returns
    -------
    table : biom.Table
        the filtered biom table
    """
    logger = logging.getLogger(__name__)
    logger.debug('filter_minreads_started. minreads=%d' % minreads)
    samp_sum = table.sum(axis='sample')
    samp_ids = table.ids(axis='sample')
    bad_samples = samp_ids[samp_sum < minreads]
    if len(bad_samples) > 0:
        logger.warn('removed %d samples with reads per sample<%d'
                    % (len(bad_samples), minreads))
        table = table.filter(bad_samples, axis='sample',
                             inplace=inplace, invert=True)
    else:
        logger.debug('all samples contain > %d reads' % minreads)
    return table


def remove_artifacts_from_biom_table(table_filename,
                                     fasta_filename,
                                     ref_fp,
                                     biom_table_dir,
                                     ref_db_fp,
                                     threads=1,
                                     verbose=False,
                                     sim_thresh=None,
                                     coverage_thresh=None):
    """Remove artifacts from a biom table using SortMeRNA

    Parameters
    ----------
    table : str
        name of the biom table file
    fasta_filename : str
        the fasta file containing all the sequences of the biom table
    """
    logger = logging.getLogger(__name__)
    logger.info('getting 16s sequences from the biom table')

    # remove artifacts from the fasta file. output is in clean_fp fasta file
    clean_fp = remove_artifacts_seqs(fasta_filename, ref_fp,
                                     working_dir=biom_table_dir,
                                     ref_db_fp=ref_db_fp,
                                     negate=False, threads=threads,
                                     verbose=verbose,
                                     sim_thresh=sim_thresh,
                                     coverage_thresh=coverage_thresh)
    logger.debug('removed artifacts from sequences input %s'
                 ' to output %s' % (fasta_filename, clean_fp))

    # read the clean fasta file
    good_seqs = {s for _, s in sequence_generator(clean_fp)}
    logger.debug('loaded %d sequences from cleaned biom table'
                 ' fasta file' % len(good_seqs))

    logger.debug('loading biom table %s' % table_filename)
    table = load_table(table_filename)

    # filter and save the artifact biom table
    artifact_table = table.filter(list(good_seqs),
                                  axis='observation', inplace=False,
                                  invert=True)
    # remove the samples with 0 reads
    filter_minreads_samples_from_table(artifact_table)
    output_artifact_fp = join(biom_table_dir, 'final.only-non16s.biom')
    write_biom_table(artifact_table, output_artifact_fp)
    logger.info('wrote artifact only filtered biom table to %s'
                % output_artifact_fp)

    # filter and save the only 16s biom table
    table.filter(list(good_seqs), axis='observation')
    filter_minreads_samples_from_table(table)
    output_fp = join(biom_table_dir, 'final.only-16s.biom')
    write_biom_table(table, output_fp)
    logger.info('wrote 16s filtered biom table to %s' % output_fp)


def remove_artifacts_seqs(seqs_fp,
                          ref_fp,
                          working_dir,
                          ref_db_fp,
                          negate=False,
                          threads=1,
                          verbose=False,
                          sim_thresh=None,
                          coverage_thresh=None):
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
    sim_thresh: float, optional
        The minimal similarity threshold (between 0 and 1)
        for keeping the sequence
        if None, the default values used are 0.65 for negate=False,
        0.95 for negate=True
    coverage_thresh: float, optional
        The minimal coverage threshold (between 0 and 1)
        for alignments for keeping the sequence
        if Nonr, the default values used are 0.3 for negate=False,
        0.95 for negate=True
    """
    logger = logging.getLogger(__name__)
    logger.info('remove_artifacts_seqs file %s' % seqs_fp)

    if stat(seqs_fp).st_size == 0:
        logger.warn('file %s has size 0, continuing' % seqs_fp)
        return

    if coverage_thresh is None:
        if negate:
            coverage_thresh = 0.95
        else:
            coverage_thresh = 0.3

    if sim_thresh is None:
        if negate:
            sim_thresh = 0.95
        else:
            sim_thresh = 0.65

    output_fp = join(working_dir,
                     "%s.no_artifacts" % basename(seqs_fp))
    blast_output = join(working_dir,
                        '%s.sortmerna' % basename(seqs_fp))
    aligned_seq_ids = set()
    for i, db in enumerate(ref_fp):
        logger.debug('running on ref_fp %s working dir %s refdb_fp %s seqs %s'
                     % (db, working_dir, ref_db_fp[i], seqs_fp))
        # run SortMeRNA
        params = ['sortmerna', '--reads', seqs_fp, '--ref', '%s,%s' %
                  (db, ref_db_fp[i]),
                  '--aligned', blast_output, '--blast', '3', '--best', '1',
                  '--print_all_reads', '-v', '-e', '10']

        sout, serr, res = _system_call(params)
        if not res == 0:
            logger.error('sortmerna error on file %s' % seqs_fp)
            logger.error('stdout : %s' % sout)
            logger.error('stderr : %s' % serr)
            return output_fp

        with open('%s.blast' % blast_output, 'r') as bfl:
            for line in bfl:
                line = line.strip().split('\t')
                # if * means no match
                if line[1] == '*':
                    continue
                # check if % identity[2] and coverage[13] are large enough
                # note e-value is [10]
                if negate:
                    if (float(line[2]) >= sim_thresh*100) and \
                       (float(line[13]) >= coverage_thresh*100):
                        aligned_seq_ids.add(line[0])
                else:
                    if float(line[10]) <= 10:
                        aligned_seq_ids.add(line[0])

    if negate:
        def op(x): return x not in aligned_seq_ids
    else:
        def op(x): return x in aligned_seq_ids

    # if negate = False, only output sequences
    # matching to at least one of the databases
    totalseqs = 0
    okseqs = 0
    badseqs = 0
    with open(output_fp, 'w') as out_f:
        for label, seq in sequence_generator(seqs_fp):
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
        number of threads to use. 0 to use all threads

    Returns
    -------
    msa_fp : str
        name of output alignment file or None if error encountered
    """
    logger = logging.getLogger(__name__)
    logger.info('multiple_sequence_alignment seqs file %s' % seqs_fp)

    # for mafft we use -1 to denote all threads and not 0
    if threads == 0:
        threads = -1

    if stat(seqs_fp).st_size == 0:
        logger.warning('msa failed. file %s has no reads' % seqs_fp)
        return None
    msa_fp = seqs_fp + '.msa'
    params = ['mafft', '--quiet', '--parttree', '--auto',
              '--thread', str(threads), seqs_fp]
    sout, serr, res = _system_call(params, stdoutfilename=msa_fp)
    if not res == 0:
        logger.info('msa failed for file %s (maybe only 1 read?)' % seqs_fp)
        logger.debug('stderr : %s' % serr)
        return None
    return msa_fp


def remove_chimeras_denovo_from_seqs(seqs_fp, working_dir, threads=1):
    """Remove chimeras de novo using UCHIME (VSEARCH implementation).

    Parameters
    ----------
    seqs_fp: string
        file path to FASTA input sequence file
    output_fp: string
        file path to store chimera-free results
    threads : int
        number of threads (0 for all cores)

    Returns
    -------
    output_fp
        the chimera removed fasta file name
    """
    logger = logging.getLogger(__name__)
    logger.info('remove_chimeras_denovo_from_seqs seqs file %s'
                'to working dir %s' % (seqs_fp, working_dir))

    output_fp = join(
        working_dir, "%s.no_chimeras" % basename(seqs_fp))

    # we use the parameters dn=0.000001, xn=1000, minh=10000000
    # so 1 mismatch in the A/B region will cancel it being labeled as chimera
    # and ~3 unique reads in each region will make it a chimera if
    # no mismatches
    params = ['vsearch', '--uchime_denovo', seqs_fp,
              '--nonchimeras', output_fp,
              '-dn', '0.000001', '-xn', '1000',
              '-minh', '10000000', '--mindiffs', '5',
              '--fasta_width', '0', '--threads', str(threads)]
    sout, serr, res = _system_call(params)
    if not res == 0:
        logger.error('problem with chimera removal for file %s' % seqs_fp)
        logger.debug('stdout : %s' % sout)
        logger.debug('stderr : %s' % serr)
    return output_fp


def sample_id_from_read_id(readid):
    """Get SampleID from the split_libraries_fastq.py output
    fasta file read header

    Parameters
    ----------
    readid : str
        the fasta file read name

    Returns
    -------
    sampleid : str
        the sample id
    """

    # get the sampleid_readid field
    sampleread = readid.split(' ')[0]

    # get the sampleid field
    sampleid = sampleread.rsplit('_', 1)[0]
    return sampleid


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

    for bits in sequence_generator(seqs):
        sample = sample_id_from_read_id(bits[0])

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
        table.to_hdf5(h5grp=f, generated_by="deblur")
        logger.debug('wrote to BIOM file %s' % biom_fp)


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
        (default '.fasta.trim.derep.no_artifacts.msa.deblur.no_chimeras')

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
    for cfile in glob(join(input_dir, "*%s" % file_end)):
        if not isfile(cfile):
            continue
        sample_id = basename(cfile)[:-len(file_end)]
        names.append((cfile, sample_id))

    logger.debug('found %d files' % len(names))
    return names


def create_otu_table(output_fp, deblurred_list,
                     outputfasta_fp=None, minreads=0):
    """Create a biom table out of all files in a directory

    Parameters
    ----------
    output_fp : string
        filepath to final BIOM table
    deblurred_list : list of (str, str)
        list of file names (including path), sampleid of all deblurred
        fasta files to add to the table
    outputfasta_fp : str, optional
        name of output fasta file (of all sequences in the table) or None
        to not write
    minreads : int, optional
        minimal number of reads per bacterial sequence in order to write
        it to the biom table and fasta file or 0 to write all
    """
    logger = logging.getLogger(__name__)
    logger.info('create_otu_table for %d samples, '
                'into output table %s' % (len(deblurred_list), output_fp))

    # the regexp for finding the number of reads of a sequence
    sizeregexp = re.compile('(?<=size=)\w+')
    seqdict = {}
    seqlist = []
    sampset = set()
    samplist = []
    # arbitrary size for the sparse results matrix so we won't run out of space
    obs = scipy.sparse.dok_matrix((1E9, len(deblurred_list)), dtype=np.int)

    # load the sequences from all samples into a sprase matrix
    for (cfilename, csampleid) in deblurred_list:
        # test if sample has already been processed
        if csampleid in sampset:
            warnings.warn('sample %s already in table!', UserWarning)
            logger.error('sample %s already in table!' % csampleid)
            continue
        sampset.add(csampleid)
        samplist.append(csampleid)
        csampidx = len(sampset)-1
        # read the fasta file and add to the matrix
        for chead, cseq in sequence_generator(cfilename):
            if cseq not in seqdict:
                seqdict[cseq] = len(seqlist)
                seqlist.append(cseq)
            cseqidx = seqdict[cseq]
            cfreq = float(sizeregexp.search(chead).group(0))
            try:
                obs[cseqidx, csampidx] = cfreq
            except IndexError:
                # exception means we ran out of space - add more OTUs
                shape = obs.shape
                obs.resize((shape[0]*2,  shape[1]))
                obs[cseqidx, csampidx] = cfreq

    logger.info('for final biom table loaded %d samples, %d unique sequences'
                % (len(samplist), len(seqlist)))

    # and now make the sparse matrix the real size
    obs.resize((len(seqlist), len(samplist)))

    # do the minimal reads per otu filtering
    if minreads > 0:
        readsperotu = obs.sum(axis=1)
        keep = np.where(readsperotu >= minreads)[0]
        logger.info('keeping %d (out of %d sequences) with >=%d reads' %
                    (len(keep), len(seqlist), minreads))
        obs = obs[keep, :]
        seqlist = list(np.array(seqlist)[keep])
        logger.debug('filtering completed')

    # convert the matrix to a biom table
    table = Table(obs, seqlist, samplist,
                  observation_metadata=None,
                  sample_metadata=None, table_id=None,
                  generated_by="deblur",
                  create_date=datetime.now().isoformat())
    logger.debug('converted to biom table')

    # remove samples with 0 reads
    filter_minreads_samples_from_table(table)

    # save the merged otu table
    write_biom_table(table, output_fp)
    logger.info('saved to biom file %s' % output_fp)

    # and save the fasta file
    if outputfasta_fp is not None:
        logger.debug('saving fasta file')
        with open(outputfasta_fp, 'w') as f:
            for cseq in seqlist:
                f.write('>%s\n%s\n' % (cseq, cseq))
        logger.info('saved sequence fasta file to %s' % outputfasta_fp)


def launch_workflow(seqs_fp, working_dir, mean_error, error_dist,
                    indel_prob, indel_max, trim_length, min_size, ref_fp,
                    ref_db_fp, negate, threads_per_sample=1,
                    sim_thresh=None, coverage_thresh=None):
    """Launch full deblur workflow for a single post split-libraries fasta file

    Parameters
    ----------
    seqs_fp: string
        a post split library fasta file for debluring
    working_dir: string
        working directory path
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
    threads_per_sample: integer, optional
        number of threads to use for SortMeRNA/mafft/vsearch
        (0 for max available)
    sim_thresh: float, optional
        the minimal similarity for a sequence to the database.
        if None, take the defaults (0.65 for negate=False,
        0.95 for negate=True)
    coverage_thresh: float, optional
        the minimal coverage for alignment of a sequence to the database.
        if None, take the defaults (0.3 for negate=False, 0.95 for negate=True)

    Return
    ------
    output_no_chimers_fp : string
        filepath to fasta file with no chimeras of None if error encountered
    """
    logger = logging.getLogger(__name__)
    logger.info('--------------------------------------------------------')
    logger.info('launch_workflow for file %s' % seqs_fp)

    # Step 1: Trim sequences to specified length
    output_trim_fp = join(working_dir, "%s.trim" % basename(seqs_fp))
    with open(output_trim_fp, 'w') as out_f:
        for label, seq in trim_seqs(
                input_seqs=sequence_generator(seqs_fp), trim_len=trim_length):
            out_f.write(">%s\n%s\n" % (label, seq))
    # Step 2: Dereplicate sequences
    output_derep_fp = join(working_dir,
                           "%s.derep" % basename(output_trim_fp))
    dereplicate_seqs(seqs_fp=output_trim_fp,
                     output_fp=output_derep_fp,
                     min_size=min_size, threads=threads_per_sample)
    # Step 3: Remove artifacts
    output_artif_fp = remove_artifacts_seqs(seqs_fp=output_derep_fp,
                                            ref_fp=ref_fp,
                                            working_dir=working_dir,
                                            ref_db_fp=ref_db_fp,
                                            negate=negate,
                                            threads=threads_per_sample,
                                            sim_thresh=sim_thresh)
    if not output_artif_fp:
        warnings.warn('Problem removing artifacts from file %s' %
                      seqs_fp, UserWarning)
        logger.warning('remove artifacts failed, aborting')
        return None
    # Step 4: Multiple sequence alignment
    output_msa_fp = join(working_dir,
                         "%s.msa" % basename(output_artif_fp))
    alignment = multiple_sequence_alignment(seqs_fp=output_artif_fp,
                                            threads=threads_per_sample)
    if not alignment:
        warnings.warn('Problem performing multiple sequence alignment '
                      'on file %s' % seqs_fp, UserWarning)
        logger.warning('msa failed. aborting')
        return None

    # Step 5: Launch deblur
    output_deblur_fp = join(working_dir,
                            "%s.deblur" % basename(output_msa_fp))
    with open(output_deblur_fp, 'w') as f:
        seqs = deblur(sequence_generator(output_msa_fp), mean_error,
                      error_dist, indel_prob, indel_max)
        if seqs is None:
            warnings.warn('multiple sequence alignment file %s contains '
                          'no sequences' % output_msa_fp, UserWarning)
            logger.warn('no sequences returned from deblur for file %s' %
                        output_msa_fp)
            return None
        for s in seqs:
            # remove '-' from aligned sequences
            s.sequence = s.sequence.replace('-', '')
            f.write(s.to_fasta())
    # Step 6: Chimera removal
    output_no_chimeras_fp = remove_chimeras_denovo_from_seqs(
        output_deblur_fp, working_dir, threads=threads_per_sample)
    logger.info('finished processing file')
    return output_no_chimeras_fp


def start_log(level=logging.DEBUG, filename=None):
    """start the logger for the run

    Parameters
    ----------
    level : int, optional
        logging.DEBUG, logging.INFO etc. for the log level (between 0-50).
    filename : str, optional
      name of the filename to save the log to or
      None (default) to use deblur.log.TIMESTAMP
    """
    if filename is None:
        tstr = time.ctime()
        tstr = tstr.replace(' ', '.')
        tstr = tstr.replace(':', '.')
        filename = 'deblur.log.%s' % tstr
    logging.basicConfig(filename=filename, level=level,
                        format='%(levelname)s(%(thread)d)'
                        '%(asctime)s:%(message)s')
    logger = logging.getLogger(__name__)
    logger.info('*************************')
    logger.info('deblurring started')


def _system_call(cmd, stdoutfilename=None):
    """Execute the command `cmd`
    Parameters
    ----------
    cmd : str
        The string containing the command to be run.
    stdoutfilename : str
        Name of the file to save stdout to or None
        (default) to not save to file
    stderrfilename : str
        Name of the file to save stderr to or None
        (default) to not save to file

    Returns
    -------
    tuple of (str, str, int)
        The standard output, standard error and exist status of the
        executed command

    Notes
    -----
    This function is ported and modified from QIIME
    (http://www.qiime.org), previously named
    qiime_system_call. QIIME is a GPL project, but we obtained permission from
    the authors of this function to port it to Qiita and keep it under BSD
    license.
    """
    logger = logging.getLogger(__name__)
    logger.debug('system call: %s' % cmd)
    if stdoutfilename:
        with open(stdoutfilename, 'w') as f:
            proc = subprocess.Popen(cmd, universal_newlines=True,
                                    shell=False, stdout=f,
                                    stderr=subprocess.PIPE)
    else:
        proc = subprocess.Popen(cmd, universal_newlines=True,
                                shell=False, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    # Communicate pulls all stdout/stderr from the PIPEs
    # This call blocks until the command is done
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return stdout, stderr, return_value
