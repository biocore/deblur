# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division
import multiprocessing as mp
import traceback
import sys
import subprocess
from functools import partial
import logging


def deblur_system_call(params, input_fp):
    """Build deblur command for subprocess.

    Parameters
    ----------
    params: list of str
        parameter settings to pass to deblur CLI
    input_fp : str
        name of the input fasta file to deblur

    Returns
    -------
    stdout: string
        process output directed to standard output
    stderr: string
        process output directed to standard error
    return_value: integer
        return code from process

    """
    logger = logging.getLogger(__name__)
    logger.debug('deblur system call params %s, input_fp %s' %
                 (params, input_fp))

    # construct command
    script_name = "deblur"
    script_subprogram = "workflow"
    command = [script_name,
               script_subprogram,
               '--seqs-fp', input_fp,
               '--is-worker-thread',
               '--keep-tmp-files']
    command.extend(params)

    logger.debug('running command %s' % command)
    return system_call(command)


def system_call(cmd):
    """Call cmd and return (stdout, stderr, return_value).

    Parameters
    ----------
    cmd: iterable
        A sequence of strings that are the tokens of the command

    Returns
    -------
    stdout: string
        process output directed to standard output
    stderr: string
        process output directed to standard error
    return_value: integer
        return code from process
    """
    logger = logging.getLogger(__name__)
    logger.info('system call cmd %s' % cmd)

    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    logger.info('system call finished return code %d' % return_value)

    return stdout, stderr, return_value


def run_functor(functor, *args, **kwargs):
    """
    Given a functor, run it and return its result. We can use this with
    multiprocessing.map and map it over a list of job functors to do them.

    Handles getting more than multiprocessing's pitiful exception output

    This function was derived from:
    http://stackoverflow.com/a/16618842/19741

    This code was adopted from the American Gut project:
    https://github.com/biocore/American-Gut/blob/master/americangut/parallel.py
    """
    try:
        # This is where you do your actual work
        return functor(*args, **kwargs)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


def parallel_deblur(inputs, params, ref_db_fp, jobs_to_start=1):
    """Dispatch execution over a pool of processors

    This code was adopted from the American Gut project:
    https://github.com/biocore/American-Gut/blob/master/americangut/parallel.py

    Parameters
    ----------
    inputs : iterable of str
        File paths to input per-sample sequence files
    params : list of str
        list of CLI parameters supplied to the deblur workflow
        (argv - first 2 are 'deblur','workflow' and are ignored)
    ref_db_fp : list of str
        the indexed sortmerna database (created in the main thread)
    jobs_to_start : int, optional
        The number of processors on the local system to use

    Returns
    -------
    all_result_paths : list
        list of expected output files
    """
    logger = logging.getLogger(__name__)
    logger.info('parallel deblur started for %d inputs' % len(inputs))

    # remove the irrelevant parameters
    remove_param_list = ['-O', '--jobs-to-start', '--seqs-fp', '--ref-db-fp']
    skipnext = False
    newparams = []
    for carg in params[2:]:
        if skipnext:
            skipnext = False
            continue
        if carg in remove_param_list:
            skipnext = True
            continue
        newparams.append(carg)

    # add the ref_db_fp (since it may be not present in the
    # original command parameters)
    new_ref_db_fp = ','.join(ref_db_fp)
    newparams.append('--ref-db-fp')
    newparams.append(new_ref_db_fp)

    logger.debug('ready for functor %s' % newparams)
    functor = partial(run_functor, deblur_system_call, newparams)
    logger.debug('ready for pool %d jobs' % jobs_to_start)
    pool = mp.Pool(processes=jobs_to_start)
    logger.debug('almost running...')
    for stdout, stderr, es in pool.map(functor, inputs):
        if es != 0:
            raise RuntimeError("stdout: %s\nstderr: %s\nexit: %d" % (stdout,
                                                                     stderr,
                                                                     es))
