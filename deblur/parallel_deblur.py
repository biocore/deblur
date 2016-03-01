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


def deblur_system_call(filetype, ref_fp_str, ref_db_fp_str, input_fp,
                       output_fp):
    """
    """
    script_name = "deblur workflow_parallel"
    command = ' '.join([script_name,
                        '--seqs-fp %s' % input_fp,
                        '--output-fp %s' % output_fp,
                        '--file-type %s' % filetype,
                        ref_fp_str,
                        ref_db_fp_str])
    #return system_call(command)
    return command


def system_call(cmd):
    """Call cmd and return (stdout, stderr, return_value).
    Parameters
    ----------
    cmd: str
        Can be either a string containing the command to be run, or a sequence
        of strings that are the tokens of the command.

    Notes
    -----
    This function is ported from QIIME (http://www.qiime.org), previously
    named qiime_system_call. QIIME is a GPL project, but we obtained permission
    from the authors of this function to port it to deblur (and keep it under
    deblur's BSD license).
    """
    proc = subprocess.Popen(cmd,
                            universal_newlines=True,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    if return_value != 0:
        raise ValueError("Failed to execute: %s\nstdout: %s\nstderr: %s" %
                         (cmd, stdout, stderr))

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


def parallel_deblur(inputs, outputs, params, jobs_to_start=None):
    """Dispatch execution over a pool of processors

    Parameters
    ----------
    inputs : iterable of str
        File paths to input per-sample sequence files
    outputs : iterable of str
        File paths to outputs
    params : dict
        A dict of parameters invariant to per-sample calls
    jobs_to_start : int, optional
        The number of processors on the local system to use.

    This code was adopted from the American Gut project:
    https://github.com/biocore/American-Gut/blob/master/americangut/parallel.py
    """
    ref_fp_str = ""
    for db in enumerate(params.get('ref_fp', [])):
        ref_fp_str = "%s --ref-fp %s" % (ref_fp_str, db[1])

    ref_db_fp_str = ""
    for db in enumerate(params.get('ref_db_fp', [])):
        ref_db_fp_str = "%s --ref-db-fp %s" % (ref_db_fp_str, db[1])

    filetype = params['file_type']

    functor = partial(run_functor, deblur_system_call, filetype, ref_fp_str,
                      ref_db_fp_str)

    args = list(zip(inputs, outputs, params))
    pool = mp.Pool(processes=jobs_to_start)
    for result in pool.map(functor, args):
        print(result)
