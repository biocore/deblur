# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division
import os
import multiprocessing as mp
import traceback
import sys
import subprocess
from functools import partial


def deblur_system_call(filetype, ref_fp_str, ref_db_fp_str, fps):
    input_fp, output_fp = fps
    script_name = "deblur workflow_parallel"
    command = [script_name,
               '--seqs-fp %s' % input_fp,
               '--output-fp %s' % output_fp,
               '--file-type %s' % filetype,
               ref_fp_str,
               ref_db_fp_str]
    return system_call(command)


def system_call(cmd):
    """Call cmd and return (stdout, stderr, return_value).

    Parameters
    ----------
    cmd: iterable
        A sequence of strings that are the tokens of the command.
    """
    proc = subprocess.Popen(cmd,
                            shell=True,  # unfortunately needed for ease
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = proc.communicate()
    return_value = proc.returncode

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


def parallel_deblur(inputs, output_dir, params, jobs_to_start=None):
    """Dispatch execution over a pool of processors

    Parameters
    ----------
    inputs : iterable of str
        File paths to input per-sample sequence files
    output_dir : filepath
        Output directory
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

    args = []
    all_result_paths = []
    for in_ in inputs:
        filename = os.path.split(in_)[-1]
        output = os.path.join(output_dir, filename)
        all_result_paths.append(output)
        args.append((in_, output))

    filetype = params['file_type']
    functor = partial(run_functor, deblur_system_call, filetype, ref_fp_str,
                      ref_db_fp_str)

    pool = mp.Pool(processes=jobs_to_start)
    for stdout, stderr, es in pool.map(functor, args):
        if es != 0:
            raise RuntimeError("stdout: %s\nstderr: %s\nexit: %d" % (stdout,
                                                                     stderr,
                                                                     es))

    return all_result_paths
