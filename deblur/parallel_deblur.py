# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division
import os
from os.path import join
import multiprocessing as mp
import traceback
import sys
import subprocess
from functools import partial


def deblur_system_call(params, fp):
    """Build deblur command for subprocess.

    Parameters
    ----------
    params: dict
        parameter settings to pass to deblur
    fp: string
        input filepath

    Returns
    -------
    stdout: string
        process output directed to standard output
    stderr: string
        process output directed to standard error
    return_value: integer
        return code from process

    """
    input_fp = fp
    script_name = "deblur"
    script_subprogram = "workflow_parallel"
    # construct command
    command = [script_name,
               script_subprogram,
               '--seqs-fp', input_fp]
    # add reference databases to command
    ref_fp_l = [db[1] for db in enumerate(params.get('ref-fp', []))]
    for ref_fp in ref_fp_l:
        command.append('--ref-fp')
        command.append(ref_fp)
    # add reference database indexes to command
    ref_db_fp_l = [db[1] for db in enumerate(params.get('ref-db-fp', []))]
    for ref_dp_fp in ref_db_fp_l:
        command.append('--ref-db-fp')
        command.append(ref_dp_fp)
    cmd_list = []
    # add the remainder of the parameters
    for key, value in params.iteritems():
        if (key != 'ref-fp' and key != 'ref-db-fp' and value is not None):
            cmd_list.append("--%s" % key)
            cmd_list.append("%s" % value)
    command.extend(cmd_list)

    print "[deblur_system_call] command = ", command

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
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode

    print(stdout)
    print(stderr)
    print(return_value)

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


def parallel_deblur(inputs, params, jobs_to_start=1):
    """Dispatch execution over a pool of processors

    This code was adopted from the American Gut project:
    https://github.com/biocore/American-Gut/blob/master/americangut/parallel.py

    Parameters
    ----------
    inputs : iterable of str
        File paths to input per-sample sequence files
    params : dict
        A dict of parameters invariant to per-sample calls
    jobs_to_start : int, optional
        The number of processors on the local system to use

    Returns
    -------
    all_result_paths : list
        list of expected output files
    """
    # Create working directory
    output_dir = params['output-dir']
    working_dir = join(output_dir, "deblur_working_dir")
    os.mkdir(working_dir)
    args = []
    all_result_paths = []
    for in_ in inputs:
        filename = os.path.split(in_)[-1]
        output = os.path.join(working_dir, "%s.biom" % filename)
        all_result_paths.append(output)
        args.append(in_)
    functor = partial(run_functor, deblur_system_call, params)
    pool = mp.Pool(processes=jobs_to_start)
    for stdout, stderr, es in pool.map(functor, args):
        if es != 0:
            raise RuntimeError("stdout: %s\nstderr: %s\nexit: %d" % (stdout,
                                                                     stderr,
                                                                     es))
    return all_result_paths
