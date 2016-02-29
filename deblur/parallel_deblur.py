# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division

from os.path import basename, join, split, splitext
from os import mkdir

from bfillings.sortmerna_v2 import build_database_sortmerna
from skbio.parse.sequences import parse_fasta

from qiime.parallel.util import ParallelWrapper


class ParallelDeblur(ParallelWrapper):
    _script_name = "deblur workflow_parallel"
    _job_prefix = 'PDeblur'
    _input_splitter = ParallelWrapper._input_existing_filepaths


    def _get_job_commands(self,
                          input_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Split multiple sample FASTA files to multiple jobs
        """
        commands = []
        result_filepaths = []
        # construct ref-fp string
        ref_fp_str = ""
        for db in enumerate(params['ref_fp']):
            ref_fp_str = "%s --ref-fp %s" % (ref_fp_str, db[1])
        # construct ref-dp-fp string
        ref_db_fp_str = ""
        for db in enumerate(params['ref_db_fp']):
            ref_db_fp_str = "%s --ref-db-fp %s" % (ref_db_fp_str, db[1])
        for i, input_fp in enumerate(input_fps):
            working_dir_t = join(working_dir, "%d" % i)
            mkdir(working_dir_t)
            input_path, input_fn = split(input_fp)
            input_basename, input_ext = splitext(input_fn)
            output_fns = ['%s.biom' % (input_basename)]
            rename_command, curr_result_fps = self._get_rename_command(
                output_fns, working_dir_t, output_dir)
            result_filepaths += curr_result_fps
            command = '%s %s --seqs-fp %s --output-fp %s --file-type %s %s %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 input_fp,
                 join(working_dir_t, output_fns[0]),
                 params['file_type'],
                 ref_fp_str,
                 ref_db_fp_str,
                 rename_command,
                 command_suffix)
            commands.append(command)
        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)
        return commands, result_filepaths
