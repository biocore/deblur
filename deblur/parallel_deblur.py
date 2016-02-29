# ----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division

from os.path import basename, join
from os import mkdir

from bfillings.sortmerna_v2 import build_database_sortmerna

from skbio.parse.sequences import parse_fasta

from qiime.parallel.util import ParallelWrapper, BufferedWriter
from qiime.parallel.poller import basic_process_run_results_f


class ParallelDeblur(ParallelWrapper):
    _script_name = "deblur"
    _job_prefix = 'PDeblur'
    _input_splitter = ParallelWrapper._input_existing_filepaths

    def _get_poller_command(self,
                            expected_files_filepath,
                            merge_map_filepath,
                            deletion_list_filepath,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate command to initiate a poller to monitior/process
           completed runs
        """

        result = '%s poller.py -f %s -p %s -m %s -d %s -t %d %s' % \
            (command_prefix,
             expected_files_filepath,
             self._process_run_results_f,
             merge_map_filepath,
             deletion_list_filepath,
             self._seconds_to_sleep,
             command_suffix)

        return result, []
