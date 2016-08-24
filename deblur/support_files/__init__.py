#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import os

__all__ = ['pos_db', 'neg_db']

pos_db = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '88_otus.fasta')
neg_db = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'artifacts.fa')
