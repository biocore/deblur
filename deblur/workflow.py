# -----------------------------------------------------------------------------
# Copyright (c) 2015, The Deblur Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from deblur.deblurring import deblur


def demultiplex_seqs(seqs_fp, mapping_fp):
	"""Step 1: demultiplex and quality filter a FASTQ file"""
	pass


def split_seqs_on_sample_ids(seqs_fp):
	"""Step 2: split demultiplexed FASTA file based on
	   Sample IDs"""
	pass


def trim_seqs(seqs_fp, trim_len):
	"""Step 3: trim FASTA sequences to specified length"""
	pass


def dereplicate_seqs(seqs_fp):
	"""Step 4a: dereplicate FASTA sequences using VSEARCH"""
	pass


def remove_singletons_seqs(seqs_fp):
	"""Step 4b: remove singletons from dereplicated FASTA
	   file using VSEARCH"""
	pass


def remove_artifacts_seqs(seqs_fp):
	"""Step 5: remove artifacts from FASTA file"""
	pass


def multiple_sequence_alignment(seqs_fp):
	"""Step 6: perform multiple sequence alignment on FASTA
	   file using MAFFT"""
	pass


def remove_chimeras_denovo_from_seqs(seqs_fp):
	"""Step 7: remove chimeras de novo using UCHIME"""
	pass


def generate_biom_table(seqs_fp):
	"""Step 8: generate BIOM table and representative
	   FASTA set"""
	pass


def assign_taxonomy(seqs_fp, biom_fp):
	"""Step 9: assign taxonomy to sequences"""
	pass


def launch_workflow(seqs_fp, mapping_fp):
	"""Launch full deblur workflow"""
	pass

