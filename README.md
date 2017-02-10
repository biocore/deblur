Deblur
======

[![Build Status](https://travis-ci.org/biocore/deblur.png?branch=master)](https://travis-ci.org/biocore/deblur)
[![Coverage Status](https://coveralls.io/repos/github/biocore/deblur/badge.svg?branch=master)](https://coveralls.io/github/biocore/deblur?branch=master)

Deblur is a greedy deconvolution algorithm for amplicon sequencing based on Illumina Miseq/Hiseq error profiles.

Install
=======
- Deblur requires Python 3.5. If Python 3.5 is not installed, you can create a [conda](http://conda.pydata.org/docs/install/quick.html) environment for Deblur using:
```
conda create -n deblurenv python=3.5 numpy
```

and activate it using:
```
source activate deblurenv
```

(note you will need to activate this environment every time you want to use Deblur)

Install Deblur dependencies and Deblur itself:
```
conda install -c bioconda -c biocore VSEARCH=2.0.3 MAFFT=7.221 biom-format SortMeRNA==2.0 deblur
```

N.B. Some dependencies are version restricted at the moment but for different reasons. SortMeRNA 2.1 has a different output format which Deblur is not compatible with yet. VSEARCH 2.0.3 and MAFFT 7.221 yield slightly different results in testing than their respective latest versions (at the time of release of Deblur 1.0.0). A review of their changelogs did not reveal any remarkable notes (e.g., bugs) about the reasons for the differences. In testing, the differences affected <0.1% of the sOTUs. As a precaution, we are advising the use of these specific versions for consistency with the manuscript. 

Example usage
=============

We recommend using Deblur via the QIIME2 plugin [q2-deblur](https://github.com/wasade/q2-deblur). Examples of its use can be found within the plugin itself. However, Deblur itself does not depend on QIIME2.

If you are running Deblur directly, we recommend focusing on the `workflow` subcommand. Detailed help can be obtained with:

```
deblur workflow --help
```

As a simple example, let's specify an input FASTA file, an output path and a sequence trim length of 150. This command will trim all sequences in `all_samples.fna` to 150nt in length; any read that is shorter will be omitted. This execution mode assumes that `all_samples.fna` is demultiplexed such that the sequence IDs are compatible with QIIME 1.9.1. On completion, a new directory `output` will be created with multiple output files (see the Input and Output Files section for more detail). 

```
deblur workflow --seqs-fp all_samples.fna --output-dir output -t 150
```

If starting from a barcode and read file, you can first use the qiime [split_libraries_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) command (we recommend using -q 19 to remove low quality reads):

```
split_libraries_fastq.py -i XXX_R1_001.fastq -m map.txt -o split -b XXXX_I1_001.fastq -q 19
```

The resulting `split/seqs.fna` file can be used as the input to the deblur workflow.

Input and Output Files
======================

The input to Deblur workflow is a directory of FASTA or FASTQ files (1 per sample) or a single demultiplexed FASTA or FASTQ file. These files can be gzip'd. The output directory will contain three BIOM tables in which the observation IDs are the Deblurred sequences. The outputs are contingent on the reference databases used and a more focused discussion on them is in the subsequent README section titled "Positive and Negative Filtering." The output files are as follows:

- reference-hit.biom : contains only Deblurred reads matching the positive filtering database. By default, a reference composed of 16S sequences is used, and this resulting table will contain only those reads which recruit at a coarse level to it will be retained. Reads are also filtered against the negative reference, which by default will remove any read which appears to be PhiX or adapter.

- reference-hit.seqs.fa : a fasta file containing all the sequences in reference-hit.biom

- reference-non-hit.biom : contains only Deblurred reads that did not align to the positive filtering database. Negative filtering is also appied to this table, so by default, PhiX and adapter are removed.

- reference-non-hit.seqs.fa : a fasta file containing all the sequences in reference-non-hit.biom

- all.biom : contains all Deblurred reads. This file represents the union of the "reference-hit.biom" and "reference-non-hit.biom" tables.

- all.seqs.fa : a fasta file containing all the sequences in all.biom

Important options
=================

Deblur cannot associate sequences with different lengths. As such, trimming reads is a required first step in the Deblur pipeline. The sequence trim length is specified by the ```-t NNN``` flag, where NNN denotes the length all sequences will be trimmed to. All reads shorter than this length will be discarded. If the input data are known to have a common length, it is possible to disable trimming by specifying a trim value of `-1`.

Deblur can operate in parallel. The number of threads can be specified by the ```-O NNN``` flag (default it 1). Running more threads than available cores is not advised. 

Positive and Negative Filtering
===============================

Deblur uses two types of filtering on the sequences:

- Negative mode - removes [known artifact sequences](deblur/support_files/artifacts.fa) (i.e. sequences aligning to PhiX or Adapter with >=95% identity and coverage).

- Positive mode - keeps only sequences similar to a reference database (by default [known 16S sequences](deblur/support_files/88_otus.fasta)). SortMeRNA is used, and any sequence with an e-value <= 10 is retained. Deblur also outputs a BIOM table without this positive filtering step (named all.biom).

The FASTA files for both of these filtering steps can be supplied via the --neg-ref-fp and --pos-ref-fp options. By default, the negative database is composed of PhiX and adapter sequence and the positive database of known 16S sequences.

Deblur uses negative mode filtering to remove known artifact (i.e. PhiX and Adapter sequences) prior to denoising. The output of Deblur contains three files: all.biom, which includes all sOTUs, reference-hit.biom, which contains the output of positive filtering of the sOTUs (default only sOTUs similar to 16S sequences), and reference-non-hit.biom, which contains only sOTUs failing the positive filtering (default only non-16S sOTUs).

Minimal Reads Filtering
=======================

Deblur runs on each sample independently. However, sometimes there is also additional information based on the total number of times an sOTU is observed in all samples (e.g. an sOTU which is observed only in one sample at low read count is more likely to be pcr/read error as opposed to an sOTU present in many samples). The --min-reads option allows to use this information by removing sOTUs with a total read count (across all samples) lower than the given threshold. The default value is set to 10, and should be useful for most cases. However, if such filtering is not wanted (e.g. if all samples in an experiment are expected not to contain the same bacteria, so no additional information is gained by combining the information from multiple samples), --min-reads can be set to 0 to skip this final filtering step.

Troubleshooting
===============
- Mac users: if you get the following error:

RuntimeError: Python is not installed as a framework. The Mac OS X backend will not be able to function correctly if Python is not installed as a framework. See the Python documentation for more information on installing Python as a framework on Mac OS X. Please either reinstall Python as a framework, or try one of the other backends. If you are Working with Matplotlib in a virtual enviroment see 'Working with Matplotlib in Virtual environments' in the Matplotlib FAQ

You can solve it by the following commands:
```
cd ~/.matplotlib
echo "backend: TkAgg" >> ~/.matplotlib/matplotlibrc
```

- "Too many open files" : This error indicates deblur is trying to split a single fasta/q file into per-sample files, and the OS does not allow so many open simultaneous open files. Current solution is to use the qiime1.9 command split_sequence_file_on_sample_ids.py or the equivalent qiime2 command to split the single fasta/q file into a directory of per sample fasta/q files and then run deblur with this directory as the input to deblur (--seqs-fp).

Code Development Note
=====================

Some of the code in the package deblur has been derived from [QIIME](http://qiime.org).
The contributors to these specific QIIME modules have granted permission
for this porting to take place and put under the BSD license.
