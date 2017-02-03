Deblur
======

[![Build Status](https://travis-ci.org/biocore/deblur.png?branch=master)](https://travis-ci.org/biocore/deblur)
[![Coverage Status](https://coveralls.io/repos/github/biocore/deblur/badge.svg?branch=master)](https://coveralls.io/github/biocore/deblur?branch=master)

Deblur is a greedy deconvolution algorithm for amplicon sequencing based on Illumina Miseq/Hiseq error profiles.

Install
=======
- Deblur requires Python 3.5. If Python 3.5 is not installed, you can create a [conda](http://conda.pydata.org/docs/install/quick.html) environment for deblur using:
```
conda create -n deblurenv python=3 numpy
```

and activate it using:
```
source activate deblurenv
```

(note you will need to activate this environment every time you want to use deblur)

At the moment, the install is a two stage process as we do not currently have deblur staged in a conda channel.

- install deblur dependencies
```
conda install -c bioconda VSEARCH MAFFT biom-format SortMeRNA==2.0
```

- Install Deblur:
```
pip install deblur
```

Example usage
=============

The input to deblur workflow is a directory of fasta files (1 per sample) or a demultiplexed FASTA or FASTQ file. The output are three biom tables (in the output directory) with sequences as the OTU ids:

- all.biom : contains all deblurred reads.

- reference-hit.biom : contains only deblurred reads matching the positive filtering database (default is only 16S reads)

- reference-non-hit.biom : contains only deblurred reads not matching the positive filtering database (default is all non-16S reads)

The simple use case just specifies the input fasta file (or directory), output directory name and the sequence trim length:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output -t 150
```

If starting from a barcode and read file, you can first use the qiime [split_libraries_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) command (we recommend using -q 19 to remove low quality reads):

```
split_libraries_fastq.py -i XXX_R1_001.fastq -m map.txt -o split -b XXXX_I1_001.fastq -q 19
```

and use the split/seqs.fna as the input to the deblur workflow.

Important options
=================
- Since deblur cannot associate together sequences with different lengths, trimming is automaticallly performed as the first step in the deblur pipeline. The sequence trim length is specified by the ```-t NNN``` flag, where NNN denotes the length all sequences will be trimmed to. Note that all reads shorter than this length will be discarded.

- In order to run in parallel, the number of threads can be specified by the ```-O NNN``` flag (default it 1). Note that running more threads than available cores will not speed up performance.

- To get a full list of options, use:
```
deblur workflow --help
```

Positive and Negative Filtering
===============================
Deblur uses two types of filtering on the sequences:

- negative mode - removes [known artifact sequences](deblur/support_files/artifacts.fa) (i.e. sequences aligning to PhiX or Adapter with >=95% identity and coverage).

- positive mode - keeps only sequences similar (sortmerna e-value<=10) to a reference database (by default [known 16S sequences](deblur/support_files/88_otus.fasta)). Note that deblur also outputs a biom table without this positive filtering step (named all.biom).

The fasta files for both of these filtering steps can be supplied via the --neg-ref-fp and --pos-ref-fp options, and by default are supplied for 16S sequences.

Deblur uses negative mode filtering to remove known artifact (i.e. PhiX and Adapter sequences) prior to denoising. The output of deblur contains three files: all.biom, which includes all sOTUs, reference-hit.biom, which contains the output of positive filtering of the sOTUs (default only sOTUs similar to 16S sequences), and reference-non-hit.biom, which contains only sOTUs failing the positive filtering (default only non-16S sOTUs).

Minimal Reads Filtering
=======================
Deblur runs on each sample independently. However, sometimes there is also additional information based on the total number of times an sOTU is observed in all samples (e.g. an sOTU which is observed only in one sample at low read count is more likely to be pcr/read error as opposed to an sOTU present in many samples). The --min-reads option allows to use this information by removing sOTUs with a total read count (across all samples) lower than the given threshold. The default value is set to 10, and should be useful for most cases. However, if such filtering is not wanted (e.g. if all samples in an experiment are expected not to contain the same bacteria, so no additional information is gained by combining the information from multiple samples), --min-reads can be set to 0 to skip this final filtering step.

Code Development Note
=====================

Some of the code in the package deblur has been derived from [QIIME](http://qiime.org).
The contributors to these specific QIIME modules have granted permission
for this porting to take place and put under the BSD license.
