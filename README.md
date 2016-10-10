Deblur
======

[![Build Status](https://travis-ci.org/biocore/deblur.png?branch=master)](https://travis-ci.org/biocore/deblur)
[![Coverage Status](https://coveralls.io/repos/github/biocore/deblur/badge.svg?branch=master)](https://coveralls.io/github/biocore/deblur?branch=master)

Deblur is a greedy deconvolution algorithm based on Illumina Miseq/Hiseq error profiles.

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
conda install -c bioconda VSEARCH MAFFT SortMeRNA==2.0 biom-format
```

- Install Deblur:
```
pip install deblur
```

Example usage
=============

The input to deblur workflow is a directory of fasta files (1 per sample) or a demultiplexed FASTA or FASTQ file. The output is a biom table with sequences as the OTU ids (final.biom in the output directory).

The simple use case just specifies the input fasta file (or directory) and output directory name:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output
```

If starting from a barcode and read file, you can first use the qiime [split_libraries_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) command (we recommend using -q 19 to remove low quality reads):

```
split_libraries_fastq.py -i XXX_R1_001.fastq -m map.txt -o split -b XXXX_I1_001.fastq -q 19
```

and use the split/seqs.fna as the input to the deblur workflow.

Important options
=================
- The sequence read length can be specified by the ```-t NNN``` flag, where NNN denotes the length all sequences will be trimmed to (default=150). Note that all reads shorter than this length will be discarded.

- In order to run in parallel, the number of threads can be specified by the ```-O NNN``` flag (default it 1). Note that running more threads than available cores will not speed up performance.

- To get a full list of options, use:
```
deblur workflow --help
```

Positive and Negative Filtering
===============================
By default, deblur uses positive filtering, keeping only 16S sequences (based on homology to [Greengenes](http://greengenes.secondgenome.com/) [88% representative set](deblur/support_files/88_otus.fasta)). For example:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output
```

Negative filtering can be selected using the '-n' flag. This causes deblurring to keep all sequences except for [known artifact sequences](deblur/support_files/artifacts.fa) (i.e. PhiX and Adapter sequences), so other non-16S sequences are retained. For example:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output -n
```

Code Development Note
=====================

Some of the code in the package deblur has been derived from [QIIME](http://qiime.org).
The contributors to these specific QIIME modules have granted permission
for this porting to take place and put under the BSD license.
