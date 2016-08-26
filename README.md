Deblur
======

[![Build Status](https://travis-ci.org/biocore/deblur.png?branch=master)](https://travis-ci.org/biocore/deblur)
[![Coverage Status](https://coveralls.io/repos/github/biocore/deblur/badge.svg?branch=master)](https://coveralls.io/github/biocore/deblur?branch=master)

Deblur is a greedy deconvolution algorithm based on Illumina Miseq/Hiseq error profiles.

Install
=======
- Deblur requires Python 3.5. If Python 3.5 is not installed, you can create a conda environment for deblur using:

```
conda create -n debenv python=3
```

and activate it using:

```
source activate debenv
```

- Install Deblur dependencies:
```
conda install -c biocore VSEARCH MAFFT SortMeRNA numpy
```

- Install Deblur:
```
pip install git+git://github.com/biocore/deblur
```

Example usage
=============

The input to deblur workflow is a directory of fasta files (1 per sample) or a demultiplexed FASTA or FASTQ file. The output is a biom table with sequences as the OTU ids (final.biom in the output directory).

The simple use case just specifies the input fasta file (or directory) and output directory name:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output
```

If starting from a barcode and read file, you can first use the qiime split_libraries_fastq.py command (we recommend using -q 19 to remove low quality reads):

```
split_libraries_fastq.py -i XXX_R1_001.fastq -m map.txt -o split -b XXXX_I1_001.fastq -q 19
```

and use the split/seqs.fna as the input to the deblur workflow.

Positive and Negative Filtering
===============================
By default, deblur uses positive filtering, keeping only 16S sequences (based on homology to greengenes 88% representative set). for example:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output
```

Negative filtering can be selected using the '-n' flag. This causes deblurring to keep all sequences except known artifact sequences (i.e. PhiX and Adapter sequences), so other non-16S sequences are retained. For example:

```
deblur workflow --seqs-fp all_samples.fna --output-dir output -n
```

Code Development Note
=====================

Some of the code in the package deblur has been derived from [QIIME](http://qiime.org).
The contributors to these specific QIIME modules have granted permission
for this porting to take place and put under the BSD license.

