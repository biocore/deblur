Deblur
======

[![Build Status](https://travis-ci.org/biocore/deblur.png?branch=master)](https://travis-ci.org/biocore/deblur)
[![Coverage Status](https://coveralls.io/repos/biocore/deblur/badge.png?branch=master)](https://coveralls.io/r/biocore/deblur)

Deblur is a greedy deconvolution algorithm based on known read error profiles.

Code Development Note
=====================

Some of the code in the package deblur has been derived from [QIIME](http://qiime.org).
The contributors to these specific QIIME modules have granted permission
for this porting to take place and put under the BSD license.

Dependencies
============

QIIME 1.9.1
burrito-fillings development version

Install
=======

Once both QIIME 1.9.1 and the development version of burrito-fillings is installed,
make sure your $PYTHONPATH is updated to have burrito-fillings appear before
qiime.

Example usage
=============

The input to deblur workflow is a demultiplexed FASTA file. The deblur algorithm is
designed to work on individual samples.

```
deblur workflow s1.fna db.fna output.biom
```

or

```
parallel_deblur.py workflow all_samples.fna db.fna output.biom
```