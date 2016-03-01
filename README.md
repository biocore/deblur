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

burrito-fillings development version

Example usage
=============

The input to deblur workflow is a demultiplexed FASTA file. The deblur algorithm is
designed to work on individual samples.

```
deblur workflow all_samples.fna db.fna output.biom
```

or

```
parallel_deblur workflow all_samples.fna db.fna output.biom
```