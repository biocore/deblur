# deblur changelog

## Version 1.1.0-dev

### Performance enhancements

* Updated to Python 3.8.
* Update MAFFT to >=7.39.

## Version 1.1.0

Official version 1.1.0 Released on 12 September 2018.

### Bug fixes
* Fixed problem causing deblur to ignore the indel error uper bound and use the mismatch error upper bound instead. This caused deblur to use the mismatch depenent default error rate instead of the default indel error rate (which is constant for up to 3 indels). See [issue #178](https://github.com/biocore/deblur/issues/178) for more details.
* Fixed problem running deblur on new python versions (due to python treating the sparse matrix size (1E9) as float rather than int).

## Version 1.0.4

Official version 1.0.4 Released on 24 January 2018.

### Bug fixes
* During checking for chimeras via vsearch, lower quality nucleotides were indicated by lower case bases. This is not informative for Deblur, but led to issues on downstream analyses, e.g. importing the resulting table into QIIME2 or a disconnect of fragment sequences between the table and SEPP generated phylogenies. It also posed the risk of distributing read count to two different versions of the same sequence, when one version was all uppercase, the other with one or more lower case characters. With this bugfix-release: Deblur count tables only contain features with all uppercase nucleotide sequences.

## Version 1.0.3

Official version 1.0.3 Released on 19 October 2017.

### Features

* Added `--left-trim-length` to allow for trimming nucleotides on the 5' end of each sequence. Please see [issue #154](https://github.com/biocore/deblur/issues/154) for more information.

### Performance enhancements
* Updated MAFFT version requirement from 7.221 to 7.310.
* Updated vsearch version requirement from 2.0.3 to >= 2.0.3

## Version 1.0.2

### Bug fixes
E-value thresholding removed in favor of a bitscore threshold which is scaled based on the length of the query sequence. This is to make each query sequence independent of each other on assessment against the positive filtering database. Please see [PR #146](https://github.com/biocore/deblur/pull/146) for further details.

## Version 1.0.1

Official version 1.0.1 Released on 28 February 2017.

### Features

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes
* Filtering thresholds for the reference database were not in use. What this means is that, in some cases, highly artifactual sequence could filter through. Specifically, no coverage and similarity thresholds were used with SortMeRNA for the postive reference filter. These are now set to 60 and 65 respectively. In addition, the e-value threshold used is now the SortMeRNA default of 1. Please see [issue #139](https://github.com/biocore/deblur/issues/139) for further discussion in addition for the rational behind the parameter values picked.

## Version 1.0

Official version 1.0. Released on 9 February 2017.

### Features
* Updated documentation to improve clarity on examples, file inputs and file outputs.
* Delete tmp files and split directory in deblur output directory

### Backward-incompatible changes [stable]
* Removed the --skip-trimming flag and use --trim-length=-1 instead
* The trim length parameter is now required.
* Negative mode is always run. What that means is that an output is always generated with reads matching the negative database removed (an unfiltered output is still generated). By default, this means that any read which appears to be PhiX or Adapter will be filtered out.
* Change output file names to all.biom, all.seqs.fa, reference-non-hit.biom, reference-non-hit.seqs.fa, reference-hit.biom, reference-hit.seqs.fa. The term "reference" in this case refers to sequences which recruited (hit) to the positive database or failed to recruit (non-hit). Recruitment is performed using SortMeRNA and a sequence is retained if it aligns with an e-value of 10 or less (i.e., a coarse filter).
* Changed --min-reads default value to 10 (was 0) in accordance with the manuscript.
* Renamed --threads shortcut to -a (was -t) in remove_artifacts CLI

### Performance enhancements

### Bug fixes

## Version 0.1.7-dev
### Features
* Adding command line option in deblur workflow to skip trimming [#124](https://github.com/biocore/deblur/pull/124)

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes


## Version 0.1.4-dev
### Features
* support for .gz FASTA and FASTQ files [#113](https://github.com/biocore/deblur/pull/113)

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes


## Version 0.1.3-dev
### Features

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes
* Positive mode using multiple threads failed [#111](https://github.com/biocore/deblur/pull/111)
* Default negative database not set correctly when not supplied [#98](https://github.com/biocore/deblur/pull/98)

### Miscellaneous


## Version 0.1.1-dev
### Features
* Create 16S only and non-16S only biom tables in negative mode [#87](https://github.com/biocore/deblur/pull/87)
* Added different flags for the positive and negative databases [#87](https://github.com/biocore/deblur/pull/87)

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes
* Fix total pipeline failure when a sample contains 1 read after dereplcation+singleton removal. Now issues a warning and continues [#75](https://github.com/biocore/deblur/issues/75)

### Miscellaneous

## Version 0.1.0-dev
* First release

### Features

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes

### Miscellaneous
