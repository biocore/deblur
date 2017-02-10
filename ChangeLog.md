# deblur changelog

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
