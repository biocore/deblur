# deblur changelog


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
First release

### Features

### Backward-incompatible changes [stable]

### Performance enhancements

### Bug fixes

### Miscellaneous
