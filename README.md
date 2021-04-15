# README

The `fuc` package is my attempt to wrap some of the most frequently used commands in the field of bioinformatics with pure Python 3 code, only using standard libraries. That means no installation for external Python packages, not even the popular ones like `numpy` and `pandas`. This also includes many famous bioinformatics tools such as `samtools` and `bedtools`. The motivation for not relying on external packages or programs is quite simple: I just got tired of managing different Python environments for doing simple things.

You can use `fuc` for both command line interface (CLI) and application programming interface (API). Click [here](doc/CLI.md) to see the CLI documentation and [here](doc/API.md) to see the API documentation.

To install `fuc`, enter the following in your terminal:

```
$ git clone https://github.com/sbslee/fuc
$ cd fuc
$ pip install .
```

For getting help:

```
$ fuc -h
usage: fuc [-h] [-v] COMMAND ...

positional arguments:
  COMMAND        name of the command
    bfintxn      [BED] find intersection of two or more BED files
    bfsum        [BED] summarize a BED file
    dfmerge      [TABLE] merge two text files
    dfsum        [TABLE] summarize a text file
    dfsumcol     [TABLE] summarize a column in a text file
    fuccheckf    [FUC] check whether files exist
    fuccompf     [FUC] compare any two files
    qfcount      [FASTQ] count sequence reads in a FASTQ file
    qfreadlen    [FASTQ] compute read lengths for a FASTQ file
    vfmerge      [VCF] merge two or more VCF files

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show the version number and exit
```

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

