# README

## Introduction

The main goal of the `fuc` package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use `fuc` for both command line interface (CLI) and application programming interface (API). Click [here](doc/CLI.md) to see the CLI documentation and [here](doc/API.md) to see the API documentation.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

## Required Packages

The following packages are required to run `fuc`:

```
numpy
pandas
pyranges
```

## Getting Started

To install `fuc`, enter the following in your terminal:

```
$ git clone https://github.com/sbslee/fuc
$ cd fuc
$ pip install .
```

For getting help on CLI:

```
$ fuc -h
usage: fuc [-h] [-v] COMMAND ...

positional arguments:
  COMMAND        name of the command
    bfintxn      [BED] find intersection of two or more BED files
    bfsum        [BED] summarize a BED file
    dfmerge      [TABLE] merge two text files
    dfsum        [TABLE] summarize a text file
    fuccompf     [FUC] compare two files
    fucexist     [FUC] check whether files/dirs exist
    qfcount      [FASTQ] count sequence reads in a FASTQ file
    qfsum        [FASTQ] summarize a FASTQ file
    vfmerge      [VCF] merge two or more VCF files

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show the version number and exit
```

For getting help on a specific command (e.g. `qfcount`):

```
$ fuc qfcount -h
usage: fuc qfcount [-h] fastq_file

This command will count sequence reads in a FASTQ file (both zipped and
unzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

Below is the list of modules available in API:

- **BedFrame** : The BedFrame module is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **FastqFrame** : The FastqFrame module is designed for working with FASTQ files (both zipped and unzipped).
- **VcfFrame** : The VcfFrame module is designed for working with VCF files (both zipped and unzipped).
- **common** : The common module is used by other fuc modules such as VcfFrame and BedFrame. It also provides many useful methods.

For getting help on a specific module (e.g. `VcfFrame`):

```
from fuc import VcfFrame
help(VcfFrame)
```

To give:

```
Python Library Documentation: module fuc.api.VcfFrame in fuc.api

NAME
    fuc.api.VcfFrame

DESCRIPTION
    The VcfFrame module is designed for working with VCF files (both zipped
    and unzipped).

CLASSES
    builtins.object
        VcfFrame
    
    class VcfFrame(builtins.object)
     |  VcfFrame(meta, df)
     |  
     |  Class for storing VCF data.
     |  
     |  This class strictly sticks to the standard Variant Call Format
     |  specification (https://samtools.github.io/hts-specs/VCFv4.3.pdf).
     |  
     |  VCF lines have nine required fields for storing variant data and
     |  variable-length fields for storing sample genotype data. In all cases,
     |  missing values are specified with a dot ('.'). The required fields are:
     |      1. CHROM - An identifier from the reference genome.
     |      2. POS - The 1-based reference position.
     |      3. ID - Semicolon-separated list of unique identifiers.
     |      4. REF - Reference base(s).
     |      5. ALT - Comma-separated list of alternate non-reference alleles.
     |      6. QUAL - Phred-scaled quality score for the assertion made in ALT.
     |      7. FILTER - PASS or a semicolon-separated list of filters that fail.
     |      8. INFO - Semicolon-separated series of additional information fields.
     |      9. FORMAT - Colon-separated series of genotype fields.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, meta, df)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  add_dp(self)
     |      Compute and add the DP subfield of the FORMAT field.
     |  
     |  compare(self, n1, n2)
     |      Compare two samples within the VcfFrame.
     |      
     |      Parameters
     |      ----------
     |      n1 : str
     |          Test sample.
     |      n2 : str
     |          Truth sample.
     |      
     |      Returns
     |      -------
     |      result : tuple
     |          Comparison result (tp, fp, fn, tn).
     |  
     |  filter_af(self, threshold=0.1)
     |      Filter rows based on the AF subfield of the FORMAT field.
     |  
     |  filter_bed(self, bed)
     |      Filter rows based on BED data.
     |      
     |      Parameters
     |      ----------
     |      bed : BedFrame or str
     |          BedFrame or path to a BED file.
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Filtered VcfFrame.
     |  
     |  filter_dp(self, threshold=200)
     |      Filter rows based on the DP subfield of the FORMAT field.
     |  
     |  filter_empty(self)
     |      Filter out rows that have no genotype calls.
     |  
     |  merge(self, other, how='inner', format='GT')
     |      Merge with the other VcfFrame.
     |      
     |      This method essentially wraps the `pandas.DataFrame.merge` method.
     |      
     |      Parameters
     |      ----------
     |      other : VcfFrame
     |          Other VcfFrame.
     |      how : str, default: 'inner'
     |          Type of merge to be performed. ['left', 'right', 'outer',
     |          'inner', 'cross']
     |      format : str, default: 'GT'
     |          FORMAT subfields to be retained (e.g. 'GT:AD:DP').
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Merged VcfFrame.
     |  
     |  strip(self, format='GT')
     |      Remove unnecessary data from the VcfFrame.
     |      
     |      Parameters
     |      ----------
     |      format : str, default: 'GT'
     |          FORMAT subfields to be retained (e.g. 'GT:AD:DP').
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Stripped VcfFrame.
     |  
     |  to_file(self, file_path)
     |      Write the VcfFrame to a VCF file.
     |  
     |  to_string(self)
     |      Render the VcfFrame to a console-friendly tabular output.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a VcfFrame from a VCF file.
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  samples
     |      Return a list of the sample IDs.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/VcfFrame.py

```
