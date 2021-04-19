# README

## Introduction

The main goal of the `fuc` package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use `fuc` for both command line interface (CLI) and application programming interface (API). Click [here](doc/CLI.md) to see the CLI documentation and [here](doc/API.md) to see the API documentation.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

## Requirements

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

- **BedFrame** : 
- **BedFrameOLD** : The BedFrame module is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **FastqFrame** : The FastqFrame module is designed for working with FASTQ files (both zipped and unzipped).
- **VcfFrame** : The VcfFrame module is designed for working with VCF files (both zipped and unzipped).
- **VcfFrameOLD** : The VcfFrame module is designed for working with VCF files (both zipped and unzipped).
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
     |  VcfFrame(meta, data)
     |  
     |  Class for storing VCF data.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, meta, data)
     |      Initialize self.  See help(type(self)) for accurate signature.
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
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a VcfFrame from a VCF file.
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
