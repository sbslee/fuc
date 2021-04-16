# README

The `fuc` package is my attempt to wrap some of the most frequently used commands in the field of bioinformatics with pure Python 3 code, only using standard libraries. That means no installation for external Python packages, not even the popular ones like `numpy` and `pandas`. This also includes many famous bioinformatics tools such as `samtools` and `bedtools`. The motivation for not relying on external packages or programs is quite simple: I just got tired of managing different Python environments for doing simple things.

You can use `fuc` for both command line interface (CLI) and application programming interface (API). Click [here](doc/CLI.md) to see the CLI documentation and [here](doc/API.md) to see the API documentation.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

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
    dfsumcol     [TABLE] summarize a column in a text file
    fuccompf     [FUC] compare two files
    fucexist     [FUC] check whether files/dirs exist
    qfcount      [FASTQ] count sequence reads in a FASTQ file
    qfreadlen    [FASTQ] compute read lengths for a FASTQ file
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

For getting help on API (e.g. `VcfFrame`):

```
from fuc import VcfFrame
help(VcfFrame)
```

To give:

```
Python Library Documentation: module fuc.api.VcfFrame in fuc.api

NAME
    fuc.api.VcfFrame

CLASSES
    builtins.object
        VcfFrame
        VcfRecord
    
    class VcfFrame(builtins.object)
     |  VcfFrame(meta: List[str], head: List[str], data: List[fuc.api.VcfFrame.VcfRecord]) -> None
     |  
     |  VcfFrame(meta: List[str], head: List[str], data: List[fuc.api.VcfFrame.VcfRecord])
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __init__(self, meta: List[str], head: List[str], data: List[fuc.api.VcfFrame.VcfRecord]) -> None
     |  
     |  __repr__(self)
     |  
     |  add_dp(self)
     |      Compute and add the DP subfield of the FORMAT field.
     |  
     |  compare(self, n1, n2)
     |      Compare two samples within the VcfFrame.
     |      
     |      Parameters
     |      ----------
     |      n1 : string or int
     |          Test sample or its index in the header row.
     |      n2 : string or int
     |          Truth sample or its index in the header row.
     |      
     |      Returns
     |      -------
     |      result : tuple
     |          Comparison result (tp, fp, fn, and tn).
     |  
     |  describe(self)
     |      Generate descriptive statistics.
     |  
     |  filter_af(self, threshold=0.1)
     |      Filter based on the AF subfield of the FORMAT field.
     |  
     |  filter_bed(self, bed)
     |      Filter VcfRecords in the VcfFrame using BED data.
     |      
     |      Parameters
     |      ----------
     |      bed : BedFrame or string
     |          BedFrame or path to a BED file.
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Filtered VcfFrame.
     |  
     |  filter_dp(self, threshold=200)
     |      Filter based on the DP subfield of the FORMAT field.
     |  
     |  filter_empty(self)
     |      Filter out VcfRecords that are empty.
     |  
     |  index(self, name)
     |      Return the sample index.
     |  
     |  merge(self, other, format_subfields=None)
     |      Merge with the other VcfFrame.
     |      
     |      Parameters
     |      ----------
     |      other : VcfFrame
     |          Other VcfFrame.
     |      format_subfields : list, optional
     |          Additional FORMAT subfields (e.g. DP and AD) to be retained
     |          other than GT, which is included as default.
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Stripped VcfFrame.
     |  
     |  multiallelic_sites(self)
     |      Return the indicies of multiallelic sites.
     |  
     |  reset_samples(self, samples)
     |      Reset the sample list.
     |  
     |  strip(self, format_subfields=None)
     |      Remove unnecessary data from the VcfFrame.
     |      
     |      Parameters
     |      ----------
     |      format_subfields : list, optional
     |          Additional FORMAT subfields (e.g. DP and AD) to be retained
     |          other than GT, which is included as default.
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Stripped VcfFrame.
     |  
     |  to_file(self, file_path)
     |      Write the VcfFrame to a file.
     |  
     |  to_string(self)
     |      Render the VcfFrame to a console-friendly tabular output.
     |  
     |  update(self, other, query_fields, missing_only=False)
     |      Copy data from another VcfFrame.
     |      
     |      This method will copy requested data from another VcfFrame for
     |      overlapping records. You can only request data from the following
     |      VCF headers: ID, QUAL, FILTER, INFO, and FORMAT. Any other
     |      requested VCF headers will be ignored.
     |      
     |      Parameters
     |      ----------
     |      other : VcfFrame
     |          Target VcfFrame.
     |      names : list
     |          List of VCF headers.
     |      missing_only : boolean, optional
     |          If True, only fields with the missing value will be updated.
     |      
     |      Returns
     |      -------
     |      vcf_result : VcfFrame
     |          Updated VcfFrame.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a VcfFrame from a file.
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  samples
     |      Return a list of the sample IDs.
     |  
     |  shape
     |      Return a tuple representing the dimensionality of the VcfFrame.
     |  
     |  vdata
     |      Return a view (copy) of the data.
     |  
     |  vhead
     |      Return a view (copy) of the headers.
     |  
     |  vmeta
     |      Return a view (copy) of the metadata.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  HEADERS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INF...
     |  
     |  __annotations__ = {'data': typing.List[fuc.api.VcfFrame.VcfRecord], 'h...
     |  
     |  __dataclass_fields__ = {'data': Field(name='data',type=typing.List[fuc...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...
     |  
     |  __hash__ = None
    
    class VcfRecord(builtins.object)
     |  VcfRecord(chrom: str, pos: int, id: str, ref: str, alt: List[str], qual: str, filter: List[str], info: List[str], format: List[str], gt: List[str]) -> None
     |  
     |  Class for storing the information of single VCF record.
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |      Test whether two VcfRecords are equal.
     |  
     |  __init__(self, chrom: str, pos: int, id: str, ref: str, alt: List[str], qual: str, filter: List[str], info: List[str], format: List[str], gt: List[str]) -> None
     |  
     |  __repr__(self)
     |  
     |  to_list(self)
     |      Convert the VcfRecord to a list of strings.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_list(l) from builtins.type
     |      Create a VcfRecord from a list of strings.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __annotations__ = {'alt': typing.List[str], 'chrom': <class 'str'>, 'f...
     |  
     |  __dataclass_fields__ = {'alt': Field(name='alt',type=typing.List[str],...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=False,o...
     |  
     |  __hash__ = None

FUNCTIONS
    has_var(x)
        Return if the GT field has a variant (e.g. 0/1).

DATA
    List = typing.List

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/VcfFrame.py

```
