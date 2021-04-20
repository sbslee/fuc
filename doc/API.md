# API

## Table of contents

* [BedFrame](#BedFrame) 
* [FastqFrame](#FastqFrame) 
* [VcfFrame](#VcfFrame) 
* [common](#common) 

## BedFrame <a name="BedFrame"></a>

```
Python Library Documentation: module fuc.api.BedFrame in fuc.api

NAME
    fuc.api.BedFrame

DESCRIPTION
    The BedFrame module is designed for working with BED files. For example,
    it can be used to find the intersection between multiple BED files.

CLASSES
    builtins.object
        BedFrame
    
    class BedFrame(builtins.object)
     |  BedFrame(meta, gr)
     |  
     |  Class for storing BED data.
     |  
     |  This class is essentially a wrapper for the `pyranges` package
     |  (https://github.com/biocore-ntnu/pyranges).
     |  
     |  BED lines have three required fields and nine additional optional fields:
     |       1. chrom (required) - The name of the chromosome.
     |       2. chromStart (required) - The starting position of the feature.
     |       3. chromEnd (required) - The ending position of the feature.
     |       4. name (optional) - Defines the name of the BED line.
     |       5. score (optional) - A score between 0 and 1000 for color density.
     |       6. strand (optional) - Either "." (=no strand) or "+" or "-".
     |       7. thickStart (optional) - The starting position for thick drawing.
     |       8. thickEnd (optional) - The ending position for thick drawing.
     |       9. itemRgb (optional) - An RGB value (e.g. 255,0,0).
     |      10. blockCount (optional) - The number of blocks (exons).
     |      11. blockSizes (optional) - A comma-separated list of the block sizes.
     |      12. blockStarts (optional) - A comma-separated list of block starts.
     |  
     |  For more information about the BED format, visit the UCSC Genome Browser
     |  FAQ (https://genome.ucsc.edu/FAQ/FAQformat.html).
     |  
     |  Methods defined here:
     |  
     |  __init__(self, meta, gr)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  intersect(self, other)
     |      Find intersection between the BedFrames.
     |  
     |  to_file(self, file_path)
     |      Write the BedFrame to a BED file.
     |  
     |  to_string(self)
     |      Render the BedFrame to a console-friendly tabular output.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a BedFrame from a BED file.
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)

DATA
    HEADERS = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'T...

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/BedFrame.py

```

## FastqFrame <a name="FastqFrame"></a>

```
Python Library Documentation: module fuc.api.FastqFrame in fuc.api

NAME
    fuc.api.FastqFrame

DESCRIPTION
    The FastqFrame module is designed for working with FASTQ files (both zipped
    and unzipped).

CLASSES
    builtins.object
        FastqFrame
        FastqRecord
    
    class FastqFrame(builtins.object)
     |  FastqFrame(data: List[fuc.api.FastqFrame.FastqRecord]) -> None
     |  
     |  FastqFrame(data: List[fuc.api.FastqFrame.FastqRecord])
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __init__(self, data: List[fuc.api.FastqFrame.FastqRecord]) -> None
     |  
     |  __repr__(self)
     |  
     |  readlen(self)
     |      Return a dictionary of read lengths and their counts.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a FastqFrame from a file.
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  shape
     |      Return the size of the FastqFrame.
     |  
     |  vdata
     |      Return a view (copy) of the data.
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
     |  __annotations__ = {'data': typing.List[fuc.api.FastqFrame.FastqRecord]...
     |  
     |  __dataclass_fields__ = {'data': Field(name='data',type=typing.List[fuc...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...
     |  
     |  __hash__ = None
    
    class FastqRecord(builtins.object)
     |  FastqRecord(id: str, seq: str, ext: str, qual: str) -> None
     |  
     |  FastqRecord(id: str, seq: str, ext: str, qual: str)
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __hash__(self)
     |  
     |  __init__(self, id: str, seq: str, ext: str, qual: str) -> None
     |  
     |  __repr__(self)
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
     |  __annotations__ = {'ext': <class 'str'>, 'id': <class 'str'>, 'qual': ...
     |  
     |  __dataclass_fields__ = {'ext': Field(name='ext',type=<class 'str'>,def...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...

DATA
    List = typing.List
        A generic version of list.

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/FastqFrame.py

```

## VcfFrame <a name="VcfFrame"></a>

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
     |      n1 : str or int
     |          Name of index of the test sample.
     |      n2 : str or int
     |          Name of index of the truth sample.
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
     |  filter_multiallelic(self)
     |      Filter out rows that have multiple alternative alleles.
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
     |  parse_snpeff(self, idx, sep=' | ')
     |      Parse SnpEff annotations.
     |      
     |      SnpEff provides the following functional annotations:
     |           1. Allele
     |           2. Annotation
     |           3. Annotation_Impact
     |           4. Gene_Name
     |           5. Gene_ID
     |           6. Feature_Type
     |           7. Feature_ID
     |           8. Transcript_BioType
     |           9. Rank
     |          10. HGVS.c
     |          11. HGVS.p
     |          12. cDNA.pos / cDNA.length
     |          13. CDS.pos / CDS.length
     |          14. AA.pos / AA.length
     |          15. Distance
     |          16. ERRORS / WARNINGS
     |          17. INFO
     |      
     |      Parameters
     |      ----------
     |      i : list
     |          List of annotation indicies.
     |      sep : str, default: ' | '
     |          Separator for joining requested annotations.
     |      
     |      Returns
     |      -------
     |      s : pandas.Series
     |          Parsed annotations.
     |  
     |  reset_samples(self, names)
     |      Reset the sample list.
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
     |  update(self, other, headers=None, missing=True)
     |      Copy data from the other VcfFrame.
     |      
     |      This method will copy and paste data from the other VcfFrame for
     |      overlapping records. By default, the following VCF headers are
     |      used: ID, QUAL, FILTER, and, INFO.
     |      
     |      Parameters
     |      ----------
     |      other : VcfFrame
     |          Other VcfFrame.
     |      headers : list, optional
     |          List of VCF headers to exclude.
     |      missing : bool, default: True
     |          If True, only fields with the missing value ('.') will be updated.
     |      
     |      Returns
     |      -------
     |      vf : VcfFrame
     |          Updated VcfFrame.
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

## common <a name="common"></a>

```
Python Library Documentation: module fuc.api.common in fuc.api

NAME
    fuc.api.common

DESCRIPTION
    The common module is used by other `fuc` modules such as VcfFrame and
    BedFrame. It also provides many useful methods.

FUNCTIONS
    fuc_dir()
        Return the path to the fuc directory.
    
    get_most_similar(a, l)
        Return the most similar string in a list.
    
    get_script_name(script)
        Return the script name.
    
    get_similarity(a, b)
        Return a value from 0 to 1 representing how similar two strings are.
    
    is_numeric(s)
        Return True if the string is numeric.
    
    is_similar(a, b, threshold=0.9)
        Return True if the similarity is equal to or greater than threshold.
    
    parse_condition(condition)
        Parse one condition in the SQLite WHERE clause.
    
    parse_where(where)
        Parse the SQLite WHERE clause.

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/common.py

```

