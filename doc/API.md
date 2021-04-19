# API

## Table of contents

* [BedFrame](#BedFrame) 
* [BedFrameOLD](#BedFrameOLD) 
* [FastqFrame](#FastqFrame) 
* [VcfFrame](#VcfFrame) 
* [VcfFrameOLD](#VcfFrameOLD) 
* [common](#common) 

## BedFrame <a name="BedFrame"></a>

```
Python Library Documentation: module fuc.api.BedFrame in fuc.api

NAME
    fuc.api.BedFrame

CLASSES
    builtins.object
        BedFrame
    
    class BedFrame(builtins.object)
     |  BedFrame(meta, data)
     |  
     |  Class for storing BED data.
     |  
     |  Methods defined here:
     |  
     |  __init__(self, meta, data)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  intersect(self, other)
     |      Find intersection between the BedFrames.
     |  
     |  to_file(self, file_path)
     |      Write the BedFrame to a BED file.
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
    BF_HEADERS = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'str...
    PR_HEADERS = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand',...

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/BedFrame.py

```

## BedFrameOLD <a name="BedFrameOLD"></a>

```
Python Library Documentation: module fuc.api.BedFrameOLD in fuc.api

NAME
    fuc.api.BedFrameOLD

DESCRIPTION
    The BedFrame module is designed for working with BED files. For example,
    it can be used to find the intersection between multiple BED files.

CLASSES
    builtins.object
        BedFrame
        BedRecord
    
    class BedFrame(builtins.object)
     |  BedFrame(meta: List[str], data: List[fuc.api.BedFrameOLD.BedRecord]) -> None
     |  
     |  Class for storing the information of single BED file.
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __init__(self, meta: List[str], data: List[fuc.api.BedFrameOLD.BedRecord]) -> None
     |  
     |  __repr__(self)
     |  
     |  intersect(self, other)
     |      Find the intersection between the two BedFrames.
     |      
     |      Metadata and optional BED fields in the other BedFrame
     |      will be ignored.
     |      
     |      Parameters
     |      ----------
     |      other : BedFrame
     |          Other BedFrame.
     |      
     |      Returns
     |      -------
     |      bf: BedFrame
     |          New BedFrame.
     |  
     |  to_file(self, file_path)
     |      Write the BedFrame to a file.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_file(file_path) from builtins.type
     |      Create a BedFrame from a file.
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  shape
     |      Return the size of the BedFrame.
     |  
     |  vdata
     |      Return a view (copy) of the data.
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
     |  __annotations__ = {'data': typing.List[fuc.api.BedFrameOLD.BedRecord],...
     |  
     |  __dataclass_fields__ = {'data': Field(name='data',type=typing.List[fuc...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...
     |  
     |  __hash__ = None
    
    class BedRecord(builtins.object)
     |  BedRecord(chrom: str, start: int, end: int, name: Optional[str] = None, score: Optional[int] = None, strand: Optional[str] = None, tstart: Optional[int] = None, tend: Optional[int] = None, itemrgb: Optional[List[int]] = None, bcount: Optional[int] = None, bsizes: Optional[List[int]] = None, bstarts: Optional[List[int]] = None) -> None
     |  
     |  Class for storing the information of single BED record.
     |  
     |  This class strictly sticks to the BED format described in the UCSC
     |  Genome Browser (https://genome.ucsc.edu/FAQ/FAQformat.html).
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
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |      Test whether two BedRecords are equal.
     |  
     |  __hash__(self)
     |  
     |  __init__(self, chrom: str, start: int, end: int, name: Optional[str] = None, score: Optional[int] = None, strand: Optional[str] = None, tstart: Optional[int] = None, tend: Optional[int] = None, itemrgb: Optional[List[int]] = None, bcount: Optional[int] = None, bsizes: Optional[List[int]] = None, bstarts: Optional[List[int]] = None) -> None
     |  
     |  __repr__(self)
     |  
     |  to_list(self)
     |      Convert the BedRecord to a list of strings.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  from_list(l) from builtins.type
     |      Create a BedRecord from a list of strings.
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
     |  __annotations__ = {'bcount': typing.Optional[int], 'bsizes': typing.Op...
     |  
     |  __dataclass_fields__ = {'bcount': Field(name='bcount',type=typing.Opti...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...
     |  
     |  bcount = None
     |  
     |  bsizes = None
     |  
     |  bstarts = None
     |  
     |  itemrgb = None
     |  
     |  name = None
     |  
     |  score = None
     |  
     |  strand = None
     |  
     |  tend = None
     |  
     |  tstart = None

DATA
    List = typing.List
        A generic version of list.
    
    Optional = typing.Optional
        Optional type.
        
        Optional[X] is equivalent to Union[X, None].

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/BedFrameOLD.py

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

## VcfFrameOLD <a name="VcfFrameOLD"></a>

```
Python Library Documentation: module fuc.api.VcfFrameOLD in fuc.api

NAME
    fuc.api.VcfFrameOLD

DESCRIPTION
    The VcfFrame module is designed for working with VCF files (both zipped
    and unzipped).

CLASSES
    builtins.object
        VcfFrame
        VcfRecord
    
    class VcfFrame(builtins.object)
     |  VcfFrame(meta: List[str], head: List[str], data: List[fuc.api.VcfFrameOLD.VcfRecord]) -> None
     |  
     |  VcfFrame(meta: List[str], head: List[str], data: List[fuc.api.VcfFrameOLD.VcfRecord])
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __init__(self, meta: List[str], head: List[str], data: List[fuc.api.VcfFrameOLD.VcfRecord]) -> None
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
     |  __annotations__ = {'data': typing.List[fuc.api.VcfFrameOLD.VcfRecord],...
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
     |  __eq__(self, other)
     |  
     |  __hash__(self)
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
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...

DATA
    List = typing.List
        A generic version of list.

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/VcfFrameOLD.py

```

## common <a name="common"></a>

```
Python Library Documentation: module fuc.api.common in fuc.api

NAME
    fuc.api.common

DESCRIPTION
    The common module is used by other fuc modules such as VcfFrame and
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

