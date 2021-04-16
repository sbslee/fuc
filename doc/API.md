# API

## Table of contents

* [BedFrame](#BedFrame) 
* [DataFrame](#DataFrame) 
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
        BedRecord
    
    class BedFrame(builtins.object)
     |  BedFrame(meta: List[str], data: List[fuc.api.BedFrame.BedRecord]) -> None
     |  
     |  Class for storing the information of single BED file.
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |  
     |  __init__(self, meta: List[str], data: List[fuc.api.BedFrame.BedRecord]) -> None
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
     |  __annotations__ = {'data': typing.List[fuc.api.BedFrame.BedRecord], 'm...
     |  
     |  __dataclass_fields__ = {'data': Field(name='data',type=typing.List[fuc...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=True,or...
     |  
     |  __hash__ = None
    
    class BedRecord(builtins.object)
     |  BedRecord(chrom: str, start: int, end: int, name: Union[str, NoneType] = None, score: Union[int, NoneType] = None, strand: Union[str, NoneType] = None, tstart: Union[int, NoneType] = None, tend: Union[int, NoneType] = None, itemrgb: Union[List[int], NoneType] = None, bcount: Union[int, NoneType] = None, bsizes: Union[List[int], NoneType] = None, bstarts: Union[List[int], NoneType] = None) -> None
     |  
     |  Class for storing the information of single BED record.
     |  
     |  Methods defined here:
     |  
     |  __eq__(self, other)
     |      Test whether two BedRecords are equal.
     |  
     |  __init__(self, chrom: str, start: int, end: int, name: Union[str, NoneType] = None, score: Union[int, NoneType] = None, strand: Union[str, NoneType] = None, tstart: Union[int, NoneType] = None, tend: Union[int, NoneType] = None, itemrgb: Union[List[int], NoneType] = None, bcount: Union[int, NoneType] = None, bsizes: Union[List[int], NoneType] = None, bstarts: Union[List[int], NoneType] = None) -> None
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
     |  __annotations__ = {'bcount': typing.Union[int, NoneType], 'bsizes': ty...
     |  
     |  __dataclass_fields__ = {'bcount': Field(name='bcount',type=typing.Unio...
     |  
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=False,o...
     |  
     |  __hash__ = None
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
    Optional = typing.Optional

FILE
    /Users/sbslee/Desktop/fuc/fuc/api/BedFrame.py

```

## DataFrame <a name="DataFrame"></a>

```
Python Library Documentation: module fuc.api.DataFrame in fuc.api

NAME
    fuc.api.DataFrame

DESCRIPTION
    The DataFrame module is designed for working with table-like text files.
    It provides many useful methods for manipulating tables such as merging and
    filtering.

CLASSES
    builtins.object
        DataFrame
    
    class DataFrame(builtins.object)
     |  Methods defined here:
     |  
     |  __init__(self)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |  
     |  filter_columns(self, headers, exclude=False)
     |      Return a new table after filtering the columns.
     |  
     |  filter_rows(self, key, values, exclude=False)
     |      Return a new table after filtering the rows.
     |  
     |  get_col(self, i)
     |      Return the column as a list.
     |  
     |  get_data(self)
     |      Return a copy of the data which can be modified safely.
     |  
     |  get_head(self)
     |      Return a copy of the headers which can be modified safely.
     |  
     |  get_index(self, name)
     |      Return the column index.
     |  
     |  get_unique(self, i)
     |      Return unique values of the column.
     |  
     |  merge(self, other, on, missing='.')
     |      Return a merged DataFrame.
     |  
     |  summarize_col(self, i)
     |      Return summary of the column.
     |  
     |  write(self, file_path, delimiter='\t')
     |      Write the DataFrame to a file.
     |  
     |  ----------------------------------------------------------------------
     |  Class methods defined here:
     |  
     |  read(file_path, delimiter='\t', header=True, skiprows=None) from builtins.type
     |      Create a DataFrame from a file.
     |  
     |  ----------------------------------------------------------------------
     |  Readonly properties defined here:
     |  
     |  shape
     |      Return a tuple representing the dimensionality of the DataFrame.
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
    /Users/sbslee/Desktop/fuc/fuc/api/DataFrame.py

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
     |      Test whether two FastqRecords are equal.
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
     |  __dataclass_params__ = _DataclassParams(init=True,repr=True,eq=False,o...
     |  
     |  __hash__ = None

DATA
    List = typing.List

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

