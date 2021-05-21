"""
The pyvcf submodule is designed for working with Variant Call Format (VCF)
files (both zipped and unzipped). It implements ``pyvcf.VcfFrame`` which
stores VCF data as ``pandas.DataFrame`` to allow fast computation and easy
manipulation. The submodule strictly adheres to the standard `VCF
specification <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_.

A VCF file contains metadata lines (prefixed with '##'), a header line
(prefixed with '#'), and genotype lines that start with a chromosome
identifier (e.g. 'chr1'). See the VCF specification above for an example
VCF file.

Genotype lines have nine required fields for storing variant information
and variable-length fields for storing sample genotype data. For some
fields, missing values are tolerated and specified with a dot ('.').
The nine required fields are:

+-----+--------+------------------------------------+------------------------+
| No. | Name   | Description                        | Examples               |
+=====+========+====================================+========================+
| 1   | CHROM  | Chromosome or contig identifier    | 'chr2', '2', 'chrM'    |
+-----+--------+------------------------------------+------------------------+
| 2   | POS    | 1-based reference position         | 10041, 23042           |
+-----+--------+------------------------------------+------------------------+
| 3   | ID     | ';'-separated variant identifiers  | '.', 'rs35', 'rs9;rs53'|
+-----+--------+------------------------------------+------------------------+
| 4   | REF    | Reference allele                   | 'A', 'GT'              |
+-----+--------+------------------------------------+------------------------+
| 5   | ALT    | ','-separated alternate alleles    | 'T', 'ACT', 'C,T'      |
+-----+--------+------------------------------------+------------------------+
| 6   | QUAL   | Phred-scaled quality score for ALT | '.', 67, 12            |
+-----+--------+------------------------------------+------------------------+
| 7   | FILTER | ';'-separated filters that failed  | '.', 'PASS', 'q10;s50' |
+-----+--------+------------------------------------+------------------------+
| 8   | INFO   | ';'-separated information fields   | '.', 'DP=14;AF=0.5;DB' |
+-----+--------+------------------------------------+------------------------+
| 9   | FORMAT | ':'-separated genotype fields      | 'GT', 'GT:AD:DP'       |
+-----+--------+------------------------------------+------------------------+

There are several common, reserved genotype keywords that are standards
across the community. Currently, the module is aware of the
following:

* AD - Total read depth for each allele (R, Integer).
* AF - Allele fraction of the event in the tumor (1, Float)
* DP - Read depth (1, Integer)
"""

import pandas as pd
import numpy as np
import gzip
from copy import deepcopy
from Bio import bgzf
from . import pybed

CONTIGS = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M',
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM'
]

# -- Private methods ---------------------------------------------------------

# -- Public methods ----------------------------------------------------------

def gt_haspolyp(g):
    """Return True if the sample genotype has a polyploid call.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if the genotype has a polyploid call.

    Examples
    --------
    >>> pyvcf.gt_haspolyp('0/1')
    False
    >>> pyvcf.gt_haspolyp('1')
    False
    >>> pyvcf.gt_haspolyp('0/1/1')
    True
    >>> pyvcf.gt_haspolyp('1|0|1')
    True
    >>> pyvcf.gt_haspolyp('0/./1/1')
    True
    """
    gt = g.split(':')[0]
    if '/' in gt:
        return gt.count('/') > 1
    else:
        return gt.count('|') > 1

def gt_hasvar(g):
    """Return True if the sample genotype has a variant call.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if the genotype has a variant call.

    Examples
    --------
    Below are some simple examples:

    >>> pyvcf.gt_hasvar('0/1:35:4')
    True
    >>> pyvcf.gt_hasvar('0/0:61:2')
    False
    >>> pyvcf.gt_hasvar('1|2:21:6:23,27')
    True
    >>> pyvcf.gt_hasvar('0|0:48:1:51,51')
    False
    """
    if g.split(':')[0].replace('/', '').replace(
        '|', '').replace('.', '').replace('0', ''):
        return True
    else:
        return False

def gt_missing(g):
    """Return True if the sample genotype has a missing call.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if the genotype has a missing call.

    Examples
    --------
    Below are some simple examples:

    >>> pyvcf.gt_missing('0|0:48:1:51,51')
    False
    >>> pyvcf.gt_missing('./.:.:.')
    True
    >>> pyvcf.gt_missing('.:.')
    True
    >>> pyvcf.gt_missing('.')
    True
    """
    return '.' in g.split(':')[0]

def gt_unphase(g):
    """Unphase the sample genotype.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    str
        Unphased genotype.

    Examples
    --------
    Below are some simple examples:

    >>> pyvcf.gt_unphase('1|2:21:6:23,27')
    '1/2:21:6:23,27'
    >>> pyvcf.gt_unphase('2|1:2:0:18,2')
    '1/2:2:0:18,2'
    >>> pyvcf.gt_unphase('0/1:35:4')
    '0/1:35:4'
    """
    l = g.split(':')
    gt = l[0]
    if '|' not in gt:
        return g
    l[0] = '/'.join([str(b) for b in sorted([int(a) for a in gt.split('|')])])
    return ':'.join(l)

def merge(vfs, how='inner', format='GT', sort=True, collapse=False):
    """Merge VcfFrame objects.

    Parameters
    ----------
    vfs : list
        List of VcfFrames to be merged.
    how : str, default: 'inner'
        Type of merge as defined in pandas.DataFrame.merge.
    format : str, default: 'GT'
        FORMAT subfields to be retained (e.g. 'GT:AD:DP').
    sort : bool, default: True
        If True, sort the VcfFrame before returning.
    collapse : bool, default: False
        If True, collapse duplicate records.

    Returns
    -------
    VcfFrame
        Merged VcfFrame.

    See Also
    --------
    VcfFrame.merge
        Merge self with another VcfFrame.

    Examples
    --------
    Assume we have the following data:

    >>> data1 = {
    ...     'CHROM': ['chr1', 'chr1'],
    ...     'POS': [100, 101],
    ...     'ID': ['.', '.'],
    ...     'REF': ['G', 'T'],
    ...     'ALT': ['A', 'C'],
    ...     'QUAL': ['.', '.'],
    ...     'FILTER': ['.', '.'],
    ...     'INFO': ['.', '.'],
    ...     'FORMAT': ['GT:DP', 'GT:DP'],
    ...     'Steven': ['0/0:32', '0/1:29'],
    ...     'Sara': ['0/1:24', '1/1:30'],
    ... }
    >>> data2 = {
    ...     'CHROM': ['chr1', 'chr1', 'chr2'],
    ...     'POS': [100, 101, 200],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['G', 'T', 'A'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['.', '.', '.'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP'],
    ...     'Dona': ['./.:.', '0/0:24', '0/0:26'],
    ...     'Michel': ['0/1:24', '0/1:31', '0/1:26'],
    ... }
    >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
    >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
    >>> vf1.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara
    0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:32  0/1:24
    1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  1/1:30
    >>> vf2.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Dona  Michel
    0  chr1  100  .   G   A    .      .    .  GT:DP   ./.:.  0/1:24
    1  chr1  101  .   T   C    .      .    .  GT:DP  0/0:24  0/1:31
    2  chr2  200  .   A   T    .      .    .  GT:DP  0/0:26  0/1:26

    We can merge the two VcfFrames with ``how='inner'`` (default):

    >>> pyvcf.merge([vf1, vf2]).df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara Dona Michel
    0  chr1  100  .   G   A    .      .    .     GT    0/0  0/1  ./.    0/1
    1  chr1  101  .   T   C    .      .    .     GT    0/1  1/1  0/0    0/1

    We can also merge with ``how='outer'``:

    >>> pyvcf.merge([vf1, vf2], how='outer').df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara Dona Michel
    0  chr1  100  .   G   A    .      .    .     GT    0/0  0/1  ./.    0/1
    1  chr1  101  .   T   C    .      .    .     GT    0/1  1/1  0/0    0/1
    2  chr2  200  .   A   T    .      .    .     GT    ./.  ./.  0/0    0/1

    Since both VcfFrames have the DP subfield, we can use ``format='GT:DP'``:

    >>> pyvcf.merge([vf1, vf2], how='outer', format='GT:DP').df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara    Dona  Michel
    0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:32  0/1:24   ./.:.  0/1:24
    1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  1/1:30  0/0:24  0/1:31
    2  chr2  200  .   A   T    .      .    .  GT:DP   ./.:.   ./.:.  0/0:26  0/1:26
    """
    merged_vf = vfs[0]
    for vf in vfs[1:]:
        merged_vf = merged_vf.merge(vf, how=how, format=format, sort=sort,
            collapse=collapse)
    return merged_vf

def read_file(fn, compression=False):
    """Read a VCF file into VcfFrame.

    If the file name ends with '.gz', the method will automatically use the
    BGZF decompression when reading the file.

    Parameters
    ----------
    fn : str
        Path to a zipped or unzipped VCF file.
    compression : bool, default: False
        If True, use the BGZF decompression.

    Returns
    -------
    VcfFrame
        VcfFrame.

    Examples
    --------
    >>> pyvcf.read_file('example.vcf')
    >>> pyvcf.read_file('example.vcf.gz')
    >>> pyvcf.read_file('example.vcf.gz', compression=True)
    """
    meta = []
    skip_rows = 0
    if fn.endswith('.gz') or compression:
        f = bgzf.open(fn, 'rt')
    else:
        f = open(fn)
    for line in f:
        if line.startswith('##'):
            meta.append(line.strip())
            skip_rows += 1
        else:
            break
    df = pd.read_table(fn, skiprows=skip_rows)
    df = df.rename(columns={'#CHROM': 'CHROM'})
    vf = VcfFrame(meta, df)
    f.close()
    return vf

def row_hasindel(r):
    """Return True if the row has an indel.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    bool
        True if the row has an indel.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102, 103],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'CT', 'A', 'C'],
    ...     'ALT': ['A', 'C', 'C,AT', 'A'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': ['.', '.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '1/2', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF   ALT QUAL FILTER INFO FORMAT Steven
    0  chr1  100  .   G     A    .      .    .     GT    0/1
    1  chr1  101  .  CT     C    .      .    .     GT    0/1
    2  chr1  102  .   A  C,AT    .      .    .     GT    1/2
    3  chr1  103  .   C     A    .      .    .     GT    0/1
    >>> vf.df.apply(pyvcf.row_hasindel, axis=1)
    0    False
    1     True
    2     True
    3    False
    dtype: bool
    """
    ref_has = len(r['REF']) > 1
    alt_has = max([len(x) for x in r['ALT'].split(',')]) > 1
    return ref_has or alt_has

def row_parseinfo(r, key):
    """Return INFO data in the row that match the given key.

    Parameters
    ----------
    r : pandas.Series
        VCF row.
    key : str
        INFO key.

    Returns
    -------
    str
        Requested INFO data. Empty string if the key is not found.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102, 103],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'T', 'A', 'C'],
    ...     'ALT': ['A', 'C', 'T', 'A'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': ['DB;AC=0', 'DB;H2;AC=1', 'DB;H2;AC=1', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/0', '0/1', '0/1', '0/0'],
    ... }
    >>>
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER        INFO FORMAT Steven
    0  chr1  100  .   G   A    .      .     DB;AC=0     GT    0/0
    1  chr1  101  .   T   C    .      .  DB;H2;AC=1     GT    0/1
    2  chr1  102  .   A   T    .      .  DB;H2;AC=1     GT    0/1
    3  chr1  103  .   C   A    .      .           .     GT    0/0
    >>> vf.df.apply(pyvcf.row_parseinfo, args=('AC',), axis=1)
    0    0
    1    1
    2    1
    3
    dtype: object
    """
    result = ''
    for field in r.INFO.split(';'):
        if field.startswith(f'{key}='):
            result = field[len(key)+1:]
    return result

def row_missval(r):
    """Return the correctly formatted missing value for the row.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    str
        Missing value.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chrX'],
    ...     'POS': [100, 101, 102, 100],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'T', 'A', 'C'],
    ...     'ALT': ['A', 'C', 'T', 'A'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': ['.', '.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT:AD', 'GT:AD:DP', 'GT'],
    ...     'Steven': ['0/1', '0/1:14,15', '0/1:13,19:32', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT        Steven
    0  chr1  100  .   G   A    .      .    .        GT           0/1
    1  chr1  101  .   T   C    .      .    .     GT:AD     0/1:14,15
    2  chr1  102  .   A   T    .      .    .  GT:AD:DP  0/1:13,19:32
    3  chrX  100  .   C   A    .      .    .        GT           0/1
    >>> vf.df.apply(pyvcf.row_missval, axis=1)
    0        ./.
    1      ./.:.
    2    ./.:.:.
    3          .
    dtype: object
    """
    if 'X' in r.CHROM or 'Y' in r.CHROM:
        m = '.'
    else:
        m = './.'
    for i in range(1, len(r.FORMAT.split(':'))):
        m += ':.'
    return m

# -- VcfFrame ----------------------------------------------------------------

class VcfFrame:
    """Class for storing VCF data.

    Parameters
    ----------
    meta : list
        List of metadata lines.
    df : pandas.DataFrame
        DataFrame containing VCF data.

    See Also
    --------
    VcfFrame.from_dict
        Construct VcfFrame from dict of array-like or dicts.
    read_file
        Read a VCF file into VcfFrame.

    Examples
    --------
    Constructing VcfFrame from pandas DataFrame:

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102],
    ...     'ID': ['.', '.', '.',],
    ...     'REF': ['G', 'T', 'A'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['.', '.', '.'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '0/1'],
    ... }
    >>>
    >>> df = pd.DataFrame(data)
    >>> vf = pyvcf.VcfFrame(['##fileformat=VCFv4.3'], df)
    >>> vf.meta
    ['##fileformat=VCFv4.3']
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
    0  chr1  100  .   G   A    .      .    .     GT    0/1
    1  chr1  101  .   T   C    .      .    .     GT    0/1
    2  chr1  102  .   A   T    .      .    .     GT    0/1
    """
    def __init__(self, meta, df):
        self._meta = meta
        self._df = df.reset_index(drop=True)

    @property
    def meta(self):
        """list : List of metadata lines."""
        return self._meta

    @meta.setter
    def meta(self, value):
        self._meta = value

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing VCF data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @property
    def samples(self):
        """list : List of the sample IDs."""
        return self.df.columns[9:].to_list()

    @property
    def shape(self):
        """tuple : Tuple representing the dimensionality of the VcfFrame."""
        return (self.df.shape[0], len(self.samples))

    @classmethod
    def from_dict(cls, meta, data):
        """Construct VcfFrame from dict of array-like or dicts.

        Parameters
        ----------
        meta : list
            List of the metadata lines.
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        VcfFrame
            VcfFrame

        See Also
        --------
        VcfFrame
            VcfFrame object creation using constructor.

        Examples
        --------
        Below is a simple example:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr2'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT', 'GT'],
        ...     'Steven': ['0/1', '1/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr2  101  .   T   C    .      .    .     GT    1/1
        """
        return cls(meta, pd.DataFrame(data))

    def compare(self, a, b, c=None):
        """Compare genotype data between Samples A & B or Samples A, B, & C.

        This method will return (Ab, aB, AB, ab) for two samples and
        (Abc, aBc, ABc, abC, AbC, aBC, ABC, abc) for three samples. Note that
        the former is equivalent to (FP, FN, TP, TN) if we assume A is
        the test sample and B is the truth sample.

        Parameters
        ----------
        a : str or int
            Name or index of Sample A (or test).
        b : str or int
            Name or index of Sample B (or truth).
        c : str or int, optional
            Name or index of Sample C.

        Returns
        -------
        tuple
            Four- or eight-element tuple depending on the number of samples.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'C'],
        ...     'ALT': ['A', 'C', 'A', 'T'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '0/0', '0/0', '0/1'],
        ...     'Sara': ['0/1', '0/1', './.', '0/1'],
        ...     'James': ['0/1', '0/1', '1/1', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/1   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/1
        2  chr1  102  .   T   A    .      .    .     GT    0/0  ./.   1/1
        3  chr1  103  .   C   T    .      .    .     GT    0/1  0/1   ./.

        We compare Steven and Sara:

        >>> vf.compare('Steven', 'Sara')
        (0, 1, 2, 1)

        Next, we compare all three:

        >>> vf.compare('Steven', 'Sara', 'James')
        (0, 0, 1, 1, 0, 1, 1, 0)
        """
        if c is None:
            result = self._compare_two(a, b)
        else:
            result = self._compare_three(a, b, c)
        return result

    def _compare_two(self, a, b):
        a = a if isinstance(a, str) else self.samples[a]
        b = b if isinstance(b, str) else self.samples[b]
        def func(r):
            a_has = gt_hasvar(r[a])
            b_has = gt_hasvar(r[b])
            if a_has and not b_has:
                return 'Ab'
            elif not a_has and b_has:
                return 'aB'
            elif a_has and b_has:
                return 'AB'
            else:
                return 'ab'
        d = self.df.apply(func, axis=1).value_counts().to_dict()
        Ab = d['Ab'] if 'Ab' in d else 0
        aB = d['aB'] if 'aB' in d else 0
        AB = d['AB'] if 'AB' in d else 0
        ab = d['ab'] if 'ab' in d else 0
        return (Ab, aB, AB, ab)

    def _compare_three(self, a, b, c):
        a = a if isinstance(a, str) else self.samples[a]
        b = b if isinstance(b, str) else self.samples[b]
        c = c if isinstance(c, str) else self.samples[c]
        def func(r):
            a_has = gt_hasvar(r[a])
            b_has = gt_hasvar(r[b])
            c_has = gt_hasvar(r[c])
            if a_has and not b_has and not c_has:
                return 'Abc'
            elif not a_has and b_has and not c_has:
                return 'aBc'
            elif a_has and b_has and not c_has:
                return 'ABc'
            elif not a_has and not b_has and c_has:
                return 'abC'
            elif a_has and not b_has and c_has:
                return 'AbC'
            elif not a_has and b_has and c_has:
                return 'aBC'
            elif a_has and b_has and c_has:
                return 'ABC'
            else:
                return 'abc'
        d = self.df.apply(func, axis=1).value_counts().to_dict()
        Abc = d['Abc'] if 'Abc' in d else 0
        aBc = d['aBc'] if 'aBc' in d else 0
        ABc = d['ABc'] if 'ABc' in d else 0
        abC = d['abC'] if 'abC' in d else 0
        AbC = d['AbC'] if 'AbC' in d else 0
        aBC = d['aBC'] if 'aBC' in d else 0
        ABC = d['ABC'] if 'ABC' in d else 0
        abc = d['abc'] if 'abc' in d else 0
        return (Abc, aBc, ABc, abC, AbC, aBC, ABC, abc)

    def collapse(self):
        """Collapse duplicate records in the VcfFrame.

        Duplicate records have the identical values for CHROM, POS, and REF.
        They can result from merging two VCF files.

        .. note::
           The method will sort the order of ALT alleles.

        Returns
        -------
        VcfFrame
            Collapsed VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
        ...     'POS': [100, 100, 200, 200],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:AD', 'GT:AD', 'GT:AD', 'GT:AD'],
        ...     'Steven': ['0/1:12,15', './.:.', '0/1:16,12', './.:.'],
        ...     'Sara': ['./.:.', '0/1:14,15', './.:.', '1/2:0,11,17'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT     Steven         Sara
        0  chr1  100  .   A    C    .      .    .  GT:AD  0/1:12,15        ./.:.
        1  chr1  100  .   A    T    .      .    .  GT:AD      ./.:.    0/1:14,15
        2  chr2  200  .   C    G    .      .    .  GT:AD  0/1:16,12        ./.:.
        3  chr2  200  .   C  G,A    .      .    .  GT:AD      ./.:.  1/2:0,11,17

        We collapse the VcfFrame:

        >>> vf.collapse().df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT       Steven         Sara
        0  chr1  100  .   A  C,T    .      .    .  GT:AD  0/1:12,15,0  0/2:14,0,15
        2  chr2  200  .   C  A,G    .      .    .  GT:AD  0/2:16,0,12  1/2:0,17,11
        """
        df = self.df.copy()
        dup_idx = df.duplicated(['CHROM', 'POS', 'REF'], keep=False)
        dups = {}
        for i, r in df[dup_idx].iterrows():
            name = f'{r.CHROM}:{r.POS}:{r.REF}'
            if name not in dups:
                dups[name] = []
            dups[name].append(i)

        def collapse_one(df):
            ref_allele = df.REF.unique()[0]
            alt_alleles = []
            for i, r in df.iterrows():
                alt_alleles += r.ALT.split(',')
            alt_alleles = sorted(list(set(alt_alleles)),
                                 key=lambda x: (len(x), x))
            all_alleles = [ref_allele] + alt_alleles

            def infunc(x, r_all_alleles, index_map):
                if gt_missing(x):
                    return ''
                old_fields = x.split(':')
                old_gt = old_fields[0]
                new_gt = '/'.join([str(x) for x in
                    sorted([index_map[int(i)] for i in old_gt.split('/')])
                ])
                new_fields = [new_gt]
                for old_field in old_fields[1:]:
                    old_subfields = old_field.split(',')
                    new_subfields = ['0' for x in all_alleles]
                    if len(old_subfields) == len(r_all_alleles):
                        for i, old_subfield in enumerate(old_subfields):
                            new_subfields[index_map[i]] = old_subfield
                        new_fields.append(','.join(new_subfields))
                    else:
                        new_fields.append(old_field)
                return ':'.join(new_fields)

            def outfunc(r):
                r_alt_alleles = r.ALT.split(',')
                r_all_alleles = [r.REF] + r_alt_alleles
                old_indicies = [i for i in range(len(r_all_alleles))]
                new_indicies = [all_alleles.index(x) for x in r_all_alleles]
                index_map = dict(zip(old_indicies, new_indicies))
                r[9:] = r[9:].apply(infunc, args=(r_all_alleles, index_map))
                return r

            df2 = df.apply(outfunc, axis=1)

            def raise_error(c):
                if sum(c.values != '') > 1:
                    message = ('cannot collapse following '
                               f'records:\n{df.loc[c.index]}')
                    raise ValueError(message)

            df2.iloc[:, 9:].apply(raise_error)
            df2 = df2.groupby(['CHROM', 'POS', 'REF']).agg(''.join)
            df2 = df2.reset_index()
            cols = list(df2)
            cols[2], cols[3] = cols[3], cols[2]
            df2 = df2[cols]
            df2.ID = df.ID.unique()[0]
            df2.ALT = ','.join(alt_alleles)
            df2.QUAL = df.QUAL.unique()[0]
            df2.FILTER = df.FILTER.unique()[0]
            df2.INFO = df.INFO.unique()[0]
            df2.FORMAT = df.FORMAT.unique()[0]
            s = df2.squeeze()
            s = s.replace('', row_missval(s))
            return s

        for name, i in dups.items():
            df.iloc[i] = collapse_one(df.iloc[i])
        df.drop_duplicates(subset=['CHROM', 'POS', 'REF'], inplace=True)

        vf = self.__class__(self.copy_meta(), df)
        return vf

    def combine(self, a, b):
        """Combine the genotype data of Sample A and Sample B.

        This method is useful when, for example, you are trying to
        consolidate data from multiple replicate samples. When the same
        variant is found (or not found) in both samples, the method will
        use the genotype data of the first sample.

        Parameters
        ----------
        a : str or int
            Name or index of the first sample (or original).
        b : str or int
            Name or index of the second sample (or replicate).

        Returns
        -------
        pandas.Series
            Resulting VCF column.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'T'],
        ...     'ALT': ['A', 'C', 'A'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Original': ['./.:.', '0/0:29', '0/1:28'],
        ...     'Replicate': ['0/1:24', '0/1:30', './.:.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Original Replicate
        0  chr1  100  .   G   A    .      .    .  GT:DP    ./.:.    0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP   0/0:29    0/1:30
        2  chr1  102  .   T   A    .      .    .  GT:DP   0/1:28     ./.:.

        We combine the two samples to get consolidated genotype data:

        >>> vf.df['Combined'] = vf.combine('Original', 'Replicate')
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Original Replicate Combined
        0  chr1  100  .   G   A    .      .    .  GT:DP    ./.:.    0/1:24   0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP   0/0:29    0/1:30   0/1:30
        2  chr1  102  .   T   A    .      .    .  GT:DP   0/1:28     ./.:.   0/1:28
        """
        a = a if isinstance(a, str) else self.samples[a]
        b = b if isinstance(b, str) else self.samples[b]
        def func(r):
            a_has = gt_hasvar(r[a])
            b_has = gt_hasvar(r[b])
            if a_has and b_has:
                return r[a]
            elif a_has and not b_has:
                return r[a]
            elif not a_has and b_has:
                return r[b]
            else:
                return r[a]
        s = self.df.apply(func, axis=1)
        return s

    def copy_meta(self):
        """Return a copy of the metadata."""
        return deepcopy(self.meta)

    def copy_df(self):
        """Return a copy of the dataframe."""
        return self.df.copy()

    def copy(self):
        """Return a copy of the VcfFrame."""
        return self.__class__(self.copy_meta(), self.copy_df())

    def to_file(self, fn, compression=False):
        """Write the VcfFrame to a VCF file.

        If the file name ends with '.gz', the method will automatically
        use the BGZF compression when writing the file.

        Parameters
        ----------
        fn : str
            VCF file path.
        compression : bool, default: False
            If True, use the BGZF compression.

        Examples
        --------

        >>> data = {
        ...     'CHROM': ['chr1', 'chr2'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT', 'GT'],
        ...     'Steven': ['0/1', '1/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict(['##fileformat=VCFv4.3'], data)
        >>> vf.to_file('unzipped.vcf')
        >>> vf.to_file('zipped.vcf.gz')
        >>> vf.to_file('zipped.vcf.gz', compression=True)
        """
        if fn.endswith('.gz') or compression:
            f = bgzf.open(fn, 'w')
        else:
            f = open(fn, 'w')
        f.write(self.to_string())
        f.close()

    def to_string(self):
        """Render the VcfFrame to a console-friendly tabular output.

        Returns
        -------
        str
            String representation of the VcfFrame.

        Examples
        --------

        >>> data = {
        ...     'CHROM': ['chr1', 'chr2'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT', 'GT'],
        ...     'Steven': ['0/1', '1/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict(['##fileformat=VCFv4.3'], data)
        >>> print(vf.to_string())
        ##fileformat=VCFv4.3
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Steven
        chr1	100	.	G	A	.	.	.	GT	0/1
        chr2	101	.	T	C	.	.	.	GT	1/1
        """
        s = ''
        if self.meta:
            s += '\n'.join(self.meta) + '\n'
        s += self.df.rename(columns={'CHROM': '#CHROM'}
            ).to_csv(index=False, sep='\t')
        return s

    def strip(self, format='GT'):
        """Remove unnecessary data from the VcfFrame.

        Parameters
        ----------
        format : str, default: 'GT'
            FORMAT subfields to be retained (e.g. 'GT:AD:DP').

        Returns
        -------
        VcfFrame
            Stripped VcfFrame.
        """
        def outfunc(r):
            idx = [r.FORMAT.split(':').index(x) for x in format.split(':')]
            infunc = lambda x: ':'.join([x.split(':')[i] for i in idx])
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        df = self.df.copy()
        df[['ID', 'QUAL', 'FILTER', 'INFO']] = '.'
        df = df.apply(outfunc, axis=1)
        df.FORMAT = format
        vf = self.__class__([], df)
        return vf

    def merge(self, other, how='inner', format='GT', sort=True,
              collapse=False):
        """Merge with the other VcfFrame.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        how : str, default: 'inner'
            Type of merge as defined in `pandas.DataFrame.merge`.
        format : str, default: 'GT'
            FORMAT subfields to be retained (e.g. 'GT:AD:DP').
        sort : bool, default: True
            If True, sort the VcfFrame before returning.
        collapse : bool, default: False
            If True, collapse duplicate records.

        Returns
        -------
        VcfFrame
            Merged VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data1 = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP'],
        ...     'Steven': ['0/0:32', '0/1:29'],
        ...     'Sara': ['0/1:24', '1/1:30'],
        ... }
        >>> data2 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2'],
        ...     'POS': [100, 101, 200],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Dona': ['./.:.', '0/0:24', '0/0:26'],
        ...     'Michel': ['0/1:24', '0/1:31', '0/1:26'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> vf1.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:32  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  1/1:30
        >>> vf2.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Dona  Michel
        0  chr1  100  .   G   A    .      .    .  GT:DP   ./.:.  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/0:24  0/1:31
        2  chr2  200  .   A   T    .      .    .  GT:DP  0/0:26  0/1:26

        We can merge the two VcfFrames with ``how='inner'`` (default):

        >>> vf1.merge(vf2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara Dona Michel
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/1  ./.    0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/1  1/1  0/0    0/1

        We can also merge with ``how='outer'``:

        >>> vf1.merge(vf2, how='outer').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara Dona Michel
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/1  ./.    0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/1  1/1  0/0    0/1
        2  chr2  200  .   A   T    .      .    .     GT    ./.  ./.  0/0    0/1

        Since both VcfFrames have the DP subfield, we can use ``format='GT:DP'``:

        >>> vf1.merge(vf2, how='outer', format='GT:DP').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara    Dona  Michel
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:32  0/1:24   ./.:.  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  1/1:30  0/0:24  0/1:31
        2  chr2  200  .   A   T    .      .    .  GT:DP   ./.:.   ./.:.  0/0:26  0/1:26
        """
        vf1 = self.strip(format=format)
        vf2 = other.strip(format=format)
        dropped = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        shared = ['CHROM', 'POS', 'REF', 'ALT']
        df = vf1.df.merge(vf2.df.drop(columns=dropped), on=shared, how=how)

        # This ensures that the column order is intact when either of the
        # dataframes is empty.
        cols = vf1.df.columns.to_list() + vf2.df.columns[9:].to_list()
        df = df[cols]

        df[dropped] = df[dropped].fillna('.')
        df.FORMAT = format
        def func(r):
            n = len(r.FORMAT.split(':'))
            x = './.'
            for i in range(1, n):
                x += ':.'
            r = r.fillna(x)
            return r
        df = df.apply(func, axis=1)
        vf3 = self.__class__([], df)
        if collapse:
            vf3 = vf3.collapse()
        if sort:
            vf3 = vf3.sort()
        return vf3

    def meta_keys(self):
        """Print metadata lines with a key."""
        for line in self.meta:
            if '=<ID=' in line:
                print(line)

    def add_dp(self):
        """Compute DP using AD and add it to the FORMAT field.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
        ...     'POS': [100, 100, 200, 200],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:AD', 'GT:AD', 'GT:AD', 'GT:AD'],
        ...     'Steven': ['0/1:12,15', '0/0:32,1', '0/1:16,12', './.:.'],
        ...     'Sara': ['0/1:13,17', '0/1:14,15', './.:.', '1/2:0,11,17'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT     Steven         Sara
        0  chr1  100  .   A    C    .      .    .  GT:AD  0/1:12,15    0/1:13,17
        1  chr1  100  .   A    T    .      .    .  GT:AD   0/0:32,1    0/1:14,15
        2  chr2  200  .   C    G    .      .    .  GT:AD  0/1:16,12        ./.:.
        3  chr2  200  .   C  G,A    .      .    .  GT:AD      ./.:.  1/2:0,11,17

        We can add the DP subfield to our genotype data:

        >>> vf.add_dp().df
          CHROM  POS ID REF  ALT QUAL FILTER INFO    FORMAT        Steven            Sara
        0  chr1  100  .   A    C    .      .    .  GT:AD:DP  0/1:12,15:27    0/1:13,17:30
        1  chr1  100  .   A    T    .      .    .  GT:AD:DP   0/0:32,1:33    0/1:14,15:29
        2  chr2  200  .   C    G    .      .    .  GT:AD:DP  0/1:16,12:28         ./.:.:.
        3  chr2  200  .   C  G,A    .      .    .  GT:AD:DP       ./.:.:.  1/2:0,11,17:28
        """
        def outfunc(r):
            i = r.FORMAT.split(':').index('AD')
            def infunc(x):
                ad = x.split(':')[i].split(',')
                dp = 0
                for depth in ad:
                    if depth == '.':
                        return f'{x}:.'
                    dp += int(depth)
                return f'{x}:{dp}'
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            r.FORMAT += ':DP'
            return r
        df = self.df.apply(outfunc, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def add_flag(self, flag, order='last', index=None):
        """Add the given flag to the INFO field.

        The default behavior is to add the flag to all rows in the VcfFrame.

        Parameters
        ----------
        flag : str
            INFO flag.
        order : {'last', 'first', False}, default: 'last'
            Determines the order in which the flag will be added.

            - ``last`` : Add to the end of the list.
            - ``first`` : Add to the beginning of the list.
            - ``False`` : Overwrite the existing field.

        index : list or pandas.Series, optional
            Boolean index array indicating which rows should be updated.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', 'DB', 'DB', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', '0/1', '0/1', '1/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/0
        1  chr1  101  .   T   C    .      .   DB     GT    0/1
        2  chr1  102  .   A   T    .      .   DB     GT    0/1
        3  chr1  103  .   C   A    .      .    .     GT    1/1

        We can add the SOMATIC flag to the INFO field:

        >>> vf.add_flag('SOMATIC').df
          CHROM  POS ID REF ALT QUAL FILTER        INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .     SOMATIC     GT    0/0
        1  chr1  101  .   T   C    .      .  DB;SOMATIC     GT    0/1
        2  chr1  102  .   A   T    .      .  DB;SOMATIC     GT    0/1
        3  chr1  103  .   C   A    .      .     SOMATIC     GT    1/1

        Setting ``order='first'`` will append the flag at the beginning:

        >>> vf.add_flag('SOMATIC', order='first').df
          CHROM  POS ID REF ALT QUAL FILTER        INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .     SOMATIC     GT    0/0
        1  chr1  101  .   T   C    .      .  SOMATIC;DB     GT    0/1
        2  chr1  102  .   A   T    .      .  SOMATIC;DB     GT    0/1
        3  chr1  103  .   C   A    .      .     SOMATIC     GT    1/1

        Setting ``order=False`` will overwrite the INFO field:

        >>> vf.add_flag('SOMATIC', order=False).df
          CHROM  POS ID REF ALT QUAL FILTER     INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .  SOMATIC     GT    0/0
        1  chr1  101  .   T   C    .      .  SOMATIC     GT    0/1
        2  chr1  102  .   A   T    .      .  SOMATIC     GT    0/1
        3  chr1  103  .   C   A    .      .  SOMATIC     GT    1/1

        We can also specify which rows should be updated:

        >>> vf.add_flag('SOMATIC', index=[True, True, False, False]).df
          CHROM  POS ID REF ALT QUAL FILTER        INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .     SOMATIC     GT    0/0
        1  chr1  101  .   T   C    .      .  DB;SOMATIC     GT    0/1
        2  chr1  102  .   A   T    .      .          DB     GT    0/1
        3  chr1  103  .   C   A    .      .           .     GT    1/1
        """
        if index is None:
            index = [True for i in range(self.shape[0])]
        def f(r):
            if not index[r.name]:
                return r
            if r.INFO == '.':
                r.INFO = flag
            elif not order:
                r.INFO = flag
            elif order == 'first':
                r.INFO = f'{flag};{r.INFO}'
            else:
                r.INFO += f';{flag}'
            return r
        df = self.df.apply(f, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def markmiss_ad(self, threshold, samples=None, full=False, as_nan=False):
        """Mark genotypes whose AD is below threshold as missing.

        By default, the marking is applied only to genotypes with ALT allele.

        Parameters
        ----------
        threshold : int
            Minimum allele depth.
        sampels : list, optional
            If provided, only these samples will be marked.
        full : bool, default: False
            If True, apply the marking to all genotypes.
        as_nan : bool, default: False
            If True, mark as ``NaN`` instead of missing value.

        Returns
        -------
        VcfFrame
            Filtered VcfFrame.

        See Also
        --------
        markmiss_af
            Similar method using AF.
        markmiss_dp
            Similar method using DP.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:AD', 'GT:AD'],
        ...     'Steven': ['0/1:15,13', '0/0:28,1'],
        ...     'Sara': ['0/1:15,15', '0/1:14,18'],
        ...     'James': ['0/0:32,0', '0/1:17,5'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT     Steven       Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AD  0/1:15,13  0/1:15,15  0/0:32,0
        1  chr1  101  .   T   C    .      .    .  GT:AD   0/0:28,1  0/1:14,18  0/1:17,5

        We mark all the genotypes whose AD is below 15 as missing:

        >>> vf.markmiss_ad(15).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven       Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AD     ./.:.  0/1:15,15  0/0:32,0
        1  chr1  101  .   T   C    .      .    .  GT:AD  0/0:28,1  0/1:14,18     ./.:.

        We can apply the marking only to a subset of the samples:

        >>> vf.markmiss_ad(15, samples=['Steven', 'Sara']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven       Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AD     ./.:.  0/1:15,15  0/0:32,0
        1  chr1  101  .   T   C    .      .    .  GT:AD  0/0:28,1  0/1:14,18  0/1:17,5

        We can mark all genotypes including those without ALT allele:

        >>> vf.markmiss_ad(15, full=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven       Sara  James
        0  chr1  100  .   G   A    .      .    .  GT:AD  ./.:.  0/1:15,15  ./.:.
        1  chr1  101  .   T   C    .      .    .  GT:AD  ./.:.  0/1:14,18  ./.:.

        We can mark as ``NaN`` instead of missing value:

        >>> vf.markmiss_ad(15, as_nan=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven       Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AD       NaN  0/1:15,15  0/0:32,0
        1  chr1  101  .   T   C    .      .    .  GT:AD  0/0:28,1  0/1:14,18       NaN
        """
        if samples is None:
            samples = self.samples
        def one_row(r):
            i = r.FORMAT.split(':').index('AD')
            if as_nan:
                m = np.nan
            else:
                m = row_missval(r)
            def one_gt(g):
                if not full and not gt_hasvar(g):
                    return g
                s = g.split(':')[i].split(',')[1]
                if s == '.' or int(s) < threshold:
                    return m
                return g
            r[samples] = r[samples].apply(one_gt)
            return r
        df = self.df.apply(one_row, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def markmiss_af(self, threshold, samples=None, full=False, as_nan=False):
        """Mark genotypes whose AF is below threshold as missing.

        By default, the marking is applied only to genotypes with ALT allele.

        Parameters
        ----------
        threshold : int
            Minimum allele fraction.
        sampels : list, optional
            If provided, only these samples will be marked.
        full : bool, default: False
            If True, apply the marking to all genotypes.
        as_nan : bool, default: False
            If True, mark as ``NaN`` instead of missing value.

        Returns
        -------
        VcfFrame
            Filtered VcfFrame.

        See Also
        --------
        markmiss_ad
            Similar method using AD.
        markmiss_dp
            Similar method using DP.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:AF', 'GT:AF'],
        ...     'Steven': ['0/0:0.01', '0/1:0.31'],
        ...     'Sara': ['0/1:0.12', '0/1:0.25'],
        ...     'James': ['0/1:0.11', '0/0:0.09'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven      Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AF  0/0:0.01  0/1:0.12  0/1:0.11
        1  chr1  101  .   T   C    .      .    .  GT:AF  0/1:0.31  0/1:0.25  0/0:0.09

        We mark all the genotypes whose AF is below 0.3 as missing:

        >>> vf.markmiss_af(0.3).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven   Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AF  0/0:0.01  ./.:.     ./.:.
        1  chr1  101  .   T   C    .      .    .  GT:AF  0/1:0.31  ./.:.  0/0:0.09

        We can apply the marking only to a subset of the samples:

        >>> vf.markmiss_af(0.3, samples=['Steven', 'Sara']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven   Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AF  0/0:0.01  ./.:.  0/1:0.11
        1  chr1  101  .   T   C    .      .    .  GT:AF  0/1:0.31  ./.:.  0/0:0.09

        We can mark all genotypes including those without ALT allele:

        >>> vf.markmiss_af(0.3, full=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven   Sara  James
        0  chr1  100  .   G   A    .      .    .  GT:AF     ./.:.  ./.:.  ./.:.
        1  chr1  101  .   T   C    .      .    .  GT:AF  0/1:0.31  ./.:.  ./.:.

        We can mark as ``NaN`` instead of missing value:

        >>> vf.markmiss_af(0.3, as_nan=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Steven  Sara     James
        0  chr1  100  .   G   A    .      .    .  GT:AF  0/0:0.01   NaN       NaN
        1  chr1  101  .   T   C    .      .    .  GT:AF  0/1:0.31   NaN  0/0:0.09
        """
        if samples is None:
            samples = self.samples
        def one_row(r):
            i = r.FORMAT.split(':').index('AF')
            if as_nan:
                m = np.nan
            else:
                m = row_missval(r)
            def one_gt(g):
                if not full and not gt_hasvar(g):
                    return g
                s = g.split(':')[i]
                if s == '.' or float(s) < threshold:
                    return m
                return g
            r[samples] = r[samples].apply(one_gt)
            return r
        df = self.df.apply(one_row, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def markmiss_dp(self, threshold, samples=None, full=False, as_nan=False):
        """Mark genotypes whose DP is below threshold as missing.

        By default, the marking is applied only to genotypes with ALT allele.

        Parameters
        ----------
        threshold : int
            Minimum read depth.
        sampels : list, optional
            If provided, only these samples will be marked.
        full : bool, default: False
            If True, apply the marking to all genotypes.
        as_nan : bool, default: False
            If True, mark as ``NaN`` instead of missing value.

        Returns
        -------
        VcfFrame
            Filtered VcfFrame.

        See Also
        --------
        markmiss_ad
            Similar method using allele depth.
        markmiss_af
            Similar method using allele fraction.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP'],
        ...     'Steven': ['0/0:26', '0/1:29'],
        ...     'Sara': ['0/1:24', '0/1:30'],
        ...     'James': ['0/1:18', '0/0:24'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara   James
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:26  0/1:24  0/1:18
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  0/1:30  0/0:24

        We mark all the genotypes whose DP is below 30 as missing:

        >>> vf.markmiss_dp(30).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara   James
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:26   ./.:.   ./.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP   ./.:.  0/1:30  0/0:24

        We can apply the marking only to a subset of the samples:

        >>> vf.markmiss_dp(30, samples=['Steven', 'Sara']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara   James
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:26   ./.:.  0/1:18
        1  chr1  101  .   T   C    .      .    .  GT:DP   ./.:.  0/1:30  0/0:24

        We can mark all genotypes including those without ALT allele:

        >>> vf.markmiss_dp(30, full=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven    Sara  James
        0  chr1  100  .   G   A    .      .    .  GT:DP  ./.:.   ./.:.  ./.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP  ./.:.  0/1:30  ./.:.

        We can mark as ``NaN`` instead of missing value:

        >>> vf.markmiss_dp(30, as_nan=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara   James
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:26     NaN     NaN
        1  chr1  101  .   T   C    .      .    .  GT:DP     NaN  0/1:30  0/0:24
        """
        if samples is None:
            samples = self.samples
        def one_row(r):
            i = r.FORMAT.split(':').index('DP')
            if as_nan:
                m = np.nan
            else:
                m = row_missval(r)
            def one_gt(g):
                if not full and not gt_hasvar(g):
                    return g
                s = g.split(':')[i]
                if s == '.' or int(s) < threshold:
                    return m
                return g
            r[samples] = r[samples].apply(one_gt)
            return r
        df = self.df.apply(one_row, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_bed(self, bed, opposite=False, index=False):
        """Select rows that overlap with the given BED data.

        Parameters
        ----------
        bed : pybed.BedFrame or str
            BedFrame or path to a BED file.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'Chromosome': ['chr1', 'chr2', 'chr3'],
        ...     'Start': [100, 400, 100],
        ...     'End': [200, 500, 200]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1    100  200
        1       chr2    400  500
        2       chr3    100  200
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr3'],
        ...     'POS': [100, 201, 450, 99],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'CT', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'AT', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '1/1', '0/1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  201  .  CT   C    .      .    .     GT    1/1
        2  chr2  450  .   A  AT    .      .    .     GT    0/1
        3  chr3   99  .   C   A    .      .    .     GT    0/1

        We can select rows that overlap with the BED data:

        >>> vf.filter_bed(bf).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr2  450  .   A  AT    .      .    .     GT    0/1

        We can also remove those rows:

        >>> vf.filter_bed(bf, opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  201  .  CT   C    .      .    .     GT    1/1
        1  chr3   99  .   C   A    .      .    .     GT    0/1

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_bed(bf, index=True)
        0     True
        1    False
        2     True
        3    False
        dtype: bool
        >>>
        """
        if isinstance(bed, pybed.BedFrame):
            bf = bed
        else:
            bf = pybed.read_file(bed)
        f = lambda r: not bf.gr[r.CHROM, r.POS:r.POS+1].empty
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        vf = self.__class__(self.copy_meta(), self.df[i])
        return vf

    def filter_empty(self, opposite=False, index=False):
        """Remove rows with no genotype calls at all.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', './.', '0/1', './.'],
        ...     'Sara': ['0/0', './.', './.', './.'],
        ...     'James': ['0/0', './.', '0/1', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.filter_indel().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/0   0/0
        1  chr1  101  .   T   C    .      .    .     GT    ./.  ./.   ./.
        2  chr1  102  .   A   T    .      .    .     GT    0/1  ./.   0/1
        3  chr1  103  .   C   A    .      .    .     GT    ./.  ./.   ./.

        We can remove empty rows:

        >>> vf.filter_empty().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/0   0/0
        1  chr1  102  .   A   T    .      .    .     GT    0/1  ./.   0/1

        We can also select those rows:

        >>> vf.filter_empty(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  101  .   T   C    .      .    .     GT    ./.  ./.   ./.
        1  chr1  103  .   C   A    .      .    .     GT    ./.  ./.   ./.

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_empty(index=True)
        0     True
        1    False
        2     True
        3    False
        dtype: bool
        """
        f = lambda r: not all(r.iloc[9:].apply(gt_missing))
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_indel(self, opposite=False, index=False):
        """Remove rows with an indel.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'CT', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'C,AT', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '0/1', '1/2', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF   ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G     A    .      .    .     GT    0/1
        1  chr1  101  .  CT     C    .      .    .     GT    0/1
        2  chr1  102  .   A  C,AT    .      .    .     GT    1/2
        3  chr1  103  .   C     A    .      .    .     GT    0/1

        We can remove rows with an indel:

        >>> vf.filter_indel().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  103  .   C   A    .      .    .     GT    0/1

        We can also select those rows:

        >>> vf.filter_indel(opposite=True).df
          CHROM  POS ID REF   ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  101  .  CT     C    .      .    .     GT    0/1
        1  chr1  102  .   A  C,AT    .      .    .     GT    1/2

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_indel(index=True)
        0     True
        1    False
        2    False
        3     True
        dtype: bool
        """
        i = ~self.df.apply(row_hasindel, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagall(self, flags, opposite=False, index=False):
        """Select rows if all of the given INFO flags are present.

        Parameters
        ----------
        flags : list
            List of INFO flags.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        filter_flagany
            Similar method that selects rows if any one of the given
            INFO flags is present.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['DB', 'DB;H2', 'DB;H2', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', '0/1', '0/1', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER   INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .     DB     GT    0/0
        1  chr1  101  .   T   C    .      .  DB;H2     GT    0/1
        2  chr1  102  .   A   T    .      .  DB;H2     GT    0/1
        3  chr1  103  .   C   A    .      .      .     GT    0/0

        We can select rows with both the H2 and DB tags:

        >>> vf.filter_flagall(['H2', 'DB']).df
          CHROM  POS ID REF ALT QUAL FILTER   INFO FORMAT Steven
        0  chr1  101  .   T   C    .      .  DB;H2     GT    0/1
        1  chr1  102  .   A   T    .      .  DB;H2     GT    0/1

        We can also remove those rows:

        >>> vf.filter_flagall(['H2', 'DB'], opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .   DB     GT    0/0
        1  chr1  103  .   C   A    .      .    .     GT    0/0

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_flagall(['H2', 'DB'], index=True)
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        def f(r):
            l = r.INFO.split(';')
            return all([x in l for x in flags])
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagany(self, flags, opposite=False, index=False):
        """Select rows if any one of the given INFO flags is present.

        Parameters
        ----------
        flags : list
            List of INFO flags.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        filter_flagall
            Similar method that selects rows if all of the given INFO
            flags are present.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['DB', 'DB;H2', 'DB;H2', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', '0/1', '0/1', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER   INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .     DB     GT    0/0
        1  chr1  101  .   T   C    .      .  DB;H2     GT    0/1
        2  chr1  102  .   A   T    .      .  DB;H2     GT    0/1
        3  chr1  103  .   C   A    .      .      .     GT    0/0

        We can select rows with the H2 tag:

        >>> vf.filter_flagany(['H2']).df
          CHROM  POS ID REF ALT QUAL FILTER   INFO FORMAT Steven
        0  chr1  101  .   T   C    .      .  DB;H2     GT    0/1
        1  chr1  102  .   A   T    .      .  DB;H2     GT    0/1

        We can also remove those rows:

        >>> vf.filter_flagany(['H2'], opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .   DB     GT    0/0
        1  chr1  103  .   C   A    .      .    .     GT    0/0

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_flagany(['H2'], index=True)
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        def f(r):
            l = r.INFO.split(';')
            return any([x in l for x in flags])
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_multialt(self, opposite=False, index=False):
        """Remove rows with multiple ALT alleles.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------

        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
        ...     'POS': [100, 100, 200, 200],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C,T', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/2', '0/0', '0/1', './.'],
        ...     'Sara': ['0/1', '0/1', './.', '1/2'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT Steven Sara
        0  chr1  100  .   A  C,T    .      .    .     GT    0/2  0/1
        1  chr1  100  .   A    T    .      .    .     GT    0/0  0/1
        2  chr2  200  .   C    G    .      .    .     GT    0/1  ./.
        3  chr2  200  .   C  G,A    .      .    .     GT    ./.  1/2

        We can remove rows with multiple ALT alleles:

        >>> vf.filter_multialt().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara
        0  chr1  100  .   A   T    .      .    .     GT    0/0  0/1
        1  chr2  200  .   C   G    .      .    .     GT    0/1  ./.

        We can also select those rows:

        >>> vf.filter_multialt(opposite=True).df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT Steven Sara
        0  chr1  100  .   A  C,T    .      .    .     GT    0/2  0/1
        1  chr2  200  .   C  G,A    .      .    .     GT    ./.  1/2

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_multialt(index=True)
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        i = self.df.apply(lambda r: ',' not in r.ALT, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_pass(self, opposite=False, index=False):
        """Select rows with PASS in the FILTER field.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['PASS', 'FAIL', 'PASS', 'FAIL'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', './.', '0/1', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .   PASS    .     GT    0/0
        1  chr1  101  .   T   C    .   FAIL    .     GT    ./.
        2  chr1  102  .   A   T    .   PASS    .     GT    0/1
        3  chr1  103  .   C   A    .   FAIL    .     GT    ./.

        We can select rows with PASS:

        >>> vf.filter_pass().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .   PASS    .     GT    0/0
        1  chr1  102  .   A   T    .   PASS    .     GT    0/1

        We can also remove those rows:

        >>> vf.filter_pass(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  101  .   T   C    .   FAIL    .     GT    ./.
        1  chr1  103  .   C   A    .   FAIL    .     GT    ./.

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_pass(index=True)
        0     True
        1    False
        2     True
        3    False
        dtype: bool
        """
        f = lambda r: r.FILTER == 'PASS'
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_phased(self, opposite=False, index=False):
        """Remove rows with phased genotypes.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'CT', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'C', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['1|0', '0/1', '0/1', '0|1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    1|0
        1  chr1  101  .  CT   C    .      .    .     GT    0/1
        2  chr1  102  .   A   C    .      .    .     GT    0/1
        3  chr1  103  .   C   A    .      .    .     GT    0|1

        We can remove rows with a phased genotype:

        >>> vf.filter_phased().df
          CHROM  POS ID REF   ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  101  .  CT     C    .      .    .     GT    0/1
        1  chr1  102  .   A  C,AT    .      .    .     GT    0/1

        We can also select those rows:

        >>> vf.filter_phased(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    1|0
        1  chr1  103  .   C   A    .      .    .     GT    0|1

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_phased(index=True)
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        f = lambda r: not any(['|' in gt.split(':')[0] for gt in r[9:]])
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_polyp(self, opposite=False, index=False):
        """Remove rows with a polyploid genotype call.

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
        ...     'POS': [100, 100, 200, 200],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0/1', '0/0', '1/1/1', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   A   C    .      .    .     GT  0/0/1
        1  chr1  100  .   A   T    .      .    .     GT    0/0
        2  chr2  200  .   C   G    .      .    .     GT  1/1/1
        3  chr2  200  .   C   G    .      .    .     GT    ./.

        We can remove rows with a polyploid genotype call:

        >>> vf.filter_polyp().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   A   T    .      .    .     GT    0/0
        1  chr2  200  .   C   G    .      .    .     GT    ./.

        We can also select those rows:

        >>> vf.filter_polyp(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   A   C    .      .    .     GT  0/0/1
        1  chr2  200  .   C   G    .      .    .     GT  1/1/1

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_polyp(index=True)
        0    False
        1     True
        2    False
        3     True
        dtype: bool
        """
        f = lambda r: not any([gt_haspolyp(x) for x in r[9:]])
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampall(self, samples=None, opposite=False, index=False):
        """Select rows if all of the given samples have the variant.

        The default behavior is to use all samples in the VcfFrame.

        Parameters
        ----------
        samples : list, optional
            List of sample names or indicies.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        filter_sampany
            Similar method that selects rows if any one of the given
            samples has the variant.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'T'],
        ...     'ALT': ['A', 'C', 'A', 'C'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '0/0', '0/1', '0/1'],
        ...     'Sara': ['0/1', '0/1', '0/0', '0/1'],
        ...     'James': ['0/1', '0/1', '0/1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/1   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/1
        2  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/1
        3  chr1  103  .   T   C    .      .    .     GT    0/1  0/1   0/1

        We can select rows where all three samples have the variant:

        >>> vf.filter_sampall().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/1   0/1
        1  chr1  103  .   T   C    .      .    .     GT    0/1  0/1   0/1

        We can also remove those rows:

        >>> vf.filter_sampall(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/1
        1  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/1

        We can select rows where both Sara and James have the variant:

        >>> vf.filter_sampall(samples=['Sara', 'James']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/1   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/1
        2  chr1  103  .   T   C    .      .    .     GT    0/1  0/1   0/1

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_sampall(index=True)
        0     True
        1    False
        2    False
        3     True
        dtype: bool
        """
        if samples is None:
            samples = self.samples
        else:
            samples = [x if isinstance(x, str) else self.samples[x]
                       for x in samples]
        f = lambda r: all(r[samples].apply(gt_hasvar))
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampany(self, samples=None, opposite=False, index=False):
        """Select rows if any one of the given samples has the variant.

        The default behavior is to use all samples in the VcfFrame.

        Parameters
        ----------
        samples : list, optional
            List of sample names or indicies.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        filter_sampall
            Similar method that selects rows if all of the given
            samples have the variant.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'T'],
        ...     'ALT': ['A', 'C', 'A', 'C'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/0', '0/0', '0/1', '0/0'],
        ...     'Sara': ['0/0', '0/1', '0/0', '0/0'],
        ...     'James': ['0/1', '0/0', '0/0', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/0
        2  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/0
        3  chr1  103  .   T   C    .      .    .     GT    0/0  0/0   0/0

        We can select rows where at least one sample has the variant:

        >>> vf.filter_sampany().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/0
        2  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/0

        We can also remove those rows:

        >>> vf.filter_sampany(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  103  .   T   C    .      .    .     GT    0/0  0/0   0/0

        We can select rows where either Sara or James has the variant:

        >>> vf.filter_sampany(samples=['Sara', 'James']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/0  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/0  0/1   0/0

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_sampany(index=True)
        0     True
        1     True
        2     True
        3    False
        dtype: bool
        """
        if samples is None:
            samples = self.samples
        else:
            samples = [x if isinstance(x, str) else self.samples[x]
                       for x in samples]
        f = lambda r: any(r[samples].apply(gt_hasvar))
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampnum(self, threshold, opposite=False, index=False):
        """Select rows if the variant is prevalent enough.

        Parameters
        ----------
        threshold : int or float
            Minimum number or fraction of samples with the variant.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'T'],
        ...     'ALT': ['A', 'C', 'A'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '0/1', '0/1'],
        ...     'Sara': ['0/0', '0/1', '0/0'],
        ...     'James': ['0/1', '0/1', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/1  0/1   0/1
        2  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/0

        We can select rows where at least two samples have the variant:

        >>> vf.filter_sampnum(2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/1  0/1   0/1

        Similarly, we can select rows where at least 50% of the samples
        have the variant:

        >>> vf.filter_sampnum(0.5).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/0   0/1
        1  chr1  101  .   T   C    .      .    .     GT    0/1  0/1   0/1

        We can also remove those rows:

        >>> vf.filter_sampnum(0.5, opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara James
        0  chr1  102  .   T   A    .      .    .     GT    0/1  0/0   0/0

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_sampnum(2, index=True)
        0     True
        1     True
        2    False
        dtype: bool
        """
        def f(r):
            n = r[9:].apply(gt_hasvar).sum()
            if isinstance(threshold, int):
                return n >= threshold
            else:
                return n / self.shape[1] >= threshold
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        vf = self.__class__(self.copy_meta(), self.df[i])
        return vf

    def subtract(self, a, b):
        """Subtract the genotype data of Sample B from Sample A.

        This method is useful when, for example, you want to distinguish
        between somatic mutations and germline variants from an individual.

        Parameters
        ----------
        a : str or int
            Name or index of Sample A (e.g. somatic).
        b : str or int
            Name or index of Sample B (e.g. germline).

        Returns
        -------
        pandas.Series
            Resulting VCF column.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'A', 'C'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Somatic': ['./.', '0/1', '0/1', '0/0'],
        ...     'Germline': ['0/1', '0/1', './.', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Somatic Germline
        0  chr1  100  .   G   A    .      .    .     GT     ./.      0/1
        1  chr1  101  .   T   C    .      .    .     GT     0/1      0/1
        2  chr1  102  .   T   A    .      .    .     GT     0/1      ./.
        3  chr1  103  .   A   C    .      .    .     GT     0/0      0/1

        We subtract the two samples to get the true somatic mutations:

        >>> vf.df['TruelySomatic'] = vf.subtract('Somatic', 'Germline')
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Somatic Germline TruelySomatic
        0  chr1  100  .   G   A    .      .    .     GT     ./.      0/1           ./.
        1  chr1  101  .   T   C    .      .    .     GT     0/1      0/1           ./.
        2  chr1  102  .   T   A    .      .    .     GT     0/1      ./.           0/1
        3  chr1  103  .   A   C    .      .    .     GT     0/0      0/1           0/0
        """
        a = a if isinstance(a, str) else self.samples[a]
        b = b if isinstance(b, str) else self.samples[b]
        def func(r):
            m = row_missval(r)
            a_has = gt_hasvar(r[a])
            b_bas = gt_hasvar(r[b])
            if a_has and b_bas:
                return m
            elif a_has and not b_bas:
                return r[a]
            elif not a_has and b_bas:
                return r[a]
            else:
                return r[a]
        s = self.df.apply(func, axis=1)
        return s

    def sort(self):
        """Sort the VcfFrame by chromosome and position.

        Returns
        -------
        VcfFrame
            Sorted VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr10', 'chr2', 'chr1', 'chr2'],
        ...     'POS': [100, 101, 102, 90],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'A', 'T'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Steven': ['./.:.', '0/0:29', '0/0:28', '0/1:17']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
           CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven
        0  chr10  100  .   G   A    .      .    .  GT:DP   ./.:.
        1   chr2  101  .   T   C    .      .    .  GT:DP  0/0:29
        2   chr1  102  .   T   A    .      .    .  GT:DP  0/0:28
        3   chr2   90  .   A   T    .      .    .  GT:DP  0/1:17

        We can sort the VcfFrame by:

        >>> vf.sort().df
           CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven
        0   chr1  102  .   T   A    .      .    .  GT:DP  0/0:28
        1   chr2   90  .   A   T    .      .    .  GT:DP  0/1:17
        2   chr2  101  .   T   C    .      .    .  GT:DP  0/0:29
        3  chr10  100  .   G   A    .      .    .  GT:DP   ./.:.
        """
        df = self.df.sort_values(by=['CHROM', 'POS'], ignore_index=True,
            key=lambda col: [CONTIGS.index(x) if isinstance(x, str)
                             else x for x in col])
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def subset(self, samples, exclude=False):
        """Subset the VcfFrame for the selected samples.

        The order of the samples matters.

        Parameters
        ----------
        samples : list
            List of samples.
        exclude : bool, default: False
            If True, exclude the selected samples.

        Returns
        -------
        VcfFrame
            Subsetted VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP'],
        ...     'Steven': ['0/1:30', '0/1:29'],
        ...     'Sara': ['0/1:24', '0/1:30'],
        ...     'James': ['0/1:18', '0/1:24'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT  Steven    Sara   James
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/1:30  0/1:24  0/1:18
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  0/1:30  0/1:24

        We can subset the VcfFrame for James and Steven (in that order too):

        >>> vf.subset(['James', 'Steven']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT   James  Steven
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/1:18  0/1:30
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:24  0/1:29

        Alternatively, we can exclude James and Steven:

        >>> vf.subset(['James', 'Steven'], exclude=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    Sara
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:30
        """
        if exclude:
            samples = [x for x in self.samples if x not in samples]
        cols = self.df.columns[:9].to_list() + samples
        df = self.df[cols]
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def unphase(self):
        """Unphase all the sample genotypes.

        Returns
        -------
        VcfFrame
            Unphased VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['1|0', './.', '0|1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    1|0
        1  chr1  101  .   T   C    .      .    .     GT    ./.
        2  chr1  102  .   A   T    .      .    .     GT    0|1
        3  chr1  103  .   C   A    .      .    .     GT    0/1

        We can unphase the samples genotypes:

        >>> vf.unphase().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  101  .   T   C    .      .    .     GT    ./.
        2  chr1  102  .   A   T    .      .    .     GT    0/1
        3  chr1  103  .   C   A    .      .    .     GT    0/1
        """
        def func(r):
            r[9:] = r[9:].apply(gt_unphase)
            return r
        df = self.df.apply(func, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def update(self, other, headers=None, missing=True):
        """Copy data from the other VcfFrame.

        This method will copy and paste data from the other VcfFrame for
        overlapping records. By default, the following VCF headers are
        used: ID, QUAL, FILTER, and, INFO.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        headers : list, optional
            List of VCF headers to exclude.
        missing : bool, default: True
            If True, only fields with the missing value ('.') will be updated.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.
        """
        targets = ['ID', 'QUAL', 'FILTER', 'INFO']
        if headers is not None:
            for header in headers:
                targets.remove(header)
        def func(r1):
            r2 = other.df[(other.df['CHROM'] == r1['CHROM']) &
                          (other.df['POS'] == r1['POS']) &
                          (other.df['REF'] == r1['REF']) &
                          (other.df['ALT'] == r1['ALT'])]
            if r2.empty:
                return r1
            for target in targets:
                if missing:
                    if r1[target] == '.':
                        r1[target] = r2.iloc[0][target]
                    else:
                        pass
                else:
                    r1[target] = r2.iloc[0][target]
            return r1
        df = self.df.apply(func, axis=1)
        return self.__class__(self.copy_meta(), df)

    def slice(self, chrom, start=None, end=None):
        """Slice the VcfFrame for the given region.

        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int, optional
            Start position.
        end : int, optional
            End position.

        Returns
        -------
        VcfFrame
            Sliced VcfFrame.

        Examples
        --------

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr2'],
        ...     'POS': [100, 205, 297, 101],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '1/1', '0/1', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  205  .   T   C    .      .    .     GT    1/1
        2  chr1  297  .   A   T    .      .    .     GT    0/1
        3  chr2  101  .   C   A    .      .    .     GT    0/1
        >>> vf.slice('chr1').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  205  .   T   C    .      .    .     GT    1/1
        2  chr1  297  .   A   T    .      .    .     GT    0/1
        >>> vf.slice('chr1', start=101, end=300).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  205  .   T   C    .      .    .     GT    1/1
        1  chr1  297  .   A   T    .      .    .     GT    0/1
        >>> vf.slice('chr1', end=296).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  205  .   T   C    .      .    .     GT    1/1
        """
        i = (self.df.CHROM == chrom)
        if start:
            i = i & (start <= self.df.POS)
        if end:
            i = i & (self.df.POS <= end)
        df = self.df[i]
        return self.__class__(self.copy_meta(), df)
