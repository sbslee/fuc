"""
The pyvcf submodule is designed for working with VCF files. It implements
``pyvcf.VcfFrame`` which stores VCF data as ``pandas.DataFrame`` to allow
fast computation and easy manipulation. The ``pyvcf.VcfFrame`` class also
contains many useful plotting methods such as ``VcfFrame.plot_comparison``
and ``VcfFrame.plot_tmb``. The submodule strictly adheres to the
standard `VCF specification
<https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_.

A typical VCF file contains metadata lines (prefixed with '##'), a header
line (prefixed with '#'), and genotype lines that begin with a chromosome
identifier (e.g. 'chr1'). See the VCF specification above for an example
VCF file.

Genotype lines have nine required fields for storing variant information
and variable-length fields for storing sample genotype data. For some
fields, missing values are tolerated and can be specified with a dot ('.').
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
across the community. Currently, the pyvcf submodule is aware of the
following:

* AD - Total read depth for each allele (R, Integer)
* AF - Allele fraction of the event in the tumor (1, Float)
* DP - Read depth (1, Integer)

If sample annotation data are available for a given VCF file, use
the :class:`AnnFrame` class to import the data.
"""

import os
import re
import gzip
from copy import deepcopy

from . import pybed, common, pymaf

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf
from Bio import bgzf
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns

HEADERS = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
           'QUAL', 'FILTER', 'INFO', 'FORMAT']

CONTIGS = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M',
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM'
]

# Below are reserved genotype keys copied from Table 2 of the VCF
# specification: https://samtools.github.io/hts-specs/VCFv4.3.pdf

RESERVED_GENOTYPE_KEYS = {
    'AD':  {'number': 'R', 'type': int},   # Read depth for each allele
    'ADF': {'number': 'R', 'type': int},   # Read depth for each allele on the forward strand
    'ADR': {'number': 'R', 'type': int},   # Read depth for each allele on the reverse strand
    'DP':  {'number': 1,   'type': int},   # Read depth
    'EC':  {'number': 'A', 'type': int},   # Expected alternate allele counts
    'FT':  {'number': 1,   'type': str},   # Filter indicating if this genotype was “called”
    'GL':  {'number': 'G', 'type': float}, # Genotype likelihoods
    'GP':  {'number': 'G', 'type': float}, # Genotype posterior probabilities
    'GQ':  {'number': 1,   'type': int},   # Conditional genotype quality
    'GT':  {'number': 1,   'type': str},   # Genotype
    'HQ':  {'number': 2,   'type': int},   # Haplotype quality
    'MQ':  {'number': 1,   'type': int},   # RMS mapping quality
    'PL':  {'number': 'G', 'type': int},   # Phred-scaled genotype likelihoods rounded to the closest integer
    'PP':  {'number': 'G', 'type': int},   # Phred-scaled genotype posterior probabilities rounded to the closest integer
    'PQ':  {'number': 1,   'type': int},   # Phasing quality
    'PS':  {'number': 1,   'type': int},   # Phase set
}

CUSTOM_GENOTYPE_KEYS = {
    'AF':  {'number': 1,   'type': float},  # Allele fraction of the event in the tumor
}

def rescue_filtered_variants(vfs, format='GT'):
    """
    Rescue filtered variants if they are PASS in at least one of the input
    VCF files.

    Parameters
    ----------
    vfs : list
        List of VcfFrame objects.

    Returns
    -------
    VcfFrame
        VcfFrame object.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> data1 = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['G', 'T', 'C'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['PASS', 'weak_evidence', 'PASS'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'A': ['0/1', '0/1', '0/1']
    ... }
    >>> data2 = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['G', 'T', 'C'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['orientation', 'weak_evidence', 'PASS'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'B': ['0/1', '0/1', '0/1']
    ... }
    >>> data3 = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [102, 103, 104],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['C', 'T', 'A'],
    ...     'ALT': ['T', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['PASS', 'weak_evidence', 'PASS'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'C': ['0/1', '0/1', '0/1']
    ... }
    >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
    >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
    >>> vf3 = pyvcf.VcfFrame.from_dict([], data3)
    >>> vf1.df
      CHROM  POS ID REF ALT QUAL         FILTER INFO FORMAT    A
    0  chr1  100  .   G   A    .           PASS    .     GT  0/1
    1  chr1  101  .   T   C    .  weak_evidence    .     GT  0/1
    2  chr1  102  .   C   T    .           PASS    .     GT  0/1
    >>> vf2.df
      CHROM  POS ID REF ALT QUAL         FILTER INFO FORMAT    B
    0  chr1  100  .   G   A    .    orientation    .     GT  0/1
    1  chr1  101  .   T   C    .  weak_evidence    .     GT  0/1
    2  chr1  102  .   C   T    .           PASS    .     GT  0/1
    >>> vf3.df
      CHROM  POS ID REF ALT QUAL         FILTER INFO FORMAT    C
    0  chr1  102  .   C   T    .           PASS    .     GT  0/1
    1  chr1  103  .   T   C    .  weak_evidence    .     GT  0/1
    2  chr1  104  .   A   T    .           PASS    .     GT  0/1
    >>> rescued_vf = pyvcf.rescue_filtered_variants([vf1, vf2, vf3])
    >>> rescued_vf.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
    0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  ./.
    1  chr1  102  .   C   T    .      .    .     GT  0/1  0/1  0/1
    2  chr1  104  .   A   T    .      .    .     GT  ./.  ./.  0/1
    """
    # Check for duplicate samples.
    samples = []
    for vf in vfs:
        samples += vf.samples
    s = pd.Series(samples)
    duplicates = s[s.duplicated()].values
    if duplicates:
        raise ValueError(f'Duplicate samples found: {duplicates}.')

    dfs = []
    for vf in vfs:
        df = vf.df[vf.df.FILTER == 'PASS']
        dfs.append(df[['CHROM', 'POS', 'REF', 'ALT']])
    df = pd.concat(dfs).drop_duplicates()
    s = df.apply(lambda r: common.Variant(r.CHROM, r.POS, r.REF, r.ALT), axis=1)
    filtered_vfs = []
    for vf in vfs:
        i = vf.df.apply(lambda r: common.Variant(r.CHROM, r.POS, r.REF, r.ALT) in s.values, axis=1)
        filtered_vf = vf.copy()
        filtered_vf.df = vf.df[i]
        filtered_vfs.append(filtered_vf)
    merged_vf = merge(filtered_vfs, how='outer', format=format)
    return merged_vf

def gt_miss(g):
    """
    Return True if sample genotype is missing.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if sample genotype is missing.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.gt_miss('0')
    False
    >>> pyvcf.gt_miss('0/0')
    False
    >>> pyvcf.gt_miss('0/1')
    False
    >>> pyvcf.gt_miss('0|0:48:1:51,51')
    False
    >>> pyvcf.gt_miss('./.:.:.')
    True
    >>> pyvcf.gt_miss('.:.')
    True
    >>> pyvcf.gt_miss('.')
    True
    >>> pyvcf.gt_miss('./.:13,3:16:41:41,0,402')
    True
    """
    return '.' in g.split(':')[0]

def gt_polyp(g):
    """Return True if sample genotype has a polyploid call.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if sample genotype has a polyploid call.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.gt_polyp('1')
    False
    >>> pyvcf.gt_polyp('0/1')
    False
    >>> pyvcf.gt_polyp('0/1/1')
    True
    >>> pyvcf.gt_polyp('1|0|1')
    True
    >>> pyvcf.gt_polyp('0/./1/1')
    True
    """
    gt = g.split(':')[0]
    if '/' in gt:
        return gt.count('/') > 1
    else:
        return gt.count('|') > 1

def gt_hasvar(g):
    """
    Return True if sample genotype has at least one variant call.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if sample genotype has a variant call.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.gt_hasvar('0')
    False
    >>> pyvcf.gt_hasvar('0/0')
    False
    >>> pyvcf.gt_hasvar('./.')
    False
    >>> pyvcf.gt_hasvar('1')
    True
    >>> pyvcf.gt_hasvar('0/1')
    True
    >>> pyvcf.gt_hasvar('1/2')
    True
    >>> pyvcf.gt_hasvar('1|0')
    True
    >>> pyvcf.gt_hasvar('1|2:21:6:23,27')
    True
    """
    if g.split(':')[0].replace('/', '').replace(
        '|', '').replace('.', '').replace('0', ''):
        return True
    else:
        return False

def gt_unphase(g):
    """
    Return unphased sample genotype.

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

    >>> from fuc import pyvcf
    >>> pyvcf.gt_unphase('1')
    '1'
    >>> pyvcf.gt_unphase('0/0')
    '0/0'
    >>> pyvcf.gt_unphase('0/1')
    '0/1'
    >>> pyvcf.gt_unphase('0/1:35:4')
    '0/1:35:4'
    >>> pyvcf.gt_unphase('0|1')
    '0/1'
    >>> pyvcf.gt_unphase('1|0')
    '0/1'
    >>> pyvcf.gt_unphase('2|1:2:0:18,2')
    '1/2:2:0:18,2'
    >>> pyvcf.gt_unphase('.')
    '.'
    >>> pyvcf.gt_unphase('./.')
    './.'
    >>> pyvcf.gt_unphase('.|.')
    './.'
    """
    l = g.split(':')
    gt = l[0]
    if '|' not in gt:
        return g
    if '.' in gt:
        return g.replace('|', '/')
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

    >>> from fuc import pyvcf
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

    >>> from fuc import pyvcf
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

    >>> from fuc import pyvcf
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

def row_updateinfo(r, key, value):
    """Update INFO data in the row that match the given key.

    Parameters
    ----------
    r : pandas.Series
        VCF row.
    key : str
        INFO key.
    value : str
        New value to be assigned.

    Returns
    -------
    str
        New INFO field.

    Examples
    --------

    >>> from fuc import pyvcf
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
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER        INFO FORMAT Steven
    0  chr1  100  .   G   A    .      .     DB;AC=0     GT    0/0
    1  chr1  101  .   T   C    .      .  DB;H2;AC=1     GT    0/1
    2  chr1  102  .   A   T    .      .  DB;H2;AC=1     GT    0/1
    3  chr1  103  .   C   A    .      .           .     GT    0/0
    >>> vf.df.apply(pyvcf.row_updateinfo, args=('AC', '4'), axis=1)
    0       DB;AC=4
    1    DB;H2;AC=4
    2    DB;H2;AC=4
    3             .
    dtype: object
    """
    fields = r.INFO.split(';')
    for i, field in enumerate(fields):
        if field.startswith(f'{key}='):
            fields[i] = field[:len(key)+1] + value
            break
    return ';'.join(fields)

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

    >>> from fuc import pyvcf
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

def simulate_genotype(
    p=0.5, noise_scale=0.05, dp_show=True, dp_loc=30, dp_scale=10,
    ad_show=True, ad_loc=0.5, ad_scale=0.05, af_show=True
):
    lower, upper = 0, 1
    mu, sigma = p, noise_scale
    X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    has_var = np.random.binomial(1, X.rvs(1)) == 1

    if has_var:
        result = '0/1'
    else:
        result = '0/0'

    dp = np.random.normal(loc=dp_loc, scale=dp_scale)
    dp = round(abs(dp))

    if has_var:
        alt = round(abs(np.random.normal(loc=ad_loc, scale=ad_scale)) * dp)
        ref = dp - alt
        ad = f'{ref},{alt}'
    else:
        ad = f'{dp},0'

    if has_var:
        af = round(alt / (ref + alt), 3)
    else:
        af = 0

    if dp_show:
        result += f':{dp}'

    if ad_show:
        result += f':{ad}'

    if af_show:
        result += f':{af}'

    return result

def simulate_sample(
    n, p=0.5, noise_scale=0.1, dp_show=True, dp_loc=30, dp_scale=10,
    ad_show=True, ad_loc=0.5, ad_scale=0.05, af_show=True
):
    l = []
    for i in range(n):
        genotype = simulate_genotype(
            p=p, noise_scale=noise_scale, dp_show=dp_show,
            dp_loc=dp_loc, dp_scale=dp_scale, ad_show=ad_show,
            ad_loc=ad_loc, ad_scale=ad_scale
        )
        l.append(genotype)
    return l

class VcfFrame:
    """
    Class for storing VCF data.

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
    VcfFrame.from_file
        Construct VcfFrame from a VCF file.

    Examples
    --------
    Constructing VcfFrame from pandas DataFrame:

    >>> from fuc import pyvcf
    >>> import pandas as pd
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
        """list : List of the sample names."""
        return self.df.columns[9:].to_list()

    @property
    def shape(self):
        """tuple : Dimensionality of VcfFrame (variants, samples)."""
        return (self.df.shape[0], len(self.samples))

    def add_af(self, decimals=3):
        """Compute AF using AD and add it to the FORMAT field.

        Parameters
        ----------
        decimals : int, default: 3
            Number of decimals to display.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        We can add the AF subfield to our genotype data:

        >>> vf.add_af().df
          CHROM  POS ID REF  ... INFO    FORMAT           Steven               Sara
        0  chr1  100  .   A  ...    .  GT:AD:AF  0/1:12,15:0.556    0/1:13,17:0.567
        1  chr1  100  .   A  ...    .  GT:AD:AF   0/0:32,1:0.030    0/1:14,15:0.517
        2  chr2  200  .   C  ...    .  GT:AD:AF  0/1:16,12:0.429            ./.:.:.
        3  chr2  200  .   C  ...    .  GT:AD:AF          ./.:.:.  1/2:0,11,17:0.607
        """
        def one_row(r):
            i = r.FORMAT.split(':').index('AD')
            def one_gt(g):
                ad = g.split(':')[i]
                if ad == '.':
                    return f'{g}:.'
                depths = [int(x) for x in ad.split(',')]
                total = sum(depths)
                af = max(depths[1:]) / total
                return f'{g}:{af:.{decimals}f}'
            r.iloc[9:] = r.iloc[9:].apply(one_gt)
            r.FORMAT += ':AF'
            return r
        df = self.df.apply(one_row, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def add_dp(self):
        """Compute DP using AD and add it to the FORMAT field.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> from fuc import pyvcf
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

    def cfilter_empty(self, opposite=False, as_list=False):
        """Remove samples whose genotype calls are all missing.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return samples that don't meet the said criteria.
        as_list : bool, default: False
             If True, return a list of sample names instead of a VcfFrame.

        Returns
        -------
        VcfFrame
            Filtered VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'G', 'T'],
        ...     'ALT': ['A', 'C', 'C', 'C'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '1/1', '1/1', '1/1'],
        ...     'Rachel': ['./.', './.', './.', './.'],
        ...     'John': ['0/0', './.', '0/0', '0/0'],
        ...     'Sara': ['./.', './.', './.', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Rachel John Sara
        0  chr1  100  .   G   A    .      .    .     GT    0/1    ./.  0/0  ./.
        1  chr1  101  .   T   C    .      .    .     GT    1/1    ./.  ./.  ./.
        2  chr1  102  .   G   C    .      .    .     GT    1/1    ./.  0/0  ./.
        3  chr1  103  .   T   C    .      .    .     GT    1/1    ./.  0/0  ./.

        We can remove samples whose genotypes are all missing:

        >>> vf.cfilter_empty().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven John
        0  chr1  100  .   G   A    .      .    .     GT    0/1  0/0
        1  chr1  101  .   T   C    .      .    .     GT    1/1  ./.
        2  chr1  102  .   G   C    .      .    .     GT    1/1  0/0
        3  chr1  103  .   T   C    .      .    .     GT    1/1  0/0

        We can also select those samples:

        >>> vf.cfilter_empty(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Rachel Sara
        0  chr1  100  .   G   A    .      .    .     GT    ./.  ./.
        1  chr1  101  .   T   C    .      .    .     GT    ./.  ./.
        2  chr1  102  .   G   C    .      .    .     GT    ./.  ./.
        3  chr1  103  .   T   C    .      .    .     GT    ./.  ./.

        Finally, we can return a list of sample names from the filtering:

        >>> vf.cfilter_empty(as_list=True)
        ['Steven', 'John']
        """
        f = lambda r: r[9:].apply(gt_miss)
        s = self.df.apply(f, axis=1).all()
        if opposite:
            s = s[s == True]
        else:
            s = s[s == False]
        l = s.index.to_list()
        if as_list:
            return l
        return self.subset(l)

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

        >>> from fuc import pyvcf
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
                if gt_miss(x):
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

    @classmethod
    def from_dict(cls, meta, data):
        """
        Construct VcfFrame from dict of array-like or dicts.

        Parameters
        ----------
        meta : list
            List of the metadata lines.
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        VcfFrame
            VcfFrame.

        See Also
        --------
        VcfFrame
            VcfFrame object creation using constructor.
        VcfFrame.from_file
            Construct VcfFrame from a VCF file.

        Examples
        --------
        Below is a simple example:

        >>> from fuc import pyvcf
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
        ...     'A': ['0/1', '1/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0/1
        1  chr2  101  .   T   C    .      .    .     GT  1/1
        """
        return cls(meta, pd.DataFrame(data))

    @classmethod
    def from_file(cls, fn, compression=False):
        """Construct VcfFrame from a VCF file.

        If the file name ends with '.gz', the method will automatically
        use the BGZF decompression when reading the file.

        Parameters
        ----------
        fn : str
            VCF file path (zipped or unzipped).
        compression : bool, default: False
            If True, use the BGZF decompression regardless of file name.

        Returns
        -------
        VcfFrame
            VcfFrame.

        See Also
        --------
        VcfFrame
            VcfFrame object creation using constructor.
        VcfFrame.from_dict
            Construct VcfFrame from dict of array-like or dicts.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> vf = pyvcf.VcfFrame.from_file('unzipped.vcf')
        >>> vf = pyvcf.VcfFrame.from_file('zipped.vcf.gz')
        >>> vf = pyvcf.VcfFrame.from_file('zipped.vcf', compression=True)
        """
        meta = []
        skip_rows = 0
        if fn.startswith('~'):
            fn = os.path.expanduser(fn)
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
        f.close()
        return cls(meta, df)

    def compare(self, a, b, c=None, mode='all'):
        """
        Compare genotype data between two (A, B) or three (A, B, C) samples.

        This method will return (Ab, aB, AB, ab) for comparison between two
        samples and (Abc, aBc, ABc, abC, AbC, aBC, ABC, abc) for three
        samples. Note that the former is equivalent to (FP, FN, TP, TN) if
        we assume A is the test sample and B is the truth sample.

        Only biallelic sites will be used for comparison. Additionally, the
        method will only consider presence or absence of variant calls (i.e.
        zygosity is ignored).

        Parameters
        ----------
        a, b : str or int
            Name or index of Sample A or B (test or truth).
        c : str or int, optional
            Name or index of Sample C.
        mode : {'all', 'snv', 'indel'}, default: 'all'
            Determines which variant types should be analyzed:

            - 'all': Include both SNVs and INDELs.
            - 'snv': Include SNVs only.
            - 'indel': Include INDELs only.

        Returns
        -------
        tuple
            Four- or eight-element tuple depending on the number of samples.

        See Also
        --------
        fuc.api.common.sumstat
            Return various summary statistics from (FP, FN, TP, TN).

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'CT', 'T', 'C', 'A'],
        ...     'ALT': ['A', 'C', 'A', 'T', 'G,C'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', '0/0', '0/0', '0/1', '0/0'],
        ...     'B': ['1/1', '0/1', './.', '0/1', '0/0'],
        ...     'C': ['0/1', '0/1', '1/1', './.', '1/2'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  100  .   G    A    .      .    .     GT  0/1  1/1  0/1
        1  chr1  101  .  CT    C    .      .    .     GT  0/0  0/1  0/1
        2  chr1  102  .   T    A    .      .    .     GT  0/0  ./.  1/1
        3  chr1  103  .   C    T    .      .    .     GT  0/1  0/1  ./.
        4  chr1  104  .   A  G,C    .      .    .     GT  0/0  0/0  1/2

        We can first compare the samples A and B:

        >>> vf.compare('A', 'B', mode='all')
        (0, 1, 2, 1)
        >>> vf.compare('A', 'B', mode='snv')
        (0, 0, 2, 1)
        >>> vf.compare('A', 'B', mode='indel')
        (0, 1, 0, 0)

        We can also compare all three samples at once:

        >>> vf.compare('A', 'B', 'C')
        (0, 0, 1, 1, 0, 1, 1, 0)
        """
        vf = self.filter_multialt()

        if mode == 'all':
            pass
        elif mode == 'snv':
            vf = vf.filter_indel()
        elif mode == 'indel':
            vf = vf.filter_indel(opposite=True)
        else:
            raise ValueError(f'Incorrect mode: {mode}.')

        if c is None:
            result = self._compare_two(vf, a, b)
        else:
            result = self._compare_three(vf, a, b, c)

        return result

    def _compare_two(self, vf, a, b):
        a = a if isinstance(a, str) else vf.samples[a]
        b = b if isinstance(b, str) else vf.samples[b]
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
        d = vf.df.apply(func, axis=1).value_counts().to_dict()
        Ab = d['Ab'] if 'Ab' in d else 0
        aB = d['aB'] if 'aB' in d else 0
        AB = d['AB'] if 'AB' in d else 0
        ab = d['ab'] if 'ab' in d else 0
        return (Ab, aB, AB, ab)

    def _compare_three(self, vf, a, b, c):
        a = a if isinstance(a, str) else vf.samples[a]
        b = b if isinstance(b, str) else vf.samples[b]
        c = c if isinstance(c, str) else vf.samples[c]
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
        d = vf.df.apply(func, axis=1).value_counts().to_dict()
        Abc = d['Abc'] if 'Abc' in d else 0
        aBc = d['aBc'] if 'aBc' in d else 0
        ABc = d['ABc'] if 'ABc' in d else 0
        abC = d['abC'] if 'abC' in d else 0
        AbC = d['AbC'] if 'AbC' in d else 0
        aBC = d['aBC'] if 'aBC' in d else 0
        ABC = d['ABC'] if 'ABC' in d else 0
        abc = d['abc'] if 'abc' in d else 0
        return (Abc, aBc, ABc, abC, AbC, aBC, ABC, abc)

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

        >>> from fuc import pyvcf
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

        >>> from fuc import pyvcf
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
        """
        Render the VcfFrame to a console-friendly tabular output.

        Returns
        -------
        str
            String representation of the VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
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
        """
        Merge with the other VcfFrame.

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

        >>> from fuc import pyvcf
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
        ...     'A': ['0/0:32', '0/1:29'],
        ...     'B': ['0/1:24', '1/1:30'],
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
        ...     'C': ['./.:.', '0/0:24', '0/0:26'],
        ...     'D': ['0/1:24', '0/1:31', '0/1:26'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> vf1.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT       A       B
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/0:32  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/1:29  1/1:30
        >>> vf2.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT       C       D
        0  chr1  100  .   G   A    .      .    .  GT:DP   ./.:.  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP  0/0:24  0/1:31
        2  chr2  200  .   A   T    .      .    .  GT:DP  0/0:26  0/1:26

        We can merge the two VcfFrames with ``how='inner'`` (default):

        >>> vf1.merge(vf2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/0  0/1  ./.  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1  1/1  0/0  0/1

        We can also merge with ``how='outer'``:

        >>> vf1.merge(vf2, how='outer').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/0  0/1  ./.  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1  1/1  0/0  0/1
        2  chr2  200  .   A   T    .      .    .     GT  ./.  ./.  0/0  0/1

        Since both VcfFrames have the DP subfield, we can use ``format='GT:DP'``:

        >>> vf1.merge(vf2, how='outer', format='GT:DP').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT       A       B       C       D
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

    def miss2ref(self):
        """
        Convert missing genotype (./.) to homozygous REF (0/0).

        Returns
        -------
        VcfFrame
            VcfFrame object.

        Examples
        --------

        >>> from fuc import pyvcf
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
        ...     'A': ['./.', '1/1'],
        ...     'B': ['./.', './.']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  ./.  ./.
        1  chr2  101  .   T   C    .      .    .     GT  1/1  ./.
        >>> new_vf = vf.miss2ref()
        >>> new_vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/0  0/0
        1  chr2  101  .   T   C    .      .    .     GT  1/1  0/0
        """
        df = self.copy_df()
        def one_gt(g):
            l = [g.split(':')[0].replace('.', '0')] + g.split(':')[1:]
            return ':'.join(l)
        df.iloc[:, 9:] = df.iloc[:, 9:].applymap(one_gt)
        return self.__class__(self.copy_meta(), df)

    def plot_comparison(
        self, a, b, c=None, labels=None, ax=None, figsize=None
    ):
        """Create a Venn diagram showing genotype concordance between groups.

        This method supports comparison between two groups (Groups A & B)
        as well as three groups (Groups A, B, & C).

        Parameters
        ----------
        a, b : list
            Sample names. The lists must have the same shape.
        c : list, optional
            Same as above.
        labels : list, optional
            List of labels to be displayed.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        matplotlib_venn._common.VennDiagram
            VennDiagram object.

        Examples
        --------

        .. plot::
            :context: close-figs

            >>> from fuc import pyvcf, common
            >>> common.load_dataset('pyvcf')
            >>> f = '~/fuc-data/pyvcf/plot_comparison.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(f)
            >>> a = ['Steven_A', 'John_A', 'Sara_A']
            >>> b = ['Steven_B', 'John_B', 'Sara_B']
            >>> c = ['Steven_C', 'John_C', 'Sara_C']
            >>> vf.plot_comparison(a, b)

        .. plot::
            :context: close-figs

            >>> vf.plot_comparison(a, b, c)
        """
        if len(a) != len(b):
            raise ValueError('Groups A and B have different length.')
        if c is not None and len(a) != len(c):
            raise ValueError('Group C has unmatched length.')
        if labels is None:
            if c is None:
                labels = ('A', 'B')
            else:
                labels = ('A', 'B', 'C')

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        venn_kws = dict(ax=ax, alpha=0.5, set_labels=labels)
        if c is None:
            out = self._plot_comparison_two(a, b, venn_kws)
        else:
            out = self._plot_comparison_three(a, b, c, venn_kws)
        return ax, out

    def _plot_comparison_two(self, a, b, venn_kws):
        n = [0, 0, 0, 0]
        for i in range(len(a)):
            n = [x + y for x, y in zip(n, self.compare(a[i], b[i]))]
        out = venn2(subsets=n[:-1], **venn_kws)
        return out

    def _plot_comparison_three(self, a, b, c, venn_kws):
        n = [0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(len(a)):
            n = [x + y for x, y in zip(n, self.compare(a[i], b[i], c[i]))]
        out = venn3(subsets=n[:-1], **venn_kws)
        return out

    def plot_hist(
        self, k, af=None, group_col=None, group_order=None, kde=True,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a histogram showing AD/AF/DP distribution.

        Parameters
        ----------
        k : {'AD', 'AF', 'DP'}
            Genotype key.
        af : common.AnnFrame
            AnnFrame containing sample annotation data.
        group_col : list, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        kde : bool, default: True
            Compute a kernel density estimate to smooth the distribution.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.histplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> from fuc import common, pyvcf
            >>> common.load_dataset('pyvcf')
            >>> vcf_file = '~/fuc-data/pyvcf/normal-tumor.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_hist('DP')

        We can draw multiple histograms with hue mapping:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/pyvcf/normal-tumor-annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col='Sample')
            >>> vf.plot_hist('DP', af=af, group_col='Tissue')

        We can show AF instead of DP:

        .. plot::
            :context: close-figs

            >>> vf.plot_hist('AF')
        """
        d = {
            'AD': lambda x: float(x.split(',')[1]),
            'AF': lambda x: float(x),
            'DP': lambda x: int(x),
        }
        df = self.extract(k, as_nan=True, func=d[k])
        df = df.T
        id_vars = ['index']
        if group_col is not None:
            df = pd.concat([df, af.df[group_col]], axis=1, join='inner')
            id_vars.append(group_col)
        df = df.reset_index()
        df = pd.melt(df, id_vars=id_vars)
        df = df.dropna()
        df = df.rename(columns={'value': k})

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.histplot(
            data=df, x=k, hue=group_col, hue_order=group_order, kde=kde,
            ax=ax, **kwargs
        )

        return ax

    def plot_tmb(
        self, af=None, group_col=None, group_order=None, kde=True, ax=None, figsize=None, **kwargs
    ):
        """
        Create a histogram showing TMB distribution.

        Parameters
        ----------
        af : common.AnnFrame
            AnnFrame containing sample annotation data (requires ``hue``).
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        kde : bool, default: True
            Compute a kernel density estimate to smooth the distribution.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.histplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> from fuc import common, pyvcf
            >>> common.load_dataset('pyvcf')
            >>> vcf_file = '~/fuc-data/pyvcf/normal-tumor.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_tmb()

        We can draw multiple histograms with hue mapping:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/pyvcf/normal-tumor-annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col='Sample')
            >>> vf.plot_tmb(af=af, group_col='Tissue')
        """
        s = self.df.iloc[:, 9:].applymap(gt_hasvar).sum()
        s.name = 'TMB'
        if af is None:
            df = s.to_frame()
        else:
            df = pd.concat([af.df, s], axis=1, join='inner')

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.histplot(
            data=df, x='TMB', ax=ax, hue=group_col, hue_order=group_order,
            kde=kde, **kwargs
        )

        return ax

    def plot_regplot(self, a, b, ax=None, figsize=None, **kwargs):
        """
        Create a scatter plot with a linear regression model fit visualizing
        correlation between TMB in two sample groups.

        The method will automatically calculate and print summary statistics
        including R-squared and p-value.

        Parameters
        ----------
        a, b : array-like
            Lists of sample names. The lists must have the same shape.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.regplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> from fuc import common, pyvcf
            >>> common.load_dataset('pyvcf')
            >>> vcf_file = '~/fuc-data/pyvcf/normal-tumor.vcf'
            >>> annot_file = '~/fuc-data/pyvcf/normal-tumor-annot.tsv'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> af = common.AnnFrame.from_file(annot_file, sample_col='Sample')
            >>> normal = af.df[af.df.Tissue == 'Normal'].index
            >>> normal.name = 'Normal'
            >>> tumor = af.df[af.df.Tissue == 'Tumor'].index
            >>> tumor.name = 'Tumor'
            >>> vf.plot_regplot(normal, tumor)
            Results for B ~ A:
            R^2 = 0.01
            P = 7.17e-01
            >>> plt.tight_layout()
        """
        s = self.df.iloc[:, 9:].applymap(gt_hasvar).sum()
        df = pd.concat([s[a].reset_index(), s[b].reset_index()], axis=1)
        df.columns = ['A_Label', 'A_TMB', 'B_Label', 'B_TMB']

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.regplot(x='A_TMB', y='B_TMB', data=df, ax=ax, **kwargs)

        try:
            ax.set_xlabel(a.name)
        except AttributeError:
            ax.set_xlabel('A')

        try:
            ax.set_ylabel(b.name)
        except AttributeError:
            ax.set_ylabel('B')

        # Print summary statistics including R-squared and p-value.
        results = smf.ols(f'B_TMB ~ A_TMB', data=df).fit()
        print(f'Results for B ~ A:')
        print(f'R^2 = {results.rsquared:.2f}')
        print(f'  P = {results.f_pvalue:.2e}')

        return ax

    def markmiss(
        self, expr, greedy=False, opposite=False, samples=None, as_nan=False
    ):
        """
        Mark all genotypes that satisfy the query expression as missing.

        Parameters
        ----------
        expr : str
            The expression to evaluate. See the examples below for details.
        greedy : bool, default: False
            If True, mark even ambiguous genotypes as missing.
        opposite : bool, default: False
            If True, mark all genotypes that do not satisfy the query
            expression as missing and leave those that do intact.
        sampels : list, optional
            If provided, apply the marking only to these samples.
        as_nan : bool, default: False
            If True, mark genotypes as ``NaN`` instead of as missing.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'T'],
        ...     'ALT': ['A', 'C', 'G'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT:DP:AD', 'GT:DP:AD', 'GT:DP:AD'],
        ...     'A': ['0/0:26:0,26', '0/1:32:16,16', '0/0:.:.'],
        ...     'B': ['./.:.:.', '0/0:31:29,2', './.:.:.'],
        ...     'C': ['0/1:18:12,6', '0/0:24:24,0', '1/1:8:0,8'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26      ./.:.:.  0/1:18:12,6
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       0/0:.:.      ./.:.:.    1/1:8:0,8

        To mark as missing all genotypes with ``0/0``:

        >>> vf.markmiss('GT == "0/0"').df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A        B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD       ./.:.:.  ./.:.:.  0/1:18:12,6
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  ./.:.:.      ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.  ./.:.:.    1/1:8:0,8

        To mark as missing all genotypes that do not have ``0/0``:

        >>> vf.markmiss('GT != "0/0"').df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT            A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD  0/0:26:0,26      ./.:.:.      ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD      ./.:.:.  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD      0/0:.:.      ./.:.:.      ./.:.:.

        To mark as missing all genotypes whose ``DP`` is below 30:

        >>> vf.markmiss('DP < 30').df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B        C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.  ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       0/0:.:.      ./.:.:.  ./.:.:.

        Note that the genotype ``0/0:.:.`` was not marked as missing because
        its ``DP`` is missing and therefore it could not be evaluated
        properly. To mark even ambiguous genotypes like this one as missing,
        you can set ``greedy`` as True:

        >>> vf.markmiss('DP < 30', greedy=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B        C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.  ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.  ./.:.:.

        To mark as missing all genotypes whose ALT allele has read depth
        below 10:

        >>> vf.markmiss('AD[1] < 10', greedy=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A        B        C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26  ./.:.:.  ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  ./.:.:.  ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.  ./.:.:.  ./.:.:.

        To mark as missing all genotypes whose ALT allele has read depth
        below 10 and ``DP`` is below 30:

        >>> vf.markmiss('AD[1] < 10 and DP < 30', greedy=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B        C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26      ./.:.:.  ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.  ./.:.:.

        To mark as missing all genotypes whose ALT allele has read depth
        below 10 or ``DP`` is below 30:

        >>> vf.markmiss('AD[1] < 10 or DP < 30', greedy=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A        B        C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD       ./.:.:.  ./.:.:.  ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  ./.:.:.  ./.:.:.
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.  ./.:.:.  ./.:.:.

        To only retain genotypes whose ALT allele has read depth below 10 or
        ``DP`` is below 30:

        >>> vf.markmiss('AD[1] < 10 or DP < 30', opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT            A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD  0/0:26:0,26      ./.:.:.  0/1:18:12,6
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD      ./.:.:.  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD      ./.:.:.      ./.:.:.    1/1:8:0,8

        To mark as missing all genotypes whose mean of ``AD`` is below 10:

        >>> vf.markmiss('np.mean(AD) < 10', greedy=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26      ./.:.:.      ./.:.:.
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.      ./.:.:.

        To do the same as above, but only for the samples A and B:

        >>> vf.markmiss('np.mean(AD) < 10', greedy=True, samples=['A', 'B']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26      ./.:.:.  0/1:18:12,6
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       ./.:.:.      ./.:.:.    1/1:8:0,8

        To mark as ``NaN`` all genotypes whose sum of ``AD`` is below 10:

        >>> vf.markmiss('sum(AD) < 10', as_nan=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B            C
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD   0/0:26:0,26      ./.:.:.  0/1:18:12,6
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD  0/1:32:16,16  0/0:31:29,2  0/0:24:24,0
        2  chr1  102  .   T   G    .      .    .  GT:DP:AD       0/0:.:.      ./.:.:.          NaN

        Marking as ``NaN`` is useful when, for example, it is necessary to
        count how many genotypes are marked:

        >>> vf.markmiss('sum(AD) < 10', as_nan=True).df.isna().sum().sum()
        1
        """
        genotype_keys = {**RESERVED_GENOTYPE_KEYS, **CUSTOM_GENOTYPE_KEYS}
        types = {}
        for k, v in genotype_keys.items():
            if v['number'] == 1:
                types[k] = v['type']
            else:
                types[k] = lambda x: [v['type'](x) for x in x.split(',')]
        # Extract unique genotype keys from the expression.
        target_keys = re.findall('[a-z]+', expr, flags=re.IGNORECASE)
        target_keys = list(set(target_keys))
        target_keys = [x for x in target_keys if x in genotype_keys]
        # Define the marking algorithm for each row.
        def one_row(r):
            row_keys = r.FORMAT.split(':')
            # Infer the most appropriate missing value.
            if as_nan:
                missing = np.nan
            else:
                missing = row_missval(r)
            # Define the marking algorithm for each genotype.
            def one_gt(g):
                if opposite:
                    matched, unmatched = g, missing
                else:
                    matched, unmatched = missing, g
                subfields = g.split(':')
                for target_key in target_keys:
                    try:
                        i = row_keys.index(target_key)
                        x = subfields[i]
                        ambiguous = not x.replace('.', '').replace('/', ''
                            ).replace('|', '')
                    except ValueError:
                        ambiguous = True
                    if ambiguous:
                        if greedy:
                            return matched
                        else:
                            return unmatched
                    locals()[target_key] = types[target_key](x)
                if eval(expr):
                    return matched
                else:
                    return unmatched
            # Apply the marking to each genotype.
            if samples is None:
                r[9:] = r[9:].apply(one_gt)
            else:
                r[samples] = r[samples].apply(one_gt)
            return r
        # Apply the marking to each row.
        df = self.df.apply(one_row, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def expand(self):
        """Expand each multiallelic locus to multiple rows.

        Only the GT subfield of FORMAT will be retained.

        Returns
        -------
        VcfFrame
            Expanded VcfFrame.

        See Also
        --------
        VcfFrame.collapse
            Collapse duplicate records in the VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C', 'T,G', 'G', 'A,G,CT'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Steven': ['0/1:32', './.:.', '0/1:27', '0/2:34'],
        ...     'Sara': ['0/0:28', '1/2:30', '1/1:29', '1/2:38'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF     ALT QUAL FILTER INFO FORMAT  Steven    Sara
        0  chr1  100  .   A       C    .      .    .  GT:DP  0/1:32  0/0:28
        1  chr1  101  .   A     T,G    .      .    .  GT:DP   ./.:.  1/2:30
        2  chr1  102  .   C       G    .      .    .  GT:DP  0/1:27  1/1:29
        3  chr1  103  .   C  A,G,CT    .      .    .  GT:DP  0/2:34  1/2:38

        We can expand each of the multiallelic loci:

        >>> vf.expand().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven Sara
        0  chr1  100  .   A   C    .      .    .     GT    0/1  0/0
        1  chr1  101  .   A   T    .      .    .     GT    ./.  0/1
        2  chr1  101  .   A   G    .      .    .     GT    ./.  0/1
        3  chr1  102  .   C   G    .      .    .     GT    0/1  1/1
        4  chr1  103  .   C   A    .      .    .     GT    0/0  0/1
        5  chr1  103  .   C   G    .      .    .     GT    0/1  0/1
        6  chr1  103  .   C  CT    .      .    .     GT    0/0  0/0
        """
        data = []
        def one_gt(g, i):
            if gt_miss(g):
                return g
            l = g.split('/')
            l = ['1' if x == str(i+1) else '0' for x in l]
            l = sorted(l)
            return '/'.join(l)
        for i, r in self.df.iterrows():
            r.FORMAT = 'GT'
            r[9:] = r[9:].apply(lambda x: x.split(':')[0])
            alt_alleles = r.ALT.split(',')
            if len(alt_alleles) == 1:
                data.append(r)
                continue
            for i, alt_allele in enumerate(alt_alleles):
                s = r.copy()
                s.ALT = alt_allele
                s[9:] = s[9:].apply(one_gt, args=(i,))
                data.append(s)
        return self.__class__(self.copy_meta(), pd.concat(data, axis=1).T)

    def filter_bed(self, bed, opposite=False, as_index=False):
        """
        Select rows that overlap with the given BED data.

        Parameters
        ----------
        bed : pybed.BedFrame or str
            BedFrame or path to a BED file.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pybed, pyvcf
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

        >>> vf.filter_bed(bf, as_index=True)
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
            bf = pybed.BedFrame.from_file(bed)
        f = lambda r: not bf.gr[r.CHROM, r.POS:r.POS+1].empty
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_empty(self, opposite=False, as_index=False):
        """Remove rows with no genotype calls at all.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_empty(as_index=True)
        0     True
        1    False
        2     True
        3    False
        dtype: bool
        """
        f = lambda r: not all(r.iloc[9:].apply(gt_miss))
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_indel(self, opposite=False, as_index=False):
        """Remove rows with an indel.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_indel(as_index=True)
        0     True
        1    False
        2    False
        3     True
        dtype: bool
        """
        i = ~self.df.apply(row_hasindel, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagall(self, flags, opposite=False, as_index=False):
        """Select rows if all of the given INFO flags are present.

        Parameters
        ----------
        flags : list
            List of INFO flags.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        VcfFrame.filter_flagany
            Similar method that selects rows if any one of the given
            INFO flags is present.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_flagall(['H2', 'DB'], as_index=True)
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
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagany(self, flags, opposite=False, as_index=False):
        """Select rows if any one of the given INFO flags is present.

        Parameters
        ----------
        flags : list
            List of INFO flags.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        VcfFrame.filter_flagall
            Similar method that selects rows if all of the given INFO
            flags are present.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_flagany(['H2'], as_index=True)
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
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_multialt(self, opposite=False, as_index=False):
        """
        Remove rows with multiple ALT alleles.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------

        Assume we have the following data:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C,T', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/2', '0/0', '0/1', './.'],
        ...     'B': ['0/1', '0/1', './.', '1/2'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   A  C,T    .      .    .     GT  0/2  0/1
        1  chr1  101  .   A    T    .      .    .     GT  0/0  0/1
        2  chr1  102  .   C    G    .      .    .     GT  0/1  ./.
        3  chr1  103  .   C  G,A    .      .    .     GT  ./.  1/2

        We can remove rows with multiple ALT alleles:

        >>> vf.filter_multialt().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  101  .   A   T    .      .    .     GT  0/0  0/1
        1  chr1  102  .   C   G    .      .    .     GT  0/1  ./.

        We can also select those rows:

        >>> vf.filter_multialt(opposite=True).df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   A  C,T    .      .    .     GT  0/2  0/1
        1  chr1  103  .   C  G,A    .      .    .     GT  ./.  1/2

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_multialt(as_index=True)
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        i = self.df.apply(lambda r: ',' not in r.ALT, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_pass(self, opposite=False, as_index=False):
        """Select rows with PASS in the FILTER field.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_pass(as_index=True)
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
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_phased(self, opposite=False, as_index=False):
        """Remove rows with phased genotypes.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_phased(as_index=True)
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
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_polyp(self, opposite=False, as_index=False):
        """Remove rows with a polyploid genotype call.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_polyp(as_index=True)
        0    False
        1     True
        2    False
        3     True
        dtype: bool
        """
        f = lambda r: not any([gt_polyp(x) for x in r[9:]])
        i = self.df.apply(f, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_qual(self, threshold, opposite=False, as_index=False):
        """Select rows with minimum QUAL value.

        Parameters
        ----------
        threshold : float
            Minimum QUAL value.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A', 'T'],
        ...     'QUAL': ['.', 30, 19, 41, 29],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'Steven': ['0/1', '1/1', '0/1', '0/1', '1/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  101  .   T   C   30      .    .     GT    1/1
        2  chr1  102  .   A   T   19      .    .     GT    0/1
        3  chr1  103  .   C   A   41      .    .     GT    0/1
        4  chr1  104  .   C   T   29      .    .     GT    1/1

        We can select rows with minimum QUAL value of 30:

        >>> vf.filter_qual(30).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  101  .   T   C   30      .    .     GT    1/1
        1  chr1  103  .   C   A   41      .    .     GT    0/1

        We can also remove those rows:

        >>> vf.filter_qual(30, opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  102  .   A   T   19      .    .     GT    0/1
        2  chr1  104  .   C   T   29      .    .     GT    1/1

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_qual(30, as_index=True)
        0    False
        1     True
        2    False
        3     True
        4    False
        dtype: bool
        """
        def one_row(r):
            if isinstance(r.QUAL, str):
                return False
            return r.QUAL >= threshold
        i = self.df.apply(one_row, axis=1)
        if opposite:
            i = ~i
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampall(self, samples=None, opposite=False, as_index=False):
        """Select rows if all of the given samples have the variant.

        The default behavior is to use all samples in the VcfFrame.

        Parameters
        ----------
        samples : list, optional
            List of sample names or indicies.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        VcfFrame.filter_sampany
            Similar method that selects rows if any one of the given
            samples has the variant.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_sampall(as_index=True)
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
        if as_index:
            return i
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampany(self, samples=None, opposite=False, as_index=False):
        """Select rows if any one of the given samples has the variant.

        The default behavior is to use all samples in the VcfFrame.

        Parameters
        ----------
        samples : list, optional
            List of sample names or indicies.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        See Also
        --------
        VcfFrame.filter_sampall
            Similar method that selects rows if all of the given
            samples have the variant.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_sampany(as_index=True)
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

    def filter_sampnum(self, threshold, opposite=False, as_index=False):
        """Select rows if the variant is prevalent enough.

        Parameters
        ----------
        threshold : int or float
            Minimum number or fraction of samples with the variant.
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        as_index : bool, default: False
            If True, return boolean index array instead of VcfFrame.

        Returns
        -------
        VcfFrame or pandas.Series
            Filtered VcfFrame or boolean index array.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

        >>> vf.filter_sampnum(2, as_index=True)
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
        if as_index:
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

        >>> from fuc import pyvcf
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

        >>> from fuc import pyvcf
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
        """
        Subset the VcfFrame for specified samples.

        The order of input samples matters.

        Parameters
        ----------
        samples : str or list
            Name or list of sample names.
        exclude : bool, default: False
            If True, exclude the selected samples.

        Returns
        -------
        VcfFrame
            Subsetted VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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
        if isinstance(samples, str):
            samples = [samples]
        if exclude:
            samples = [x for x in self.samples if x not in samples]
        cols = self.df.columns[:9].to_list() + samples
        df = self.df[cols]
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def unphase(self):
        """
        Unphase all the sample genotypes.

        Returns
        -------
        VcfFrame
            Unphased VcfFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
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

    def slice(self, region):
        """
        Slice the VcfFrame for the region.

        Parameters
        ----------
        region : str
            Region ('chrom:start-end').

        Returns
        -------
        VcfFrame
            Sliced VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
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
        >>> vf.slice('chr1:101-300').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  205  .   T   C    .      .    .     GT    1/1
        1  chr1  297  .   A   T    .      .    .     GT    0/1
        >>> vf.slice('chr1').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  205  .   T   C    .      .    .     GT    1/1
        2  chr1  297  .   A   T    .      .    .     GT    0/1
        >>> vf.slice('chr1:-296').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  100  .   G   A    .      .    .     GT    0/1
        1  chr1  205  .   T   C    .      .    .     GT    1/1
        >>> vf.slice('chr1:101').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  205  .   T   C    .      .    .     GT    1/1
        1  chr1  297  .   A   T    .      .    .     GT    0/1
        >>> vf.slice('chr1:101-').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Steven
        0  chr1  205  .   T   C    .      .    .     GT    1/1
        1  chr1  297  .   A   T    .      .    .     GT    0/1
        """
        chrom, start, end = common.parse_region(region)
        df = self.df[self.df.CHROM == chrom]
        if start:
            df = df[df.POS >= start]
        if end:
            df = df[df.POS <= end]
        return self.__class__(self.copy_meta(), df)

    def to_bed(self):
        """Write BedFrame from the VcfFrame.

        Returns
        -------
        BedFrame
            BedFrame.

        Examples
        --------
        Assume we have the following data:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C', 'ACGT'],
        ...     'ALT': ['C', 'T,G', 'G', 'A,G,CT', 'A'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Steven': ['0/1:32', './.:.', '0/1:27', '0/2:34', '0/0:31'],
        ...     'Sara': ['0/0:28', '1/2:30', '1/1:29', '1/2:38', '0/1:27'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID   REF     ALT QUAL FILTER INFO FORMAT  Steven    Sara
        0  chr1  100  .     A       C    .      .    .  GT:DP  0/1:32  0/0:28
        1  chr1  101  .     A     T,G    .      .    .  GT:DP   ./.:.  1/2:30
        2  chr1  102  .     C       G    .      .    .  GT:DP  0/1:27  1/1:29
        3  chr1  103  .     C  A,G,CT    .      .    .  GT:DP  0/2:34  1/2:38
        4  chr1  104  .  ACGT       A    .      .    .  GT:DP  0/0:31  0/1:27

        We can construct BedFrame from the VcfFrame:

        >>> bf = vf.to_bed()
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1    100  100
        1       chr1    101  101
        2       chr1    102  102
        3       chr1    103  103
        4       chr1    103  104
        5       chr1    105  107
        """
        def one_row(r):
            if len(r.REF) == len(r.ALT) == 1:
                start = r.POS
                end = r.POS
            elif len(r.REF) > len(r.ALT):
                start = r.POS + 1
                end = r.POS + len(r.REF) - len(r.ALT)
            else:
                start = r.POS
                end = r.POS + 1
            return pd.Series([r.CHROM, start, end])
        df = self.expand().df.apply(one_row, axis=1)
        df.columns = ['Chromosome', 'Start', 'End']
        df = df.drop_duplicates()
        bf = pybed.BedFrame.from_frame([], df)
        return bf

    def extract(self, k, func=None, as_nan=False):
        """
        Create a DataFrame containing analysis-ready data for the given
        genotype key.

        Parameters
        ----------
        k : str
            Genotype key such as 'DP' and 'AD'.
        func : function, optional
            Function to apply to each of the returned values. For example, if
            the goal is to extract AD only for reference alleles, one can use
            ``func=lambda x: x.split(',')[0]``
        as_nan : bool, default: False
            If True, missing values will be returned as ``NaN``.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing requested data.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['A', 'T'],
        ...     'ALT': ['A', 'T'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT:AD:DP', 'GT:AD:DP'],
        ...     'A': ['0/1:15,13:28', '0/0:28,1:29'],
        ...     'B': ['./.:.:.', '1/1:0,26:26'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A            B
        0  chr1  100  .   A   A    .      .    .  GT:AD:DP  0/1:15,13:28      ./.:.:.
        1  chr1  101  .   T   T    .      .    .  GT:AD:DP   0/0:28,1:29  1/1:0,26:26
        >>> vf.extract('GT')
             A    B
        0  0/1  ./.
        1  0/0  1/1
        >>> vf.extract('GT', as_nan=True)
             A    B
        0  0/1  NaN
        1  0/0  1/1
        >>> vf.extract('AD')
               A     B
        0  15,13     .
        1   28,1  0,26
        >>> vf.extract('AD', as_nan=True, func=lambda x: float(x.split(',')[1]))
              A     B
        0  13.0   NaN
        1   1.0  26.0
        """
        def one_row(r):
            i = r.FORMAT.split(':').index(k)
            def one_gt(g):
                v = g.split(':')[i]
                if as_nan and not i and gt_miss(v):
                    return np.nan
                if as_nan and v == '.':
                    return np.nan
                if func is not None:
                    return func(v)
                return v
            return r[9:].apply(one_gt)
        df = self.df.apply(one_row, axis=1)
        return df

    def rename(self, names, indicies=None):
        """
        Rename the samples.

        Parameters
        ----------
        names : dict or list
            Dict of old names to new names or list of new names.
        indicies : list or tuple, optional
            List of 0-based sample indicies. Alternatively, a tuple
            (int, int) can be used to specify an index range.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
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
        ...     'A': ['0/1', '0/1'],
        ...     'B': ['0/1', '0/1'],
        ...     'C': ['0/1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1
        >>> vf.rename(['X', 'Y', 'Z']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    X    Y    Z
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1
        >>> vf.rename({'B': 'X', 'C': 'Y'}).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    X    Y
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1
        >>> vf.rename(['X'], indicies=[1]).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    X    C
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1
        >>> vf.rename(['X', 'Y'], indicies=(1, 3)).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    X    Y
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1
        """
        samples = self.samples

        if not isinstance(names, list) and not isinstance(names, dict):
            raise TypeError("Argument 'names' must be dict or list.")

        if len(names) > len(samples):
            raise ValueError("There are too many names.")

        if isinstance(names, list) and indicies is not None:
            if isinstance(indicies, tuple):
                if len(indicies) != 2:
                    raise ValueError("Index range must be two integers.")
                l = len(range(indicies[0], indicies[1]))
            elif isinstance(indicies, list):
                l = len(indicies)
            else:
                raise TypeError("Argument 'indicies' must be list or tuple.")

            if len(names) != l:
                raise ValueError("Names and indicies have different lengths.")

        if isinstance(names, list):
            if len(names) == len(samples):
                names = dict(zip(samples, names))
            else:
                if indicies is None:
                    message = ("There are too few names. If this was "
                        "intended, use the 'indicies' argument.")
                    raise ValueError(message)
                elif isinstance(indicies, tuple):
                    names = dict(zip(samples[indicies[0]:indicies[1]], names))
                else:
                    names = dict(zip([samples[i] for i in indicies], names))

        for old, new in names.items():
            i = samples.index(old)
            samples[i] = new

        if len(samples) > len(set(samples)):
            raise ValueError('There are more than one duplicate names.')

        columns = self.df.columns[:9].to_list() + samples
        vf = self.copy()
        vf.df.columns = columns

        return vf

    def drop_duplicates(self, subset=None, keep='first'):
        """
        Return VcfFrame with duplicate rows removed.

        This method essentially wraps the
        :meth:`pandas.DataFrame.drop_duplicates` method.

        Considering certain columns is optional.

        Parameters
        ----------
        subset : column label or sequence of labels, optional
            Only consider certain columns for identifying duplicates, by
            default use all of the columns.
        keep : {'first', 'last', False}, default 'first'
            Determines which duplicates (if any) to keep.

            - ``first`` : Drop duplicates except for the first occurrence.
            - ``last`` : Drop duplicates except for the last occurrence.
            - False : Drop all duplicates.

        Returns
        -------
        VcfFrame
            VcfFrame with duplicates removed.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
        ...     'POS': [100, 100, 200, 200],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'A', 'C', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', './.', '0/1', './.'],
        ...     'B': ['./.', '0/1', './.', '1/2'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   A    C    .      .    .     GT  0/1  ./.
        1  chr1  100  .   A    T    .      .    .     GT  ./.  0/1
        2  chr2  200  .   C    G    .      .    .     GT  0/1  ./.
        3  chr2  200  .   C  G,A    .      .    .     GT  ./.  1/2
        >>> vf.drop_duplicates(['CHROM', 'POS', 'REF']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   A   C    .      .    .     GT  0/1  ./.
        1  chr2  200  .   C   G    .      .    .     GT  0/1  ./.
        >>> vf.drop_duplicates(['CHROM', 'POS', 'REF'], keep='last').df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   A    T    .      .    .     GT  ./.  0/1
        1  chr2  200  .   C  G,A    .      .    .     GT  ./.  1/2
        """
        df = self.df.drop_duplicates(subset=subset, keep=keep)
        return self.__class__(self.copy_meta(), df)

    def plot_titv(
        self, af=None, group_col=None, group_order=None, flip=False, ax=None,
        figsize=None, **kwargs
    ):
        """
        Create a box plot showing the :ref:`Ti/Tv <glossary:Transitions (Ti)
        and transversions (Tv)>` proportions of samples.

        Under the hood, this method simply converts the VcfFrame to the
        :class:`pymaf.MafFrame` class and then applies the
        :meth:`pymaf.MafFrame.plot_titv` method.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.boxplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pymaf.MafFrame.plot_titv
            Similar method for the :class:`fuc.api.pymaf.MafFrame` class.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pyvcf
            >>> common.load_dataset('tcga-laml')
            >>> vcf_file = '~/fuc-data/tcga-laml/tcga_laml.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_titv()
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> vf.plot_titv(af=af,
            ...              group_col='FAB_classification',
            ...              group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        mf = pymaf.MafFrame.from_vcf(self)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_titv(
            af=af, group_col=group_col, group_order=group_order, flip=flip,
            ax=ax, figsize=figsize, **kwargs
        )

        return ax

    def plot_rainfall(
        self, sample, palette=None, ax=None, figsize=None, legend='auto',
        **kwargs
    ):
        """
        Create a rainfall plot visualizing inter-variant distance on a linear
        genomic scale for single sample.

        Under the hood, this method simply converts the VcfFrame to the
        :class:`fuc.api.pymaf.MafFrame` class and then applies the
        :meth:`fuc.api.pymaf.MafFrame.plot_rainfall` method.

        Parameters
        ----------
        sample : str
            Name of the sample.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.scatterplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pymaf.MafFrame.plot_rainfall
            Similar method for the :meth:`fuc.api.pymaf.MafFrame` class.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pyvcf
            >>> common.load_dataset('brca')
            >>> vcf_file = '~/fuc-data/brca/brca.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_rainfall('TCGA-A8-A08B',
            ...                  figsize=(14, 7),
            ...                  palette=sns.color_palette('Set2')[:6])
            >>> plt.tight_layout()
        """
        mf = pymaf.MafFrame.from_vcf(self)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_rainfall(
            sample, palette=palette, ax=ax, legend=legend, **kwargs
        )

        return ax

    def plot_snvclsc(
        self, af=None, group_col=None, group_order=None, palette=None,
        flip=False, ax=None, figsize=None, **kwargs
    ):
        """
        Create a bar plot summarizing the count distrubtions of the six
        :ref:`glossary:SNV classes` for all samples.

        A grouped bar plot can be created with ``group_col`` (requires an AnnFrame).

        Under the hood, this method simply converts the VcfFrame to the
        :class:`fuc.api.pymaf.MafFrame` class and then applies the
        :meth:`fuc.api.pymaf.MafFrame.plot_snvclsc` method.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pymaf.MafFrame.plot_snvclsc
            Similar method for the :meth:`fuc.api.pymaf.MafFrame` class.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pyvcf
            >>> common.load_dataset('tcga-laml')
            >>> vcf_file = '~/fuc-data/tcga-laml/tcga_laml.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_snvclsc(palette=sns.color_palette('Pastel1'))
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> vf.plot_snvclsc(af=af,
            ...                 group_col='FAB_classification',
            ...                 group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        mf = pymaf.MafFrame.from_vcf(self)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_snvclsc(
            af=af, group_col=group_col, group_order=group_order,
            palette=palette, flip=flip, ax=ax, **kwargs
        )

        return ax

    def plot_snvclsp(
        self, af=None, group_col=None, group_order=None, palette=None, flip=False,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a box plot summarizing the proportion distrubtions of the six
        :ref:`glossary:SNV classes` for all sample.

        Under the hood, this method simply converts the VcfFrame to the
        :class:`fuc.api.pymaf.MafFrame` class and then applies the
        :meth:`fuc.api.pymaf.MafFrame.plot_snvclsp` method.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.boxplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pymaf.MafFrame.plot_snvclsp
            Similar method for the :meth:`fuc.api.pymaf.MafFrame` class.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pyvcf
            >>> common.load_dataset('tcga-laml')
            >>> vcf_file = '~/fuc-data/tcga-laml/tcga_laml.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_snvclsp(palette=sns.color_palette('Pastel1'))
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> vf.plot_snvclsp(af=af,
            ...                 group_col='FAB_classification',
            ...                 group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        mf = pymaf.MafFrame.from_vcf(self)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_snvclsp(
            af=af, group_col=group_col, group_order=group_order,
            palette=palette, flip=flip, ax=ax, **kwargs
        )

        return ax

    def plot_snvclss(
        self, color=None, colormap=None, width=0.8, legend=True, flip=False,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a bar plot showing the proportions of the six
        :ref:`glossary:SNV classes` for individual samples.

        Under the hood, this method simply converts the VcfFrame to the
        :class:`fuc.api.pymaf.MafFrame` class and then applies the
        :meth:`fuc.api.pymaf.MafFrame.plot_snvclss` method.

        Parameters
        ----------
        color : list, optional
            List of color tuples. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        colormap : str or matplotlib colormap object, optional
            Colormap to select colors from. See the :ref:`tutorials:Control
            plot colors` tutorial for details.
        width : float, default: 0.8
            The width of the bars.
        legend : bool, default: True
            Place legend on axis subplots.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.bar` or
            :meth:`pandas.DataFrame.plot.barh`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pymaf.MafFrame.plot_snvclss
            Similar method for the :meth:`fuc.api.pymaf.MafFrame` class.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> ax = mf.plot_snvclss(width=1, color=plt.get_cmap('Pastel1').colors)
            >>> ax.legend(loc='upper right')
            >>> plt.tight_layout()
        """
        mf = pymaf.MafFrame.from_vcf(self)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_snvclss(
            color=color, colormap=colormap, width=width, legend=legend,
            flip=flip, af=af, **kwargs
        )

        return ax

    def chr_prefix(self, mode='remove'):
        """
        Add or remove the (annoying) 'chr' string from the CHROM column.

        Parameters
        ----------
        mode : {'add', 'remove'}, default: 'remove'
            Whether to add or remove the 'chr' string.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
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
        ...     'A': ['0/1', '1/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0/1
        1  chr2  101  .   T   C    .      .    .     GT  1/1
        >>> vf.chr_prefix().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0     1  100  .   G   A    .      .    .     GT  0/1
        1     2  101  .   T   C    .      .    .     GT  1/1
        """
        if mode == 'remove':
            def one_row(r):
                r.CHROM = r.CHROM.replace('chr', '')
                return r
        elif mode == 'add':
            def one_row(r):
                r.CHROM = 'chr' + r.CHROM
                return r
        else:
            raise ValueError(f'Incorrect mode: {mode}')
        df = self.df.apply(one_row, axis=1)
        return self.__class__(self.copy_meta(), df)
