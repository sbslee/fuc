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

Genotype lines usually consist of nine columns for storing variant
information (all fixed and mandatory except for the FORMAT column) plus
additional sample-specific columns for expressing individual genotype calls
(e.g. '0/1'). Missing values are allowed in some cases and can be specified
with a dot ('.'). The first nine columns are:

.. list-table::
   :header-rows: 1

   * - No.
     - Column
     - Description
     - Required
     - Missing
     - Examples
   * - 1
     - CHROM
     - Chromosome or contig identifier
     - ✅
     - ❌
     - 'chr2', '2', 'chrM'
   * - 2
     - POS
     - 1-based reference position
     - ✅
     - ❌
     - 10041, 23042
   * - 3
     - ID
     - ';'-separated variant identifiers
     - ✅
     - ✅
     - '.', 'rs35', 'rs9;rs53'
   * - 4
     - REF
     - Reference allele
     - ✅
     - ❌
     - 'A', 'GT'
   * - 5
     - ALT
     - ','-separated alternate alleles
     - ✅
     - ❌
     - 'T', 'ACT', 'C,T'
   * - 6
     - QUAL
     - Phred-scaled quality score for ALT
     - ✅
     - ✅
     - '.', 67, 12
   * - 7
     - FILTER
     - ';'-separated filters that failed
     - ✅
     - ✅
     - '.', 'PASS', 'q10;s50'
   * - 8
     - INFO
     - ';'-separated information fields
     - ✅
     - ✅
     - '.', 'DP=14;AF=0.5;DB'
   * - 9
     - FORMAT
     - ':'-separated genotype fields
     - ❌
     - ❌
     - 'GT', 'GT:AD:DP'

You will sometimes come across VCF files that have only eight columns, and
do not contain the FORMAT column or sample-specific information. These are
called "sites-only" VCF files, and normally represent genetic variation that
has been observed in a large population. Generally, information about the
population of origin should be included in the header.

There are several reserved keywords in the INFO and FORMAT columns that are
standards across the community. Popular keywords are listed below:

.. list-table::
   :header-rows: 1

   * - Column
     - Key
     - Number
     - Type
     - Description
   * - INFO
     - AC
     - A
     - Integer
     - Allele count in genotypes, for each ALT allele, in the same order as listed
   * - INFO
     - AN
     - 1
     - Integer
     - Total number of alleles in called genotypes
   * - INFO
     - AF
     - A
     - Float
     - Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)
   * - FORMAT
     - AD
     - R
     - Integer
     - Total read depth for each allele
   * - FORMAT
     - AF
     - 1
     - Float
     - Allele fraction of the event in the tumor
   * - FORMAT
     - DP
     - 1
     - Integer
     - Read depth

If sample annotation data are available for a given VCF file, use
the :class:`common.AnnFrame` class to import the data.
"""

import os
import re
import sys
import gzip
from copy import deepcopy
import warnings
import tempfile

from . import pybed, common, pymaf, pybam

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.formula.api as smf
from Bio import bgzf
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
from pysam import VariantFile, bcftools
from io import StringIO

HEADERS = {
    'CHROM': str,
    'POS': int,
    'ID': str,
    'REF': str,
    'ALT': str,
    'QUAL': str,
    'FILTER': str,
    'INFO': str,
    'FORMAT': str,
}

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

INFO_SPECIAL_KEYS = {
    '#AC': ['AC', lambda x: sum([int(x) for x in x.split(',')]), True],
    '#AF': ['AF', lambda x: sum([float(x) for x in x.split(',')]), True],
}

FORMAT_SPECIAL_KEYS = {
    '#DP': ['DP', lambda x: int(x), True],
    '#AD_REF': ['AD', lambda x: float(x.split(',')[0]), True],
    '#AD_ALT': ['AD', lambda x: sum([int(y) for y in x.split(',')[1:]]), True],
    '#AD_FRAC_REF': ['AD', lambda x: np.nan if sum([int(y) for y in x.split(',')]) == 0 else int(x.split(',')[0]) / sum([int(y) for y in x.split(',')]), True],
    '#AD_FRAC_ALT': ['AD', lambda x: np.nan if sum([int(y) for y in x.split(',')]) == 0 else sum([int(y) for y in x.split(',')[1:]]) / sum([int(y) for y in x.split(',')]), True],
}

def call(
    fasta, bams, regions=None, path=None, min_mq=1, max_depth=250,
    dir_path=None, gap_frac=0.002, group_samples=None
):
    """
    Call SNVs and indels from BAM files.

    Under the hood, the method utilizes the bcftool program to call variants.

    Parameters
    ----------
    fasta : str
        Reference FASTA file.
    bams : str or list
        One or more input BAM files. Alternatively, you can provide a text
        file (.txt, .tsv, .csv, or .list) containing one BAM file per line.
    regions : str, list, or pybed.BedFrame, optional
        By default (``regions=None``), the method looks at each genomic
        position with coverage in BAM files, which can be excruciatingly slow
        for large files (e.g. whole genome sequencing). Therefore, use this
        argument to only call variants in given regions. Each region must
        have the format chrom:start-end and be a half-open interval with
        (start, end]. This means, for example, chr1:100-103 will extract
        positions 101, 102, and 103. Alternatively, you can provide a BED
        file (compressed or uncompressed) or a :class:`pybed.BedFrame` object
        to specify regions. Note that the 'chr' prefix in contig names (e.g.
        'chr1' vs. '1') will be automatically added or removed as necessary
        to match the input VCF's contig names.
    path : str, optional
        Output VCF file. Writes to stdout when ``path='-'``. If None is
        provided the result is returned as a string.
    min_mq : int, default: 1
        Minimum mapping quality for an alignment to be used.
    max_depth : int, default: 250
        At a position, read maximally this number of reads per input file.
    dir_path : str, optional
        By default, intermediate files (likelihoods.bcf, calls.bcf, and
        calls.normalized.bcf) will be stored in a temporary directory, which
        is automatically deleted after creating final VCF. If you provide a
        directory path, intermediate files will be stored there.
    gap_frac : float, default: 0.002
        Minimum fraction of gapped reads for calling indels.
    group_samples : str, optional
        By default, all samples are assumed to come from a single population.
        This option allows to group samples into populations and apply the
        HWE assumption within but not across the populations. To use this
        option, provide a tab-delimited text file with sample names in the
        first column and group names in the second column. If
        ``group_samples='-'`` is given instead, no HWE assumption is made at
        all and single-sample calling is performed. Note that in low coverage
        data this inflates the rate of false positives. Therefore, make sure
        you know what you are doing.

    Returns
    -------
    str
        VcfFrame object.
    """
    # Parse input BAM files.
    bams = common.parse_list_or_file(bams)

    # Check the 'chr' prefix.
    if all([pybam.has_chr_prefix(x) for x in bams]):
        chr_prefix = 'chr'
    else:
        chr_prefix = ''

    # Parse target regions, if provided.
    if regions is not None:
        if isinstance(regions, str):
            regions = [regions]
        elif isinstance(regions, list):
            pass
        elif isinstance(regions, pybed.BedFrame):
            regions = regions.to_regions()
        else:
            raise TypeError("Incorrect type of argument 'regions'")
        if '.bed' in regions[0]:
            bf = pybed.BedFrame.from_file(regions[0])
            regions = bf.to_regions()
        else:
            regions = common.sort_regions(regions)
        regions = [chr_prefix + x.replace('chr', '') for x in regions]

    if dir_path is None:
        t = tempfile.TemporaryDirectory()
        temp_dir = t.name
    else:
        temp_dir = dir_path

    # Step 1: Get genotype likelihoods.
    args = ['-Ou', '-a', 'AD']
    args += ['-q', str(min_mq)]
    args += ['--max-depth', str(max_depth)]
    args += ['-f', fasta]
    args += ['-F', str(gap_frac)]
    args += ['-o', f'{temp_dir}/likelihoods.ubcf']
    if regions is not None:
        args += ['-r', ','.join(regions)]
    bcftools.mpileup(*(args + bams), catch_stdout=False)

    # Step 2: Call variants.
    args = [f'{temp_dir}/likelihoods.ubcf', '-Ou', '-mv']
    args += ['-o', f'{temp_dir}/calls.ubcf']
    if group_samples is not None:
        args += ['-G', group_samples]
    bcftools.call(*args, catch_stdout=False)

    # Step 3: Normalize indels.
    args = [f'{temp_dir}/calls.ubcf', '-Ou', '-f', fasta]
    args += ['-o', f'{temp_dir}/calls.normalized.ubcf']
    bcftools.norm(*args, catch_stdout=False)

    # Step 4: Filter variant.
    args = [f'{temp_dir}/calls.normalized.ubcf', '-Ov', '--IndelGap', '5']
    results = bcftools.filter(*args)

    if path is None:
        return results
    elif path == '-':
        sys.stdout.write(results)
    else:
        with open(path, 'w') as f:
            f.write(results)

    if dir_path is None:
        t.cleanup()

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
    For given genotype, return True if it has missing value.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if genotype is missing.

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

def gt_ploidy(g):
    """
    For given genotype, return its ploidy number.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    int
        Ploidy number.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.gt_ploidy('1')
    1
    >>> pyvcf.gt_ploidy('.')
    1
    >>> pyvcf.gt_ploidy('0/1')
    2
    >>> pyvcf.gt_ploidy('./.')
    2
    >>> pyvcf.gt_ploidy('0|1')
    2
    >>> pyvcf.gt_ploidy('1|0|1')
    3
    >>> pyvcf.gt_ploidy('0/./1/1')
    4
    """
    gt = g.split(':')[0]
    if '/' in gt:
        return gt.count('/') + 1
    elif '|' in gt:
        return gt.count('|') + 1
    else:
        return 1

def gt_polyp(g):
    """
    For given genotype, return True if it is polyploid.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if genotype is polyploid.

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
    return gt_ploidy(g) > 2

def gt_hasvar(g):
    """
    For given genotype, return True if it has variation.

    Parameters
    ----------
    g : str
        Sample genotype.

    Returns
    -------
    bool
        True if genotype has variation.

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
    For given genotype, return its unphased form.

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

def gt_het(g):
    """
    For given genotype, return True if it is heterozygous.

    Parameters
    ----------
    g : str
        Genotype call.

    Returns
    -------
    bool
        True if genotype is heterozygous.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.gt_het('0/1')
    True
    >>> pyvcf.gt_het('0/0')
    False
    >>> pyvcf.gt_het('0|0')
    False
    >>> pyvcf.gt_het('1|0')
    True
    """
    l = g.split(':')
    gt = l[0]
    if '/' in gt:
        gt = gt.split('/')
    elif '|' in gt:
        gt = gt.split('|')
    else:
        return False
    return gt[0] != gt[1]

def gt_pseudophase(g):
    """
    For given genotype, return its pseudophased form.

    Parameters
    ----------
    g : str
        Genotype call.

    Returns
    -------
    str
        Pseudophased genotype call.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> pyvcf.pseudophase('0/1')
    '0|1'
    >>> pyvcf.pseudophase('0/0:34:10,24')
    '0|0:34:10,24'
    """
    l = g.split(':')
    l[0] = l[0].replace('/', '|')
    return ':'.join(l)

def has_chr_prefix(file, size=1000):
    """
    Return True if all of the sampled contigs from a VCF file have the
    (annoying) 'chr' string.

    Parameters
    ----------
    file : str
        VCF file (compressed or uncompressed).
    size : int, default: 1000
        Sampling size.

    Returns
    -------
    bool
        True if the 'chr' string is found.
    """
    n = 0
    vcf = VariantFile(file)
    for record in vcf.fetch():
        n += 1
        if 'chr' not in record.chrom:
            return False
        if n > size:
            break
    vcf.close()
    return True

def merge(
    vfs, how='inner', format='GT', sort=True, collapse=False
):
    """
    Merge VcfFrame objects.

    Parameters
    ----------
    vfs : list
        List of VcfFrames to be merged. Note that the 'chr' prefix in contig
        names (e.g. 'chr1' vs. '1') will be automatically added or removed as
        necessary to match the contig names of the first VCF.
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
    """
    For given row, return True if it has indel.

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


def row_computeinfo(r, key, decimals=3):
    """
    For given row, return AC/AN/AF calculation for INFO column.

    Parameters
    ----------
    r : pandas.Series
        VCF row.
    key : {'AC', 'AN', 'AF'}
        INFO key.
    decimals : int, default: 3
        Number of decimals to display for AF.

    Returns
    -------
    str
        Requested INFO data.

    Example
    -------

    >>> from fuc import pyvcf
    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chrX'],
    ...     'POS': [100, 101, 102, 100],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'T', 'A', 'A'],
    ...     'ALT': ['A', 'C', 'T,G', 'G'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': ['.', '.', '.', '.'],
    ...     'FORMAT': ['GT:DP', 'GT', 'GT', 'GT'],
    ...     'A': ['1|0:29', '0|0', '1|0', '0'],
    ...     'B': ['1/1:34', './.', '0/0', '0/1'],
    ...     'C': ['0/0:23', '0/0', '1/2', '1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT       A       B       C
    0  chr1  100  .   G    A    .      .    .  GT:DP  1|0:29  1/1:34  0/0:23
    1  chr1  101  .   T    C    .      .    .     GT     0|0     ./.     0/0
    2  chr1  102  .   A  T,G    .      .    .     GT     1|0     0/0     1/2
    3  chrX  100  .   A    G    .      .    .     GT       0     0/1       1
    >>> pyvcf.row_computeinfo(vf.df.iloc[0, :], 'AC')
    '3'
    >>> pyvcf.row_computeinfo(vf.df.iloc[0, :], 'AN')
    '6'
    >>> pyvcf.row_computeinfo(vf.df.iloc[0, :], 'AF')
    '0.500'
    >>> pyvcf.row_computeinfo(vf.df.iloc[1, :], 'AC')
    '0'
    >>> pyvcf.row_computeinfo(vf.df.iloc[1, :], 'AN')
    '4'
    >>> pyvcf.row_computeinfo(vf.df.iloc[1, :], 'AF')
    '0.000'
    >>> pyvcf.row_computeinfo(vf.df.iloc[2, :], 'AC')
    '2,1'
    >>> pyvcf.row_computeinfo(vf.df.iloc[2, :], 'AN')
    '6'
    >>> pyvcf.row_computeinfo(vf.df.iloc[2, :], 'AF')
    '0.333,0.167'
    >>> pyvcf.row_computeinfo(vf.df.iloc[3, :], 'AC')
    '2'
    >>> pyvcf.row_computeinfo(vf.df.iloc[3, :], 'AN')
    '4'
    >>> pyvcf.row_computeinfo(vf.df.iloc[3, :], 'AF')
    '0.500'
    """
    def get_ac(r):
        def one_gt(g, i):
            gt = g.split(':')[0]
            if '/' in gt:
                l = gt.split('/')
            else:
                l = gt.split('|')
            return l.count(str(i + 1))
        counts = []
        for i, allele in enumerate(r.ALT.split(',')):
            count = r[9:].apply(one_gt, args=(i, ))
            counts.append(sum(count))
        return counts

    def get_an(r):
        def one_gt(g):
            if '.' in g:
                return 0
            else:
                return gt_ploidy(g)
        return r[9:].apply(one_gt).sum()

    def get_af(r):
        return [x / get_an(r) for x in get_ac(r)]

    methods = {'AC': get_ac, 'AN': get_an, 'AF': get_af}

    results = methods[key](r)

    if key == 'AC':
        results = ','.join([str(x) for x in results])
    elif key == 'AN':
        results = str(results)
    else:
        results = ','.join([f'{x:.{decimals}f}' for x in results])

    return results

def row_parseinfo(r, key):
    """
    For given row, return requested data from INFO column.

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

def row_phased(r):
    """
    For given row, return True if all genotypes are phased.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    bool
        True if the row is fully phased.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['G', 'T', 'A'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['.', '.', '.'],
    ...     'INFO': ['.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'A': ['1|1', '0/0', '1|0'],
    ...     'B': ['1|0', '0/1', '1/0'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
    0  chr1  100  .   G   A    .      .    .     GT  1|1  1|0
    1  chr1  101  .   T   C    .      .    .     GT  0/0  0/1
    2  chr1  102  .   A   T    .      .    .     GT  1|0  1/0
    >>> vf.df.apply(pyvcf.row_phased, axis=1)
    0     True
    1    False
    2    False
    dtype: bool
    """
    def one_gt(g):
        return '|' in g.split(':')[0]
    return r[9:].apply(one_gt).all()

def row_updateinfo(r, key, value, force=True, missing=False):
    """
    For given row, return updated data from INFO column.

    Parameters
    ----------
    r : pandas.Series
        VCF row.
    key : str
        INFO key.
    value : str
        New value to be assigned.
    force : bool, default: True
        If True, overwrite any existing data.
    missing : bool, default: False
        By default (``missing=False``), the method will not update INFO field
        if there is missing value ('.'). Setting ``missing=True`` will stop
        this behavior.

    Returns
    -------
    str
        New INFO field.

    Examples
    --------

    >>> from fuc import pyvcf
    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102],
    ...     'ID': ['.', '.', '.'],
    ...     'REF': ['G', 'T', 'A'],
    ...     'ALT': ['A', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.'],
    ...     'FILTER': ['.', '.', '.'],
    ...     'INFO': ['DB;AC=1', 'DB', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT'],
    ...     'A': ['0/1', '1/1', '0/0'],
    ...     'B': ['0/0', '0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict([], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER     INFO FORMAT    A    B
    0  chr1  100  .   G   A    .      .  DB;AC=1     GT  0/1  0/0
    1  chr1  101  .   T   C    .      .       DB     GT  1/1  0/1
    2  chr1  102  .   A   T    .      .        .     GT  0/0  0/1
    >>> pyvcf.row_updateinfo(vf.df.iloc[0, :], 'AC', '100')
    'DB;AC=100'
    >>> pyvcf.row_updateinfo(vf.df.iloc[0, :], 'AC', '100', force=False)
    'DB;AC=1'
    >>> pyvcf.row_updateinfo(vf.df.iloc[1, :], 'AC', '100')
    'DB;AC=100'
    >>> pyvcf.row_updateinfo(vf.df.iloc[1, :], 'AC', '100', force=False)
    'DB;AC=100'
    >>> pyvcf.row_updateinfo(vf.df.iloc[2, :], 'AC', '100')
    '.'
    >>> pyvcf.row_updateinfo(vf.df.iloc[2, :], 'AC', '100', force=False)
    '.'
    >>> pyvcf.row_updateinfo(vf.df.iloc[2, :], 'AC', '100', missing=True)
    'AC=100'
    """
    data = f'{key}={value}'

    if not missing and r.INFO == '.':
        return r.INFO

    fields = r.INFO.split(';')
    found = False

    for i, field in enumerate(fields):
        if field.startswith(f'{key}='):
            found = True
            if force:
                fields[i] = data
            else:
                pass
            break

    if not found:
        if r.INFO == '.':
            fields = [data]
        else:
            fields.append(data)

    return ';'.join(fields)

def row_missval(r):
    """
    For given row, return formatted missing genotype.

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

def slice(file, regions, path=None):
    """
    Slice a VCF file for specified regions.

    Parameters
    ----------
    file : str
        Input VCF file must be already BGZF compressed (.gz) and indexed
        (.tbi) to allow random access.
    regions : str, list, or pybed.BedFrame
        One or more regions to be sliced. Each region must have the format
        chrom:start-end and be a half-open interval with (start, end]. This
        means, for example, chr1:100-103 will extract positions 101, 102, and
        103. Alternatively, you can provide a BED file (compressed or
        uncompressed) to specify regions. Note that the 'chr' prefix in
        contig names (e.g. 'chr1' vs. '1') will be automatically added or
        removed as necessary to match the input VCF's contig names.
    path : str, optional
        Output VCF file. Writes to stdout when ``path='-'``. If None is
        provided the result is returned as a string.

    Returns
    -------
    None or str
        If ``path`` is None, returns the resulting VCF format as a string.
        Otherwise returns None.
    """
    if isinstance(regions, str):
        regions = [regions]

    if isinstance(regions, pybed.BedFrame):
        regions = regions.to_regions()
    elif isinstance(regions, list):
        if '.bed' in regions[0]:
            regions = pybed.BedFrame.from_file(regions[0]).to_regions()
        else:
            regions = common.sort_regions(regions)
    else:
        raise TypeError('Incorrect regions type')

    if has_chr_prefix(file):
        regions = common.update_chr_prefix(regions, mode='add')
    else:
        regions = common.update_chr_prefix(regions, mode='remove')

    vcf = VariantFile(file)

    if path is None:
        data = ''
        data += str(vcf.header)
        for region in regions:
            chrom, start, end = common.parse_region(region)
            if np.isnan(start):
                start = None
            if np.isnan(end):
                end = None
            for record in vcf.fetch(chrom, start, end):
                data += str(record)
    else:
        data = None
        output = VariantFile(path, 'w', header=vcf.header)
        for region in regions:
            chrom, start, end = common.parse_region(region)
            for record in vcf.fetch(chrom, start, end):
                output.write(record)
        output.close()

    vcf.close()

    return data

def split(vcf, clean=True):
    """
    Split VcfFrame by individual.

    Parameters
    ----------
    vcf : str or VcfFrame
        VCF file or VcfFrame to be split.
    clean : bool, default: True
        If True, only return variants present in each individual. In other
        words, setting this as False will make sure that all individuals
        have the same number of variants.

    Returns
    -------
    list
        List of individual VcfFrames.

    Examples
    --------

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
    ...     'FORMAT': ['GT', 'GT'],
    ...     'A': ['0/1', '0/0'],
    ...     'B': ['0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict(['##fileformat=VCFv4.3'], data)
    >>> vf.df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
    0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1
    1  chr1  101  .   T   C    .      .    .     GT  0/0  0/1
    >>> vfs = pyvcf.split(vf)
    >>> vfs[0].df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
    0  chr1  100  .   G   A    .      .    .     GT  0/1
    >>> vfs = pyvcf.split(vf, clean=False)
    >>> vfs[0].df
      CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
    0  chr1  100  .   G   A    .      .    .     GT  0/1
    1  chr1  101  .   T   C    .      .    .     GT  0/0
    """
    # Parse the input VCF.
    if isinstance(vcf, str):
        vf = VcfFrame.from_file(vcf)
    else:
        vf = vcf

    vfs = []

    for sample in vf.samples:
        temp_vf = vf.subset(sample)
        if clean:
            temp_vf = temp_vf.filter_sampall()
        vfs.append(temp_vf)

    return vfs

def plot_af_correlation(vf1, vf2, ax=None, figsize=None):
    """
    Create a scatter plot showing the correlation of allele frequency between
    two VCF files.

    This method will exclude the following sites:

        - non-onverlapping sites
        - multiallelic sites
        - sites with one or more missing genotypes

    Parameters
    ----------
    vf1, vf2 : VcfFrame
        VcfFrame objects to be compared.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, crete a new one.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        The matplotlib axes containing the plot.

    Examples
    --------
    .. plot::
        :context: close-figs

        >>> from fuc import pyvcf, common
        >>> import matplotlib.pyplot as plt
        >>> data1 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104, 105],
        ...     'ID': ['.', '.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'C', 'G,A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1:30', '0/0', '1/1', '0/1', '1/1', '0/1'],
        ...     'B': ['0/0:30', '0/0', '0/1', '0/1', '1/1', '0/1'],
        ...     'C': ['1/1:30', '0/0', '1/1', '0/1', '1/1', '0/1'],
        ...     'D': ['0/0:30', '0/0', '0/0', '0/0', '1/1', '0/1'],
        ...     'E': ['0/0:30', '0/0', '0/0', '1/2', '1/1', '0/1'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> data2 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [101, 102, 103, 104, 105],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['T', 'G', 'T', 'A', 'C'],
        ...     'ALT': ['C', 'C', 'G,A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'F': ['0/0', '0/1', '0/1', '1/1', '0/0'],
        ...     'G': ['0/0', '0/1', '0/1', '1/1', './.'],
        ...     'H': ['0/0', '0/1', '0/1', '1/1', '1/1'],
        ...     'I': ['0/0', '0/1', '0/0', '1/1', '1/1'],
        ...     'J': ['0/0', '0/1', '1/2', '1/1', '0/1'],
        ... }
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> pyvcf.plot_af_correlation(vf1, vf2)
        >>> plt.tight_layout()
    """
    def one_gt(g):
        alleles = g.split(':')[0].split('/')
        alleles = [x for x in alleles if x != '0']
        return len(alleles)

    def one_row(r):
        locus = f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}'
        ac = r[9:].apply(one_gt).sum()
        if 'X' in r.CHROM or 'Y' in r.CHROM:
            total = len(r[9:])
        else:
            total = len(r[9:]) * 2
        af = ac / total
        return pd.Series([locus, af])

    s1 = vf1.filter_multialt().filter_empty(threshold=1).df.apply(one_row, axis=1)
    s2 = vf2.filter_multialt().filter_empty(threshold=1).df.apply(one_row, axis=1)

    s1.columns = ['Locus', 'First']
    s2.columns = ['Locus', 'Second']

    s1 = s1.set_index('Locus')
    s2 = s2.set_index('Locus')

    df = pd.concat([s1, s2], axis=1).dropna()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(data=df, x='First', y='Second', ax=ax)

    return ax

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
        Construct VcfFrame from a dict of array-like or dicts.
    VcfFrame.from_file
        Construct VcfFrame from a VCF file.
    VcfFrame.from_string
        Construct VcfFrame from a string.

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

    def _check_df(self, df):
        df = df.reset_index(drop=True)
        df = df.astype(HEADERS)
        return df

    def __init__(self, meta, df):
        self._meta = meta
        self._df = self._check_df(df)

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
        self._df = self._check_df(value)

    @property
    def samples(self):
        """list : List of sample names."""
        return self.df.columns[9:].to_list()

    @property
    def shape(self):
        """tuple : Dimensionality of VcfFrame (variants, samples)."""
        return (self.df.shape[0], len(self.samples))

    @property
    def contigs(self):
        """list : List of contig names."""
        return list(self.df.CHROM.unique())

    @property
    def has_chr_prefix(self):
        """bool : Whether the (annoying) 'chr' string is found."""
        for contig in self.contigs:
            if 'chr' in contig:
                return True
        return False

    @property
    def sites_only(self):
        """bool : Whether the VCF is sites-only."""
        return not self.samples or 'FORMAT' not in self.df.columns

    @property
    def phased(self):
        """
        Return True if every genotype in VcfFrame is haplotype phased.

        Returns
        -------
        bool
            If VcfFrame is fully phased, return True, if not return False.
            Also return False if VcfFrame is empty.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data1 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT'],
        ...     'A': ['1|1', '0|0', '1|0'],
        ...     'B': ['1|0', '0|1', '1|0'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> vf1.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  1|1  1|0
        1  chr1  101  .   T   C    .      .    .     GT  0|0  0|1
        2  chr1  102  .   A   T    .      .    .     GT  1|0  1|0
        >>> vf1.phased
        True
        >>> data2 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT'],
        ...     'C': ['1|1', '0/0', '1|0'],
        ...     'D': ['1|0', '0/1', '1|0'],
        ... }
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> vf2.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    C    D
        0  chr1  100  .   G   A    .      .    .     GT  1|1  1|0
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/1
        2  chr1  102  .   A   T    .      .    .     GT  1|0  1|0
        >>> vf2.phased
        False
        """
        if self.empty:
            return False
        s = self.df.apply(row_phased, axis=1)
        return s.all()

    @property
    def empty(self):
        """
        Indicator whether VcfFrame is empty.

        Returns
        -------
        bool
            If VcfFrame is empty, return True, if not return False.

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
        >>> vf.df = vf.df[0:0]
        >>> vf.df
        Empty DataFrame
        Columns: [CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, A]
        Index: []
        >>> vf.empty
        True
        """
        return self.df.empty

    def add_af(self, decimals=3):
        """
        Compute AF from AD and then add it to the FORMAT field.

        This method will compute allele fraction for each ALT allele in the
        same order as listed.

        Parameters
        ----------
        decimals : int, default: 3
            Number of decimals to display.

        Returns
        -------
        VcfFrame
            Updated VcfFrame object.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'G', 'A', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:AD', 'GT:AD', 'GT', 'GT:AD'],
        ...     'A': ['0/1:12,15', '0/0:32,1', '0/1', './.:.'],
        ...     'B': ['0/1:13,17', '0/1:14,15', './.', '1/2:0,11,17'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO FORMAT          A            B
        0  chr1  100  .   A    C    .      .    .  GT:AD  0/1:12,15    0/1:13,17
        1  chr1  101  .   G    T    .      .    .  GT:AD   0/0:32,1    0/1:14,15
        2  chr1  102  .   A    G    .      .    .     GT        0/1          ./.
        3  chr1  103  .   C  G,A    .      .    .  GT:AD      ./.:.  1/2:0,11,17
        >>> vf.add_af().df
          CHROM  POS ID REF  ALT QUAL FILTER INFO    FORMAT                      A                              B
        0  chr1  100  .   A    C    .      .    .  GT:AD:AF  0/1:12,15:0.444,0.556          0/1:13,17:0.433,0.567
        1  chr1  101  .   G    T    .      .    .  GT:AD:AF   0/0:32,1:0.970,0.030          0/1:14,15:0.483,0.517
        2  chr1  102  .   A    G    .      .    .     GT:AF                  0/1:.                          ./.:.
        3  chr1  103  .   C  G,A    .      .    .  GT:AD:AF                ./.:.:.  1/2:0,11,17:0.000,0.393,0.607
        """
        def one_row(r):
            try:
                i = r.FORMAT.split(':').index('AD')
            except ValueError:
                i = None

            def one_gt(g):
                if i is None:
                    ad = None
                else:
                    ad = g.split(':')[i]

                if ad is None or ad == '.':
                    af = '.'
                else:
                    depths = [int(x) for x in ad.split(',')]
                    total = sum(depths)
                    if total == 0:
                        af = '.'
                    else:
                        af = ','.join([f'{x/total:.{decimals}f}' for x in depths])

                return f'{g}:{af}'

            r.iloc[9:] = r.iloc[9:].apply(one_gt)
            r.FORMAT += ':AF'

            return r

        df = self.df.apply(one_row, axis=1)

        return self.__class__(self.copy_meta(), df)

    def add_dp(self):
        """
        Compute DP using AD and add it to the FORMAT field.

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
        """
        Add the given flag to the INFO field.

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

    def collapse(self):
        """
        Collapse duplicate records in the VcfFrame.

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

    def compute_info(self, key):
        """
        Compute AC/AN/AF for INFO column.

        The method will ignore and overwrite any existing data for selected
        key.

        Returns
        -------
        VcfFrame
            Updated VcfFrame.
        key : {'AC', 'AN', 'AF'}
            INFO key.

        Example
        -------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chrX'],
        ...     'POS': [100, 101, 102, 100],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'T,G', 'A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['AC=100', 'MQ=59', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT', 'GT', 'GT'],
        ...     'A': ['1|0:34', '0|0', '1|0', '0'],
        ...     'B': ['1/1:23', '0/1', '0/0', '0/0'],
        ...     'C': ['0/0:28', './.', '1/2', '1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER    INFO FORMAT       A       B       C
        0  chr1  100  .   G    A    .      .  AC=100  GT:DP  1|0:34  1/1:23  0/0:28
        1  chr1  101  .   T    C    .      .   MQ=59     GT     0|0     0/1     ./.
        2  chr1  102  .   A  T,G    .      .       .     GT     1|0     0/0     1/2
        3  chrX  100  .   C    A    .      .       .     GT       0     0/0       1
        >>> vf = vf.compute_info('AC')
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER        INFO FORMAT       A       B       C
        0  chr1  100  .   G    A    .      .        AC=1  GT:DP  1|0:34  1/1:23  0/0:28
        1  chr1  101  .   T    C    .      .  MQ=59;AC=1     GT     0|0     0/1     ./.
        2  chr1  102  .   A  T,G    .      .      AC=1,1     GT     1|0     0/0     1/2
        3  chrX  100  .   C    A    .      .        AC=1     GT       0     0/0       1
        >>> vf = vf.compute_info('AN')
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER             INFO FORMAT       A       B       C
        0  chr1  100  .   G    A    .      .        AC=1;AN=6  GT:DP  1|0:34  1/1:23  0/0:28
        1  chr1  101  .   T    C    .      .  MQ=59;AC=1;AN=4     GT     0|0     0/1     ./.
        2  chr1  102  .   A  T,G    .      .      AC=1,1;AN=6     GT     1|0     0/0     1/2
        3  chrX  100  .   C    A    .      .        AC=1;AN=4     GT       0     0/0       1
        >>> vf = vf.compute_info('AF')
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER                        INFO FORMAT       A       B       C
        0  chr1  100  .   G    A    .      .          AC=1;AN=6;AF=0.167  GT:DP  1|0:34  1/1:23  0/0:28
        1  chr1  101  .   T    C    .      .    MQ=59;AC=1;AN=4;AF=0.250     GT     0|0     0/1     ./.
        2  chr1  102  .   A  T,G    .      .  AC=1,1;AN=6;AF=0.167,0.167     GT     1|0     0/0     1/2
        3  chrX  100  .   C    A    .      .          AC=1;AN=4;AF=0.250     GT       0     0/0       1
        """
        def one_row(r, key):
            data = row_computeinfo(r, key)
            r.INFO = row_updateinfo(r, key, data, missing=True)
            return r

        df = self.df.apply(one_row, args=(key,), axis=1)

        return self.__class__(self.copy_meta(), df)

    @classmethod
    def from_dict(cls, meta, data):
        """
        Construct VcfFrame from a dict of array-like or dicts.

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
        VcfFrame.from_string
            Construct VcfFrame from a string.

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
    def from_file(
        cls, fn, compression=False, meta_only=False, regions=None
    ):
        """
        Construct VcfFrame from a VCF file.

        The method will automatically use BGZF decompression if the filename
        ends with '.gz'.

        If the file is large you can speicfy regions of interest to speed up
        data processing. Note that this requires the file be BGZF compressed
        and indexed (.tbi) for random access. Each region to be sliced must
        have the format chrom:start-end and be a half-open interval with
        (start, end]. This means, for example, 'chr1:100-103' will extract
        positions 101, 102, and 103. Alternatively, you can provide BED data
        to specify regions.

        Parameters
        ----------
        fn : str or file-like object
            VCF file (compressed or uncompressed). By file-like object, we
            refer to objects with a :meth:`read()` method, such as a file
            handle.
        compression : bool, default: False
            If True, use BGZF decompression regardless of the filename.
        meta_only : bool, default: False
            If True, only read metadata and header lines.
        regions : str, list, or pybed.BedFrame, optional
            Region or list of regions to be sliced. Also accepts a BED file
            or a BedFrame.

        Returns
        -------
        VcfFrame
            VcfFrame object.

        See Also
        --------
        VcfFrame
            VcfFrame object creation using constructor.
        VcfFrame.from_dict
            Construct VcfFrame from a dict of array-like or dicts.
        VcfFrame.from_string
            Construct VcfFrame from a string.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> vf = pyvcf.VcfFrame.from_file('unzipped.vcf')
        >>> vf = pyvcf.VcfFrame.from_file('zipped.vcf.gz')
        >>> vf = pyvcf.VcfFrame.from_file('zipped.vcf', compression=True)
        """
        if isinstance(fn, str):
            if regions is None:
                s = ''
                if fn.startswith('~'):
                    fn = os.path.expanduser(fn)
                if fn.endswith('.gz') or compression:
                    f = bgzf.open(fn, 'rt')
                else:
                    f = open(fn)
                for line in f:
                    s += line
                f.close()
            else:
                s = slice(fn, regions)
        elif hasattr(fn, 'read'):
            s = fn.read()
            try:
                s = s.decode('utf-8')
            except AttributeError:
                pass
        else:
            raise TypeError(f'Incorrect input type: {type(fn).__name__}')
        vf = cls.from_string(s)
        return vf

    @classmethod
    def from_string(cls, s, meta_only=False):
        """
        Construct VcfFrame from a string.

        Parameters
        ----------
        s : str
            String representation of a VCF file.

        Returns
        -------
        VcfFrame
            VcfFrame object.

        See Also
        --------
        VcfFrame
            VcfFrame object creation using constructor.
        VcfFrame.from_file
            Construct VcfFrame from a VCF file.
        VcfFrame.from_dict
            Construct VcfFrame from a dict of array-like or dicts.

        Examples
        --------

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
        ...     'FORMAT': ['GT', 'GT'],
        ...     'A': ['0/1', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict(['##fileformat=VCFv4.3'], data)
        >>> s = vf.to_string()
        >>> print(s[:20])
        ##fileformat=VCFv4.3
        >>> vf = pyvcf.VcfFrame.from_string(s)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1
        """
        skiprows = 0
        meta = []
        for line in s.split('\n'):
            if line.startswith('##'):
                meta.append(line.strip())
                skiprows += 1
            elif line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                skiprows += 1
            else:
                break
        columns[0] = 'CHROM'
        for header in HEADERS:
            if header not in columns and header != 'FORMAT':
                raise ValueError(f"Required VCF column missing: '{header}'")
        if meta_only:
            df = pd.DataFrame(columns=columns)
        else:
            dtype = {**HEADERS, **{x: str for x in columns[9:]}}
            df = pd.read_table(StringIO(s), skiprows=skiprows,
                names=columns, dtype=dtype)
        return cls(meta, df)

    def calculate_concordance(self, a, b, c=None, mode='all'):
        """
        Calculate genotype concordance between two (A, B) or three (A, B, C)
        samples.

        This method will return (Ab, aB, AB, ab) for comparison between two
        samples and (Abc, aBc, ABc, abC, AbC, aBC, ABC, abc) for three
        samples. Note that the former is equivalent to (FP, FN, TP, TN) if
        we assume A is the test sample and B is the truth sample.

        Only biallelic sites will be used for calculation. Additionally, the
        method will ignore zygosity and only consider presence or absence of
        variant calls (e.g. ``0/1`` and ``1/1`` will be treated the same).

        Parameters
        ----------
        a, b : str or int
            Name or index of Samples A and B.
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

        >>> vf.calculate_concordance('A', 'B', mode='all')
        (0, 1, 2, 1)
        >>> vf.calculate_concordance('A', 'B', mode='snv')
        (0, 0, 2, 1)
        >>> vf.calculate_concordance('A', 'B', mode='indel')
        (0, 1, 0, 0)

        We can also compare all three samples at once:

        >>> vf.calculate_concordance('A', 'B', 'C')
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
        """
        Combine genotype data from two samples (A, B).

        This method can be especially useful when you want to consolidate
        genotype data from replicate samples. See examples below for more
        details.

        Parameters
        ----------
        a, b : str or int
            Name or index of Samples A and B.

        Returns
        -------
        pandas.Series
            Resulting VCF column.

        See Also
        --------
        VcfFrame.subtract
            Subtract genotype data between two samples (A, B).

        Examples
        --------
        Assume we have following data where a cancer patient's tissue sample
        has been sequenced twice:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'A', 'C'],
        ...     'ALT': ['A', 'C', 'A', 'C', 'G'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT:DP', 'GT:DP', 'GT:DP', 'GT:DP', 'GT:DP'],
        ...     'Tissue1': ['./.:.', '0/0:7', '0/1:28', '0/1:4', '0/1:32'],
        ...     'Tissue2': ['0/1:24', '0/1:42', './.:.', './.:.', '0/1:19'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Tissue1 Tissue2
        0  chr1  100  .   G   A    .      .    .  GT:DP   ./.:.  0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP   0/0:7  0/1:42
        2  chr1  102  .   T   A    .      .    .  GT:DP  0/1:28   ./.:.
        3  chr1  103  .   A   C    .      .    .  GT:DP   0/1:4   ./.:.
        4  chr1  104  .   C   G    .      .    .  GT:DP  0/1:32  0/1:19

        We can combine genotype data from 'Tissue1' and 'Tissue2' to get a
        more comprehensive variant profile:

        >>> vf.df['Combined'] = vf.combine('Tissue1', 'Tissue2')
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT Tissue1 Tissue2 Combined
        0  chr1  100  .   G   A    .      .    .  GT:DP   ./.:.  0/1:24   0/1:24
        1  chr1  101  .   T   C    .      .    .  GT:DP   0/0:7  0/1:42   0/1:42
        2  chr1  102  .   T   A    .      .    .  GT:DP  0/1:28   ./.:.   0/1:28
        3  chr1  103  .   A   C    .      .    .  GT:DP   0/1:4   ./.:.    0/1:4
        4  chr1  104  .   C   G    .      .    .  GT:DP  0/1:32  0/1:19   0/1:32
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

    def to_bed(self):
        """
        Convert VcfFrame to BedFrame.

        Returns
        -------
        BedFrame
            BedFrame obejct.

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

    def to_file(self, fn, compression=False):
        """
        Write VcfFrame to a VCF file.

        If the filename ends with '.gz', the method will automatically
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
        ...     'CHROM': ['chr1', 'chr1'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT', 'GT'],
        ...     'A': ['0/1', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict(['##fileformat=VCFv4.3'], data)
        >>> print(vf.to_string())
        ##fileformat=VCFv4.3
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A
        chr1	100	.	G	A	.	.	.	GT	0/1
        chr1	101	.	T	C	.	.	.	GT	0/1
        """
        s = ''
        if self.meta:
            s += '\n'.join(self.meta) + '\n'
        s += self.df.rename(columns={'CHROM': '#CHROM'}
            ).to_csv(index=False, sep='\t')
        return s

    def to_variants(self):
        """
        List unique variants in VcfFrame.

        Returns
        -------
        list
            List of unique variants.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr2'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'A,C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ...     'FORMAT': ['GT', 'GT'],
        ...     'A': ['0/1', '1/2']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.variants()
        ['chr1-100-G-A', 'chr2-101-T-A', 'chr2-101-T-C']
        """
        def one_row(r):
            result = []
            for alt in r.ALT.split(','):
                result.append(f'{r.CHROM}-{r.POS}-{r.REF}-{alt}')
            return ','.join(result)
        s = self.df.apply(one_row, axis=1)
        s = ','.join(s)
        return s.split(',')

    def strip(self, format='GT', metadata=False):
        """
        Remove any unnecessary data.

        Parameters
        ----------
        format : str, default: 'GT'
            FORMAT keys to retain (e.g. 'GT:AD:DP').
        metadata : bool, default: False
            If True, keep the metadata.

        Returns
        -------
        VcfFrame
            Stripped VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'T', 'A'],
        ...     'ALT': ['A', 'C', 'T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT:DP:AD', 'GT:DP:AD', 'GT'],
        ...     'A': ['0/1:30:15,15', '1/1:28:0,28', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO    FORMAT             A
        0  chr1  100  .   G   A    .      .    .  GT:DP:AD  0/1:30:15,15
        1  chr1  101  .   T   C    .      .    .  GT:DP:AD   1/1:28:0,28
        2  chr1  102  .   A   T    .      .    .        GT           0/1
        >>> vf.strip('GT:DP').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT       A
        0  chr1  100  .   G   A    .      .    .  GT:DP  0/1:30
        1  chr1  101  .   T   C    .      .    .  GT:DP  1/1:28
        2  chr1  102  .   A   T    .      .    .  GT:DP   0/1:.
        """
        new_keys = format.split(':')

        def one_row(r):
            old_keys = r.FORMAT.split(':')
            indicies = [
                old_keys.index(x) if x in old_keys else None for x in new_keys
            ]

            def one_gt(g):
                old_fields = g.split(':')
                new_fields = [
                    '.' if x is None
                    else '.' if x >= len(old_fields)
                    else old_fields[x]
                    for x in indicies
                ]
                return ':'.join(new_fields)

            r.iloc[9:] = r.iloc[9:].apply(one_gt)

            return r

        df = self.df.copy()
        df[['ID', 'QUAL', 'FILTER', 'INFO']] = '.'
        df = df.apply(one_row, axis=1)
        df.FORMAT = format

        if metadata:
            meta = self.copy_meta()
        else:
            meta = []
        vf = self.__class__(meta, df)
        return vf

    def merge(
        self, other, how='inner', format='GT', sort=True, collapse=False
    ):
        """
        Merge with the other VcfFrame.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame. Note that the 'chr' prefix in contig names (e.g.
            'chr1' vs. '1') will be automatically added or removed as
            necessary to match the contig names of ``self``.
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
        if self.has_chr_prefix and other.has_chr_prefix:
            pass
        elif self.has_chr_prefix and not other.has_chr_prefix:
            other = other.update_chr_prefix('add')
        elif not self.has_chr_prefix and other.has_chr_prefix:
            other = other.update_chr_prefix('remove')
        else:
            pass

        if self.sites_only and other.sites_only:
            df = pd.concat([self.df, other.df])
            merged = self.__class__([], df)
            merged = merged.drop_duplicates()
        else:
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
            merged = self.__class__([], df)

        if collapse:
            merged = merged.collapse()
        if sort:
            merged = merged.sort()

        return merged

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

    def plot_region(
        self, sample, k='#DP', color=None, region=None, label=None, ax=None,
        figsize=None, **kwargs
    ):
        """
        Create a scatter plot showing read depth profile of a sample for
        the specified region.

        Parameters
        ----------
        sample : str or int
            Name or index of target sample.
        k : str, default: '#DP'
            Genotype key to use for extracting data:

            - '#DP': Return read depth.
            - '#AD_REF': Return REF allele depth.
            - '#AD_ALT': Return ALT allele depth.
            - '#AD_FRAC_REF': Return REF allele fraction.
            - '#AD_FRAC_ALT': Return ALT allele fraction.

        color : str, optional
            Marker color.
        region : str, optional
            Target region ('chrom:start-end').
        label : str, optional
            Label to use for the data points.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.scatter`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> from fuc import pyvcf, common
            >>> import matplotlib.pyplot as plt
            >>> common.load_dataset('pyvcf')
            >>> vcf_file = '~/fuc-data/pyvcf/getrm-cyp2d6-vdr.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_region('NA18973')
            >>> plt.tight_layout()

        We can display allele fraction of REF and ALT instead of DP:

        .. plot::
            :context: close-figs

            >>> ax = vf.plot_region('NA18973', k='#AD_FRAC_REF', label='REF')
            >>> vf.plot_region('NA18973', k='#AD_FRAC_ALT', label='ALT', ax=ax)
            >>> plt.tight_layout()
        """
        if self.df.empty:
            raise ValueError('VcfFrame is empty')

        sample = sample if isinstance(sample, str) else self.samples[sample]

        if region is None:
            if len(self.contigs) == 1:
                vf = self.copy()
            else:
                raise ValueError('Multiple contigs found.')
        else:
            vf = self.slice(region)

        df = vf.extract_format(k)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            x=vf.df.POS, y=df[sample], c=color, label=label, **kwargs
        )

        ax.set_xlabel('Position')
        ax.set_ylabel('Depth')

        return ax

    def plot_comparison(
        self, a, b, c=None, labels=None, ax=None, figsize=None
    ):
        """
        Create a Venn diagram showing genotype concordance between groups.

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
            n = [x + y for x, y in zip(n, self.calculate_concordance(a[i], b[i]))]
        out = venn2(subsets=n[:-1], **venn_kws)
        return out

    def _plot_comparison_three(self, a, b, c, venn_kws):
        n = [0, 0, 0, 0, 0, 0, 0, 0]
        for i in range(len(a)):
            n = [x + y for x, y in zip(n, self.calculate_concordance(a[i], b[i], c[i]))]
        out = venn3(subsets=n[:-1], **venn_kws)
        return out

    def plot_hist_format(
        self, k, af=None, group_col=None, group_order=None, kde=True,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a histogram showing the distribution of data for the
        specified FORMAT key.

        Parameters
        ----------
        k : str
            One of the special FORMAT keys as defined in
            :meth:`VcfFrame.extract_format`.
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
            >>> vf.plot_hist_format('#DP')

        We can draw multiple histograms with hue mapping:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/pyvcf/normal-tumor-annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col='Sample')
            >>> vf.plot_hist_format('#DP', af=af, group_col='Tissue')

        We can show AF instead of DP:

        .. plot::
            :context: close-figs

            >>> vf.plot_hist_format('#AD_FRAC_REF')
        """
        if k not in FORMAT_SPECIAL_KEYS:
            raise ValueError('Incorrect FORMAT key.')

        df = self.extract_format(k)

        df = df.T
        id_vars = ['index']
        if group_col is not None:
            df = pd.concat([df, af.df[group_col]], axis=1, join='inner')
            id_vars.append(group_col)
        df = df.reset_index()
        df = pd.melt(df, id_vars=id_vars)
        df = df.dropna()
        df = df.rename(columns={'value': k})

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.histplot(
            data=df, x=k, hue=group_col, hue_order=group_order, kde=kde,
            ax=ax, **kwargs
        )

        return ax

    def plot_hist_info(
        self, k, kde=True, ax=None, figsize=None, **kwargs
    ):
        """
        Create a histogram showing the distribution of data for the
        specified INFO key.

        Parameters
        ----------
        k : str
            One of the special INFO keys as defined in
            :meth:`VcfFrame.extract_info`.
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
            >>> vcf_file = '~/fuc-data/pyvcf/getrm-cyp2d6-vdr.vcf'
            >>> vf = pyvcf.VcfFrame.from_file(vcf_file)
            >>> vf.plot_hist_info('#AC')

        We can show AF instead of AC:

        .. plot::
            :context: close-figs

            >>> vf.plot_hist_info('#AF')
        """
        if k not in INFO_SPECIAL_KEYS:
            raise ValueError('Incorrect INFO key.')

        s = self.extract_info(k)

        if s.isnull().all():
            raise ValueError('Dataset is empty.')

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.histplot(
            data=s, kde=kde, ax=ax, **kwargs
        )

        ax.set_xlabel(k)

        return ax

    def plot_tmb(
        self, af=None, group_col=None, group_order=None, kde=True, ax=None,
        figsize=None, **kwargs
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

    def empty_samples(self, threshold=0, opposite=False, as_list=False):
        """
        Remove samples with high missingness.

        Samples with missingness >= threshold will be removed.

        Parameters
        ----------
        threshold : int or float, default: 0
            Number or fraction of missing variants. By default
            (``threshold=0``), only samples with 100% missingness will be
            removed.
        opposite : bool, default: False
            If True, return samples that don't meet the said criteria.
        as_list : bool, default: False
             If True, return a list of sample names instead of a VcfFrame.

        Returns
        -------
        VcfFrame
            Subsetted VcfFrame.

        Examples
        --------

        >>> from fuc import pyvcf
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
        ...     'A': ['0/0', '0/0', '0/0', '0/0'],
        ...     'B': ['./.', '0/0', '0/0', '0/0'],
        ...     'C': ['./.', './.', '0/0', '0/0'],
        ...     'D': ['./.', './.', './.', '0/0'],
        ...     'E': ['./.', './.', './.', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D    E
        0  chr1  100  .   G   A    .      .    .     GT  0/0  ./.  ./.  ./.  ./.
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/0  ./.  ./.  ./.
        2  chr1  102  .   G   C    .      .    .     GT  0/0  0/0  0/0  ./.  ./.
        3  chr1  103  .   T   C    .      .    .     GT  0/0  0/0  0/0  0/0  ./.
        >>> vf.empty_samples().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/0  ./.  ./.  ./.
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/0  ./.  ./.
        2  chr1  102  .   G   C    .      .    .     GT  0/0  0/0  0/0  ./.
        3  chr1  103  .   T   C    .      .    .     GT  0/0  0/0  0/0  0/0
        >>> vf.empty_samples(threshold=2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/0  ./.
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/0
        2  chr1  102  .   G   C    .      .    .     GT  0/0  0/0
        3  chr1  103  .   T   C    .      .    .     GT  0/0  0/0
        >>> vf.empty_samples(threshold=0.5).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/0  ./.
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/0
        2  chr1  102  .   G   C    .      .    .     GT  0/0  0/0
        3  chr1  103  .   T   C    .      .    .     GT  0/0  0/0
        >>> vf.empty_samples(threshold=0.5, opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    C    D    E
        0  chr1  100  .   G   A    .      .    .     GT  ./.  ./.  ./.
        1  chr1  101  .   T   C    .      .    .     GT  ./.  ./.  ./.
        2  chr1  102  .   G   C    .      .    .     GT  0/0  ./.  ./.
        3  chr1  103  .   T   C    .      .    .     GT  0/0  0/0  ./.
        >>> vf.empty_samples(threshold=0.5, opposite=True, as_list=True)
        ['C', 'D', 'E']
        """
        total = self.shape[0]

        if not threshold:
            threshold = total

        if isinstance(threshold, int):
            use_fraction = False
        elif isinstance(threshold, float):
            use_fraction = True
        else:
            raise TypeError("Incorrect argument type 'threshold'")

        def one_col(c):
            data = sum(c.apply(gt_miss))
            if use_fraction:
                data /= total
            return data

        s = self.df.iloc[:, 9:].apply(one_col, axis=0) < threshold

        if opposite:
            s = s[s == False]
        else:
            s = s[s == True]

        l = s.index.to_list()

        if as_list:
            return l

        return self.subset(l)

    def expand(self):
        """
        Expand each multiallelic locus to multiple rows.

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
        Filter rows intersecting with given BED.

        Only variants intersecting with given BED data will remain.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_empty(self, threshold=0, opposite=False, as_index=False):
        """
        Filter rows with high missingness.

        Variants with missingness >= threshold will be removed.

        Parameters
        ----------
        threshold: int, default: 0
            Exclude the row if it has a number of missing genotypes that is
            greater than or equal to this number. When 0 (default), exclude
            rows where all of the samples have a missing genotype.
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
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'A', 'C', 'C'],
        ...     'ALT': ['A', 'C', 'T', 'A', 'T'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', './.', './.', './.', './.'],
        ...     'B': ['0/0', '0/1', './.', './.', './.'],
        ...     'C': ['0/0', '0/0', '0/1', './.', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0  0/0
        1  chr1  101  .   T   C    .      .    .     GT  ./.  0/1  0/0
        2  chr1  102  .   A   T    .      .    .     GT  ./.  ./.  0/1
        3  chr1  103  .   C   A    .      .    .     GT  ./.  ./.  ./.
        4  chr1  104  .   C   T    .      .    .     GT  ./.  ./.  ./.

        We can remove rows that are completely empty:

        >>> vf.filter_empty().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0  0/0
        1  chr1  101  .   T   C    .      .    .     GT  ./.  0/1  0/0
        2  chr1  102  .   A   T    .      .    .     GT  ./.  ./.  0/1

        We can remove rows where at least two samples have missing genotype:

        >>> vf.filter_empty(threshold=2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0  0/0
        1  chr1  101  .   T   C    .      .    .     GT  ./.  0/1  0/0

        We can show rows that are completely empty:

        >>> vf.filter_empty(opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C
        0  chr1  103  .   C   A    .      .    .     GT  ./.  ./.  ./.
        1  chr1  104  .   C   T    .      .    .     GT  ./.  ./.  ./.

        Finally, we can return boolean index array from the filtering:

        >>> vf.filter_empty(as_index=True)
        0     True
        1     True
        2     True
        3    False
        4    False
        dtype: bool
        """
        if not threshold:
            threshold = self.shape[1]

        def one_row(r):
            s = r[9:].apply(gt_miss)
            return s.sum() < threshold

        i = self.df.apply(one_row, axis=1)

        if opposite:
            i = ~i
        if as_index:
            return i
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagall(self, flags, opposite=False, as_index=False):
        """
        Filter rows with all given INFO flags.

        Only variants with all given INFO flags will remain.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_flagany(self, flags, opposite=False, as_index=False):
        """
        Filter rows with any given INFO flags.

        Only variants with any given INFO flags will remain.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_indel(self, opposite=False, as_index=False):
        """
        Filter rows with indel.

        Variants with indel will be removed.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_multialt(self, opposite=False, as_index=False):
        """
        Filter rows with multiple ALT alleles.

        Variants with multiple ALT alleles will be removed.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_pass(self, opposite=False, as_index=False):
        """
        Filter rows with PASS in FILTER column.

        Only variants with PASS in the FILTER column will remain.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_phased(self, opposite=False, as_index=False):
        """
        Filter rows with phased genotypes.

        Variants with phased genotypes will be removed.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_polyp(self, opposite=False, as_index=False):
        """
        Filter rows with polyploid genotypes.

        Variants with polyploid genotypes will be removed.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_qual(self, threshold, opposite=False, as_index=False):
        """
        Filter rows with low QUAL values.

        Only variants with QUAL >= threashold will remain.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampall(self, samples=None, opposite=False, as_index=False):
        """
        Filter rows where all given samples have variant.

        Only variants where all given samples have variant. The default
        behavior is to use all samples in the VcfFrame.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampany(self, samples=None, opposite=False, as_index=False):
        """
        Filter rows where any given samples have variant.

        Only variants where any given samples have variant will remain. The
        default behavior is to use all samples in the VcfFrame.

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
        if as_index:
            return i
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_sampnum(self, threshold, opposite=False, as_index=False):
        """
        Filter rows with high variant prevalence.

        Only variants with variant prevalence >= threshold will remian.

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
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def filter_vcf(self, vcf, opposite=False, as_index=False):
        """
        Filter rows intersecting with given VCF.

        Only variants intersecting with given VCF data will remain.

        Parameters
        ----------
        vcf : VcfFrame or str
            VcfFrame or VCF file.
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
        >>> data1 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr4', 'chr8', 'chr8'],
        ...     'POS': [100, 203, 192, 52, 788],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['A', 'C', 'T', 'T', 'GA'],
        ...     'ALT': ['C', 'G', 'A', 'G', 'G'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', '0/1', '0/1', '0/1', '0/1'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> vf1.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   A   C    .      .    .     GT  0/1
        1  chr1  203  .   C   G    .      .    .     GT  0/1
        2  chr4  192  .   T   A    .      .    .     GT  0/1
        3  chr8   52  .   T   G    .      .    .     GT  0/1
        4  chr8  788  .  GA   G    .      .    .     GT  0/1
        >>> data2 = {
        ...     'CHROM': ['chr1', 'chr8'],
        ...     'POS': [100, 788],
        ...     'ID': ['.', '.'],
        ...     'REF': ['A', 'GA'],
        ...     'ALT': ['C', 'G'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['.', '.'],
        ... }
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> vf2.df
          CHROM  POS ID REF ALT QUAL FILTER INFO
        0  chr1  100  .   A   C    .      .    .
        1  chr8  788  .  GA   G    .      .    .

        We can select rows that overlap with the VCF data:

        >>> vf1.filter_vcf(vf2).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   A   C    .      .    .     GT  0/1
        1  chr8  788  .  GA   G    .      .    .     GT  0/1

        We can also remove those rows:

        >>> vf1.filter_vcf(vf2, opposite=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  203  .   C   G    .      .    .     GT  0/1
        1  chr4  192  .   T   A    .      .    .     GT  0/1
        2  chr8   52  .   T   G    .      .    .     GT  0/1

        Finally, we can return boolean index array from the filtering:

        >>> vf1.filter_vcf(vf2, as_index=True)
        0     True
        1    False
        2    False
        3    False
        4     True
        dtype: bool
        """
        if isinstance(vcf, VcfFrame):
            vf = vcf
        else:
            vf = VcfFrame.from_file(vcf)
        df1 = self.df[['CHROM', 'POS', 'REF', 'ALT']]
        df2 = vf.df[['CHROM', 'POS', 'REF', 'ALT', 'ID']]
        df3 = df1.merge(df2, how='left')
        i = ~pd.isna(df3.ID)
        i.name = None
        if opposite:
            i = ~i
        if as_index:
            return i
        if i.empty:
            return self.__class__(self.copy_meta(), self.copy_df())
        return self.__class__(self.copy_meta(), self.df[i])

    def subtract(self, a, b):
        """
        Subtract genotype data between two samples (A, B).

        This method can be especially useful when you want to distinguish
        between somatic and germline variants for an individual. See examples
        below for more details.

        Parameters
        ----------
        a, b : str or int
            Name or index of Samples A and B.

        Returns
        -------
        pandas.Series
            Resulting VCF column.

        See Also
        --------
        VcfFrame.combine
            Combine genotype data from two samples (A, B).

        Examples
        --------
        Assume we have following data for a cancer patient:

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
        ...     'REF': ['G', 'T', 'C', 'A', 'G'],
        ...     'ALT': ['A', 'G', 'T', 'C', 'C'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'Tissue': ['./.', '0/1', '0/1', '0/0', '0/1'],
        ...     'Blood': ['0/1', '0/1', './.', '0/1', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS   ID REF ALT QUAL FILTER INFO FORMAT Tissue Blood
        0  chr1  100  rs1   G   A    .      .    .     GT    ./.   0/1
        1  chr1  101  rs2   T   G    .      .    .     GT    0/1   0/1
        2  chr1  102  rs3   C   T    .      .    .     GT    0/1   ./.
        3  chr1  103  rs4   A   C    .      .    .     GT    0/0   0/1
        4  chr1  104  rs5   G   C    .      .    .     GT    0/1   0/0

        We can compare genotype data between 'Tissue' and 'Blood' to identify
        somatic variants (i.e. rs3 and rs5; rs2 is most likely germline):

        >>> vf.df['Somatic'] = vf.subtract('Tissue', 'Blood')
        >>> vf.df
          CHROM  POS   ID REF ALT QUAL FILTER INFO FORMAT Tissue Blood Somatic
        0  chr1  100  rs1   G   A    .      .    .     GT    ./.   0/1     ./.
        1  chr1  101  rs2   T   G    .      .    .     GT    0/1   0/1     ./.
        2  chr1  102  rs3   C   T    .      .    .     GT    0/1   ./.     0/1
        3  chr1  103  rs4   A   C    .      .    .     GT    0/0   0/1     0/0
        4  chr1  104  rs5   G   C    .      .    .     GT    0/1   0/0     0/1
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
        """
        Sort the VcfFrame by chromosome and position.

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
        Subset VcfFrame for specified samples.

        Parameters
        ----------
        samples : str or list
            Sample name or list of names (the order matters).
        exclude : bool, default: False
            If True, exclude specified samples.

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
        ...     'FORMAT': ['GT', 'GT'],
        ...     'A': ['0/1', '0/1'],
        ...     'B': ['0/0', '0/1'],
        ...     'C': ['0/0', '0/0'],
        ...     'D': ['0/1', '0/0'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0  0/0  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1  0/1  0/0  0/0

        We can subset the VcfFrame for the samples A and B:

        >>> vf.subset(['A', 'B']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0
        1  chr1  101  .   T   C    .      .    .     GT  0/1  0/1

        Alternatively, we can exclude those samples:

        >>> vf.subset(['A', 'B'], exclude=True).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/0  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/0  0/0
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
        Slice VcfFrame for specified region.

        Parameters
        ----------
        region : str
            Region to slice ('chrom:start-end'). The 'chr' string will be
            handled automatically. For example, 'chr1' and '1' will have the
            same effect in slicing.

        Returns
        -------
        VcfFrame
            Sliced VcfFrame.

        Examples
        --------

        Assume we have following data:

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
        ...     'A': ['0/1', '1/1', '0/0', '0/0'],
        ...     'B': ['0/0', '0/0', '0/1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0
        1  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        2  chr1  297  .   A   T    .      .    .     GT  0/0  0/1
        3  chr2  101  .   C   A    .      .    .     GT  0/0  0/1

        We can specify a full region:

        >>> vf.slice('chr1:101-300').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        1  chr1  297  .   A   T    .      .    .     GT  0/0  0/1

        We can specify a region without the 'chr' string:

        >>> vf.slice('1:101-300').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        1  chr1  297  .   A   T    .      .    .     GT  0/0  0/1

        We can specify the contig only:

        >>> vf.slice('chr1').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0
        1  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        2  chr1  297  .   A   T    .      .    .     GT  0/0  0/1

        We can omit the start position:

        >>> vf.slice('chr1:-296').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/0
        1  chr1  205  .   T   C    .      .    .     GT  1/1  0/0

        We can omit the end position as well:

        >>> vf.slice('chr1:101').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        1  chr1  297  .   A   T    .      .    .     GT  0/0  0/1

        You can also omit the end position this way:

        >>> vf.slice('chr1:101-').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B
        0  chr1  205  .   T   C    .      .    .     GT  1/1  0/0
        1  chr1  297  .   A   T    .      .    .     GT  0/0  0/1
        """
        if self.has_chr_prefix:
            region = common.update_chr_prefix(region, mode='add')
        else:
            region = common.update_chr_prefix(region, mode='remove')
        chrom, start, end = common.parse_region(region)
        df = self.df[self.df.CHROM == chrom]
        if not pd.isna(start):
            df = df[df.POS >= start]
        if not pd.isna(end):
            df = df[df.POS <= end]
        return self.__class__(self.copy_meta(), df)

    def extract_format(self, k, func=None, as_nan=False):
        """
        Extract data for the specified FORMAT key.

        By default, this method will return string data. Use ``func`` and
        ``as_nan`` to output numbers. Alternatvely, select one of the special
        keys for ``k``, which have predetermined values of ``func`` and
        ``as_nan`` for convenience.

        Parameters
        ----------
        k : str
            FORMAT key to use when extracting data. In addition to regular
            FORMAT keys (e.g. 'DP', 'AD'), the method also accepts the
            special keys listed below:

            - '#DP': Return numeric DP.
            - '#AD_REF': Return numeric AD for REF.
            - '#AD_ALT': Return numeric AD for ALT. If multiple values are
              available (i.e. multiallelic site), return the sum.
            - '#AD_FRAC_REF': Return allele fraction for REF.
            - '#AD_FRAC_ALT': Return allele fraction for ALT. If multiple
              values are available (i.e. multiallelic site), return the sum.

        func : function, optional
            Function to apply to each of the extracted results.
        as_nan : bool, default: False
            If True, return missing values as ``NaN``.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing requested data.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['A', 'C', 'A'],
        ...     'ALT': ['G', 'T', 'C,T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT:AD:DP', 'GT', 'GT:AD:DP'],
        ...     'A': ['0/1:15,13:28', '0/0', '0/1:9,14,0:23'],
        ...     'B': ['./.:.:.', '1/1', '1/2:0,11,15:26'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO    FORMAT              A               B
        0  chr1  100  .   A    G    .      .    .  GT:AD:DP   0/1:15,13:28         ./.:.:.
        1  chr1  101  .   C    T    .      .    .        GT            0/0             1/1
        2  chr1  102  .   A  C,T    .      .    .  GT:AD:DP  0/1:9,14,0:23  1/2:0,11,15:26
        >>> vf.extract_format('GT')
             A    B
        0  0/1  ./.
        1  0/0  1/1
        2  0/1  1/2
        >>> vf.extract_format('GT', as_nan=True)
             A    B
        0  0/1  NaN
        1  0/0  1/1
        2  0/1  1/2
        >>> vf.extract_format('AD')
                A        B
        0   15,13        .
        1     NaN      NaN
        2  9,14,0  0,11,15
        >>> vf.extract_format('DP', func=lambda x: int(x), as_nan=True)
              A     B
        0  28.0   NaN
        1   NaN   NaN
        2  23.0  26.0
        >>> vf.extract_format('#DP') # Same as above
              A     B
        0  28.0   NaN
        1   NaN   NaN
        2  23.0  26.0
        >>> vf.extract_format('AD', func=lambda x: float(x.split(',')[0]), as_nan=True)
              A    B
        0  15.0  NaN
        1   NaN  NaN
        2   9.0  0.0
        >>> vf.extract_format('#AD_REF') # Same as above
              A    B
        0  15.0  NaN
        1   NaN  NaN
        2   9.0  0.0
        """
        def one_row(r, k, func, as_nan):
            try:
                i = r.FORMAT.split(':').index(k)
            except ValueError:
                return pd.Series([np.nan] * len(self.samples), index=self.samples)
            def one_gt(g):
                if k == 'GT' and gt_miss(g) and as_nan:
                    return np.nan
                v = g.split(':')[i]
                if v == '.' and as_nan:
                    return np.nan
                if func is not None:
                    return func(v)
                return v
            return r[9:].apply(one_gt)

        if k in FORMAT_SPECIAL_KEYS:
            k, func, as_nan = FORMAT_SPECIAL_KEYS[k]

        df = self.df.apply(one_row, args=(k, func, as_nan), axis=1)

        return df

    def extract_info(vf, k, func=None, as_nan=False):
        """
        Extract data for the specified INFO key.

        By default, this method will return string data. Use ``func`` and
        ``as_nan`` to output numbers. Alternatvely, select one of the special
        keys for ``k``, which have predetermined values of ``func`` and
        ``as_nan`` for convenience.

        Parameters
        ----------
        k : str
            INFO key to use when extracting data. In addition to regular
            INFO keys (e.g. 'AC', 'AF'), the method also accepts the
            special keys listed below:

            - '#AC': Return numeric AC. If multiple values are available
              (i.e. multiallelic site), return the sum.
            - '#AF': Similar to '#AC'.

        func : function, optional
            Function to apply to each of the extracted results.
        as_nan : bool, default: False
            If True, return missing values as ``NaN``.

        Returns
        -------
        pandas.Series
            Requested data.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'C', 'A', 'A'],
        ...     'ALT': ['G', 'T', 'C,T', 'T'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['AC=1;AF=0.167;H2', 'AC=2;AF=0.333', 'AC=1,2;AF=0.167,0.333;H2', 'AC=.;AF=.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', '0/0', '0/1', './.'],
        ...     'B': ['0/0', '1/1', '0/2', './.'],
        ...     'C': ['0/0', '0/0', '0/2', './.'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER                      INFO FORMAT    A    B    C
        0  chr1  100  .   A    G    .      .          AC=1;AF=0.167;H2     GT  0/1  0/0  0/0
        1  chr1  101  .   C    T    .      .             AC=2;AF=0.333     GT  0/0  1/1  0/0
        2  chr1  102  .   A  C,T    .      .  AC=1,2;AF=0.167,0.333;H2     GT  0/1  0/2  0/2
        3  chr1  103  .   A    T    .      .                 AC=.;AF=.     GT  ./.  ./.  ./.
        >>> vf.extract_info('H2')
        0     H2
        1    NaN
        2     H2
        3    NaN
        dtype: object
        >>> vf.extract_info('AC')
        0      1
        1      2
        2    1,2
        3      .
        dtype: object
        >>> vf.extract_info('AC', as_nan=True)
        0      1
        1      2
        2    1,2
        3    NaN
        dtype: object
        >>> vf.extract_info('AC', func=lambda x: sum([int(x) for x in x.split(',')]), as_nan=True)
        0    1.0
        1    2.0
        2    3.0
        3    NaN
        dtype: float64
        >>> vf.extract_info('#AC') # Same as above
        0    1.0
        1    2.0
        2    3.0
        3    NaN
        dtype: float64
        """
        def one_row(r, k, func, as_nan):
            d = {}

            for field in r.INFO.split(';'):
                if '=' not in field:
                    d[field] = field
                else:
                    d[field.split('=')[0]] = field.split('=')[1]

            try:
                result = d[k]
            except KeyError:
                return np.nan

            if result == '.' and as_nan:
                return np.nan

            if func is not None:
                return func(result)

            return result

        if k in INFO_SPECIAL_KEYS:
            k, func, as_nan = INFO_SPECIAL_KEYS[k]

        s = vf.df.apply(one_row, args=(k, func, as_nan), axis=1)

        return s

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
        ...     'D': ['0/1', '0/1'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    B    C    D
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1  0/1
        >>> vf.rename(['1', '2', '3', '4']).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    1    2    3    4
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1  0/1
        >>> vf.rename({'B': '2', 'C': '3'}).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    2    3    D
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1  0/1
        >>> vf.rename(['2', '4'], indicies=[1, 3]).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    2    C    4
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1  0/1
        >>> vf.rename(['2', '3'], indicies=(1, 3)).df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A    2    3    D
        0  chr1  100  .   G   A    .      .    .     GT  0/1  0/1  0/1  0/1
        1  chr2  101  .   T   C    .      .    .     GT  0/1  0/1  0/1  0/1
        """
        samples = common.rename(self.samples, names, indicies=indicies)
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

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        mf.plot_snvclss(
            color=color, colormap=colormap, width=width, legend=legend,
            flip=flip, af=af, **kwargs
        )

        return ax

    def update_chr_prefix(self, mode='remove'):
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
        ...     'CHROM': ['chr1', 'chr1', '2', '2'],
        ...     'POS': [100, 101, 100, 101],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['G', 'T', 'T', 'C'],
        ...     'ALT': ['A', 'C', 'C', 'G'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
        ...     'A': ['0/1', '0/1', '0/1', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1
        2     2  100  .   T   C    .      .    .     GT  0/1
        3     2  101  .   C   G    .      .    .     GT  0/1
        >>> vf.update_chr_prefix(mode='remove').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0     1  100  .   G   A    .      .    .     GT  0/1
        1     1  101  .   T   C    .      .    .     GT  0/1
        2     2  100  .   T   C    .      .    .     GT  0/1
        3     2  101  .   C   G    .      .    .     GT  0/1
        >>> vf.update_chr_prefix(mode='add').df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0/1
        1  chr1  101  .   T   C    .      .    .     GT  0/1
        2  chr2  100  .   T   C    .      .    .     GT  0/1
        3  chr2  101  .   C   G    .      .    .     GT  0/1
        """
        if mode == 'remove':
            def one_row(r):
                r.CHROM = r.CHROM.replace('chr', '')
                return r
        elif mode == 'add':
            def one_row(r):
                if 'chr' not in r.CHROM:
                    r.CHROM = 'chr' + r.CHROM
                return r
        else:
            raise ValueError(f'Incorrect mode: {mode}')
        df = self.df.apply(one_row, axis=1)
        return self.__class__(self.copy_meta(), df)

    def compare(self, other):
        """
        Compare to another VcfFrame and show the differences in genotype
        calling.

        Parameters
        ----------
        other : VcfFrame
            VcfFrame to compare with.

        Returns
        -------
        pandas.DataFrame
            DataFrame comtaining genotype differences.

        Examples
        --------

        >>> from fuc import pyvcf
        >>> data1 = {
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
        >>> data2 = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103, 104],
        ...     'ID': ['.', '.', '.', '.', '.'],
        ...     'REF': ['G', 'CT', 'T', 'C', 'A'],
        ...     'ALT': ['A', 'C', 'A', 'T', 'G,C'],
        ...     'QUAL': ['.', '.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
        ...     'A': ['./.', '0/0', '0/0', '0/1', '0/0'],
        ...     'B': ['1/1', '0/1', './.', '1/1', '0/0'],
        ...     'C': ['0/1', '0/1', '0/1', './.', '1/2'],
        ... }
        >>> vf1 = pyvcf.VcfFrame.from_dict([], data1)
        >>> vf2 = pyvcf.VcfFrame.from_dict([], data2)
        >>> vf1.compare(vf2)
                  Locus Sample Self Other
        0  chr1-100-G-A      A  0/1   ./.
        1  chr1-102-T-A      C  1/1   0/1
        2  chr1-103-C-T      B  0/1   1/1
        """
        df1 = self.df.copy()
        df2 = other.df.copy()
        i1 = df1['CHROM'] + '-' + df1['POS'].astype(str) + '-' + df1['REF'] + '-' + df1['ALT']
        i2 = df2['CHROM'] + '-' + df2['POS'].astype(str) + '-' + df2['REF'] + '-' + df2['ALT']
        df1 = df1.set_index(i1)
        df2 = df2.set_index(i2)
        df1 = df1.iloc[:, 9:]
        df2 = df2.iloc[:, 9:]
        if df1.equals(df2):
            return pd.DataFrame(columns=['Locus', 'Sample', 'Self', 'Other'])
        df = df1.compare(df2, align_axis=0)
        df = df.stack().to_frame().reset_index().pivot(
            columns='level_1', index=['level_0', 'level_2'], values=0)
        df = df.reset_index()
        df.columns = ['Locus', 'Sample', 'Other', 'Self']
        df = df[['Locus', 'Sample', 'Self', 'Other']]
        return df

    def fetch(self, variant):
        """
        Fetch the VCF row that matches specified variant.

        Parameters
        ----------
        variant : str
            Target variant.

        Returns
        -------
        pandas.Series
            VCF row.

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
        >>> vf.fetch('chr1-100-G-A')
        CHROM     chr1
        POS        100
        ID           .
        REF          G
        ALT          A
        QUAL         .
        FILTER       .
        INFO         .
        FORMAT      GT
        A          0/1
        Name: 0, dtype: object
        """
        def one_row(r):
            return f'{r.CHROM}-{r.POS}-{r.REF}-{r.ALT}' == variant
        i = self.df.apply(one_row, axis=1)
        return self.df[i].squeeze()

    def pseudophase(self):
        """
        Pseudophase VcfFrame.

        Returns
        -------
        VcfFrame
            Pseudophased VcfFrame.

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
        >>> vf.pseudophase().df
          CHROM  POS ID REF ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .   G   A    .      .    .     GT  0|1
        1  chr2  101  .   T   C    .      .    .     GT  1|1
        """
        def one_row(r):
            r[9:] = r[9:].apply(pseudophase)
            return r
        df = self.df.apply(one_row, axis=1)
        return self.__class__(self.copy_meta(), df)

    def get_af(self, sample, variant):
        """
        Get allele fraction for a pair of sample and variant.

        The method will return ``numpy.nan`` if the value is missing.

        Parameters
        ----------
        sample : str
            Sample name.
        variant : str
            Variant name.

        Returns
        -------
        float
            Allele fraction.

        Examples
        --------

        >>> from fuc import pyvcf, common
        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 101, 102, 103],
        ...     'ID': ['.', '.', '.', '.'],
        ...     'REF': ['A', 'G', 'A', 'C'],
        ...     'ALT': ['C', 'T', 'G', 'G,A'],
        ...     'QUAL': ['.', '.', '.', '.'],
        ...     'FILTER': ['.', '.', '.', '.'],
        ...     'INFO': ['.', '.', '.', '.'],
        ...     'FORMAT': ['GT:AD:AF', 'GT:AD:AF', 'GT:AF', 'GT:AD:AF'],
        ...     'A': ['0/1:12,15:0.444,0.556', '0/0:32,1:0.970,0.030', '0/1:.', './.:.:.'],
        ...     'B': ['0/1:13,17:0.433,0.567', '0/1:14,15:0.483,0.517', './.:.', '1/2:0,11,17:0.000,0.393,0.607'],
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF  ALT QUAL FILTER INFO    FORMAT                      A                              B
        0  chr1  100  .   A    C    .      .    .  GT:AD:AF  0/1:12,15:0.444,0.556          0/1:13,17:0.433,0.567
        1  chr1  101  .   G    T    .      .    .  GT:AD:AF   0/0:32,1:0.970,0.030          0/1:14,15:0.483,0.517
        2  chr1  102  .   A    G    .      .    .     GT:AF                  0/1:.                          ./.:.
        3  chr1  103  .   C  G,A    .      .    .  GT:AD:AF                ./.:.:.  1/2:0,11,17:0.000,0.393,0.607
        >>> vf.get_af('A', 'chr1-100-A-C')
        0.556
        >>> vf.get_af('B', 'chr1-102-A-G')
        nan
        """
        chrom, pos, ref, alt = common.parse_variant(variant)
        r = self.df[(self.df.CHROM == chrom) & (self.df.POS == pos) & (self.df.REF == ref)]

        try:
            i = r.FORMAT.values[0].split(':').index('AF')
        except ValueError:
            return np.nan

        alts = r.ALT.values[0].split(',')

        if alt in alts:
            j = r.ALT.values[0].split(',').index(alt)
        else:
            raise ValueError(f'ALT allele not found, possible choices: {alts}')

        field = r[sample].values[0].split(':')[i]

        if field == '.':
            af = np.nan
        else:
            af = float(field.split(',')[j+1])

        return af
