"""
The pybam submodule is designed for working with sequence alignment files
(SAM/BAM/CRAM). It essentially wraps the `pysam
<https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast
computation and easy manipulation. If you are mainly interested in working
with depth of coverage data, please check out the pycov submodule which is
specifically designed for the task.
"""

import sys

from . import common, pybed

import pysam
import pandas as pd

def slice(bam, regions, format='BAM', path=None, fasta=None):
    """
    Slice a BAM file for specified regions.

    Parameters
    ----------
    bam : str
        Input BAM file. It must be already indexed to allow random access.
        You can index a BAM file with the :meth:`pybam.index` method.
    regions : str, list, or pybed.BedFrame
        One or more regions to be sliced. Each region must have the format
        chrom:start-end and be a half-open interval with (start, end]. This
        means, for example, chr1:100-103 will extract positions 101, 102, and
        103. Alternatively, you can provide a BED file (compressed or
        uncompressed) to specify regions. Note that the 'chr' prefix in
        contig names (e.g. 'chr1' vs. '1') will be automatically added or
        removed as necessary to match the input BED's contig names.
    path : str, optional
        Output BAM file. Writes to stdout when ``path='-'``. If None is
        provided the result is returned as a string.
    format : {'BAM', 'SAM', 'CRAM'}, default: 'BAM'
        Output file format.
    fasta
        FASTA file. Required when ``format`` is 'CRAM'.

    Returns
    -------
    None or str
        If ``path`` is None, returns the resulting BAM format as a string.
        Otherwise returns None.
    """
    formats = ['BAM', 'SAM', 'CRAM']

    if format not in formats:
        raise ValueError(
            f"Output format must be BAM, SAM, or CRAM, not {format}"
        )

    options = ['-h', '--no-PG']

    # Determine the output format.
    if format == 'BAM':
        options += ['-b']
    elif format == 'CRAM':
        if fasta is None:
            raise ValueError(
                "A FASTA file must be specified with '--fasta' "
                "when '--format' is 'CRAM'."
            )
        options += ['-C', '-T', fasta]
    else:
        pass

    # Parse the regions.
    if isinstance(regions, pybed.BedFrame):
        regions = regions.to_regions()
    elif isinstance(regions, list):
        if '.bed' in regions[0]:
            regions = pybed.BedFrame.from_file(regions[0]).to_regions()
        else:
            regions = common.sort_regions(regions)
    else:
        raise TypeError('Incorrect regions type')

    # Handle the 'chr' prefix.
    if has_chr_prefix(bam):
        regions = common.update_chr_prefix(regions, mode='add')
    else:
        regions = common.update_chr_prefix(regions, mode='remove')

    if path is None:
        data = pysam.view(bam, *regions, *options)
    else:
        data = None
        if path == '-':
            pysam.view(bam, *regions, *options, catch_stdout=False)
        else:
            pysam.view(bam, *regions, *options, '-o', path, catch_stdout=False)

    return data

def index(fn):
    """
    Index a BAM file.

    This simply wraps the :meth:`pysam.index` method.

    Parameters
    ----------
    fn : str
        BAM file.
    """
    pysam.index(fn)

def tag_sm(fn):
    """
    Extract SM tags (sample names) from a BAM file.

    Parameters
    ----------
    fn : str
        BAM file.

    Returns
    -------
    list
        List of SM tags.

    Examples
    --------

    >>> from fuc import pybam
    >>> pybam.tag_sm('NA19920.bam')
    ['NA19920']
    """
    lines = pysam.view('-H', fn, '--no-PG').strip().split('\n')
    tags = []
    for line in lines:
        fields = line.split('\t')
        if fields[0] == '@RG':
            for field in fields:
                if 'SM:' in field:
                    tags.append(field.replace('SM:', ''))
    return list(set(tags))

def tag_sn(fn):
    """
    Extract SN tags (contig names) from a BAM file.

    Parameters
    ----------
    fn : str
        BAM file.

    Returns
    -------
    list
        List of SN tags.

    Examples
    --------

    >>> from fuc import pybam
    >>> pybam.tag_sn('NA19920.bam')
    ['chr3', 'chr15', 'chrY', 'chr19', 'chr22', 'chr5', 'chr18', 'chr14', 'chr11', 'chr20', 'chr21', 'chr16', 'chr10', 'chr13', 'chr9', 'chr2', 'chr17', 'chr12', 'chr6', 'chrM', 'chrX', 'chr4', 'chr8', 'chr1', 'chr7']
    """
    lines = pysam.view('-H', fn, '--no-PG').strip().split('\n')
    tags = []
    for line in lines:
        fields = line.split('\t')
        if fields[0] == '@SQ':
            for field in fields:
                if 'SN:' in field:
                    tags.append(field.replace('SN:', ''))
    return list(set(tags))

def has_chr_prefix(fn):
    """
    Return True if contigs have the (annoying) 'chr' string.

    Parameters
    ----------
    fn : str
        BAM file.

    Returns
    -------
    bool
        Whether the 'chr' string is found.
    """
    contigs = tag_sn(fn)
    for contig in contigs:
        if 'chr' in contig:
            return True
    return False

def count_allelic_depth(bam, sites):
    """
    Count allelic depth for specified sites.

    Parameters
    ----------
    bam : str
        BAM file.
    sites : str or list
        Genomic site or list of sites. Each site should consist of chromosome
        and 1-based position in the format that can be recognized by
        :meth:`common.parse_variant` (e.g. '22-42127941').

    Returns
    -------
    pandas.DataFrame
        DataFrame containing allelic depth.

    Examples
    --------

    >>> from fuc import pybam
    >>> pybam.count_allelic_depth('in.bam', ['19-41510048', '19-41510053', '19-41510062'])
      Chromosome  Position  Total    A  C    G    T  N  DEL  INS
    0         19  41510048    119  106  7    4    0  0    2    0
    1         19  41510053    120    1  2    0  116  0    0    1
    2         19  41510062    115    0  0  115    0  0    0    0
    """
    if isinstance(sites, str):
        sites = [sites]

    has_prefix = has_chr_prefix(bam)

    f = pysam.AlignmentFile(bam, 'rb')

    rows = []

    for site in sites:
        chrom, pos, _, _ = common.parse_variant(site)

        if has_prefix:
            if 'chr' not in chrom:
                formatted_chrom = 'chr' + chrom
            else:
                formatted_chrom = chrom
        else:
            if 'chr' in chrom:
                formatted_chrom = chrom.replace('chr', '')
            else:
                formatted_chrom = chrom
        row = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0 , 'DEL': 0, 'INS': 0}
        kwargs = dict(
            min_base_quality=0,
            ignore_overlaps=False,
            ignore_orphans=False,
            truncate=True
        )
        for i, col in enumerate(f.pileup(formatted_chrom, pos-2, pos, **kwargs)):
            for read in col.pileups:
                if i == 0:
                    if read.indel < 0:
                        row['DEL'] += 1
                    elif read.indel > 0:
                        row['INS'] += 1
                    else:
                        continue
                else:
                    if read.is_del:
                        continue
                    allele = read.alignment.query_sequence[read.query_position]
                    row[allele] += 1

        data = list(row.values())
        rows.append([chrom, pos, sum(data)] + data)

    f.close()

    return pd.DataFrame(rows, columns=['Chromosome', 'Position', 'Total']+list(row))
