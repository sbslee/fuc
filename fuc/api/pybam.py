"""
The pybam submodule is designed for working with sequence alignment files
(SAM/BAM/CRAM). It essentially wraps the `pysam
<https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast
computation and easy manipulation.
"""

import pysam

def tag_sm(fn):
    """
    Extract the SM tags (sample names) from a BAM file.

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
    Extract the SN tags (contig names) from a BAM file.

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

def has_chr(fn):
    """
    Return True if the 'chr' string is present in the contig names.

    Parameters
    ----------
    fn : str
        BAM file.

    Returns
    -------
    bool
        Whether or not the 'chr' string is present.
    """
    contigs = tag_sn(fn)
    for contig in contigs:
        if 'chr' in contig:
            return True
    return False

def count_allelic_depth(bam, chrom, pos):
    """
    Count depth of coverage for every possible allele.

    Parameters
    ----------
    bam : str
        BAM file.
    chrom : str
        Contig name.
    pos : int
        Genomic position (1-based).

    Returns
    -------
    dict
        Dictionary containing allelic depth counts.

    Examples
    --------
    >>> from fuc import pybam
    >>> pybam.count_allelic_depth('in.bam', '19', 41510062)
    {'A': 0, 'C': 0, 'G': 115, 'T': 0, 'N': 0, 'D': 0, 'I': 0}
    >>> pybam.count_allelic_depth('in.bam', '19', 41510048)
    {'A': 106, 'C': 7, 'G': 4, 'T': 0, 'N': 0, 'D': 2, 'I': 0}
    >>> pybam.count_allelic_depth('in.bam', '19', 41510053)
    {'A': 1, 'C': 2, 'G': 0, 'T': 116, 'N': 0, 'D': 0, 'I': 1}
    """
    d = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0 , 'D': 0, 'I': 0}
    alignment_file = pysam.AlignmentFile(bam, 'rb')
    kwargs = dict(
        min_base_quality=0,
        ignore_overlaps=False,
        ignore_orphans=False,
        truncate=True
    )
    for pileupcolumn in alignment_file.pileup(chrom, pos-1, pos, **kwargs):
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                continue
            allele = pileupread.alignment.query_sequence[pileupread.query_position]
            d[allele] += 1
    for pileupcolumn in alignment_file.pileup(chrom, pos-2, pos-1, **kwargs):
        for pileupread in pileupcolumn.pileups:
            if pileupread.indel > 0:
                d['I'] += 1
            elif pileupread.indel < 0:
                d['D'] += 1
            else:
                continue
    alignment_file.close()
    return d
