"""
The pybam submodule is designed for working with sequence alignment files
(SAM/BAM/CRAM). It essentially wraps the `pysam
<https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast
computation and easy manipulation.
"""

import pysam

def header(fn):
    """Return the header of a BAM file."""
    return pysam.view('-H', fn).strip()

def tag_sm(fn):
    """
    Extract the SM tags (sample names) from a BAM file.

    Parameters
    ----------
    fn : str
        Path to the BAM file.

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
    lines = pysam.view('-H', fn).strip().split('\n')
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
        Path to the BAM file.

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
    lines = pysam.view('-H', fn).strip().split('\n')
    tags = []
    for line in lines:
        fields = line.split('\t')
        if fields[0] == '@SQ':
            for field in fields:
                if 'SN:' in field:
                    tags.append(field.replace('SN:', ''))
    return list(set(tags))
