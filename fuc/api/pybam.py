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
    """Extract the SM tags (sample names) from a BAM file.

    Parameters
    ----------
    fn : str
        Path to the BAM file.

    Returns
    -------
    list
        SM tags.
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
    """Extract the SN tags (contig names) from a BAM file.

    Parameters
    ----------
    fn : str
        Path to the BAM file.

    Returns
    -------
    list
        SN tags.
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
