"""
The pybam submodule is designed for working with BAM files.
"""

import pysam

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
    header = pysam.view('-H', fn).strip().split('\n')
    tags = []
    for line in header:
        fields = line.split('\t')
        if '@RG' == fields[0]:
            for field in fields:
                if 'SM:' in field:
                    tags.append(field.replace('SM:', ''))
    return list(set(tags))

def tag_sn(fn):
    """Extract the SN tags (sequence or contig names) from a BAM file.

    Parameters
    ----------
    fn : str
        Path to the BAM file.

    Returns
    -------
    list
        SN tags.
    """
    header = pysam.view('-H', fn).strip().split('\n')
    tags = []
    for line in header:
        fields = line.split('\t')
        if '@SQ' == fields[0]:
            for field in fields:
                if 'SN:' in field:
                    tags.append(field.replace('SN:', ''))
    return list(set(tags))
