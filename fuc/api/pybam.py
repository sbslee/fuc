"""
The pybam submodule is designed for working with sequence alignment files
(SAM/BAM/CRAM). It essentially wraps the `pysam
<https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast
computation and easy manipulation.
"""

import pysam
import tempfile

def header(fn):
    """Return the header of a BAM file."""
    return pysam.view('-H', fn).strip()

def rename(fn, name, out, index=True):
    """Update the SM tags in a BAM file.

    This method will write a new BAM file with updated SM tags. By default,
    the method will also create an index file.

    .. warning::
        This method will only update the SM tags and nothing else;
        therefore, the old sample name(s) may still be present in other
        tags.

    Parameters
    ----------
    fn : str
        Path to the input BAM file.
    name : str
        New sample name.
    out : str
        Path to the output BAM file.
    index : bool, default: True
        If False, don't index the output BAM file.
    """
    lines1 = pysam.view('-H', fn).strip().split('\n')
    lines2 = []
    for line in lines1:
        fields = line.split('\t')
        if fields[0] != '@RG':
            lines2.append(line)
            continue
        t = []
        for field in fields:
            if 'SM:' not in field:
                t.append(field)
                continue
            t.append(f'SM:{name}')
        lines2.append('\t'.join(t))
    with tempfile.TemporaryDirectory() as t:
        with open(f'{t}/header.sam', 'w') as f1:
            f1.write('\n'.join(lines2))
        with open(out, 'wb') as f2:
            f2.write(pysam.reheader(f'{t}/header.sam', fn))
    if index:
        pysam.index(out)

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
