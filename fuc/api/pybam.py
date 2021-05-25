"""
The pybam submodule is designed for working with BAM files.
"""

import pysam
import tempfile

def rename(fn, name, out):
    """Extract the SM tags (sample names) from a BAM file.

    .. warning::
        This method will only change the SM tags and nothing else;
        therefore, the old sample name(s) may still be present in other
        tags.

    Parameters
    ----------
    fn : str
        Path to the BAM file.
    name : str
        New sample name.

    Returns
    -------
    list
        SM tags.
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
