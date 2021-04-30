"""
The pybed submodule is designed for working with BED files. It
implements ``pybed.BedFrame`` which stores BED data as
``pyranges.PyRanges`` to allow fast computation and easy manipulation.

This class is essentially a wrapper for the pyranges package
(https://github.com/biocore-ntnu/pyranges).

BED lines have three required fields and nine additional optional fields:

1. chrom (required) - The name of the chromosome.
2. chromStart (required) - The starting position of the feature.
3. chromEnd (required) - The ending position of the feature.
4. name (optional) - Defines the name of the BED line.
5. score (optional) - A score between 0 and 1000 for color density.
6. strand (optional) - Either '.' (=no strand) or '+' or '-'.
7. thickStart (optional) - The starting position for thick drawing.
8. thickEnd (optional) - The ending position for thick drawing.
9. itemRgb (optional) - An RGB value (e.g. 255,0,0).
10. blockCount (optional) - The number of blocks (exons).
11. blockSizes (optional) - A comma-separated list of the block sizes.
12. blockStarts (optional) - A comma-separated list of block starts.

For more information about the BED format, visit the UCSC Genome Browser
FAQ (https://genome.ucsc.edu/FAQ/FAQformat.html).
"""

import pandas as pd
import pyranges as pr
from copy import deepcopy

HEADERS = ['Chromosome', 'Start', 'End', 'Name',
           'Score', 'Strand', 'ThickStart', 'ThickEnd',
           'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts']

def read_file(fn):
    """Read a BED file into BedFrame."""
    meta = []
    skip_rows = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'browser' in line or 'track' in line:
                meta.append(line.strip())
                skip_rows += 1
            else:
                headers = HEADERS[:len(line.strip().split())]
                break
    df = pd.read_table(fn, header=None, names=headers, skiprows=skip_rows)
    bf = BedFrame(meta, pr.PyRanges(df))
    return bf

class BedFrame:
    """Class for storing BED data."""
    def __init__(self, meta, gr):
        self._meta = meta
        self._gr = gr

    @property
    def meta(self):
        """list : List of metadata lines."""
        return self._meta

    @meta.setter
    def meta(self, value):
        self._meta = value

    @property
    def gr(self):
        """pyranges.PyRanges : Two-dimensional representation of genomic
        intervals and their annotations."""
        return self._gr

    @gr.setter
    def gr(self, value):
        self._gr = value

    def to_file(self, fn):
        """Write the BedFrame to a BED file."""
        with open(fn, 'w') as f:
            if self.meta:
                f.write('\n'.join(self.meta) + '\n')
            self.gr.df.to_csv(f, sep='\t', index=False, header=False)

    def to_string(self):
        """Render the BedFrame to a console-friendly tabular output."""
        s = ''
        if self.meta:
            s += '\n'.join(self.meta) + '\n'
        s += self.gr.df.to_csv(header=False, index=False, sep='\t')
        return s

    def intersect(self, other):
        """Find intersection between the BedFrames."""
        bf = self.__class__(deepcopy(self.meta), self.gr.intersect(other.gr))
        return bf

    @classmethod
    def from_dict(cls, meta, data):
        """Construct BedFrame from dict of array-like or dicts.

        Parameters
        ----------
        meta : list
            List of the metadata lines.
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        BedFrame
            BedFrame

        See Also
        --------
        BedFrame : BedFrame object creation using constructor.

        Examples
        --------
        Below is a simple example:

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
        """
        return cls(meta, pr.PyRanges(pd.DataFrame(data)))
