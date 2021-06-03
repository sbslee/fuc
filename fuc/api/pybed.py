"""
The pybed submodule is designed for working with BED files. It
implements ``pybed.BedFrame`` which stores BED data as ``pandas.DataFrame``
via the `pyranges <https://github.com/biocore-ntnu/pyranges>`_ package to
allow fast computation and easy manipulation. The submodule strictly adheres
to the standard `BED specification
<https://genome.ucsc.edu/FAQ/FAQformat.html>`_.

BED lines can have the following fields (the first three are required):

+-----+-------------+-----------------------------------+--------------+
| No. | Name        | Description                       | Examples     |
+=====+=============+===================================+==============+
| 1   | Chromosome  | Chromosome                        | 'chr2', '2'  |
+-----+-------------+-----------------------------------+--------------+
| 2   | Start       | Start position                    | 10041, 23042 |
+-----+-------------+-----------------------------------+--------------+
| 3   | End         | End position                      | 10041, 23042 |
+-----+-------------+-----------------------------------+--------------+
| 4   | Name        | Feature name                      | 'TP53'       |
+-----+-------------+-----------------------------------+--------------+
| 5   | Score       | Score for color density (0, 1000) | 342, 544     |
+-----+-------------+-----------------------------------+--------------+
| 6   | Strand      | '+' or '-' ('.' for no strand)    | '+', '-'     |
+-----+-------------+-----------------------------------+--------------+
| 7   | ThickStart  | Start position for thick drawing  | 10041, 23042 |
+-----+-------------+-----------------------------------+--------------+
| 8   | ThickEnd    | End position for thick drawing    | 10041, 23042 |
+-----+-------------+-----------------------------------+--------------+
| 9   | ItemRGB     | RGB value                         | '255,0,0'    |
+-----+-------------+-----------------------------------+--------------+
| 10  | BlockCount  | Number of blocks (e.g. exons)     | 12, 8        |
+-----+-------------+-----------------------------------+--------------+
| 11  | BlockSizes  | ','-separated block sizes         | '224,423'    |
+-----+-------------+-----------------------------------+--------------+
| 12  | BlockStarts | ','-separated block starts        | '2345,5245'  |
+-----+-------------+-----------------------------------+--------------+
"""

import pandas as pd
import pyranges as pr
from copy import deepcopy

HEADERS = [
    'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ThickStart',
    'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts'
]

class BedFrame:
    """Class for storing BED data.

    Parameters
    ----------
    meta : list
        Metadata lines.
    gr : pyranges.PyRanges
        PyRanges object containing BED data.

    See Also
    --------
    BedFrame.from_dict
        Construct BedFrame from dict of array-like or dicts.
    BedFrame.from_file
        Construct BedFrame from a BED file.
    BedFrame.from_frame
        Construct BedFrame from DataFrame.

    Examples
    --------

    >>> import pandas as pd
    >>> import pyranges as pr
    >>> from fuc import pybed
    >>> data = {
    ...     'Chromosome': ['chr1', 'chr2', 'chr3'],
    ...     'Start': [100, 400, 100],
    ...     'End': [200, 500, 200]
    ... }
    >>> df = pd.DataFrame(data)
    >>> gr = pr.PyRanges(df)
    >>> bf = pybed.BedFrame([], gr)
    >>> bf.gr.df
      Chromosome  Start  End
    0       chr1    100  200
    1       chr2    400  500
    2       chr3    100  200
    """
    def __init__(self, meta, gr):
        self._meta = meta
        self._gr = gr

    @property
    def meta(self):
        """list : Metadata lines."""
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
            Metadata lines.
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        BedFrame
            BedFrame.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_file
            Construct BedFrame from a BED file.
        BedFrame.from_frame
            Construct BedFrame from DataFrame.

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

    @classmethod
    def from_file(cls, fn):
        """Construct BedFrame from a BED file.

        Parameters
        ----------
        fn : str
            BED file path.

        Returns
        -------
        BedFrame
            BedFrame.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_dict
            Construct BedFrame from dict of array-like or dicts.
        BedFrame.from_frame
            Construct BedFrame from DataFrame.

        Examples
        --------

        >>> from fuc import pybed
        >>> bf = pybed.BedFrame.from_file('example.bed')
        """
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
        return cls(meta, pr.PyRanges(df))

    @classmethod
    def from_frame(cls, meta, data):
        """Construct BedFrame from DataFrame.

        Parameters
        ----------
        meta : list
            Metadata lines.
        data : pandas.DataFrame
            DataFrame containing BED data.

        Returns
        -------
        BedFrame
            BedFrame.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_dict
            Construct BedFrame from dict of array-like or dicts.
        BedFrame.from_file
            Construct BedFrame from a BED file.

        Examples
        --------

        >>> import pandas as pd
        >>> from fuc import pybed
        >>> data = {
        ...     'Chromosome': ['chr1', 'chr2', 'chr3'],
        ...     'Start': [100, 400, 100],
        ...     'End': [200, 500, 200]
        ... }
        >>> df = pd.DataFrame(data)
        >>> bf = pybed.BedFrame.from_frame([], df)
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1    100  200
        1       chr2    400  500
        2       chr3    100  200
        """
        return cls(meta, pr.PyRanges(data))
