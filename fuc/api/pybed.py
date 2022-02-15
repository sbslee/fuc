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
from . import common

HEADERS = [
    'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ThickStart',
    'ThickEnd', 'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts'
]

class BedFrame:
    """
    Class for storing BED data.

    Parameters
    ----------
    meta : list
        Metadata lines.
    gr : pyranges.PyRanges
        PyRanges object containing BED data.

    See Also
    --------
    BedFrame.from_dict
        Construct BedFrame from a dict of array-like or dicts.
    BedFrame.from_file
        Construct BedFrame from a BED file.
    BedFrame.from_frame
        Construct BedFrame from a dataframe.
    BedFrame.from_region
        Construct BedFrame from a list of regions.

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

    @property
    def contigs(self):
        """list : List of contig names."""
        return self.gr.chromosomes

    @property
    def shape(self):
        """tuple : Dimensionality of BedFrame (intervals, columns)."""
        return self.gr.df.shape

    @property
    def has_chr_prefix(self):
        """bool : Whether the (annoying) 'chr' string is found."""
        for contig in self.contigs:
            if 'chr' in contig:
                return True
        return False

    def copy_meta(self):
        """Return a copy of the metadata."""
        return deepcopy(self.meta)

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
        """
        Construct BedFrame from a dict of array-like or dicts.

        Parameters
        ----------
        meta : list
            Metadata lines.
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        BedFrame
            BedFrame object.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_file
            Construct BedFrame from a BED file.
        BedFrame.from_frame
            Construct BedFrame from a dataframe.
        BedFrame.from_region
            Construct BedFrame from a list of regions.

        Examples
        --------

        >>> from fuc import pybed
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
        """
        Construct BedFrame from a BED file.

        Parameters
        ----------
        fn : str
            BED file path.

        Returns
        -------
        BedFrame
            BedFrame object.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_dict
            Construct BedFrame from a dict of array-like or dicts.
        BedFrame.from_frame
            Construct BedFrame from a dataframe.
        BedFrame.from_region
            Construct BedFrame from a list of regions.

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
        """
        Construct BedFrame from a dataframe.

        Parameters
        ----------
        meta : list
            Metadata lines.
        data : pandas.DataFrame
            DataFrame containing BED data.

        Returns
        -------
        BedFrame
            BedFrame object.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_dict
            Construct BedFrame from a dict of array-like or dicts.
        BedFrame.from_file
            Construct BedFrame from a BED file.
        BedFrame.from_region
            Construct BedFrame from a list of regions.

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

    @classmethod
    def from_regions(cls, meta, regions):
        """
        Construct BedFrame from a list of regions.

        Parameters
        ----------
        meta : list
            Metadata lines.
        regions : str or list
            Region or list of regions.

        Returns
        -------
        BedFrame
            BedFrame object.

        See Also
        --------
        BedFrame
            BedFrame object creation using constructor.
        BedFrame.from_dict
            Construct BedFrame from a dict of array-like or dicts.
        BedFrame.from_file
            Construct BedFrame from a BED file.
        BedFrame.from_frame
            Construct BedFrame from a dataframe.

        Examples
        --------

        >>> from fuc import pybed
        >>> data = ['chr1:100-200', 'chr2:100-200', 'chr3:100-200']
        >>> bf = pybed.BedFrame.from_regions([], data)
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1    100  200
        1       chr2    100  200
        2       chr3    100  200
        """
        if isinstance(regions, str):
            regions = [regions]
        chroms = []
        starts = []
        ends = []
        for region in regions:
            chrom, start, end = common.parse_region(region)
            chroms.append(chrom)
            starts.append(start)
            ends.append(end)
        data = {'Chromosome': chroms, 'Start': starts, 'End': ends}
        return cls.from_dict([], data)

    def update_chr_prefix(self, mode='remove'):
        """
        Add or remove the (annoying) 'chr' string from the Chromosome column.

        Parameters
        ----------
        mode : {'add', 'remove'}, default: 'remove'
            Whether to add or remove the 'chr' string.

        Returns
        -------
        BedFrame
            Updated BedFrame.

        Examples
        --------

        >>> from fuc import pybed
        >>> data = {
        ...     'Chromosome': ['1', '1', 'chr2', 'chr2'],
        ...     'Start': [100, 400, 100, 200],
        ...     'End': [200, 500, 200, 300]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.gr.df
          Chromosome  Start  End
        0          1    100  200
        1          1    400  500
        2       chr2    100  200
        3       chr2    200  300
        >>> bf.update_chr_prefix(mode='remove').gr.df
          Chromosome  Start  End
        0          1    100  200
        1          1    400  500
        2          2    100  200
        3          2    200  300
        >>> bf.update_chr_prefix(mode='add').gr.df
          Chromosome  Start  End
        0       chr1    100  200
        1       chr1    400  500
        2       chr2    100  200
        3       chr2    200  300
        """
        if mode == 'remove':
            def one_row(r):
                r.Chromosome = r.Chromosome.replace('chr', '')
                return r
        elif mode == 'add':
            def one_row(r):
                if 'chr' not in r.Chromosome:
                    r.Chromosome = 'chr' + r.Chromosome
                return r
        else:
            raise ValueError(f'Incorrect mode: {mode}')
        df = self.gr.df.apply(one_row, axis=1)
        return self.__class__([], pr.PyRanges(df))

    def sort(self):
        """
        Sort the BedFrame by chromosome and position.

        Returns
        -------
        BedFrame
            Sorted BedFrame.

        Examples
        --------

        >>> from fuc import pybed
        >>> data = {
        ...     'Chromosome': ['chr1', 'chr3', 'chr1'],
        ...     'Start': [400, 100, 100],
        ...     'End': [500, 200, 200]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1    400  500
        1       chr1    100  200
        2       chr3    100  200
        >>> bf.sort().gr.df
          Chromosome  Start  End
        0       chr1    100  200
        1       chr1    400  500
        2       chr3    100  200
        """
        return self.__class__(self.copy_meta(), self.gr.sort())

    def merge(self):
        """
        Merge overlapping intervals within BedFrame.

        Returns
        -------
        BedFrame
            Merged BedFrame.

        Examples
        --------

        >>> from fuc import pybed
        >>> data = {
        ...     'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3', 'chr3'],
        ...     'Start': [10, 30, 15, 25, 50, 61],
        ...     'End': [40, 50, 25, 35, 60, 80]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.gr.df
          Chromosome  Start  End
        0       chr1     10   40
        1       chr1     30   50
        2       chr2     15   25
        3       chr2     25   35
        4       chr3     50   60
        5       chr3     61   80
        >>> bf.merge().gr.df
          Chromosome  Start  End
        0       chr1     10   50
        1       chr2     15   35
        2       chr3     50   60
        3       chr3     61   80
        """
        return self.__class__(self.copy_meta(), self.gr.merge())

    def to_regions(self, merge=True):
        """
        Return a list of regions from BedFrame.

        Parameters
        ----------
        merge : bool, default: True
            Whether to merge overlapping intervals.

        Returns
        -------
        list
            List of regions.

        Examples
        --------
        >>> from fuc import pybed
        >>> data = {
        ...     'Chromosome': ['chr1', 'chr1', 'chr2', 'chr2', 'chr3', 'chr3'],
        ...     'Start': [10, 30, 15, 25, 50, 61],
        ...     'End': [40, 50, 25, 35, 60, 80]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.to_regions()
        ['chr1:10-50', 'chr2:15-35', 'chr3:50-60', 'chr3:61-80']
        >>> bf.to_regions(merge=False)
        ['chr1:10-40', 'chr1:30-50', 'chr2:15-25', 'chr2:25-35', 'chr3:50-60', 'chr3:61-80']
        """
        if merge:
            bf = self.merge()
        else:
            bf = self
        s = bf.gr.df.apply(lambda r: f'{r.Chromosome}:{r.Start}-{r.End}', axis=1)
        return s.to_list()
