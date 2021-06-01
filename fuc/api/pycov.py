"""
The pycov submodule is designed for working with depth of coverage data
from sequence alingment files (SAM/BAM/CRAM). It implements
``pycov.CovFrame`` which stores read depth data as ``pandas.DataFrame`` via
the `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ package to
allow fast computation and easy manipulation.
"""

import pysam
import numpy as np
import pandas as pd
from io import StringIO
from . import pybam
import seaborn as sns
import matplotlib.pyplot as plt

class CovFrame:
    """Class for storing read depth data from one or more SAM/BAM/CRAM files.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing read depth data.

    See Also
    --------
    CovFrame.from_dict
        Construct CovFrame from dict of array-like or dicts.
    CovFrame.from_file
        Construct CovFrame from one or more SAM/BAM/CRAM files.

    Examples
    --------

    >>> import numpy as np
    >>> import pandas as pd
    >>> from fuc import pycov
    >>> data = {
    ...     'Chromosome': ['chr1'] * 1000,
    ...     'Position': np.arange(1000, 2000),
    ...     'Steven': np.random.normal(35, 5, 1000),
    ...     'Jane': np.random.normal(25, 7, 1000)
    ... }
    >>> df = pd.DataFrame(data)
    >>> cf = pycov.CovFrame(df)
    >>> cf.df.head()
      Chromosome  Position     Steven       Jane
    0       chr1      1000  26.476796  30.974197
    1       chr1      1001  30.869688  20.510165
    2       chr1      1002  36.039488  25.444763
    3       chr1      1003  39.002843  19.379962
    4       chr1      1004  35.641304  20.359149
    """
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @property
    def samples(self):
        """list : List of the sample names."""
        return self.df.columns[2:].to_list()

    @property
    def shape(self):
        """Return the size of CovFrame."""
        return self.df.shape

    @classmethod
    def from_dict(cls, data):
        """Construct CovFrame from dict of array-like or dicts.

        Parameters
        ----------
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        CovFrame
            CovFrame.

        See Also
        --------
        CovFrame
            BedFrame object creation using constructor.
        CovFrame.from_file
            Construct CovFrame from one or more SAM/BAM/CRAM files.

        Examples
        --------

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1'] * 1000,
        ...     'Position': np.arange(1000, 2000),
        ...     'Steven': np.random.normal(35, 5, 1000),
        ...     'Jane': np.random.normal(25, 7, 1000)
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df.head()
          Chromosome  Position     Steven       Jane
        0       chr1      1000  42.357973  29.578929
        1       chr1      1001  33.598807  32.370608
        2       chr1      1002  44.424704  20.425198
        3       chr1      1003  29.228379  22.406113
        4       chr1      1004  34.395085  29.066962
        """
        return cls(pd.DataFrame(data))

    @classmethod
    def from_file(
        cls, fn, zero=False, region=None, map_qual=None, names=None
    ):
        """Construct CovFrame from one or more SAM/BAM/CRAM files.

        This method essentially wraps the ``samtools depth`` command.

        Parameters
        ----------
        fn : str or list
            SAM/BAM/CRAM file path(s).
        zero : bool, default: False
            If True, output all positions (including those with zero depth).
        region : str, optional
            Only report depth in specified region (format: CHR:FROM-TO).
        map_qual: int, optional
            Only count reads with mapping quality greater than orequal to
            this number.
        names : list, optional
            Use these as sample names instead of the SM tags.

        Returns
        -------
        CovFrame
            CovFrame.
        CovFrame.from_dict
            Construct CovFrame from dict of array-like or dicts.

        Examples
        --------

        >>> from fuc import pycov
        >>> cf = pycov.CovFrame.from_file(bam)
        >>> cf = pycov.CovFrame.from_file([bam1, bam2])
        >>> cf = pycov.CovFrame.from_file(bam, region='19:41497204-41524301')
        """
        bam_files = []
        if isinstance(fn, str):
            bam_files.append(fn)
        else:
            bam_files += fn
        args = []
        if zero:
            args += ['-a']
        if region is not None:
            args += ['-r', region]
        if map_qual is not None:
            args += ['-Q', str(map_qual)]
        args += bam_files
        s = pysam.depth(*args)
        headers = ['Chromosome', 'Position']
        dtype = {'Chromosome': str,'Position': np.int32}
        for i, bam_file in enumerate(bam_files):
            if names:
                name = names[i]
            else:
                samples = pybam.tag_sm(bam_file)
                if len(samples) > 1:
                    m = f'multiple sample names detected: {bam_file}'
                    raise ValueError(m)
                name = samples[0]
            headers.append(name)
            dtype[name] = np.int32
        df = pd.read_csv(StringIO(s), sep='\t', header=None,
                         names=headers, dtype=dtype)
        return cls(df)

    def plot_region(
        self, chrom, start=None, end=None, names=None, ax=None,
        figsize=None, **kwargs
    ):
        """Create a read depth profile for the given region.

        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int, optional
            Start position.
        end : int, optional
            End position.
        names : str or list, optional
            Sample name or list of sample names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.plot` and :meth:`seaborn.lineplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> from fuc import pycov
            >>> data = {
            ...    'Chromosome': ['chr1'] * 1000,
            ...    'Position': np.arange(1000, 2000),
            ...    'Steven': np.random.normal(35, 5, 1000),
            ...    'Jane': np.random.normal(25, 7, 1000)
            ... }
            >>> cf = pycov.CovFrame.from_dict(data)
            >>> cf.plot_region('chr1')
            >>> plt.tight_layout()
        """
        cf = self.slice(chrom, start=start, end=end)
        if names is None:
            names = cf.samples
        if isinstance(names, str):
            names = [names]
        headers = ['Position'] + names
        df = cf.df[headers]
        df = df.set_index('Position')
        if kwargs is None:
            kwargs = {}
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(data=df, ax=ax, **kwargs)
        ax.set_ylabel('Depth')
        return ax

    def slice(self, chrom, start=None, end=None):
        """Return a sliced CovFrame for the given region.

        Parameters
        ----------
        chrom : str
            Chromosome.
        start : int, optional
            Start position.
        end : int, optional
            End position.

        Returns
        -------
        CovFrame
            Sliced CovFrame.

        Examples
        --------

        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1']*500 + ['chr2']*500,
        ...     'Position': np.arange(1000, 2000),
        ...     'Steven': np.random.normal(35, 5, 1000),
        ...     'Jane': np.random.normal(25, 7, 1000)
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.slice('chr2').df.head()
          Chromosome  Position     Steven       Jane
        0       chr2      1500  31.408382  20.409376
        1       chr2      1501  36.173578  22.327581
        2       chr2      1502  22.867945   9.962580
        3       chr2      1503  29.047738  17.589284
        4       chr2      1504  41.343270  26.612494
        """
        df = self.df[self.df.Chromosome == chrom]
        if start:
            df = df[df.Position >= start]
        if end:
            df = df[df.Position <= end]
        return self.__class__(df)
