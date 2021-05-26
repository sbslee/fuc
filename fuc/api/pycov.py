"""
The pycov submodule is designed for working with depth of coverage data
from BAM files. It implements the ``pycov.CovFrame`` class which stores read
depth data as ``pandas.DataFrame`` to allow fast computation and easy
manipulation. The ``pycov.CovFrame`` class also provides many useful plotting
methods.
"""

import pysam
import numpy as np
import pandas as pd
from io import StringIO
from . import pybam
import seaborn as sns
import matplotlib.pyplot as plt

def read_file(fn, zero=False, region=None, map_qual=None):
    """Create CovFrame from one or more BAM files.

    This method essentially wraps the ``samtools depth`` command.

    Parameters
    ----------
    fn : str or list
        Single BAM file path or list of multiple BAM file paths.
    zero : bool, default: False
        If True, output all positions (including those with zero depth).
    region : str, optional
        If provided, only report depth in specified region (format:
        CHR:FROM-TO).
    map_qual: int, optional
        If provided, only count reads with mapping quality greater than or
        equal to this number.

    Returns
    -------
    CovFrame
        CovFrame.

    Examples
    --------
    >>> from fuc import pycov
    >>> cf = pycov.read_file(bam)
    >>> cf = pycov.read_file([bam1, bam2])
    >>> cf = pycov.read_file(bam, region='19:41497204-41524301')
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
        args += ['-Q', map_qual]
    args += bam_files
    s = pysam.depth(*args)
    names = ['Chromosome', 'Position']
    dtype = {'Chromosome': str,'Position': np.int32}
    for bam_file in bam_files:
        samples = pybam.tag_sm(bam_file)
        if len(samples) > 1:
            raise ValueError(f'multiple sample names detected: {bam_file}')
        names.append(samples[0])
        dtype[samples[0]] = np.int32
    df = pd.read_csv(StringIO(s), sep='\t', header=None,
                     names=names, dtype=dtype)
    return CovFrame(df)

class CovFrame:
    """Class for storing read depth data from BAM files.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing read depth data.

    See Also
    --------
    CovFrame.from_dict
        Construct CovFrame from dict of array-like or dicts.
    read_file
        Read one or more BAM files into CovFrame.

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
        self.df = df

    @property
    def samples(self):
        """list : List of the sample names."""
        return self.df.columns[2:].to_list()

    @property
    def shape(self):
        """Return the size of the CovFrame."""
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

    def plot(
            self, names=None, ax=None, figsize=None, lineplot_kwargs=None,
            plot_kwargs=None
    ):
        """Plot read depth profile.

        Parameters
        ----------
        names : str or list
            Sample names.
        ax : matplotlib.axes.Axes, optional
            Axes object to draw the plot onto, otherwise uses the current Axes.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        lineplot_kwargs: dict, optional
            Keyword arguments passed down to the `seaborn.lineplot` method.
        plot_kwargs: dict, optional
            Keyword arguments passed down to the `matplotlib.axes.Axes.plot`
            method.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

           A plotting example:
           
           >>> import matplotlib.pyplot as plt
           >>> plt.plot([1,2,3], [4,5,6])
        """
        if names is None:
            names = self.samples
        headers = ['Position'] + names
        df = self.df[headers]
        df = df.set_index('Position')
        if lineplot_kwargs is None:
            lineplot_kwargs = {}
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        sns.lineplot(data=df, ax=ax, **lineplot_kwargs)
