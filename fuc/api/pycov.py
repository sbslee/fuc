"""
The pycov submodule is designed for working with depth of coverage data
from sequence alingment files (SAM/BAM/CRAM). It implements
``pycov.CovFrame`` which stores read depth data as ``pandas.DataFrame`` via
the `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ package to
allow fast computation and easy manipulation.
"""
from io import StringIO

from . import common, pybam

import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def simulate(mode='wgs', loc=30, scale=5, size=1000):
    """
    Simulate read depth data for single sample.

    Generated read depth will be integer and non-negative.

    Parameters
    ----------
    mode : {'wgs'}, default: 'wgs'
        Additional modes will be made available in future releases.
    loc : float, default: 30
        Mean ("centre") of the distribution.
    scale : float, default: 5
        Standard deviation (spread or "width") of the distribution. Must be
        non-negative.
    size : int, default: 1000
        Number of base pairs to return.

    Returns
    -------
    numpy.ndarray
        Numpy array object.

    Examples
    --------

    >>> from fuc import pycov
    >>> pycov.simulate(size=10)
    array([25, 32, 30, 31, 26, 25, 33, 29, 28, 35])
    """
    a = np.random.normal(loc=loc, scale=scale, size=size)

    # Read depth must be integer.
    a = a.astype(int)

    # Read dpeth must be non-negative.
    a = np.absolute(a)

    return a

class CovFrame:
    """
    Class for storing read depth data from one or more SAM/BAM/CRAM files.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing read depth data.

    See Also
    --------
    CovFrame.from_bam
        Construct CovFrame from one or more SAM/BAM/CRAM files.
    CovFrame.from_dict
        Construct CovFrame from dict of array-like or dicts.
    CovFrame.from_file
        Construct CovFrame from a text file containing read depth data.

    Examples
    --------

    >>> import numpy as np
    >>> import pandas as pd
    >>> from fuc import pycov
    >>> data = {
    ...     'Chromosome': ['chr1'] * 1000,
    ...     'Position': np.arange(1000, 2000),
    ...     'A': pycov.simulate(loc=35, scale=5),
    ...     'B': pycov.simulate(loc=25, scale=7),
    ... }
    >>> df = pd.DataFrame(data)
    >>> cf = pycov.CovFrame(df)
    >>> cf.df.head()
      Chromosome  Position   A   B
    0       chr1      1000  22  23
    1       chr1      1001  34  30
    2       chr1      1002  33  27
    3       chr1      1003  32  21
    4       chr1      1004  32  15
    """
    def __init__(self, df):
        self._df = df.reset_index(drop=True)

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing read depth data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @property
    def samples(self):
        """list : List of the sample names."""
        return self.df.columns[2:].to_list()

    @property
    def shape(self):
        """tuple : Dimensionality of CovFrame (positions, samples)."""
        return (self.df.shape[0], len(self.samples))

    @classmethod
    def from_bam(
        cls, bam=None, fn=None, bed=None, zero=False, region=None,
        map_qual=None, names=None
    ):
        """
        Construct CovFrame from one or more SAM/BAM/CRAM files.

        Either the 'bam' or 'fn' parameter must be provided, but not both.

        Under the hood, this method computes read depth from the input files
        using the :command:`samtools depth` command.

        Some parameters such as 'bed' and 'region' require that the input
        files be indexed.

        Parameters
        ----------
        bam : str or list, optional
            One or more input files.
        fn : str, optional
            File containing one input filename per line.
        bed : str, optional
            BED file.
        zero : bool, default: False
            If True, output all positions (including those with zero depth).
        region : str, optional
            Only report depth in the specified region ('chrom:start-end').
        map_qual: int, optional
            Only count reads with mapping quality greater than orequal to
            this number.
        names : list, optional
            Use these as sample names instead of the SM tags.

        Returns
        -------
        CovFrame
            CovFrame object.

        See Also
        --------
        CovFrame
            CovFrame object creation using constructor.
        CovFrame.from_dict
            Construct CovFrame from dict of array-like or dicts.
        CovFrame.from_file
            Construct CovFrame from a text file containing read depth data.

        Examples
        --------

        >>> from fuc import pycov
        >>> cf = pycov.CovFrame.from_bam(bam)
        >>> cf = pycov.CovFrame.from_bam([bam1, bam2])
        >>> cf = pycov.CovFrame.from_bam(bam, region='19:41497204-41524301')
        """
        bam_files = []

        if bam is None and fn is None:
            raise ValueError(
                "Either the 'bam' or 'fn' parameter must be provided.")
        elif bam is not None and fn is not None:
            raise ValueError(
                "The 'bam' and 'fn' parameters cannot be used together.")
        elif bam is not None and fn is None:
            if isinstance(bam, str):
                bam_files.append(bam)
            else:
                bam_files += bam
        else:
            bam_files += common.convert_file2list(fn)

        args = []

        if zero:
            args += ['-a']
        if region is not None:
            args += ['-r', region]
        if map_qual is not None:
            args += ['-Q', str(map_qual)]
        if bed is not None:
            args += ['-b', bed]

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

    @classmethod
    def from_dict(cls, data):
        """
        Construct CovFrame from dict of array-like or dicts.

        Parameters
        ----------
        data : dict
            Of the form {field : array-like} or {field : dict}.

        Returns
        -------
        CovFrame
            CovFrame object.

        See Also
        --------
        CovFrame
            CovFrame object creation using constructor.
        CovFrame.from_bam
            Construct CovFrame from one or more SAM/BAM/CRAM files.
        CovFrame.from_file
            Construct CovFrame from a text file containing read depth data.

        Examples
        --------

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1'] * 1000,
        ...     'Position': np.arange(1000, 2000),
        ...     'A': pycov.simulate(loc=35, scale=5),
        ...     'B': pycov.simulate(loc=25, scale=7),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df.head()
          Chromosome  Position   A   B
        0       chr1      1000  36  22
        1       chr1      1001  39  35
        2       chr1      1002  33  19
        3       chr1      1003  36  20
        4       chr1      1004  31  24
        """
        return cls(pd.DataFrame(data))

    @classmethod
    def from_file(cls, fn):
        """
        Construct CovFrame from a text file containing read depth data.

        Parameters
        ----------
        fn : str
            Text file containing read depth data.

        Returns
        -------
        CovFrame
            CovFrame object.

        See Also
        --------
        CovFrame
            CovFrame object creation using constructor.
        CovFrame.from_bam
            Construct CovFrame from one or more SAM/BAM/CRAM files.
        CovFrame.from_dict
            Construct CovFrame from dict of array-like or dicts.
        """
        return cls(pd.read_table(fn))

    def plot_region(
        self, region, names=None, ax=None, figsize=None, **kwargs
    ):
        """
        Create a read depth profile for the region.

        Parameters
        ----------
        region : str
            Region (‘chrom:start-end’).
        names : str or list, optional
            Sample name or list of sample names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.lineplot`.

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
            ...     'A': pycov.simulate(loc=35, scale=5),
            ...     'B': pycov.simulate(loc=25, scale=7),
            ... }
            >>> cf = pycov.CovFrame.from_dict(data)
            >>> cf.plot_region('chr1:1500-1800')
            >>> plt.tight_layout()
        """
        chrom, start, end = common.parse_region(region)
        cf = self.slice(region)
        if names is None:
            names = cf.samples
        if isinstance(names, str):
            names = [names]
        headers = ['Position'] + names
        df = cf.df[headers]
        df = df.set_index('Position')
        if kwargs is None:
            kwargs = {}

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(data=df, ax=ax, **kwargs)

        ax.set_ylabel('Depth')

        return ax

    def slice(self, region):
        """
        Slice the CovFrame for the region.

        Parameters
        ----------
        region : str
            Region (‘chrom:start-end’).

        Returns
        -------
        CovFrame
            Sliced CovFrame.

        Examples
        --------

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1']*500 + ['chr2']*500,
        ...     'Position': np.arange(1000, 2000),
        ...     'A': pycov.simulate(loc=35, scale=5),
        ...     'B': pycov.simulate(loc=25, scale=7),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.slice('chr2').df.head()
          Chromosome  Position   A   B
        0       chr2      1500  37  34
        1       chr2      1501  28  12
        2       chr2      1502  35  29
        3       chr2      1503  34  34
        4       chr2      1504  32  21
        >>> cf.slice('chr2:1500-1504').df
          Chromosome  Position   A   B
        0       chr2      1500  37  34
        1       chr2      1501  28  12
        2       chr2      1502  35  29
        3       chr2      1503  34  34
        4       chr2      1504  32  21
        >>> cf.slice('chr2:-1504').df
          Chromosome  Position   A   B
        0       chr2      1500  37  34
        1       chr2      1501  28  12
        2       chr2      1502  35  29
        3       chr2      1503  34  34
        4       chr2      1504  32  21
        """
        chrom, start, end = common.parse_region(region)
        df = self.df[self.df.Chromosome == chrom]
        if start:
            df = df[df.Position >= start]
        if end:
            df = df[df.Position <= end]
        return self.__class__(df)

    def to_string(self):
        """
        Render the CovFrame to a console-friendly tabular output.

        Returns
        -------
        str
            String representation of the CovFrame.
        """
        return self.df.to_csv(index=False, sep='\t')

    def plot_uniformity(
        self, mode='aggregated', frac=0.1, n=20, m=None, marker=None,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a line plot visualizing the uniformity in read depth.

        Parameters
        ----------
        mode : {'aggregated', 'individual'}, default: 'aggregated'
            Determines how to display the lines:

            - 'aggregated': Aggregate over repeated values to show the mean
              and 95% confidence interval.
            - 'individual': Show data for individual samples.

        frac : float, default: 0.1
            Fraction of data to be sampled (to speed up the process).
        n : int, default: 20
            Number of points to generate for the x-axis.
        m : float, optional
            Maximum point in the x-axis. By default, it will be the maximum
            depth in the entire dataset.
        marker : str, optional
            Marker style string (e.g. 'o').
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.lineplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        By default (``mode='aggregated'``), the method will aggregate over
        repeated values:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> from fuc import pycov
            >>> data = {
            ...     'Chromosome': ['chr1'] * 1000,
            ...     'Position': np.arange(1000, 2000),
            ...     'A': pycov.simulate(loc=35, scale=5),
            ...     'B': pycov.simulate(loc=25, scale=7),
            ... }
            >>> cf = pycov.CovFrame.from_dict(data)
            >>> cf.plot_uniformity(mode='aggregated')
            >>> plt.tight_layout()

        We can display data for individual samples:

        .. plot::
            :context: close-figs

            >>> cf.plot_uniformity(mode='individual')
            >>> plt.tight_layout()
        """
        # Sample positions to speed up the process.
        df = self.df.sample(frac=frac)

        # Determine x-axis points.
        if m is None:
            m = df.iloc[:, 2:].max().max()
        coverages = np.linspace(1, m, n, endpoint=True)

        data = {'Coverage': coverages}

        for sample in self.samples:
            fractions = []
            for coverage in coverages:
                count = sum(df[sample] >= coverage)
                fractions.append(count / df.shape[0])
            data[sample] = fractions
        df = pd.DataFrame(data)
        df = df.melt(id_vars=['Coverage'], var_name='Sample',
            value_name='Fraction')

        # Determines how to display the lines.
        if mode == 'aggregated':
            hue = None
        else:
            hue = 'Sample'

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(
            x='Coverage', y='Fraction', data=df, hue=hue, marker=marker,
            ax=ax, **kwargs
        )

        ax.set_xlabel('Sequencing coverage')
        ax.set_ylabel('Fraction of sampled bases')

        return ax

    def plot_distribution(
        self, mode='aggregated', frac=0.1, ax=None, figsize=None, **kwargs
    ):
        """
        Create a line plot visualizaing the distribution of per-base read depth.

        Parameters
        ----------
        mode : {'aggregated', 'individual'}, default: 'aggregated'
            Determines how to display the lines:

            - 'aggregated': Aggregate over repeated values to show the mean
              and 95% confidence interval.
            - 'individual': Show data for individual samples.

        frac : float, default: 0.1
            Fraction of data to be sampled (to speed up the process).
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.lineplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        By default (``mode='aggregated'``), the method will aggregate over
        repeated values:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import numpy as np
            >>> from fuc import pycov
            >>> data = {
            ...     'Chromosome': ['chr1'] * 1000,
            ...     'Position': np.arange(1000, 2000),
            ...     'A': pycov.simulate(loc=35, scale=5),
            ...     'B': pycov.simulate(loc=25, scale=7),
            ... }
            >>> cf = pycov.CovFrame.from_dict(data)
            >>> cf.plot_distribution(mode='aggregated', frac=0.9)
            >>> plt.tight_layout()

        We can display data for individual samples:

        .. plot::
            :context: close-figs

            >>> cf.plot_distribution(mode='individual', frac=0.9)
            >>> plt.tight_layout()
        """
        # Sample positions to speed up the process.
        df = self.df.sample(frac=frac)

        def one_col(c):
            s = c.value_counts()
            s = s / s.sum()
            return s

        df = df.iloc[:, 2:].apply(one_col)
        df = df.fillna(0)
        df = df.reset_index().rename(columns=dict(index='Coverage'))
        df = df.melt(id_vars='Coverage', var_name='Sample',
            value_name='Fraction')

        # Determines how to display the lines.
        if mode == 'aggregated':
            hue = None
        else:
            hue = 'Sample'

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(
            x='Coverage', y='Fraction', data=df, hue=hue, ax=ax, **kwargs
        )

        ax.set_xlabel('Sequencing coverage')
        ax.set_ylabel('Fraction of sampled bases')

        return ax
