"""
The pycov submodule is designed for working with depth of coverage data
from sequence alingment files (SAM/BAM/CRAM). It implements
``pycov.CovFrame`` which stores read depth data as ``pandas.DataFrame`` via
the `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ package to
allow fast computation and easy manipulation. The ``pycov.CovFrame`` class
also contains many useful plotting methods such as ``CovFrame.plot_region``
and ``CovFrame.plot_uniformity``.
"""
from io import StringIO, IOBase
import warnings
import gzip
from pathlib import Path

from . import common, pybam, pybed

import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def concat(cfs, axis=0):
    """
    Concatenate CovFrame objects along a particular axis.

    Parameters
    ----------
    cfs : list
        List of CovFrame objects.
    axis : {0/'index', 1/'columns'}, default: 0
        The axis to concatenate along.

    Returns
    -------
    CovFrame
        Concatenated CovFrame.
    """
    if axis:
        df = pd.concat([x.df.set_index(['Chromosome', 'Position'])
            for x in cfs], axis=axis).reset_index()
    else:
        df = pd.concat([x.df for x in cfs], axis=axis)
    return CovFrame(df)

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
        Construct CovFrame from BAM files.
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
    def _check_df(self, df):
        df = df.reset_index(drop=True)
        df.Chromosome = df.Chromosome.astype(str)
        return df

    def __init__(self, df):
        self._df = self._check_df(df)

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing read depth data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = self._check_df(value)

    @property
    def contigs(self):
        """list : List of contig names."""
        return list(self.df.Chromosome.unique())

    @property
    def has_chr_prefix(self):
        """bool : Whether the (annoying) 'chr' string is found."""
        for contig in self.contigs:
            if 'chr' in contig:
                return True
        return False

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
        cls, bams, regions=None, zero=False, map_qual=None, names=None
    ):
        """
        Construct CovFrame from BAM files.

        Under the hood, the method computes read depth using the
        :command:`samtools depth` command.

        Parameters
        ----------
        bams : str or list
            One or more input BAM files. Alternatively, you can provide a
            text file (.txt, .tsv, .csv, or .list) containing one BAM file
            per line.
        regions : str, list, or pybed.BedFrame, optional
            By default (``regions=None``), the method counts all reads in BAM
            files, which can be excruciatingly slow for large files (e.g.
            whole genome sequencing). Therefore, use this argument to only
            output positions in given regions. Each region must have the
            format chrom:start-end and be a half-open interval with (start,
            end]. This means, for example, chr1:100-103 will extract
            positions 101, 102, and 103. Alternatively, you can provide a BED
            file (compressed or uncompressed) or a :class:`pybed.BedFrame`
            object to specify regions. Note that the 'chr' prefix in contig
            names (e.g. 'chr1' vs. '1') will be automatically added or
            removed as necessary to match the input BAM's contig names.
        zero : bool, default: False
            If True, output all positions including those with zero depth.
        map_qual: int, optional
            Only count reads with mapping quality greater than or equal to
            this number.
        names : list, optional
            By default (``names=None``), sample name is extracted using SM
            tag in BAM files. If the tag is missing, the method will set the
            filename as sample name. Use this argument to manually provide
            sample names.

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
        # Parse input BAM files.
        bams = common.parse_list_or_file(bams)

        # Check the 'chr' prefix.
        if all([pybam.has_chr_prefix(x) for x in bams]):
            chr_prefix = 'chr'
        else:
            chr_prefix = ''

        # Run the 'samtools depth' command.
        args = []
        if zero:
            args += ['-a']
        if map_qual is not None:
            args += ['-Q', str(map_qual)]
        results = ''
        if regions is None:
            results += pysam.depth(*(bams + args))
        else:
            if isinstance(regions, str):
                regions = [regions]
            elif isinstance(regions, list):
                pass
            elif isinstance(regions, pybed.BedFrame):
                regions = bf.to_regions()
            else:
                raise TypeError("Incorrect type of argument 'regions'")
            if '.bed' in regions[0]:
                bf = pybed.BedFrame.from_file(regions[0])
                regions = bf.to_regions()
            else:
                regions = common.sort_regions(regions)
            for region in regions:
                region = chr_prefix + region.replace('chr', '')
                results += pysam.depth(*(bams + args + ['-r', region]))

        headers = ['Chromosome', 'Position']
        dtype = {'Chromosome': str,'Position': int}

        for i, bam in enumerate(bams):
            if names:
                name = names[i]
            else:
                samples = pybam.tag_sm(bam)
                if not samples:
                    basename = Path(bam).stem
                    message = (
                        f'SM tags were not found for {bam}, will use '
                        f'file name as sample name ({basename})'
                    )
                    samples = [basename]
                    warnings.warn(message)
                if len(samples) > 1:
                    m = f'multiple sample names detected: {bam}'
                    raise ValueError(m)
                name = samples[0]
            headers.append(name)
            dtype[name] = int

        df = pd.read_csv(
            StringIO(results), sep='\t', header=None,
            names=headers, dtype=dtype
        )

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
            Construct CovFrame from BAM files.
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
    def from_file(cls, fn, compression=False):
        """
        Construct CovFrame from a TSV file containing read depth data.

        Parameters
        ----------
        fn : str or file-like object
            TSV file (compressed or uncompressed). By file-like object, we
            refer to objects with a :meth:`read()` method, such as a file
            handle.
        compression : bool, default: False
            If True, use GZIP decompression regardless of filename.

        Returns
        -------
        CovFrame
            CovFrame object.

        See Also
        --------
        CovFrame
            CovFrame object creation using constructor.
        CovFrame.from_bam
            Construct CovFrame from BAM files.
        CovFrame.from_dict
            Construct CovFrame from dict of array-like or dicts.

        Examples
        --------

        >>> from fuc import pycov
        >>> cf = pycov.CovFrame.from_file('unzipped.tsv')
        >>> cf = pycov.CovFrame.from_file('zipped.tsv.gz')
        >>> cf = pycov.CovFrame.from_file('zipped.tsv', compression=True)
        """
        if isinstance(fn, IOBase):
            return cls(pd.read_table(fn))

        if fn.endswith('.gz') or compression:
            f = gzip.open(fn, 'rt')
        else:
            f = open(fn)

        headers = f.readline().strip().split('\t')

        f.close()

        if 'Chromosome' not in headers:
            raise ValueError("Input file is missing 'Chromosome' column")

        if 'Position' not in headers:
            raise ValueError("Input file is missing 'Position' column")

        dtype = {}

        for header in headers:
            if header == 'Chromosome':
                dtype[header] = str
            elif header == 'Position':
                dtype[header] = float
            else:
                dtype[header] = float

        return cls(pd.read_table(fn, dtype=dtype))

    def plot_region(
        self, sample, region=None, samples=None, label=None, ax=None, figsize=None,
        **kwargs
    ):
        """
        Create read depth profile for specified region.

        Region can be omitted if there is only one contig in the CovFrame.

        Parameters
        ----------
        region : str, optional
            Target region ('chrom:start-end').
        label : str, optional
            Label to use for the data points.
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
        Below is a simple example:

        .. plot::
            :context: close-figs

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
            >>> ax = cf.plot_region('A')
            >>> plt.tight_layout()

        We can draw multiple profiles in one plot:

        .. plot::
            :context: close-figs

            >>> ax = cf.plot_region('A', label='A')
            >>> cf.plot_region('B', label='B', ax=ax)
            >>> ax.legend()
            >>> plt.tight_layout()
        """
        if self.df.empty:
            raise ValueError('CovFrame is empty')

        if region is None:
            if len(self.contigs) == 1:
                cf = self.copy()
            else:
                raise ValueError('Multiple contigs found')
        else:
            cf = self.slice(region)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(
            data=cf.df, x='Position', y=sample, ax=ax, label=label, **kwargs
        )

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
        if not pd.isna(start):
            df = df[df.Position >= start]
        if not pd.isna(end):
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

    def matrix_uniformity(
        self, frac=0.1, n=20, m=None
    ):
        """
        Compute a matrix of fraction of sampled bases >= coverage with a
        shape of (coverages, samples).

        Parameters
        ----------
        frac : float, default: 0.1
            Fraction of data to be sampled (to speed up the process).
        n : int or list, default: 20
            Number of evenly spaced points to generate for the x-axis.
            Alternatively, positions can be manually specified by providing
            a list.
        m : float, optional
            Maximum point in the x-axis. By default, it will be the maximum
            depth in the entire dataset.

        Returns
        -------
        pandas.DataFrame
            Matrix of fraction of sampled bases >= coverage.

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
        >>> cf.matrix_uniformity()
                      A     B
        Coverage
        1.000000   1.00  1.00
        3.368421   1.00  1.00
        5.736842   1.00  1.00
        8.105263   1.00  1.00
        10.473684  1.00  1.00
        12.842105  1.00  0.98
        15.210526  1.00  0.93
        17.578947  1.00  0.87
        19.947368  1.00  0.77
        22.315789  1.00  0.64
        24.684211  1.00  0.50
        27.052632  0.97  0.35
        29.421053  0.84  0.25
        31.789474  0.70  0.16
        34.157895  0.51  0.07
        36.526316  0.37  0.07
        38.894737  0.21  0.03
        41.263158  0.09  0.02
        43.631579  0.04  0.00
        46.000000  0.02  0.00
        """
        # Sample positions to speed up the process.
        df = self.df.sample(frac=frac)

        # Determine x-axis points.
        if isinstance(n, list):
            coverages = n
        else:
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
        df = df.set_index("Coverage")

        return df

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
        n : int or list, default: 20
            Number of evenly spaced points to generate for the x-axis.
            Alternatively, positions can be manually specified by providing
            a list.
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
        df = self.matrix_uniformity(frac=frac, n=n, m=m)
        df = df.reset_index()
        df = df.melt(id_vars=['Coverage'], var_name='Sample',
            value_name='Fraction')

        # Determines how to display the lines.
        if mode == 'aggregated':
            hue = None
        else:
            hue = 'Sample'

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(
            x='Coverage', y='Fraction', data=df, hue=hue, marker=marker,
            ax=ax, **kwargs
        )

        ax.set_xlabel('Coverage')
        ax.set_ylabel('Fraction of sampled bases >= coverage')

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

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(
            x='Coverage', y='Fraction', data=df, hue=hue, ax=ax, **kwargs
        )

        ax.set_xlabel('Coverage')
        ax.set_ylabel('Fraction of sampled bases')

        return ax

    def copy_df(self):
        """Return a copy of the dataframe."""
        return self.df.copy()

    def copy(self):
        """Return a copy of the CovFrame."""
        return self.__class__(self.copy_df())

    def to_string(self):
        """
        Render the CovFrame to a console-friendly tabular output.

        Returns
        -------
        str
            String representation of the CovFrame.
        """
        return self.df.to_csv(index=False, sep='\t')

    def to_file(self, fn, compression=False):
        """
        Write the CovFrame to a TSV file.

        If the file name ends with '.gz', the method will automatically
        use the GZIP compression when writing the file.

        Parameters
        ----------
        fn : str
            TSV file (compressed or uncompressed).
        compression : bool, default: False
            If True, use the GZIP compression.
        """
        if fn.endswith('.gz') or compression:
            f = gzip.open(fn, 'wt')
        else:
            f = open(fn, 'w')
        f.write(self.to_string())
        f.close()

    def mask_bed(self, bed, opposite=False):
        """
        Mask rows that overlap with BED data.

        Parameters
        ----------
        bed : pybed.BedFrame or str
            BedFrame object or BED file.
        opposite : bool, default: False
            If True, mask rows that don't overlap with BED data.

        Returns
        -------
        CovFrame
            Masked CovFrame.

        Examples
        --------
        Assume we have the following data:

        >>> import numpy as np
        >>> from fuc import pycov, pybed
        >>> data = {
        ...     'Chromosome': ['chr1'] * 1000,
        ...     'Position': np.arange(1000, 2000),
        ...     'A': pycov.simulate(loc=35, scale=5),
        ...     'B': pycov.simulate(loc=25, scale=7),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df.head()
          Chromosome  Position   A   B
        0       chr1      1000  34  31
        1       chr1      1001  31  20
        2       chr1      1002  41  22
        3       chr1      1003  28  41
        4       chr1      1004  34  23
        >>> data = {
        ...     'Chromosome': ['chr1', 'chr1'],
        ...     'Start': [1000, 1003],
        ...     'End': [1002, 1004]
        ... }
        >>> bf = pybed.BedFrame.from_dict([], data)
        >>> bf.gr.df
          Chromosome  Start   End
        0       chr1   1000  1002
        1       chr1   1003  1004

        We can mask rows that overlap with the BED data:

        >>> cf.mask_bed(bf).df.head()
          Chromosome  Position     A     B
        0       chr1      1000   NaN   NaN
        1       chr1      1001   NaN   NaN
        2       chr1      1002  41.0  22.0
        3       chr1      1003   NaN   NaN
        4       chr1      1004  34.0  23.0

        We can also do the opposite:

        >>> cf.mask_bed(bf, opposite=True).df.head()
          Chromosome  Position     A     B
        0       chr1      1000  34.0  31.0
        1       chr1      1001  31.0  20.0
        2       chr1      1002   NaN   NaN
        3       chr1      1003  28.0  41.0
        4       chr1      1004   NaN   NaN
        """
        if isinstance(bed, pybed.BedFrame):
            bf = bed
        else:
            bf = pybed.BedFrame.from_file(bed)
        def one_row(r):
            if opposite:
                if bf.gr[r.Chromosome, r.Position:r.Position+1].empty:
                    r[2:] = r[2:].apply(lambda x: np.nan)
            else:
                if not bf.gr[r.Chromosome, r.Position:r.Position+1].empty:
                    r[2:] = r[2:].apply(lambda x: np.nan)
            return r
        df = self.df.apply(one_row, axis=1)
        return self.__class__(df)

    def update_chr_prefix(self, mode='remove'):
        """
        Add or remove the (annoying) 'chr' string from the Chromosome column.

        Parameters
        ----------
        mode : {'add', 'remove'}, default: 'remove'
            Whether to add or remove the 'chr' string.

        Returns
        -------
        CovFrame
            Updated CovFrame.

        Examples
        --------

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1'] * 3 + ['2'] * 3,
        ...     'Position': np.arange(1, 7),
        ...     'A': pycov.simulate(loc=35, scale=5, size=6),
        ...     'B': pycov.simulate(loc=25, scale=7, size=6),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df
          Chromosome  Position   A   B
        0       chr1         1  35  25
        1       chr1         2  23  14
        2       chr1         3  32  23
        3          2         4  38  25
        4          2         5  33   8
        5          2         6  21  22
        >>> cf.update_chr_prefix(mode='remove').df
          Chromosome  Position   A   B
        0          1         1  35  25
        1          1         2  23  14
        2          1         3  32  23
        3          2         4  38  25
        4          2         5  33   8
        5          2         6  21  22
        >>> cf.update_chr_prefix(mode='add').df
          Chromosome  Position   A   B
        0       chr1         1  35  25
        1       chr1         2  23  14
        2       chr1         3  32  23
        3       chr2         4  38  25
        4       chr2         5  33   8
        5       chr2         6  21  22
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
        df = self.df.apply(one_row, axis=1)
        return self.__class__(df)

    def subset(self, samples, exclude=False):
        """
        Subset CovFrame for specified samples.

        Parameters
        ----------
        samples : str or list
            Sample name or list of names (the order matters).
        exclude : bool, default: False
            If True, exclude specified samples.

        Returns
        -------
        CovFrame
            Subsetted CovFrame.

        Examples
        --------
        Assume we have the following data:

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1'] * 1000,
        ...     'Position': np.arange(1000, 2000),
        ...     'A': pycov.simulate(loc=35, scale=5),
        ...     'B': pycov.simulate(loc=25, scale=7),
        ...     'C': pycov.simulate(loc=15, scale=2),
        ...     'D': pycov.simulate(loc=45, scale=8),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df.head()
          Chromosome  Position   A   B   C   D
        0       chr1      1000  30  30  15  37
        1       chr1      1001  25  24  11  43
        2       chr1      1002  33  24  16  50
        3       chr1      1003  29  22  15  46
        4       chr1      1004  34  30  11  32

        We can subset the CovFrame for the samples A and B:

        >>> cf.subset(['A', 'B']).df.head()
          Chromosome  Position   A   B
        0       chr1      1000  30  30
        1       chr1      1001  25  24
        2       chr1      1002  33  24
        3       chr1      1003  29  22
        4       chr1      1004  34  30

        Alternatively, we can exclude those samples:

        >>> cf.subset(['A', 'B'], exclude=True).df.head()
          Chromosome  Position   C   D
        0       chr1      1000  15  37
        1       chr1      1001  11  43
        2       chr1      1002  16  50
        3       chr1      1003  15  46
        4       chr1      1004  11  32
        """
        if isinstance(samples, str):
            samples = [samples]
        if exclude:
            samples = [x for x in self.samples if x not in samples]
        cols = self.df.columns[:2].to_list() + samples
        return self.__class__(self.df[cols])

    def rename(self, names, indicies=None):
        """
        Rename the samples.

        Parameters
        ----------
        names : dict or list
            Dict of old names to new names or list of new names.
        indicies : list or tuple, optional
            List of 0-based sample indicies. Alternatively, a tuple
            (int, int) can be used to specify an index range.

        Returns
        -------
        CovFrame
            Updated CovFrame.

        Examples
        --------

        >>> import numpy as np
        >>> from fuc import pycov
        >>> data = {
        ...     'Chromosome': ['chr1'] * 2,
        ...     'Position': np.arange(1, 3),
        ...     'A': pycov.simulate(loc=35, scale=5, size=2),
        ...     'B': pycov.simulate(loc=25, scale=7, size=2),
        ...     'C': pycov.simulate(loc=25, scale=7, size=2),
        ...     'D': pycov.simulate(loc=25, scale=7, size=2),
        ... }
        >>> cf = pycov.CovFrame.from_dict(data)
        >>> cf.df
          Chromosome  Position   A   B   C   D
        0       chr1         1  31  19  28  15
        1       chr1         2  35  24  22  17
        >>> cf.rename(['1', '2', '3', '4']).df
          Chromosome  Position   1   2   3   4
        0       chr1         1  31  19  28  15
        1       chr1         2  35  24  22  17
        >>> cf.rename({'B': '2', 'C': '3'}).df
          Chromosome  Position   A   2   3   D
        0       chr1         1  31  19  28  15
        1       chr1         2  35  24  22  17
        >>> cf.rename(['2', '4'], indicies=[1, 3]).df
          Chromosome  Position   A   2   C   4
        0       chr1         1  31  19  28  15
        1       chr1         2  35  24  22  17
        >>> cf.rename(['2', '3'], indicies=(1, 3)).df
          Chromosome  Position   A   2   3   D
        0       chr1         1  31  19  28  15
        1       chr1         2  35  24  22  17
        """
        samples = common.rename(self.samples, names, indicies=indicies)
        columns = self.df.columns[:2].to_list() + samples
        cf = self.copy()
        cf.df.columns = columns
        return cf
