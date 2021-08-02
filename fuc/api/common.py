"""
The common submodule is used by other fuc submodules such as pyvcf and
pybed. It also provides many day-to-day actions used in the field of
bioinformatics.
"""

import pathlib
import re
import os
import warnings
import inspect
from argparse import RawTextHelpFormatter, SUPPRESS
from pathlib import Path
from difflib import SequenceMatcher
from urllib.request import urlretrieve

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import matplotlib.patches as mpatches
import seaborn as sns

FUC_PATH = pathlib.Path(__file__).parent.parent.parent.absolute()

class AnnFrame:
    """
    Class for storing sample annotation data.

    This class stores sample annotation data as :class:`pandas.DataFrame`
    with sample names as index.

    Note that an AnnFrame can have a different set of samples than its
    accompanying MafFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing sample annotation data. The index must be
        sample names.

    See Also
    --------
    AnnFrame.from_dict
        Construct AnnFrame from dict of array-like or dicts.
    AnnFrame.from_file
        Construct AnnFrame from a delimited text file.

    Examples
    --------

    >>> import pandas as pd
    >>> from fuc import pymaf, common
    >>> data = {
    ...     'Tumor_Sample_Barcode': ['Steven_N', 'Steven_T', 'Sara_N', 'Sara_T'],
    ...     'Subject': ['Steven', 'Steven', 'Sara', 'Sara'],
    ...     'Type': ['Normal', 'Tumor', 'Normal', 'Tumor'],
    ...     'Age': [30, 30, 57, 57]
    ... }
    >>> df = pd.DataFrame(data)
    >>> df = df.set_index('Tumor_Sample_Barcode')
    >>> af = common.AnnFrame(df)
    >>> af.df
                         Subject    Type  Age
    Tumor_Sample_Barcode
    Steven_N              Steven  Normal   30
    Steven_T              Steven   Tumor   30
    Sara_N                  Sara  Normal   57
    Sara_T                  Sara   Tumor   57
    """

    def _check_df(self, df):
        if type(df.index) == pd.RangeIndex:
            m = "Index must be sample names, not 'pandas.RangeIndex'."
            raise ValueError(m)
        if df.isin([np.inf, -np.inf]).any().any():
            raise ValueError('Found positive or negative infinity.')
        return df

    def __init__(self, df):
        self._df = self._check_df(df)

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing sample annotation data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = self._check_df(value)

    @property
    def samples(self):
        """list : List of the sample names."""
        return list(self.df.index.to_list())

    @property
    def shape(self):
        """tuple : Dimensionality of AnnFrame (samples, annotations)."""
        return self.df.shape

    def filter_mf(self, mf):
        """
        Filter the AnnFrame for the samples in the MafFrame.

        Parameters
        ----------
        mf : MafFrame
            MafFrame containing target samples.

        Returns
        -------
        AnnFrame
            Filtered AnnFrame object.
        """
        df = self.df.loc[mf.samples]
        return self.__class__(df)

    @classmethod
    def from_dict(cls, data, sample_col='Tumor_Sample_Barcode'):
        """Construct AnnFrame from dict of array-like or dicts.

        The dictionary must have at least one column that represents sample
        names which are used as index for pandas.DataFrame.

        Parameters
        ----------
        data : dict
            Of the form {field : array-like} or {field : dict}.
        sample_col : str, default: 'Tumor_Sample_Barcode'
            Column containing sample names.

        Returns
        -------
        AnnFrame
            AnnFrame object.

        See Also
        --------
        AnnFrame
            AnnFrame object creation using constructor.
        AnnFrame.from_file
            Construct AnnFrame from a delimited text file.

        Examples
        --------

        >>> from fuc import pymaf, common
        >>> data = {
        ...     'Tumor_Sample_Barcode': ['Steven_Normal', 'Steven_Tumor', 'Sara_Normal', 'Sara_Tumor'],
        ...     'Subject': ['Steven', 'Steven', 'Sara', 'Sara'],
        ...     'Type': ['Normal', 'Tumor', 'Normal', 'Tumor'],
        ...     'Age': [30, 30, 57, 57]
        ... }
        >>> af = common.AnnFrame.from_dict(data)
        >>> af.df
                             Subject    Type  Age
        Tumor_Sample_Barcode
        Steven_Normal         Steven  Normal   30
        Steven_Tumor          Steven   Tumor   30
        Sara_Normal             Sara  Normal   57
        Sara_Tumor              Sara   Tumor   57
        """
        df = pd.DataFrame(data)
        df = df.set_index(sample_col)
        return cls(df)

    @classmethod
    def from_file(cls, fn, sample_col='Tumor_Sample_Barcode', sep='\t'):
        """
        Construct AnnFrame from a delimited text file.

        The text file must have at least one column that represents
        sample names which are used as index for pandas.DataFrame.

        Parameters
        ----------
        fn : str
            Text file path (zipped or unzipped).
        sample_col : str, default: 'Tumor_Sample_Barcode'
            Column containing sample names.
        sep : str, default: '\\\\t'
            Delimiter to use.

        Returns
        -------
        AnnFrame
            AnnFrame.

        See Also
        --------
        AnnFrame
            AnnFrame object creation using constructor.
        AnnFrame.from_dict
            Construct AnnFrame from dict of array-like or dicts.

        Examples
        --------

        >>> from fuc import pymaf, common
        >>> af1 = common.AnnFrame.from_file('sample-annot-1.tsv')
        >>> af2 = common.AnnFrame.from_file('sample-annot-2.csv', sample_col='SampleID', sep=',')
        """
        df = pd.read_table(fn, sep=sep)
        df = df.set_index(sample_col)
        return cls(df)

    def legend_handles(
        self, col, samples=None, numeric=False, segments=5, decimals=0,
        cmap='Pastel1'
    ):
        """Return legend handles for the given column.

        In the case of a numeric column, use ``numeric=True`` which will
        divide the values into equal-sized intervals, with the number of
        intervals determined by the `segments` option.

        Parameters
        ----------
        col : str
            Column name.
        samples : list, optional
            If provided, these samples will be used to create legend handles.
        numeric : bool, default: False
            If True, the column will be treated as numeric.
        segments : int, default: 5
            If ``numeric`` is True, the numbers will be divided
            into this number of equal-sized segments.
        decimals : int, default: 0
            If ``numeric`` is True, the numbers will be rounded up to this
            number of decimals.
        cmap : str, default: 'Pastel1'
            Color map.

        Returns
        -------
        list
            Legend handles.

        See Also
        --------
        AnnFrame.plot_annot
            Create a 1D categorical heatmap for the given column.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(f)
            >>> fig, ax = plt.subplots()
            >>> handles1 = af.legend_handles('FAB_classification',
            ...                              cmap='Dark2')
            >>> handles2 = af.legend_handles('days_to_last_followup',
            ...                              numeric=True,
            ...                              cmap='viridis')
            >>> handles3 = af.legend_handles('Overall_Survival_Status')
            >>> leg1 = ax.legend(handles=handles1,
            ...                  title='FAB_classification',
            ...                  loc='center left')
            >>> leg2 = ax.legend(handles=handles2,
            ...                  title='days_to_last_followup',
            ...                  loc='center')
            >>> leg3 = ax.legend(handles=handles3,
            ...                  title='Overall_Survival_Status',
            ...                  loc='center right')
            >>> ax.add_artist(leg1)
            >>> ax.add_artist(leg2)
            >>> ax.add_artist(leg3)
            >>> plt.tight_layout()
        """
        s, l = self._get_col(col, numeric=numeric, segments=segments,
            samples=samples)
        colors = plt.cm.get_cmap(cmap)(np.linspace(0, 1, len(l)))
        handles = []
        for i, label in enumerate(l):
            handles.append(mpatches.Patch(color=colors[i], label=label))
        return handles

    def _get_col(
        self, col, numeric=False, segments=5, samples=None, decimals=0
    ):
        s = self.df[col]
        if numeric:
            boundaries = list(np.linspace(s.min(), s.max(),
                segments+1, endpoint=True))
            intervals = list(zip(boundaries[:-1], boundaries[1:]))
            def f(x):
                if pd.isna(x):
                    return x
                for i, interval in enumerate(intervals):
                    a, b = interval
                    if a <= x <= b:
                        return f'G{i} ({b:.{decimals}f})'
            s = s.apply(f)
        if samples is not None:
            s = s[samples]
        l = sorted([x for x in s.unique() if x == x])
        return s, l

    def plot_annot(
        self, col, samples=None, numeric=False, segments=5, decimals=0,
        cmap='Pastel1', ax=None, figsize=None, **kwargs
    ):
        """
        Create a 1D categorical heatmap of the column.

        In the case of a numeric column, set ``numeric`` as True which will
        divide the values into equal-sized intervals, with the number of
        intervals determined by ``segments``.

        See the :ref:`tutorials:Create customized oncoplots` tutorial to
        learn how to create customized oncoplots.

        Parameters
        ----------
        col : str
            Column name.
        samples : list, optional
            If provided, these samples will be used to create legend handles.
        numeric : bool, default: False
            If True, the column will be treated as numeric.
        segments : int, default: 5
            If ``numeric`` is True, the numbers will be divided
            into this number of equal-sized segments.
        decimals : int, default: 0
            If ``numeric`` is True, the numbers will be rounded up to this
            number of decimals.
        cmap : str, default: 'Pastel1'
            Color map.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        list
            Legend handles.

        See Also
        --------
        AnnFrame.legend_handles
            Return legend handles for the given column.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(f)
            >>> fig, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(10, 5))
            >>> af.plot_annot('FAB_classification', ax=ax1, linewidths=1, cmap='Dark2')
            >>> af.plot_annot('days_to_last_followup', ax=ax2, linewidths=1, cmap='viridis')
            >>> af.plot_annot('Overall_Survival_Status', ax=ax3, linewidths=1)
            >>> ax1.set_xlabel('')
            >>> ax2.set_xlabel('')
            >>> ax1.set_ylabel('')
            >>> ax2.set_ylabel('')
            >>> ax3.set_ylabel('')
            >>> plt.tight_layout()
            >>> plt.subplots_adjust(wspace=0.01, hspace=0.01)
        """
        s, l = self._get_col(col, numeric=numeric, segments=segments,
            samples=samples)
        d = {k: v for v, k in enumerate(l)}
        df = s.to_frame().applymap(lambda x: x if pd.isna(x) else d[x])

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(df.T, ax=ax, cmap=cmap, cbar=False, **kwargs)

        ax.set_xlabel('Samples')
        ax.set_ylabel(col)
        ax.set_xticks([])
        ax.set_yticks([])

        return ax

    def sorted_samples(self, by, mf=None, keep_empty=False, nonsyn=False):
        """
        Return a sorted list of sample names.

        Parameters
        ----------
        df : str or list
            Column or list of columns to sort by.
        """
        df = self.df.copy()

        if nonsyn:
            samples = mf.matrix_waterfall(keep_empty=keep_empty).columns
            df = df.loc[samples]

        df = df.sort_values(by=by)

        return df.index.to_list()

def _script_name():
    """Return the current script's filename."""
    fn = inspect.stack()[1].filename
    return pathlib.Path(fn).stem.replace('_', '-')

def _add_parser(subparsers, name, **kwargs):
    """Return the pre-formatted parser."""
    parser = subparsers.add_parser(
        name,
        add_help=False,
        formatter_class=RawTextHelpFormatter,
        **kwargs,
    )
    parser._positionals.title = 'Positional arguments'
    parser._optionals.title = 'Optional arguments'
    parser.add_argument(
        '-h',
        '--help',
        action='help',
        default=SUPPRESS,
        help='Show this help message and exit.',
    )
    return parser

def get_similarity(a, b):
    """Return a value from 0 to 1 representing how similar two strings are."""
    return SequenceMatcher(None, a, b).ratio()

def is_similar(a, b, threshold=0.9):
    """Return True if the similarity is equal to or greater than threshold."""
    return get_similarity(a, b) >= threshold

def get_most_similar(a, l):
    """Return the most similar string in a list."""
    s = [get_similarity(a, x) for x in l]
    m = max(s)
    i = [i for i, x in enumerate(s) if x == m][0]
    return l[i]

def sumstat(fp, fn, tp, tn):
    """
    Return various summary statistics from (FP, FN, TP, TN).

    This method will return the following statistics:

    +------------------------------------------------------------+-------------------------------------------------+
    | Terminology                                                | Derivation                                      |
    +============================================================+=================================================+
    | sensitivity, recall, hit rate, or true positive rate (TPR) | :math:`TPR = TP / P = TP / (TP + FN) = 1 - FNR` |
    +------------------------------------------------------------+-------------------------------------------------+
    | specificity, selectivity or true negative rate (TNR)       | :math:`TNR = TN / N = TN / (TN + FP) = 1 - FPR` |
    +------------------------------------------------------------+-------------------------------------------------+
    | precision or positive predictive value (PPV)               | :math:`PPV = TP / (TP + FP) = 1 - FDR`          |
    +------------------------------------------------------------+-------------------------------------------------+
    | negative predictive value (NPV)                            | :math:`NPV = TN / (TN + FN) = 1 - FOR`          |
    +------------------------------------------------------------+-------------------------------------------------+
    | miss rate or false negative rate (FNR)                     | :math:`FNR = FN / P = FN / (FN + TP) = 1 - TPR` |
    +------------------------------------------------------------+-------------------------------------------------+
    | fall-out or false positive rate (FPR)                      | :math:`FPR = FP / N = FP / (FP + TN) = 1 - TNR` |
    +------------------------------------------------------------+-------------------------------------------------+
    | false discovery rate (FDR)                                 | :math:`FDR = FP / (FP + TP) = 1 - PPV`          |
    +------------------------------------------------------------+-------------------------------------------------+
    | false omission rate (FOR)                                  | :math:`FOR = FN / (FN + TN) = 1 - NPV`          |
    +------------------------------------------------------------+-------------------------------------------------+
    | accuracy (ACC)                                             | :math:`ACC = (TP + TN)/(TP + TN + FP + FN)`     |
    +------------------------------------------------------------+-------------------------------------------------+

    Parameters
    ----------
    fp, fn, tp, tn : int
        Input statistics.

    Returns
    -------
    dict
        Dictionary containing summary statistics.

    Examples
    --------

    This example is taken from the Wiki page `Sensitivity and specificity <https://en.wikipedia.org/wiki/Sensitivity_and_specificity>`__.

    >>> results = common.sumstat(180, 10, 20, 1820)
    >>> for k, v in results.items():
    ...     print(k, v)
    ...
    tpr 0.6666666666666666
    tnr 0.91
    ppv 0.1
    npv 0.994535519125683
    fnr 0.3333333333333333
    fpr 0.09
    fdr 0.9
    for 0.00546448087431694
    acc 0.9064039408866995
    """
    data = {
        'tpr': tp / (tp + fn), # sensitivity, recall, hit rate
        'tnr': tn / (tn + fp), # specificity, selectivity
        'ppv': tp / (tp + fp), # precision
        'npv': tn / (tn + fn),
        'fnr': fn / (fn + tp), # miss rate
        'fpr': fp / (fp + tn), # fall-out rate
        'fdr': fp / (fp + tp),
        'for': fn / (fn + tn),
        'acc': (tp + tn) / (tp + tn + fp + fn),
    }
    return data

def load_dataset(name, force=False):
    """
    Load an example dataset from the online repository (requires internet).

    Parameters
    ----------
    name : str
        Name of the dataset in https://github.com/sbslee/fuc-data.
    force : bool, default: False
        If True, overwrite the existing files.
    """
    home_dir = str(Path.home()) + '/fuc-data'
    data_dir = f'{home_dir}/{name}'
    if not os.path.exists(home_dir):
        os.makedirs(home_dir)
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    datasets = {
        'tcga-laml': [
            'tcga_cohort.txt.gz',
            'tcga_laml.maf.gz',
            'tcga_laml_annot.tsv',
            'tcga_laml.vcf',
            'tcga_laml_vep.vcf',
        ],
        'brca': [
            'brca.maf.gz',
            'brca.vcf',
        ],
        'pyvcf': [
            'plot_comparison.vcf',
            'normal-tumor.vcf',
            'normal-tumor-annot.tsv',
        ],
        'cytoband': [
            'cytoBandIdeo.txt.gz',
            'ucsc_genes.bed.gz',
        ],
    }
    base_url = ('https://raw.githubusercontent.com/sbslee/fuc-data/main')
    for f in datasets[name]:
        file_url = f'{base_url}/{name}/{f}'
        file_path = f'{data_dir}/{f}'
        download = False
        if force:
            download = True
        else:
            if os.path.exists(file_path):
                pass
            else:
                download = True
        if download:
            urlretrieve(file_url, file_path)

def parse_region(region):
    """
    Parse the region.

    Parameters
    ----------
    region : str
        Region ('chrom:start-end').

    Returns
    -------
    tuple
        The output tuple will always have a shape of (chrom, start, end) with
        the following data types: (str, int, int).

    Examples
    --------

    >>> from fuc import common
    >>> common.parse_region('chr1:100-150')
    ('chr1', 100, 150)
    >>> common.parse_region('chr1')
    ('chr1', 0, 0)
    >>> common.parse_region('chr1:100')
    ('chr1', 100, 0)
    >>> common.parse_region('chr1:100-')
    ('chr1', 100, 0)
    >>> common.parse_region('chr1:-100')
    ('chr1', 0, 100)
    """
    chrom = region.split(':')[0]

    try:
        start = region.split(':')[1].split('-')[0]
        if not start:
            start = 0
        else:
            start = int(start)
    except IndexError:
        start = 0

    try:
        end = region.split(':')[1].split('-')[1]
        if not end:
            end = 0
        else:
            end = int(end)
    except IndexError:
        end = 0

    return (chrom, start, end)

def extract_sequence(fasta, region):
    """
    Extract the region's DNA sequence from the FASTA file.

    Parameters
    ----------
    fasta : str
        FASTA file.
    region : str
        Region ('chrom:start-end').

    Returns
    -------
    str
        DNA sequence. Empty string if there is no matching sequence.
    """
    try:
        sequence = pysam.faidx(fasta, region).split('\n')[1]
    except pysam.SamtoolsError as e:
        warnings.warn(str(e))
        sequence = ''
    return sequence

def plot_cytobands(cytoband, bed, ax=None, figsize=None):
    """
    Create chromosome ideograms along with BED data.

    The method's source code is derived from a Python script (ideograms.py)
    written by Ryan Dale. The original script can be found at:
    https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py

    Parameters
    ----------
    cytoband : str
        Text file containing cytoband ideogram information.
    bed : str
        BED file to be displayed.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, crete a new one.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Examples
    --------

    .. plot::
        :context: close-figs

        >>> import matplotlib.pyplot as plt
        >>> from fuc import common
        >>> common.load_dataset('cytoband')
        >>> cytoband_file = '~/fuc-data/cytoband/cytoBandIdeo.txt.gz'
        >>> bed_file = '~/fuc-data/cytoband/ucsc_genes.bed.gz'
        >>> common.plot_cytobands(cytoband_file, bed_file, figsize=(10, 8))
    """
    def chromosome_collections(df, y_positions, height, **kwargs):
        del_width = False
        if 'width' not in df.columns:
            del_width = True
            df['width'] = df['end'] - df['start']
        for chrom, group in df.groupby('chrom'):
            yrange = (y_positions[chrom], height)
            xranges = group[['start', 'width']].values
            yield BrokenBarHCollection(
                xranges, yrange, edgecolors=("black",), facecolors=group['colors'], **kwargs)
        if del_width:
            del df['width']

    # Height of each ideogram
    chrom_height = 1

    # Spacing between consecutive ideograms
    chrom_spacing = 1

    # Height of the gene track. Should be smaller than `chrom_spacing` in order to
    # fit correctly
    gene_height = 0.4

    # Padding between the top of a gene track and its corresponding ideogram
    gene_padding = 0.1

    # Decide which chromosomes to use
    chromosome_list = [f'chr{i}' for i in list(range(1, 23)) + ['M', 'X', 'Y']]

    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # Read in ideogram.txt, downloaded from UCSC Table Browser
    ideo = pd.read_table(
        cytoband,
        names=['chrom', 'start', 'end', 'name', 'gieStain']
    )

    # Filter out chromosomes not in our list
    ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Colors for different chromosome stains
    color_lookup = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': (.8, .4, .4),
        'gvar': (.8, .8, .8),
        'stalk': (.9, .9, .9),
    }

    # Add a new column for colors
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])

    # Same thing for genes
    genes = pd.read_table(
        bed,
        names=['chrom', 'start', 'end', 'name'],
        usecols=range(4))
    genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
    genes['width'] = genes.end - genes.start
    genes['colors'] = '#2243a8'

    # Determine which matplotlib axes to plot on.
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Now all we have to do is call our function for the ideogram data...
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    # ...and the gene data
    for collection in chromosome_collections(
        genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
    ):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')

    return ax

def file2list(fn):
    """
    Return a list of filenames from the input file.

    Parameters
    ----------
    fn : str
        File containing one filename per line.

    Returns
    -------
    list
        List of filenames.

    Examples
    --------

    >>> from fuc import common
    >>> common.file2list('bam.list')
    ['1.bam', '2.bam', '3.bam']
    """
    l = []
    with open(fn) as f:
        for line in f:
            l.append(line.strip())
    return l

def legend_handles(labels, colors='tab10'):
    """
    Create custom legend handles.

    Parameters
    ----------
    labels : list
        List of labels.
    colors : str or list, default: 'tab10'
        Colormap name or list of colors.

    Returns
    -------
    list
        List of legend handles.

    Examples
    --------

    .. plot::

        >>> import matplotlib.pyplot as plt
        >>> from fuc import common
        >>> fig, ax = plt.subplots()
        >>> handles1 = common.legend_handles(['A', 'B'], colors='tab10')
        >>> handles2 = common.legend_handles(['C', 'D'], colors=['yellow', 'green'])
        >>> legend1 = ax.legend(handles=handles1, loc='center left')
        >>> legend2 = ax.legend(handles=handles2)
        >>> ax.add_artist(legend1)
        >>> ax.add_artist(legend2)
        >>> plt.tight_layout()
    """
    if isinstance(colors, str):
        colors = plt.get_cmap(colors).colors
    elif isinstance(colors, list):
        pass
    else:
        raise TypeError(f'Incorrect type of colors: {type(colors)}')

    handles = []

    for i, label in enumerate(labels):
        handles.append(mpatches.Patch(color=colors[i], label=label))

    return handles
