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

class Variant:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def __members(self):
        return (self.chrom, self.pos, self.ref, self.alt)

    def __eq__(self, other):
        if type(other) is type(self):
            return self.__members() == other.__members()
        else:
            return False

    def __hash__(self):
        return hash(self.__members())

    def __repr__(self):
        s = ', '.join([str(x) for x in self.__members()])
        return f'Variant({s})'

class AnnFrame:
    """
    Class for storing sample annotation data.

    This class stores sample annotation data as :class:`pandas.DataFrame`
    with sample names as index.

    Note that an AnnFrame can have a different set of samples than its
    accompanying :class:`pymaf.MafFrame`, :class:`pyvcf.VcfFrame`, etc.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing sample annotation data. The index must be
        unique sample names.

    See Also
    --------
    AnnFrame.from_dict
        Construct AnnFrame from dict of array-like or dicts.
    AnnFrame.from_file
        Construct AnnFrame from a delimited text file.

    Examples
    --------

    >>> import pandas as pd
    >>> from fuc import common
    >>> data = {
    ...     'SampleID': ['A', 'B', 'C', 'D'],
    ...     'PatientID': ['P1', 'P1', 'P2', 'P2'],
    ...     'Tissue': ['Normal', 'Tissue', 'Normal', 'Tumor'],
    ...     'Age': [30, 30, 57, 57]
    ... }
    >>> df = pd.DataFrame(data)
    >>> df = df.set_index('SampleID')
    >>> af = common.AnnFrame(df)
    >>> af.df
             PatientID  Tissue  Age
    SampleID
    A               P1  Normal   30
    B               P1  Tissue   30
    C               P2  Normal   57
    D               P2   Tumor   57
    """

    def _check_df(self, df):
        if type(df.index) == pd.RangeIndex:
            raise ValueError("Index cannot be 'pandas.RangeIndex'.")
        if df.isin([np.inf, -np.inf]).any().any():
            raise ValueError('Found positive or negative infinity.')
        if df.index.has_duplicates:
            raise ValueError('Index has duplicates.')
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

    @classmethod
    def from_dict(cls, data, sample_col):
        """
        Construct AnnFrame from dict of array-like or dicts.

        The dictionary must contain a column that represents sample names.

        Parameters
        ----------
        data : dict
            Of the form {field : array-like} or {field : dict}.
        sample_col : str or int
            Column containing unique sample names, either given as string
            name or column index.

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

        >>> from fuc import common
        >>> data = {
        ...     'SampleID': ['A', 'B', 'C', 'D'],
        ...     'PatientID': ['P1', 'P1', 'P2', 'P2'],
        ...     'Tissue': ['Normal', 'Tissue', 'Normal', 'Tumor'],
        ...     'Age': [30, 30, 57, 57]
        ... }
        >>> af = common.AnnFrame.from_dict(data, sample_col='SampleID') # or sample_col=0
        >>> af.df
                 PatientID  Tissue  Age
        SampleID
        A               P1  Normal   30
        B               P1  Tissue   30
        C               P2  Normal   57
        D               P2   Tumor   57
        """
        df = pd.DataFrame(data)
        if isinstance(sample_col, int):
            sample_col = list(data)[sample_col]
        df = df.set_index(sample_col)
        return cls(df)

    @classmethod
    def from_file(cls, fn, sample_col, sep='\t'):
        """
        Construct AnnFrame from a delimited text file.

        The file must contain a column that represents sample names.

        Parameters
        ----------
        fn : str
            Text file (zipped or unzipped).
        sample_col : str or int
            Column containing unique sample names, either given as string
            name or column index.
        sep : str, default: '\\\\t'
            Delimiter to use.

        Returns
        -------
        AnnFrame
            AnnFrame object.

        See Also
        --------
        AnnFrame
            AnnFrame object creation using constructor.
        AnnFrame.from_dict
            Construct AnnFrame from dict of array-like or dicts.

        Examples
        --------

        >>> from fuc import common
        >>> af = common.AnnFrame.from_file('sample-annot.tsv', sample_col='SampleID')
        >>> af = common.AnnFrame.from_file('sample-annot.csv', sample_col=0, sep=',')
        """
        df = pd.read_table(fn, index_col=sample_col, sep=sep)
        return cls(df)

    def plot_annot(
        self, group_col, group_order=None, samples=None, colors='tab10',
        sequential=False, xticklabels=True, ax=None, figsize=None
    ):
        """
        Create a categorical heatmap for the selected column using unmatched
        samples.

        See this :ref:`tutorial <tutorials:Create customized oncoplots>` to
        learn how to create customized oncoplots.

        Parameters
        ----------
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        samples : list, optional
            Display only specified samples (in that order too).
        colors : str or list, default: 'tab10'
            Colormap name or list of colors.
        sequential : bool, default: False
            Whether the column is sequential data.
        xticklabels : bool, default: True
            If True, plot the sample names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        list
            Legend handles.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> ax, handles = af.plot_annot('FAB_classification', samples=af.samples[:10])
            >>> legend = ax.legend(handles=handles)
            >>> ax.add_artist(legend)
            >>> plt.tight_layout()

        We can display only selected groups:

        .. plot::
            :context: close-figs

            >>> ax, handles = af.plot_annot('FAB_classification', group_order=['M7', 'M6'])
            >>> legend = ax.legend(handles=handles)
            >>> ax.add_artist(legend)
            >>> plt.tight_layout()

        We can also display sequenital data in the following way:

        .. plot::
            :context: close-figs

            >>> ax, handles = af.plot_annot('FAB_classification',
            ...                             samples=af.samples[:10],
            ...                             colors='viridis',
            ...                             sequential=True)
            >>> legend = ax.legend(handles=handles)
            >>> ax.add_artist(legend)
            >>> plt.tight_layout()
        """
        # Get the selected column.
        s = self.df[group_col]

        # Subset the samples, if necessary.
        if samples is not None:
            s = s.reindex(samples)

        # Establish mapping from groups to numbers.
        if group_order is None:
            group_order = sorted([x for x in s.unique() if x == x])
        else:
            s = s[s.isin(group_order)]
        d = {k: v for v, k in enumerate(group_order)}
        df = s.to_frame().applymap(lambda x: x if pd.isna(x) else d[x])

        # Determine the colors to use.
        if isinstance(colors, str):
            if sequential:
                c = plt.get_cmap(colors).colors
                l = list(np.linspace(0, len(c)-1, len(group_order)))
                colors = [c[round(i)] for i in l]
            else:
                colors = list(plt.get_cmap(colors).colors[:len(group_order)])

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        # Plot the heatmap.
        sns.heatmap(
            df.T, ax=ax, cmap=colors, xticklabels=xticklabels, cbar=False,
            linewidths=0.5
        )
        ax.set_xlabel('')
        ax.set_ylabel(group_col)
        ax.set_yticks([])

        # Get the legend handles.
        handles = legend_handles(group_order, colors=colors)

        return ax, handles

    def plot_annot_matched(
        self, patient_col, group_col, annot_col, patient_order=None,
        group_order=None, annot_order=None, colors='tab10', sequential=False,
        xticklabels=True, ax=None, figsize=None
    ):
        """
        Create a categorical heatmap for the selected column using matched
        samples.

        See this :ref:`tutorial <tutorials:Create customized oncoplots>` to
        learn how to create customized oncoplots.

        Parameters
        ----------
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        annot_col : str
            Column to plot.
        patient_order : list, optional
            Plot only specified patients (in that order too).
        group_order : list, optional
            List of sample group names.
        annot_order : list, optional
            Plot only specified annotations (in that order too).
        colors : str or list, default: 'tab10'
            Colormap name or list of colors.
        sequential : bool, default: False
            Whether the column is sequential data.
        xticklabels : bool, default: True
            If True, plot the sample names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        list
            Legend handles.
        """

        if annot_order is None:
            annot_order = self.df[annot_col].unique()

        df = self.df.pivot(columns=group_col, index=patient_col,
            values=annot_col).T

        if patient_order is not None:
            df = df[patient_order]

        d = {k: v for v, k in enumerate(annot_order)}
        df = df.applymap(lambda x: x if pd.isna(x) else d[x])

        if group_order is not None:
            df = df.loc[group_order]

        # Determine the colors to use.
        if isinstance(colors, str):
            if sequential:
                c = plt.get_cmap(colors).colors
                l = list(np.linspace(0, len(c)-1, len(group_order)))
                colors = [c[round(i)] for i in l]
            else:
                colors = list(plt.get_cmap(colors).colors[:len(group_order)])

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        # Plot the heatmap.
        sns.heatmap(
            df, cmap=colors, cbar=False, xticklabels=xticklabels, ax=ax
        )
        ax.set_xlabel('')
        ax.set_ylabel(annot_col)
        ax.set_yticks([])

        # Add vertical lines.
        for i, sample in enumerate(df.columns, start=1):
            ax.axvline(i, color='white')

        # Get the legend handles.
        handles = legend_handles(annot_order, colors=colors)

        return ax, handles

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

def convert_file2list(fn):
    """
    Convert a text file to a list of filenames.

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
    >>> common.convert_file2list('bam.list')
    ['1.bam', '2.bam', '3.bam']
    """
    l = []
    with open(fn) as f:
        for line in f:
            l.append(line.strip())
    return l

def convert_num2cat(s, n=5, decimals=0):
    """
    Convert numeric values to categorical variables.

    Parameters
    ----------
    pandas.Series
        Series object containing numeric values.
    n : int, default: 5
        Number of variables to output.

    Returns
    -------
    pandas.Series
        Series object containing categorical variables.

    Examples
    --------

    >>> import matplotlib.pyplot as plt
    >>> from fuc import common, pymaf
    >>> common.load_dataset('tcga-laml')
    >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
    >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
    >>> s = af.df.days_to_last_followup
    >>> s[:10]
    Tumor_Sample_Barcode
    TCGA-AB-2802     365.0
    TCGA-AB-2803     792.0
    TCGA-AB-2804    2557.0
    TCGA-AB-2805     577.0
    TCGA-AB-2806     945.0
    TCGA-AB-2807     181.0
    TCGA-AB-2808    2861.0
    TCGA-AB-2809      62.0
    TCGA-AB-2810      31.0
    TCGA-AB-2811     243.0
    Name: days_to_last_followup, dtype: float64
    >>> s = common.convert_num2cat(s)
    >>> s.unique()
    array([ 572.2, 1144.4, 2861. , 2288.8, 1716.6,    nan])
    >>> s[:10]
    Tumor_Sample_Barcode
    TCGA-AB-2802     572.2
    TCGA-AB-2803    1144.4
    TCGA-AB-2804    2861.0
    TCGA-AB-2805    1144.4
    TCGA-AB-2806    1144.4
    TCGA-AB-2807     572.2
    TCGA-AB-2808    2861.0
    TCGA-AB-2809     572.2
    TCGA-AB-2810     572.2
    TCGA-AB-2811     572.2
    Name: days_to_last_followup, dtype: float64
    """
    boundaries = list(np.linspace(s.min(), s.max(), n+1, endpoint=True))
    intervals = list(zip(boundaries[:-1], boundaries[1:]))

    def f(x):
        if pd.isna(x):
            return x
        for i, interval in enumerate(intervals):
            a, b = interval
            if a <= x <= b:
                return b

    return s.apply(f).round(decimals=decimals)

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
