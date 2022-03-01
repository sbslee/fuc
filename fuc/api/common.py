"""
The common submodule is used by other fuc submodules such as pyvcf and
pybed. It also provides many day-to-day actions used in the field of
bioinformatics.
"""

import pathlib
import re
import os
import sys
import warnings
import inspect
import copy
from argparse import RawTextHelpFormatter, SUPPRESS
from pathlib import Path
from difflib import SequenceMatcher
from urllib.request import urlretrieve

from . import pyvcf

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
            Text file (compressed or uncompressed).
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
            AnnFrame column containing sample group information. If the
            column has NaN values, they will be converted to 'N/A' string.
        group_order : list, optional
            List of sample group names (in that order too). You can use this
            to subset samples belonging to specified groups only. You must
            include all relevant groups when also using ``samples``.
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
        s = s.fillna('N/A')

        # Subset the samples, if necessary.
        if samples is not None:
            s = s.reindex(samples)

        # Establish mapping from groups to numbers.
        if group_order is None:
            group_order = sorted([x for x in s.unique() if x == x])
        else:
            # Make sure all specified groups are valid.
            for group in group_order:
                groups = ', '.join([f"'{x}'" for x in sorted(s.unique())])
                if group not in s.unique():
                    raise ValueError(f"The group '{group}' does not exist. "
                        f"The following groups are available: {groups}.")

            if len(group_order) < len(s.unique()):
                if samples is None:
                    s = s[s.isin(group_order)]
                else:
                    missing = ', '.join([f"'{x}'" for x in s.unique()
                        if x not in group_order])
                    raise ValueError("The 'group_order' argumnet must "
                        "include all groups when used with the 'samples' "
                        "argument. Following groups are currently missing: "
                        f"{missing}.")

        d = {k: v for v, k in enumerate(group_order)}
        df = s.to_frame().applymap(lambda x: d[x])

        # Determine the colors to use.
        if isinstance(colors, str):
            if sequential:
                c = plt.get_cmap(colors).colors
                l = list(np.linspace(0, len(c)-1, len(group_order)))
                colors = [c[round(i)] for i in l]
            else:
                colors = list(plt.get_cmap(colors).colors[:len(group_order)])

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

    def subset(self, samples, exclude=False):
        """
        Subset AnnFrame for specified samples.

        Parameters
        ----------
        samples : str or list
            Sample name or list of names (the order matters).
        exclude : bool, default: False
            If True, exclude specified samples.

        Returns
        -------
        AnnFrame
            Subsetted AnnFrame.

        Examples
        --------

        >>> from fuc import common
        >>> data = {
        ...     'SampleID': ['A', 'B', 'C', 'D'],
        ...     'PatientID': ['P1', 'P1', 'P2', 'P2'],
        ...     'Tissue': ['Normal', 'Tumor', 'Normal', 'Tumor'],
        ...     'Age': [30, 30, 57, 57]
        ... }
        >>> af = common.AnnFrame.from_dict(data, sample_col='SampleID') # or sample_col=0
        >>> af.df
                 PatientID  Tissue  Age
        SampleID
        A               P1  Normal   30
        B               P1   Tumor   30
        C               P2  Normal   57
        D               P2   Tumor   57

        We can subset the AnnFrame for the normal samples A and C:

        >>> af.subset(['A', 'C']).df
                 PatientID  Tissue  Age
        SampleID
        A               P1  Normal   30
        C               P2  Normal   57

        Alternatively, we can exclude those samples:

        >>> af.subset(['A', 'C'], exclude=True).df
                 PatientID Tissue  Age
        SampleID
        B               P1  Tumor   30
        D               P2  Tumor   57
        """
        if isinstance(samples, str):
            samples = [samples]
        if exclude:
            samples = [x for x in self.samples if x not in samples]
        return self.__class__(self.df.loc[samples])

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

    This example is directly taken from the Wiki page `Sensitivity and specificity <https://en.wikipedia.org/wiki/Sensitivity_and_specificity>`__.

    >>> from fuc import common
    >>> results = common.sumstat(180, 10, 20, 1820)
    >>> for k, v in results.items():
    ...     print(k, f'{v:.3f}')
    ...
    tpr 0.667
    tnr 0.910
    ppv 0.100
    npv 0.995
    fnr 0.333
    fpr 0.090
    fdr 0.900
    for 0.005
    acc 0.906
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
            'getrm-cyp2d6-vdr.vcf',
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
    Parse specified genomic region.

    The method will return parsed region as a tuple with a shape of
    ``(chrom, start, end)`` which has data types of ``(str, int, int)``.

    Note that only ``chrom`` is required when specifing a region. If
    ``start`` and ``end`` are omitted, the method will return ``NaN`` in
    their respective positions in the output tuple.

    Parameters
    ----------
    region : str
        Region ('chrom:start-end').

    Returns
    -------
    tuple
        Parsed region.

    Examples
    --------

    >>> from fuc import common
    >>> common.parse_region('chr1:100-150')
    ('chr1', 100, 150)
    >>> common.parse_region('chr1')
    ('chr1', nan, nan)
    >>> common.parse_region('chr1:100')
    ('chr1', 100, nan)
    >>> common.parse_region('chr1:100-')
    ('chr1', 100, nan)
    >>> common.parse_region('chr1:-100')
    ('chr1', nan, 100)
    """
    chrom = region.split(':')[0]

    try:
        start = region.split(':')[1].split('-')[0]
        if not start:
            start = np.nan
        else:
            start = int(start)
    except IndexError:
        start = np.nan

    try:
        end = region.split(':')[1].split('-')[1]
        if not end:
            end = np.nan
        else:
            end = int(end)
    except IndexError:
        end = np.nan

    return (chrom, start, end)

def parse_variant(variant):
    """
    Parse specified genomic variant.

    Generally speaking, the input string should consist of chromosome,
    position, reference allele, and alternative allele separated by any one
    or combination of the following delimiters: ``-``, ``:``, ``>`` (e.g.
    '22-42127941-G-A'). The method will return parsed variant as a tuple with
    a shape of ``(chrom, pos, ref, alt)`` which has data types of ``(str,
    int, str, str)``.

    Note that it's possible to omit reference allele and alternative allele
    from the input string to indicate position-only data (e.g.
    '22-42127941'). In this case, the method will return empty string for
    the alleles -- i.e. ``(str, int, '', '')`` if both are omitted and
    ``(str, int, str, '')`` if only alternative allele is omitted.

    Parameters
    ----------
    variant : str
        Genomic variant.

    Returns
    -------
    tuple
        Parsed variant.

    Examples
    --------

    >>> from fuc import common
    >>> common.parse_variant('22-42127941-G-A')
    ('22', 42127941, 'G', 'A')
    >>> common.parse_variant('22:42127941-G>A')
    ('22', 42127941, 'G', 'A')
    >>> common.parse_variant('22-42127941')
    ('22', 42127941, '', '')
    >>> common.parse_variant('22-42127941-G')
    ('22', 42127941, 'G', '')
    """
    fields = re.split('-|:|>', variant)
    chrom = fields[0]
    pos = int(fields[1])

    try:
        ref = fields[2]
    except IndexError:
        ref = ''

    try:
        alt = fields[3]
    except IndexError:
        alt = ''

    return (chrom, pos, ref, alt)

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

def plot_exons(
    starts, ends, name=None, offset=1, fontsize=None, color='black', y=0,
    height=1, ax=None, figsize=None
):
    """
    Create a gene model where exons are drawn as boxes.

    Parameters
    ----------
    starts : list
        List of exon start positions.
    ends : list
        List of exon end positions.
    name : str, optional
        Gene name. Use ``name='$text$'`` to italicize the text.
    offset : float, default: 1
        How far gene name should be plotted from the gene model.
    color : str, default: 'black'
        Box color.
    y : float, default: 0
        Y position of the backbone.
    height : float, default: 1
        Height of the gene model.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, crete a new one.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        The matplotlib axes containing the plot.

    Examples
    --------

    .. plot::

        >>> import matplotlib.pyplot as plt
        >>> from fuc import common
        >>> cyp2d6_starts = [42522500, 42522852, 42523448, 42523843, 42524175, 42524785, 42525034, 42525739, 42526613]
        >>> cyp2d6_ends = [42522754, 42522994, 42523636, 42523985, 42524352, 42524946, 42525187, 42525911, 42526883]
        >>> ax = common.plot_exons(cyp2d6_starts, cyp2d6_ends, name='CYP2D6', fontsize=20)
        >>> ax.set_ylim([-2, 2])
        >>> plt.tight_layout()
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.hlines(
        y=y, xmin=starts[0], xmax=ends[-1], color=color
    )

    for i in range(len(starts)):
        ax.add_patch(mpatches.Rectangle(
            xy=(starts[i], y - height/2),
            width=ends[i] - starts[i],
            height=height,
            color="black"
        ))

    if name is not None:
        ax.text(
            x=(starts[0]+ends[-1]) / 2,
            y=y-offset,
            s=name,
            horizontalalignment='center',
            fontsize=fontsize,
        )

    return ax

def conda_env():
    """str : Name of the current conda environment."""
    return sys.executable.split('/')[-3]

def color_print(s, color='green', bold=False):
    """Print colored text."""
    colors = {
        'red' :  '31',
        'green': '32',
    }

    c = colors[color]

    if bold:
        b = 1
    else:
        b = 0

    print(f'\x1b[{b};{c};{c}m' + s + '\x1b[0m')

def rename(original, names, indicies=None):
    """
    Rename sample names flexibly.

    Parameters
    ----------
    original : list
        List of original names.
    names : dict or list
        Dict of old names to new names or list of new names.
    indicies : list or tuple, optional
        List of 0-based sample indicies. Alternatively, a tuple
        (int, int) can be used to specify an index range.

    Returns
    -------
    list
        List of updated names.

    Examples
    --------

    >>> from fuc import common
    >>> original = ['A', 'B', 'C', 'D']
    >>> common.rename(original, ['1', '2', '3', '4'])
    ['1', '2', '3', '4']
    >>> common.rename(original, {'B': '2', 'C': '3'})
    ['A', '2', '3', 'D']
    >>> common.rename(original, ['2', '4'], indicies=[1, 3])
    ['A', '2', 'C', '4']
    >>> common.rename(original, ['2', '3'], indicies=(1, 3))
    ['A', '2', '3', 'D']
    """
    samples = copy.deepcopy(original)

    if not isinstance(names, list) and not isinstance(names, dict):
        raise TypeError("Argument 'names' must be dict or list.")

    if len(names) > len(samples):
        raise ValueError("There are too many names.")

    if isinstance(names, list) and indicies is not None:
        if isinstance(indicies, tuple):
            if len(indicies) != 2:
                raise ValueError("Index range must be two integers.")
            l = len(range(indicies[0], indicies[1]))
        elif isinstance(indicies, list):
            l = len(indicies)
        else:
            raise TypeError("Argument 'indicies' must be list or tuple.")

        if len(names) != l:
            raise ValueError("Names and indicies have different lengths.")

    if isinstance(names, list):
        if len(names) == len(samples):
            names = dict(zip(samples, names))
        else:
            if indicies is None:
                message = ("There are too few names. If this was "
                    "intended, use the 'indicies' argument.")
                raise ValueError(message)
            elif isinstance(indicies, tuple):
                names = dict(zip(samples[indicies[0]:indicies[1]], names))
            else:
                names = dict(zip([samples[i] for i in indicies], names))

    for old, new in names.items():
        i = samples.index(old)
        samples[i] = new

    if len(samples) > len(set(samples)):
        raise ValueError('There are more than one duplicate names.')

    return samples

def sort_variants(variants):
    """
    Return sorted list of variants.

    Parameters
    ----------
    variants : list
        List of variants.

    Returns
    -------
    list
        Sorted list.

    Examples
    --------

    >>> from fuc import common
    >>> variants = ['5-200-G-T', '5:100:T:C', '1:100:A>C', '10-100-G-C']
    >>> sorted(variants) # Lexicographic sorting (not what we want)
    ['10-100-G-C', '1:100:A>C', '5-200-G-T', '5:100:T:C']
    >>> common.sort_variants(variants)
    ['1:100:A>C', '5:100:T:C', '5-200-G-T', '10-100-G-C']
    """
    def func(x):
        chrom, pos, ref, alt = parse_variant(x)
        if chrom in pyvcf.CONTIGS:
            chrom = pyvcf.CONTIGS.index(chrom)
        return (chrom, pos, ref, alt)
    return sorted(variants, key=func)

def sort_regions(regions):
    """
    Return sorted list of regions.

    Parameters
    ----------
    regions : list
        List of regions.

    Returns
    -------
    list
        Sorted list.

    Examples
    --------

    >>> from fuc import common
    >>> regions = ['chr22:1000-1500', 'chr16:100-200', 'chr22:200-300', 'chr16_KI270854v1_alt', 'chr3_GL000221v1_random', 'HLA-A*02:10']
    >>> sorted(regions) # Lexicographic sorting (not what we want)
    ['HLA-A*02:10', 'chr16:100-200', 'chr16_KI270854v1_alt', 'chr22:1000-1500', 'chr22:200-300', 'chr3_GL000221v1_random']
    >>> common.sort_regions(regions)
    ['chr16:100-200', 'chr22:200-300', 'chr22:1000-1500', 'chr16_KI270854v1_alt', 'chr3_GL000221v1_random', 'HLA-A*02:10']
    """
    def func(x):
        chrom, start, end = parse_region(x)
        if chrom in pyvcf.CONTIGS:
            chrom = pyvcf.CONTIGS.index(chrom)
            alt = ''
        else:
            chrom = len(pyvcf.CONTIGS)
            alt = chrom
        return (chrom, alt, start, end)
    return sorted(regions, key=func)

def update_chr_prefix(regions, mode='remove'):
    """
    Add or remove the (annoying) 'chr' string from specified regions.

    The method will automatically detect regions that don't need to be
    updated and will return them unchanged.

    Parameters
    ----------
    regions : str or list
        One or more regions to be updated.
    mode : {'add', 'remove'}, default: 'remove'
        Whether to add or remove the 'chr' string.

    Returns
    -------
    VcfFrame
        str or list.

    Example
    -------

    >>> from fuc import common
    >>> common.update_chr_prefix(['chr1:100-200', '2:300-400'], mode='remove')
    ['1:100-200', '2:300-400']
    >>> common.update_chr_prefix(['chr1:100-200', '2:300-400'], mode='add')
    ['chr1:100-200', 'chr2:300-400']
    >>> common.update_chr_prefix('chr1:100-200', mode='remove')
    '1:100-200'
    >>> common.update_chr_prefix('chr1:100-200', mode='add')
    'chr1:100-200'
    >>> common.update_chr_prefix('2:300-400', mode='add')
    'chr2:300-400'
    >>> common.update_chr_prefix('2:300-400', mode='remove')
    '2:300-400'
    """
    def remove(x):
        return x.replace('chr', '')

    def add(x):
        if 'chr' not in x:
            x = 'chr' + x
        return x

    modes = {'remove': remove, 'add': add}

    if isinstance(regions, str):
        return modes[mode](regions)

    return [modes[mode](x) for x in regions]

def parse_list_or_file(obj, extensions=['txt', 'tsv', 'csv', 'list']):
    """
    Parse the input variable and then return a list of items.

    This method is useful when parsing a command line argument that accepts
    either a list of items or a text file containing one item per line.

    Parameters
    ----------
    obj : str or list
        Object to be tested. Must be non-empty.
    extensions : list, default: ['txt', 'tsv', 'csv', 'list']
        Recognized file extensions.

    Returns
    -------
    list
        List of items.

    Examples
    --------

    >>> from fuc import common
    >>> common.parse_list_or_file(['A', 'B', 'C'])
    ['A', 'B', 'C']
    >>> common.parse_list_or_file('A')
    ['A']
    >>> common.parse_list_or_file('example.txt')
    ['A', 'B', 'C']
    >>> common.parse_list_or_file(['example.txt'])
    ['A', 'B', 'C']
    """
    if not isinstance(obj, str) and not isinstance(obj, list):
        raise TypeError(
            f'Input must be str or list, not {type(obj).__name__}')

    if not obj:
        raise ValueError('Input is empty')

    if isinstance(obj, str):
        obj = [obj]

    if len(obj) > 1:
        return obj

    for extension in extensions:
        if obj[0].endswith(f'.{extension}'):
            return convert_file2list(obj[0])

    return obj
