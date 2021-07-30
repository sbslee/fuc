"""
The common submodule is used by other fuc submodules such as pyvcf and
pybed. It also provides many day-to-day actions used in the field of
bioinformatics.
"""

import pathlib
import re
import os
from difflib import SequenceMatcher
from urllib.request import urlretrieve
from pathlib import Path
import pysam
import warnings
import inspect
from argparse import RawTextHelpFormatter, SUPPRESS
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd

FUC_PATH = pathlib.Path(__file__).parent.parent.parent.absolute()

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
