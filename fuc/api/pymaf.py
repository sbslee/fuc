"""
The pymaf submodule is designed for working with MAF files. It implements
the ``pymaf.MafFrame`` class which stores MAF data as ``pandas.DataFrame``
to allow fast computation and easy manipulation. The submodule strictly
adheres to the standard `MAF specification
<https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_.
"""

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

VARCLS_DICT = {
    "3'Flank": {'NONSYN': False, 'COLOR': None},
    "3'UTR": {'NONSYN': False, 'COLOR': None},
    "5'Flank": {'NONSYN': False, 'COLOR': None},
    "5'UTR": {'NONSYN': False, 'COLOR': None},
    'De_novo_Start_InFrame': {'NONSYN': False, 'COLOR': None},
    'De_novo_Start_OutOfFrame': {'NONSYN': False, 'COLOR': None},
    'Frame_Shift_Del': {'NONSYN': True, 'COLOR': 'tab:blue'},
    'Frame_Shift_Ins': {'NONSYN': True, 'COLOR': 'tab:purple'},
    'IGR': {'NONSYN': False, 'COLOR': None},
    'In_Frame_Del': {'NONSYN': True, 'COLOR': 'tab:olive'},
    'In_Frame_Ins': {'NONSYN': True, 'COLOR': 'tab:gray'},
    'Intron': {'NONSYN': False, 'COLOR': None},
    'Missense_Mutation': {'NONSYN': True, 'COLOR': 'tab:green'},
    'Nonsense_Mutation': {'NONSYN': True, 'COLOR': 'tab:red'},
    'Nonstop_Mutation': {'NONSYN': True, 'COLOR': 'tab:pink'},
    'RNA': {'NONSYN': False, 'COLOR': None},
    'Silent': {'NONSYN': False, 'COLOR': None},
    'Splice_Region': {'NONSYN': False, 'COLOR': None},
    'Splice_Site': {'NONSYN': True, 'COLOR': 'tab:orange'},
    'Start_Codon_Ins': {'NONSYN': False, 'COLOR': None},
    'Start_Codon_SNP': {'NONSYN': False, 'COLOR': None},
    'Stop_Codon_Del': {'NONSYN': False, 'COLOR': None},
    'Targeted_Region': {'NONSYN': False, 'COLOR': None},
    'Translation_Start_Site': {'NONSYN': True, 'COLOR': 'tab:brown'},
    'lincRNA': {'NONSYN': False, 'COLOR': None},
}

NONSYN_NAMES = sorted([k for k,v in VARCLS_DICT.items() if v['NONSYN']])
NONSYN_COLORS = [VARCLS_DICT[x]['COLOR'] for x in NONSYN_NAMES]

SNVCLS = {
    'A>C': {'REP': 'T>G'},
    'A>G': {'REP': 'T>C'},
    'A>T': {'REP': 'T>A'},
    'C>A': {'REP': 'C>A'},
    'C>G': {'REP': 'C>G'},
    'C>T': {'REP': 'C>T'},
    'G>A': {'REP': 'C>T'},
    'G>C': {'REP': 'C>G'},
    'G>T': {'REP': 'C>A'},
    'T>A': {'REP': 'T>A'},
    'T>C': {'REP': 'T>C'},
    'T>G': {'REP': 'T>G'},
}

class MafFrame:
    """Class for storing MAF data.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing MAF data.

    See Also
    --------
    MafFrame.from_file
        Construct MafFrame from a MAF file.
    """
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @property
    def shape(self):
        """tuple : Dimensionality of MafFrame (variants, samples)."""
        return (self.df.shape[0], len(self.samples))

    @property
    def samples(self):
        """list : List of the sample names."""
        return list(self.df.Tumor_Sample_Barcode.unique())

    @property
    def genes(self):
        """list : List of the genes."""
        return list(self.df.Hugo_Symbol.unique())

    def filter_nonsyn(self, opposite=False, index=False):
        """Select rows with a nonsynonymous variant.

        Parameters
        ----------
        opposite : bool, default: False
            If True, return rows that don't meet the said criteria.
        index : bool, default: False
            If True, return boolean index array instead of MafFrame.

        Returns
        -------
        MafFrame or pandas.Series
            Filtered MafFrame or boolean index array.
        """
        nonsyn_list = [k for k, v in VARCLS_DICT.items() if v['NONSYN']]
        one_row = lambda r: r.Variant_Classification in nonsyn_list
        i = self.df.apply(one_row, axis=1)
        if opposite:
            i = ~i
        if index:
            return i
        return self.__class__(self.df[i])

    @classmethod
    def from_file(cls, fn):
        """Construct MafFrame from a MAF file.

        Parameters
        ----------
        fn : str
            MAF file path (zipped or unzipped).

        Returns
        -------
        MafFrame
            MafFrame.

        See Also
        --------
        MafFrame
            MafFrame object creation using constructor.
        """
        return cls(pd.read_table(fn))

    def plot_genenum(self, count=None, ax=None, figsize=None, **kwargs):
        """Create a bar plot for mutated genes.

        Parameters
        ----------
        count : int, optional
            Number of top genes to display.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs : key, value mappings
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_genenum(count=5)
            >>> plt.tight_layout()
        """
        s = self.df.Hugo_Symbol.value_counts()
        if count is not None:
            s = s[:count]
        df = s.to_frame().reset_index()
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        if kwargs is None:
            kwargs = {}
        sns.barplot(x='Hugo_Symbol', y='index', data=df, ax=ax, **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

    def plot_sampnum(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for variants per sample.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs : key, value mappings
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_sampnum()
            >>> plt.tight_layout()
        """
        s = self.df.Tumor_Sample_Barcode.value_counts()
        df = s.to_frame().reset_index()
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        if kwargs is None:
            kwargs = {}
        df.plot.bar(
            x='index', y='Tumor_Sample_Barcode', ax=ax, width=1.0,
            legend=False, **kwargs
        )
        ax.set_xlabel('Samples')
        ax.set_ylabel('Count')
        ax.set_xticks([])
        return ax

    def plot_snvcls(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for SNV class.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs : key, value mappings
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_snvcls()
            >>> plt.tight_layout()
        """
        def one_row(r):
            ref = r.Reference_Allele
            alt = r.Tumor_Seq_Allele2
            if (ref == '-' or
                alt == '-' or
                len(alt) != 1 or
                len(ref) != 1 or
                ref == alt):
                return np.nan
            return SNVCLS[f'{ref}>{alt}']['REP']
        s = self.df.apply(one_row, axis=1).value_counts()
        i = sorted(set([v['REP'] for k, v in SNVCLS.items()]))
        s = s.reindex(index=i)
        df = s.to_frame().reset_index()
        df.columns = ['SNV', 'Count']
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        if kwargs is None:
            kwargs = {}
        sns.barplot(x='Count', y='SNV', data=df, ax=ax, **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

    def plot_varcls(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for the nonsynonymous variant classes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs : key, value mappings
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_varcls()
            >>> plt.tight_layout()
        """
        d = self.df.Variant_Classification.value_counts().to_dict()
        counts = {}
        for varcls in NONSYN_NAMES:
            if varcls in d:
                counts[varcls] = d[varcls]
            else:
                counts[varcls] = 0
        s = pd.Series(counts).reindex(index=NONSYN_NAMES)
        df = s.to_frame().reset_index()
        df.columns = ['Variant_Classification', 'Count']
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        if kwargs is None:
            kwargs = {}
        sns.barplot(x='Count', y='Variant_Classification', data=df,
                    ax=ax, palette=NONSYN_COLORS, **kwargs)
        ax.set_ylabel('')
        return ax

    def plot_vartype(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for viaration type.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs : key, value mappings
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_vartype()
            >>> plt.tight_layout()
        """
        s = self.df.Variant_Type.value_counts()
        df = s.to_frame().reset_index()
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        if kwargs is None:
            kwargs = {}
        sns.barplot(x='Variant_Type', y='index', data=df, ax=ax, **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

class AnnFrame:
    """Class for storing annotation data.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing annotation data.

    See Also
    --------
    AnnFrame.from_file
        Construct AnnFrame from an annotation file.
    """
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """Construct AnnFrame from an annotation file.

        Parameters
        ----------
        fn : str
            Annotation file path (zipped or unzipped).

        Returns
        -------
        AnnFrame
            AnnFrame.

        See Also
        --------
        AnnFrame
            AnnFrame object creation using constructor.
        """
        return cls(pd.read_table(fn))
