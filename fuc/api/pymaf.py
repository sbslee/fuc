"""
The pymaf submodule is designed for working with MAF files. It implements
``pymaf.MafFrame`` which stores MAF data as ``pandas.DataFrame`` to allow
fast computation and easy manipulation. The class also contains many useful
plotting methods such as ``MafFrame.plot_varcls`` and
``MafFrame.plot_waterfall``. The submodule strictly adheres to the
standard `MAF specification
<https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_.
"""

import pandas as pd
import seaborn as sns
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from . import pyvcf
import copy

# Below is the list of calculated variant consequences from Ensemble VEP:
# https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
# (accessed on 2021-05-31)
#
# Note that both frameshift_variant and protein_altering_variant require
# additional information to find their correct Variant_Classification.

VEP_CONSEQUENCES = {
    'transcript_ablation':                'Splice_Site',
    'splice_acceptor_variant':            'Splice_Site',
    'splice_donor_variant':               'Splice_Site',
    'stop_gained':                        'Nonsense_Mutation',
    'frameshift_variant':                 'AMBIGUOUS',
    'stop_lost':                          'Nonstop_Mutation',
    'start_lost':                         'Translation_Start_Site',
    'transcript_amplification':           'Intron',
    'inframe_insertion':                  'In_Frame_Ins',
    'inframe_deletion':                   'In_Frame_Del',
    'missense_variant':                   'Missense_Mutation',
    'protein_altering_variant':           'AMBIGUOUS',
    'splice_region_variant':              'Splice_Region',
    'incomplete_terminal_codon_variant':  'Silent',
    'start_retained_variant':             'Silent',
    'stop_retained_variant':              'Silent',
    'synonymous_variant':                 'Silent',
    'coding_sequence_variant':            'Missense_Mutation',
    'mature_miRNA_variant':               'RNA',
    '5_prime_UTR_variant':                "5'UTR",
    '3_prime_UTR_variant':                "3'UTR",
    'non_coding_transcript_exon_variant': 'RNA',
    'intron_variant':                     'Intron',
    'NMD_transcript_variant':             'Silent',
    'non_coding_transcript_variant':      'RNA',
    'upstream_gene_variant':              "5'Flank",
    'downstream_gene_variant':            "3'Flank",
    'TFBS_ablation':                      'Targeted_Region',
    'TFBS_amplification':                 'Targeted_Region',
    'TF_binding_site_variant':            'IGR',
    'regulatory_region_ablation':         'Targeted_Region',
    'regulatory_region_amplification':    'Targeted_Region',
    'feature_elongation':                 'Targeted_Region',
    'regulatory_region_variant':          'IGR',
    'feature_truncation':                 'Targeted_Region',
    'intergenic_variant':                 'IGR',
}

VARCLS_LIST = [
    "3'Flank",
    "3'UTR",
    "5'Flank",
    "5'UTR",
    'De_novo_Start_InFrame',
    'De_novo_Start_OutOfFrame',
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'IGR',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Intron',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Nonstop_Mutation',
    'RNA',
    'Silent',
    'Splice_Region',
    'Splice_Site',
    'Start_Codon_Ins',
    'Start_Codon_SNP',
    'Stop_Codon_Del',
    'Targeted_Region',
    'Translation_Start_Site',
    'lincRNA',
]

NONSYN_NAMES = [
    'Missense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
    'In_Frame_Del', 'In_Frame_Ins', 'Nonsense_Mutation',
    'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'
]

NONSYN_COLORS = [
    'tab:green', 'tab:blue', 'tab:purple', 'tab:olive', 'tab:red',
    'tab:cyan', 'tab:pink', 'tab:orange', 'tab:brown'
]

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

def plot_legend(name='regular', ax=None, figsize=None, **kwargs):
    """Create one of the pre-defined legends.

    Parameters
    ----------
    name : {'regaulr', 'waterfall'}, default: 'regular'
        Type of legend to be drawn.
    ax : matplotlib.axes.Axes, optional
        Pre-existing axes for the plot. Otherwise, crete a new one.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    kwargs
        Other keyword arguments will be passed down to
        :meth:`Axes.legend`.

    Returns
    -------
    matplotlib.axes.Axes
        The matplotlib axes containing the plot.

    Examples
    --------
    We can add legend for the :meth:`MafFrame.plot_genes` method:

    .. plot::
        :context: close-figs

        >>> import matplotlib.pyplot as plt
        >>> from fuc import common, pymaf
        >>> common.load_dataset('tcga-laml')
        >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
        >>> mf = pymaf.MafFrame.from_file(f)
        >>> fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 6),
        ...     gridspec_kw={'width_ratios': [10, 1]})
        >>> mf.plot_genes(ax=ax1)
        >>> pymaf.plot_legend(name='regular', ax=ax2, loc='center left')
        >>> plt.tight_layout()

    We can also add legend for :meth:`MafFrame.plot_waterfall` method:

    .. plot::
        :context: close-figs

        >>> fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(10, 6),
        ...     gridspec_kw={'height_ratios': [10, 1]})
        >>> mf.plot_waterfall(ax=ax1, linewidths=0.5)
        >>> pymaf.plot_legend(name='waterfall', ax=ax2, ncol=4,
        ...     loc='upper center')
        >>> plt.tight_layout()
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    h = []
    labels = copy.deepcopy(NONSYN_NAMES)
    colors = copy.deepcopy(NONSYN_COLORS)
    if name == 'regular':
        pass
    elif name == 'waterfall':
        labels += ['Multi_Hit']
        colors += ['k']
    else:
        raise ValueError(f'Found incorrect name: {name}')
    for i, label in enumerate(labels):
        h.append(mpatches.Patch(color=colors[i], label=label))
    ax.legend(handles=h, **kwargs)
    ax.axis('off')
    return ax

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

    @classmethod
    def from_vcf(cls, vcf):
        """Construct MafFrame from a VCF file or VcfFrame.

        Parameters
        ----------
        vcf : str or VcfFrame
            VCF file path or VcfFrame.
        """
        # Parse the input VCF.
        if isinstance(vcf, str):
            vf = pyvcf.VcfFrame.from_file(vcf)
        else:
            vf = vcf

        # Get the NCBI_Build data.
        for line in vf.meta:
            if line.startswith('##VEP'):
                ncbi_build = re.search(r'assembly="(.*?)"', line).group(1)

        # Define the conversion algorithm.
        def one_row(r):
            fields = r.INFO.replace('CSQ=', '').split(',')[0].split('|')

            # Get the sequence data.
            inframe = abs(len(r.REF) - len(r.ALT)) / 3 == 0
            if len(r.REF) == len(r.ALT) == 1:
                variant_type = 'SNP'
                start_position = r.POS
                end_position = r.POS
                reference_allele = r.REF
                tumor_seq_allele1 = r.ALT
                tumor_seq_allele2 = r.ALT
            elif len(r.REF) > len(r.ALT):
                variant_type = 'DEL'
                start_position = r.POS + 1
                end_position = r.POS + len(r.REF) - len(r.ALT)
                reference_allele = r.REF[1:]
                tumor_seq_allele1 = '-'
                tumor_seq_allele2 = '-'
            else:
                variant_type = 'INS'
                start_position = r.POS
                end_position = r.POS + 1
                reference_allele = '-'
                tumor_seq_allele1 = r.ALT[1:]
                tumor_seq_allele2 = r.ALT[1:]

            # Get the Strand data.
            strand = '+' if fields[19] == '1' else '-'

            # Get the Variant_Classification data.
            consequence = fields[1].split('&')[0]
            if consequence == 'frameshift_variant':
                if variant_type == 'DEL':
                    variant_classification = 'Frame_Shift_Del'
                else:
                    variant_classification = 'Frame_Shift_Ins'
            elif consequence == 'protein_altering_variant':
                if inframe:
                    if variant_type == 'DEL':
                        variant_classification = 'In_Frame_Del'
                    else:
                        variant_classification = 'In_Frame_Ins'
                else:
                    if variant_type == 'DEL':
                        variant_classification = 'Frame_Shift_Del'
                    else:
                        variant_classification = 'Frame_Shift_Ins'
            elif consequence in VEP_CONSEQUENCES:
                variant_classification = VEP_CONSEQUENCES[consequence]
            else:
                m = f'Found unknown Ensemble VEP consequence: {consequence}'
                raise ValueError(m)

            # Get the Tumor_Sample_Barcode data.
            s = r[9:].apply(pyvcf.gt_hasvar)
            tumor_sample_barcode = ','.join(s[s].index.to_list())

            # Get the Protein_Change data.
            if fields[14]:
                pos = fields[14]
                aa = fields[15].split('/')
                protein_change = f'p.{aa[0]}{pos}{aa[1]}'
            else:
                protein_change = '.'

            d = dict(
                Hugo_Symbol = fields[3],
                Entrez_Gene_Id = fields[4],
                Center = '.',
                NCBI_Build = ncbi_build,
                Chromosome = r.CHROM,
                Start_Position = start_position,
                End_Position = end_position,
                Strand = strand,
                Variant_Classification = variant_classification,
                Variant_Type = variant_type,
                Reference_Allele = reference_allele,
                Tumor_Seq_Allele1 = tumor_seq_allele1,
                Tumor_Seq_Allele2 = tumor_seq_allele2,
                Tumor_Sample_Barcode = tumor_sample_barcode,
                Protein_Change = protein_change,
            )

            return pd.Series(d)

        # Apply the conversion algorithm.
        df = vf.df.apply(one_row, axis=1)

        # Expand the Tumor_Sample_Barcode column to multiple rows.
        s = df['Tumor_Sample_Barcode'].str.split(',').apply(
            pd.Series, 1).stack()
        s.index = s.index.droplevel(-1)
        s.name = 'Tumor_Sample_Barcode'
        del df['Tumor_Sample_Barcode']
        df = df.join(s)

        return cls(df)

    def compute_genes(self, count=10, mode='variants'):
        """Compute a matrix of counts for genes and variant classifications.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing TMB data.
        """
        if mode == 'variants':
            df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]
            df = df.groupby('Hugo_Symbol')[
                'Variant_Classification'].value_counts().to_frame()
            df.columns = ['Count']
            df = df.reset_index()
            df = df.pivot(index='Hugo_Symbol', columns='Variant_Classification',
                values='Count')
            df = df.fillna(0)
            for varcls in NONSYN_NAMES:
                if varcls not in df.columns:
                    df[varcls] = 0
            i = df.sum(axis=1).sort_values(ascending=False).index
            df = df.reindex(index=i)
            df = df[NONSYN_NAMES]
            df = df[:count]
            df = df.rename_axis(None, axis=1)
        elif mode == 'samples':
            df = self.compute_waterfall(count)
            df = df.apply(lambda r: r.value_counts(), axis=1)
            for varcls in NONSYN_NAMES + ['Multi_Hit']:
                if varcls not in df.columns:
                    df[varcls] = np.nan
            df = df[NONSYN_NAMES + ['Multi_Hit']]
            df = df.fillna(0)
        else:
            raise ValueError(f'Found incorrect mode: {mode}')
        return df

    def compute_tmb(self):
        """Compute a matrix of counts for samples and variant classifications.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing TMB data.
        """
        df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]
        df = df.groupby('Tumor_Sample_Barcode')[
            'Variant_Classification'].value_counts().to_frame()
        df.columns = ['Count']
        df = df.reset_index()
        df = df.pivot(index='Tumor_Sample_Barcode',
            columns='Variant_Classification', values='Count')
        df = df.fillna(0)
        for varcls in NONSYN_NAMES:
            if varcls not in df.columns:
                df[varcls] = 0
        i = df.sum(axis=1).sort_values(ascending=False).index
        df = df.reindex(index=i)
        df = df[NONSYN_NAMES]
        df = df.rename_axis(None, axis=1)
        return df

    def compute_waterfall(self, count=10):
        """Compute a matrix of variant classifications for genes and samples.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.

        Returns
        -------
        pandas.DataFrame
            Dataframe containing waterfall data.
        """
        df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]
        f = lambda x: ''.join(x) if len(x) == 1 else 'Multi_Hit'
        df = df.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode'])[
            'Variant_Classification'].apply(f).to_frame()
        df = df.reset_index()
        df = df.pivot(index='Hugo_Symbol', columns='Tumor_Sample_Barcode',
            values='Variant_Classification')

        # Sort the rows (genes).
        i = df.isnull().sum(axis=1).sort_values(ascending=True).index
        df = df.reindex(index=i)

        # Select the top mutated genes.
        df = df[:count]

        # Remove columns (samples) with all NaN's.
        df = df.dropna(axis=1, how='all')

        # Sort the columns (samples).
        c = df.applymap(lambda x: 0 if pd.isnull(x) else 1).sort_values(
            df.index.to_list(), axis=1, ascending=False).columns
        df = df[c]
        df = df.fillna('None')
        df = df.rename_axis(None, axis=1)

        return df

    def plot_genes(
        self, count=10, mode='variants', ax=None, figsize=None, **kwargs
    ):
        """Create a bar plot for top mutated genes.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.barh`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_genes(mode='variants')
            >>> plt.tight_layout()

        .. plot::
            :context: close-figs

            >>> mf.plot_genes(mode='samples')
            >>> plt.tight_layout()
        """
        if mode == 'variants':
            colors = NONSYN_COLORS
        elif mode == 'samples':
            colors = NONSYN_COLORS + ['k']
        else:
            raise ValueError(f'Found incorrect mode: {mode}')
        df = self.compute_genes(count, mode=mode)
        df = df.iloc[::-1]
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        df.plot.barh(stacked=True, ax=ax, color=colors,
            legend=False, **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

    def plot_oncoplot(self, count=10, figsize=(15, 10), fontsize=10):
        """Create an oncoplot.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.
        figsize : tuple, default: (15, 10)
            Width, height in inches. Format: (float, float).
        fontsize : int, default: 10
            Font size.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(f)
            >>> mf.plot_oncoplot(fontsize=14)
        """
        samples = list(self.compute_waterfall(count).columns)

        plt.rc('font', size=fontsize)
        fig, axes = plt.subplots(3, 2, figsize=figsize,
            gridspec_kw={'height_ratios': [1, 10, 1],
                         'width_ratios': [10, 1]})
        [[ax1, ax2], [ax3, ax4], [ax5, ax6]] = axes

        # Create the TMB plot.
        self.plot_tmb(ax=ax1, samples=samples)
        ax1.set_xlabel('')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.set_ylabel('TMB')
        ax1.set_yticks([0, self.compute_tmb().sum(axis=1).max()])

        # Remove the top right plot.
        ax2.remove()

        # Create the waterfall plot.
        self.plot_waterfall(ax=ax3, linewidths=1)
        ax3.set_xlabel('')

        # Create the genes plot.
        self.plot_genes(ax=ax4, mode='samples', width=0.95)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax4.set_yticks([])
        ax4.set_xlabel('Samples')
        ax4.set_xticks([0, self.compute_genes(
            10, mode='samples').sum(axis=1).max()])
        ax4.set_ylim(-0.5, count-0.5)

        # Create the legend.
        plot_legend(name='waterfall', ax=ax5, ncol=4, loc='upper center')

        # Remove the bottom right plot.
        ax6.remove()

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.01, hspace=0.01)

    def plot_snvcls(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for SNV class.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
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
        sns.barplot(x='Count', y='SNV', data=df, ax=ax,
            palette='pastel', **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

    def plot_tmb(self, ax=None, figsize=None, samples=None, **kwargs):
        """Create a bar plot for tumor mutational burden (TMB) per sample.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        samples : list, optional
            Samples to be drawn (in the exact order).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.barh`.

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
            >>> mf.plot_tmb()
            >>> plt.tight_layout()
        """
        df = self.compute_tmb()
        if samples is not None:
            df = df.loc[samples]
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        df.plot.bar(stacked=True, ax=ax, width=1.0, legend=False,
            color=NONSYN_COLORS, **kwargs)
        ax.set_xlabel('Samples')
        ax.set_ylabel('Count')
        ax.set_xticks([])
        return ax

    def plot_varcls(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for the nonsynonymous variant classes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
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
        kwargs
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
        sns.barplot(x='Variant_Type', y='index', data=df, ax=ax,
            palette='pastel', **kwargs)
        ax.set_xlabel('Count')
        ax.set_ylabel('')
        return ax

    def plot_waterfall(self, count=10, ax=None, figsize=None, **kwargs):
        """Create a waterfall plot.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.pcolormesh()` and
            :meth:`seaborn.heatmap`.

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
            >>> mf.plot_waterfall(linewidths=0.5)
            >>> plt.tight_layout()
        """
        df = self.compute_waterfall(count)

        # Apply the mapping between items and integers.
        l = reversed(NONSYN_NAMES + ['Multi_Hit', 'None'])
        d = {k: v for v, k in enumerate(l)}
        df = df.applymap(lambda x: d[x])

        # Plot the heatmap.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        colors = list(reversed(NONSYN_COLORS + ['k', 'lightgray']))
        sns.heatmap(df, cmap=colors, ax=ax, xticklabels=False,
            cbar=False, **kwargs)
        ax.set_xlabel('Samples')
        ax.set_ylabel('')

        return ax

    def to_file(self, fn):
        """Write MafFrame to a MAF file.

        Parameters
        ----------
        fn : str
            VCF file path.
        """
        self.df.to_csv(fn, index=False, sep='\t')

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
