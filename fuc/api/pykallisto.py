"""
The pykallisto submodule is designed for working with RNAseq quantification
data from Kallisto. It implements ``pykallisto.KallistoFrame`` which stores
Kallisto's output data as ``pandas.DataFrame`` to allow fast computation and
easy manipulation. The ``pykallisto.KallistoFrame`` class also contains many
useful plotting methods such as ``KallistoFrame.plot_differential_abundance``.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

class KallistoFrame:
    """
    Class for working with RNAseq quantification data from Kallisto.

    Parameters
    ----------
    metadata : pandas.DataFrame
        List of metadata lines.
    tx2gene : pandas.DataFrame
        DataFrame containing transcript to gene mapping data.
    aggregation_column : str
        Column name in ``tx2gene`` to aggregate transcripts to the gene level.
    filter_func : func, optional
        Filtering function to be applied to each row (i.e. transcript). By
        default, the :meth:`pykallisto.basic_filter` method will be used.
    filter_target_id : list, optional
        Transcripts to filter using methods that can't be implemented using
        ``filter_func``. If provided, this will override ``filter_func``.
    filter_off : bool, default: False
        If True, do not apply any filtering. Useful for generating a simple
        count or tpm matrix.
    """

    def _import_data(
        self, metadata, filter_func=None, filter_target_id=None,
        filter_off=False
    ):
        dfs = {}
        for i, r in metadata.iterrows():
            df = pd.read_table(r['path'] + '/abundance.tsv', index_col=0)
            dfs[i] = df
        df_tx_count = pd.concat([v['est_counts'] for k, v in dfs.items()], axis=1)
        df_tx_count.columns = metadata.index
        df_tx_tpm = pd.concat([v['tpm'] for k, v in dfs.items()], axis=1)
        df_tx_tpm.columns = metadata.index

        if filter_off:
            filtered_ids = None
        else:
            if filter_target_id is None:
                if filter_func is None:
                    filter_func = basic_filter
                filtered_ids = df_tx_count.apply(filter_func, axis=1)
            else:
                filtered_ids = pd.Series(df_tx_count.index.isin(filter_target_id),
                    index=df_tx_count.index)

        return df_tx_count, df_tx_tpm, filtered_ids

    def __init__(
        self, metadata, tx2gene, aggregation_column, filter_func=None,
        filter_target_id=None, filter_off=False
    ):
        self.metadata = metadata
        self.tx2gene = tx2gene
        self.aggregation_column = aggregation_column
        results = self._import_data(metadata, filter_func,
            filter_target_id, filter_off)
        self.df_tx_count = results[0]
        self.df_tx_tpm = results[1]
        self.filtered_ids = results[2]
        self.df_gene_count = None
        self.df_gene_tpm = None

    def aggregate(self, filter=True):
        """
        Aggregate transcript-level data to obtain gene-level data.

        Running this method will set the attributes
        ``KallistoFrame.df_gene_count`` and ``KallistoFrame.df_gene_tpm``.

        Parameters
        ----------
        filter : bool, default: True
            If true, use filtered transcripts only. Otherwise, use all.
        """
        for unit in ['count', 'tpm']:
            tx_df = getattr(self, f'df_tx_{unit}')
            gene_s = self.tx2gene[self.aggregation_column]
            df = pd.concat([tx_df, gene_s], axis=1)
            if filter:
                if self.filtered_ids is not None:
                    df = df[self.filtered_ids]
            df = df.groupby(['gene_symbol']).sum()
            setattr(self, f'df_gene_{unit}', df)

    def plot_differential_abundance(
        self, gene, group, aggregate=True, filter=True, name='target_id',
        unit='tpm', ax=None, figsize=None
    ):
        """
        Plot differential abundance results for single gene.

        Parameters
        ----------
        gene : str
            Gene to compare.
        group : str
            Column in ``KallistoFrame.metadata`` specifying group information.
        aggregate : bool, default: True
            If true, display gene-level data (the
            :meth:`KallistoFrame.aggregate` method must be run beforehand).
            Otherwise, display transcript-level data.
        filter : bool, default: True
            If true, use filtered transcripts only. Otherwise, use all.
            Ignored when ``aggregate=True``.
        name : str, default: 'target_id'
            Column in ``KallistoFrame.tx2gene`` specifying transcript name to
            be displayed in the legend. Ignored when ``aggregate=True``.
        unit : {'tpm', 'count'}, default: 'tpm'
            Abundance unit to display.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if aggregate:
            s1 = getattr(self, f'df_gene_{unit}').loc[gene]
            s1.name = unit
            s2 = self.metadata[group]
            df = pd.concat([s1, s2], axis=1)
            sns.boxplot(x=group, y=unit, data=df, ax=ax)
        else:
            df = getattr(self, f'df_tx_{unit}')
            if filter:
                if self.filtered_ids is not None:
                    df = df[self.filtered_ids]
            s = self.tx2gene[self.tx2gene[self.aggregation_column] == gene]
            df = df[df.index.isin(s.index.to_list())]
            if name != 'target_id':
                df = df.merge(self.tx2gene[name], left_on='target_id',
                    right_index=True)
                df = df.set_index(name)
            df = df.T.melt(ignore_index=False, value_name=unit)
            df = df.merge(self.metadata[group], left_index=True,
                right_index=True)
            sns.boxplot(x=group, y=unit, data=df, hue=name, ax=ax)

        return ax

    def compute_fold_change(self, group, genes, unit='tpm', flip=False):
        """
        Compute fold change of gene expression between two groups.

        Parameters
        ----------
        group : str
            Column in ``KallistoFrame.metadata`` specifying group information.
        gene : list
            Genes to compare.
        unit : {'tpm', 'count'}, default: 'tpm'
            Abundance unit to display.
        flip : bool, default: False
            If true, flip the denominator and numerator.
        """
        df = getattr(self, f'df_gene_{unit}')
        df = df[df.index.isin(genes)].T
        df = df.merge(self.metadata[group], left_index=True, right_index=True)
        df = df.groupby(group).mean().T
        a, b = df.columns
        if flip:
            a, b = b, a
        s = np.log2(df[b] / df[a])
        print(f'# Fold change = log2( {b} / {a} )')
        return s[genes]

def basic_filter(row, min_reads=5, min_prop=0.47):
    """
    A basic filter to be used.

    By default, the method will filter out rows (i.e. transcripts) that do
    not have at least 5 estimated counts in at least 47% of the samples. Note
    that this is equivalent to the :meth:`sleuth.basic_filter` method.

    Parameters
    ----------
    row : pandas.Series
        This is a vector of numerics that will be passed in.
    min_reads : int, default: 5
        The minimum number of estimated counts.
    min_prop : float, default: 0.47
        The minimum proportion of samples.

    Returns
    -------
    pd.Series
        A pandas series of boolean.
    """
    return (row >= min_reads).mean() >= min_prop
