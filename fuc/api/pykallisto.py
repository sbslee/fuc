"""
The pykallisto submodule is designed for working with RNAseq quantification
data from Kallisto. It implements ``pykallisto.KallistoFrame`` which stores
Kallisto's output data as ``pandas.DataFrame`` to allow fast computation and
easy manipulation. The ``pykallisto.KallistoFrame`` class also contains many
useful plotting methods such as ``KallistoFrame.plot_txdiff``.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

class KallistoFrame:
    """
    Class for working with RNAseq quantification data from Kallisto.

    Parameters
    ----------
    meta : pandas.DataFrame
        List of meta lines.
    tx2gene : pandas.DataFrame, optional
        DataFrame containing VCF data.
    filter_func : func, optional
        Filtering function to be applied to each row.
    """

    def _import_data(self, meta, filter_func):
        dfs = {}
        for i, r in meta.iterrows():
            df = pd.read_table(r['path'] + '/abundance.tsv', index_col=0)
            dfs[i] = df
        count_data = pd.concat([v['est_counts'] for k, v in dfs.items()], axis=1)
        count_data.columns = meta.index
        tpm_data = pd.concat([v['tpm'] for k, v in dfs.items()], axis=1)
        tpm_data.columns = meta.index
        if filter_func is None:
            filter_func = basic_filter
        filtered_ids = count_data.apply(filter_func, axis=1)
        return count_data, tpm_data, filtered_ids

    def __init__(self, meta, tx2gene=None, filter_func=None):
        self.meta = meta
        self.tx2gene = tx2gene
        results = self._import_data(meta, filter_func)
        self.count_data, self.tpm_data, self.filtered_ids = results

    def plot_txdiff(
        self, gene, name_col, color_col, tpm=False, tx_col=None, ax=None,
        figsize=None
    ):
        """
        Create a grouped box plot showing differential expressions of a
        given gene at the transcript level.

        Parameters
        ----------
        gene : list
            List of meta lines.
        name_col : pandas.DataFrame
            DataFrame containing VCF data.
        color_col
        tpm
        tx_col
        ax
        figsize
        """
        if tpm:
            unit = 'tpm'
        else:
            unit = 'est_counts'

        target_ids = self.tx2gene[self.tx2gene[name_col] == gene].index.to_list()

        df = self.count_data[self.filtered_ids]
        df = df[df.index.isin(target_ids)].T
        df = df.melt(ignore_index=False, value_name='expression')
        df = df.merge(self.meta[color_col], left_index=True, right_index=True)

        if tx_col is None:
            hue = 'target_id'
        else:
            hue = tx_col
            df = df.merge(self.tx2gene[tx_col], left_on='target_id', right_index=True)

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.boxplot(x='lauren', y='expression', data=df, hue=hue)

        return ax

def basic_filter(row, min_reads=5, min_prop=0.47):
    """
    A basic filter to be used.

    By default, the method will filter out rows (transcripts or genes) that
    do not have at least 5 estimated counts in at least 47% of the samples.

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
