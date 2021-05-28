"""
The pymaf submodule is designed for working with MAF files. It implements
the ``pymaf.MafFrame`` class which stores MAF data as ``pandas.DataFrame``
to allow fast computation and easy manipulation.
"""

import pandas as pd

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
