"""
The pymaf submodule is designed for working with MAF files (both zipped
and unzipped). It implements the ``pymaf.MafFrame`` class which stores MAF
data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
"""

import pandas as pd

class MafFrame:
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """Construct MafFrame from a MAF file.

        Parameters
        ----------
        fn : str
            Path to the input MAF file.

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
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """Construct AnnFrame from a MAF file.

        Parameters
        ----------
        fn : str
            Path to the input annotation file.

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
