"""
The pycov submodule is designed for working with depth of coverage data
from BAM files. It implements ``pycov.CovFrame`` which stores read depth
data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
"""

import pysam
import pandas as pd
from io import StringIO
from . import pybam

def read_file(fn, region=None):
    """Create CovFrame from a BAM file.
    """
    args = ['-a', '-Q', '1']
    if region:
        args.append('-r')
        args.append(region)
    args.append(fn)
    s = pysam.depth(*args)
    df = pd.read_csv(StringIO(s), sep='\t')
    names = pybam.tag_sm(fn)
    if len(names) > 1:
        raise ValueError('multiple sample names detected')
    df.columns = ['Chromosome', 'Position', names[0]]
    return CovFrame(df)

class CovFrame:
    """Class for storing read depth data."""
    def __init__(self, df):
        self.df = df

    @property
    def shape(self):
        """Return the size of the CovFrame."""
        return self.df.shape
