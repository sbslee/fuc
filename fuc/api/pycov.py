"""
The pycov submodule is designed for working with depth of coverage data
from BAM files. It implements ``pycov.CovFrame`` which stores read depth
data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
"""

import pysam
import numpy as np
import pandas as pd
from io import StringIO
from . import pybam

def read_file(fn, zero=False, region=None, map_qual=None):
    """Create CovFrame from BAM files.

    Parameters
    ----------
    fn : str
        Path to the BAM file.
    zero : bool, default: False
        If True, output all positions (including those with zero depth).
    region : str, optional
        Only report depth in specified region (format: CHR:FROM-TO).
    map_qual: int, optional
        Only count reads with mapping quality greater than or equal to INT

    Returns
    -------
    CovFrame
        CovFrame.
    """
    args = []
    if zero:
        args += ['-a']
    if region is not None:
        args += ['-r', region]
    if map_qual is not None:
        args += ['-Q', map_qual]
    args.append(fn)
    s = pysam.depth(*args)
    names = ['Chromosome', 'Position']
    samples = pybam.tag_sm(fn)
    if len(samples) > 1:
        raise ValueError('multiple sample names detected')
    dtype = {'Chromosome': str,'Position': np.int32, samples[0]: np.int32}
    names.append(samples[0])
    df = pd.read_csv(StringIO(s), sep='\t', header=None, names=names,
                     dtype=dtype)
    return CovFrame(df)

class CovFrame:
    """Class for storing read depth data."""
    def __init__(self, df):
        self.df = df

    @property
    def shape(self):
        """Return the size of the CovFrame."""
        return self.df.shape
