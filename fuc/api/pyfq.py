"""
The ``pyfq`` submodule is designed for working with FASTQ files (both zipped
and unzipped).
"""

import gzip
import pandas as pd

def read_file(fn):
    """Create a FqFrame from a FASTQ file (both zipped and unzipped)."""
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rt')
    else:
        f = open(fn)
    data = {'ID': [], 'SEQ': [], 'EXT': [], 'QUAL': []}
    n = 0
    for line in f:
        x = line.strip()
        if n == 0:
            data['ID'].append(x)
            n += 1
        elif n == 1:
            data['SEQ'].append(x)
            n += 1
        elif n == 2:
            data['EXT'].append(x)
            n += 1
        else:
            data['QUAL'].append(x)
            n = 0
    df = pd.DataFrame.from_dict(data)
    qf = FqFrame(df)
    f.close()
    return qf

class FqFrame:
    """Class for storing FASTQ data."""
    def __init__(self, df):
        self.df = df

    @property
    def vdata(self):
        """Return a view (copy) of the data."""
        return deepcopy(self.data)

    @property
    def shape(self):
        """Return the size of the FqFrame."""
        return len(self.data)

    def readlen(self):
        """Return a dictionary of read lengths and their counts."""
        lengths = {}
        for r in self.data:
            length = len(r.seq)
            if length not in lengths:
                lengths[length] = 0
            lengths[length] += 1
        return lengths