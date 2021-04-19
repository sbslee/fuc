import pandas as pd
import pyranges as pr
from copy import deepcopy

BF_HEADERS = ['chrom', 'chromStart', 'chromEnd', 'name',
              'score', 'strand', 'thickStart', 'thickEnd',
              'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

PR_HEADERS = ['Chromosome', 'Start', 'End', 'Name',
              'Score', 'Strand', 'ThickStart', 'ThickEnd',
              'ItemRGB', 'BlockCount', 'BlockSizes', 'BlockStarts']

class BedFrame:
    """Class for storing BED data."""
    def __init__(self, meta, data):
        self.meta = meta
        self.data = data

    @classmethod
    def from_file(cls, file_path):
        """Create a BedFrame from a BED file."""
        meta = []
        skip_rows = 0
        with open(file_path, 'r') as f:
            for line in f:
                if 'browser' in line or 'track' in line:
                    meta.append(line.strip())
                    skip_rows += 1
                else:
                    headers = BF_HEADERS[:len(line.strip().split())]
                    break
        data = pd.read_table(file_path, header=None, names=headers,
            skiprows=skip_rows)
        bf = cls(meta, data)
        return bf

    def to_file(self, file_path):
        """Write the BedFrame to a BED file."""
        with open(file_path, 'w') as f:
            if self.meta:
                f.write('\n'.join(self.meta) + '\n')
            self.data.to_csv(f, sep='\t', index=False, header=False)

    def intersect(self, other):
        """Find intersection between the BedFrames."""
        n = min(self.data.shape[1], other.data.shape[1])
        d = dict(zip(BF_HEADERS[:n], PR_HEADERS[:n]))
        df1 = deepcopy(self.data).rename(columns=d)
        df2 = deepcopy(other.data).rename(columns=d)
        gr1, gr2 = pr.PyRanges(df1), pr.PyRanges(df2)
        d = dict(zip(PR_HEADERS[:n], BF_HEADERS[:n]))
        data = gr1.intersect(gr2).df.rename(columns=d)
        bf = self.__class__(deepcopy(self.meta), data)
        return bf
