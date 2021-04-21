"""
The FastqFrame module is designed for working with FASTQ files (both zipped
and unzipped).
"""

from dataclasses import dataclass, field
from typing import List
from copy import deepcopy
import gzip

@dataclass(unsafe_hash=True)
class FastqRecord:
    """Class for storing sequence read data."""
    id   : str = field(compare=False)
    seq  : str = field(compare=True)
    ext  : str = field(compare=False)
    qual : str = field(compare=False)

@dataclass
class FastqFrame:
    """Class for storing FASTQ data."""
    data : List[FastqRecord]

    @property
    def vdata(self):
        """Return a view (copy) of the data."""
        return deepcopy(self.data)

    @property
    def shape(self):
        """Return the size of the FastqFrame."""
        return len(self.data)

    @classmethod
    def from_file(cls, file_path):
        """Create a FastqFrame from a file."""
        if file_path.endswith('.gz'):
            f = gzip.open(file_path, 'rt')
        else:
            f = open(file_path)
        fields = []
        data = []
        n = 0
        for line in f:
            if n in [0, 1, 2]:
                fields.append(line.strip())
                n += 1
            else:
                fields.append(line.strip())
                data.append(FastqRecord(*fields))
                n = 0
                fields = []
        qf = cls(data)
        f.close()
        return qf

    def readlen(self):
        """Return a dictionary of read lengths and their counts."""
        lengths = {}
        for r in self.data:
            length = len(r.seq)
            if length not in lengths:
                lengths[length] = 0
            lengths[length] += 1
        return lengths
