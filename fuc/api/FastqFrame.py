from dataclasses import dataclass
from typing import List
from copy import deepcopy

@dataclass(eq=False)
class FastqRecord:
    id   : str
    seq  : str
    ext  : str
    qual : str

    def __eq__(self, other):
        """Test whether two FastqRecords are equal."""
        return self.seq == other.seq

@dataclass
class FastqFrame:
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
            f = gzip.open(file_path, 'rb')
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
