import copy
import gzip

class FASTQResult():
    def __init__(self):
        self.data = []

    @property
    def shape(self):
        """Return the size of the FASTQResult."""
        return len(self.data)

    def get_data(self):
        """Return a copy of the data which can be modified safely."""
        return copy.deepcopy(self.data)

    def get_read_length(self):
        """Return a dictionary of read lengths and their counts."""
        lengths = {}
        for fields in self.get_data():
            length = len(fields[1])
            if length not in lengths:
                lengths[length] = 0
            lengths[length] += 1
        return lengths

    @classmethod
    def read(cls, fastq_path):
        """Create a FASTQResult from a file."""
        fastq_result = cls()
        n = 0
        fields = []
        if fastq_path.endswith('.gz'):
            f = gzip.open(fastq_path, 'rb')
        else:
            f = open(fastq_path)
        for line in f:
            if n in [0, 1, 2]:
                fields.append(line)
                n += 1
            else:
                fields.append(line)
                fastq_result.data.append(fields)
                n = 0
                fields = []
        f.close()
        return fastq_result