import copy

class DataFrame():
    def __init__(self):
        self.head = []
        self.data = []

    @property
    def shape(self):
        """Return a tuple representing the dimensionality of the DataFrame."""
        return (len(self.data), len(self.head))

    def get_data(self):
        """Return a copy of the data which can be modified safely."""
        return copy.deepcopy(self.data)

    def write(self, file_path, delimiter='\t'):
        """Write the DataFrame to a file."""
        with open(file_path, 'w') as f:
            f.write(delimiter.join(self.head) + '\n')
            for fields in self.get_data():
                f.write(delimiter.join(fields) + '\n')

    def merge(self, other, on, missing='.'):
        """Return a merged DataFrame."""
        df = DataFrame()
        df.head = self.head + other.head
        on1 = self.head.index(on)
        on2 = other.head.index(on)
        other_indicies = []
        for f1 in self.get_data():
            match_found = False
            for i, f2 in enumerate(other.get_data()):
                if f1[on1] == f2[on2]:
                    f1 += f2
                    match_found = True
                    other_indicies.append(i)
                    break
            if not match_found:
                f1 += [missing] * len(f2)
            df.data.append(f1)
        for i in range(other.shape[0]):
            if i not in other_indicies:
                f1 = [missing] * self.shape[1]
                f2 = other.get_data()[i]
                f1[on1] = f2[on2]
                df.data.append(f1+f2)
        del df.head[self.shape[1] + on2]
        for fields in df.data:
            del fields[self.shape[1] + on2]
        return df

    @classmethod
    def read(cls, file_path, delimiter='\t', header=True):
        """Create a DataFrame from a file."""
        df = cls()
        with open(file_path) as f:
            if header:
                df.head = next(f).strip().split(delimiter)
            for line in f:
                fields = line.strip().split(delimiter)
                df.data.append(fields)
        return df
