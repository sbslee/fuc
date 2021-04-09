import copy
import statistics
from collections import Counter
from operator import itemgetter

from common import is_numeric

class DataFrame():
    def __init__(self):
        self.head = []
        self.data = []
        self.dtypes = []

    @property
    def shape(self):
        """Return a tuple representing the dimensionality of the DataFrame."""
        return (len(self.data), len(self.head))

    def get_index(self, name):
        """Return the column index."""
        return self.head.index(name)

    def get_col(self, i):
        """Return the column as a list."""
        return [x[i] for x in self.get_data()]

    def get_unique(self, i):
        """Return unique values of the column."""
        return list(dict.fromkeys(self.get_col(i)))

    def summarize_col(self, i):
        """Return summary of the column."""
        col = self.get_col(i)
        if self.dtypes[i] == 'categorical':
            results = dict(Counter(col))
        else:
            col = [float(x) for x in col]
            maximum = max(col)
            minimum = min(col)
            mean = '{:.4f}'.format(statistics.mean(col))
            median = statistics.median(col)
            records = len(col)
            results = dict(minimum=minimum, maximum=maximum, mean=mean,
                median=median, records=records)
        return results

    def get_head(self):
        """Return a copy of the headers which can be modified safely."""
        return copy.deepcopy(self.head)

    def get_data(self):
        """Return a copy of the data which can be modified safely."""
        return copy.deepcopy(self.data)

    def write(self, file_path, delimiter='\t'):
        """Write the DataFrame to a file."""
        with open(file_path, 'w') as f:
            f.write(delimiter.join(self.head) + '\n')
            for fields in self.get_data():
                f.write(delimiter.join(fields) + '\n')

    def filter_columns(self, headers, exclude=False):
        """Return a new table after filtering the columns."""
        # Get the final columns to be included.
        if exclude:
            cols = [x for x in self.get_head() if x not in headers]
        else:
            cols = [x for x in self.get_head() if x in headers]
        # Create a filtered table.
        df = self.__class__()
        if len(cols) == 1:
            i = self.get_index(cols[0])
            df.head = [self.head[i]]
            for fields in self.get_data():
                df.data.append([fields[i]])
            df.dtypes = [self.dtypes[i]]
        else:
            indicies = [self.get_index(x) for x in cols]
            df.head = list(itemgetter(*indicies)(self.head))
            for fields in self.get_data():
                df.data.append(list(itemgetter(*indicies)(fields)))
            df.dtypes = list(itemgetter(*indicies)(self.dtypes))
        return df

    def filter_rows(self, key, values, exclude=False):
        """Return a new table after filtering the rows."""
        df = self.__class__()
        df.head = copy.deepcopy(self.head)
        df.data = copy.deepcopy(self.data)
        df.dtypes = copy.deepcopy(self.dtypes)
        i = self.get_index(key)
        if exclude:
            df.data = [x for x in df.data if x[i] not in values]
        else:
            df.data = [x for x in df.data if x[i] in values]
        return df

    def merge(self, other, on, missing='.'):
        """Return a merged DataFrame."""
        df = DataFrame()
        df.head = self.head + other.head
        on1 = [self.get_index(x) for x in on]
        on2 = [other.get_index(x) for x in on]
        other_indicies = []
        for f1 in self.get_data():
            match_found = False
            for i, f2 in enumerate(other.get_data()):
                if [f1[x] for x in on1] == [f2[x] for x in on2]:
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
                for j in range(len(on1)):
                    f1[on1[j]] = f2[on2[j]]
                df.data.append(f1+f2)
        # Delete the redundant headers.
        for i in sorted(on2, reverse=True):
            del df.head[self.shape[1] + i]
        # Delete the redundant columns.
        for fields in df.data:
            for i in sorted(on2, reverse=True):
                del fields[self.shape[1] + i]
        return df

    @classmethod
    def read(cls, file_path, delimiter='\t', header=True, skiprows=None):
        """Create a DataFrame from a file."""
        df = cls()
        with open(file_path, encoding='utf-8-sig') as f:
            if header:
                df.head = next(f).strip().split(delimiter)
            for line in f:
                fields = line.strip().split(delimiter)
                df.data.append(fields)
        # If provided, remove the specified rows.
        if isinstance(skiprows, list):
            pass
        else:
            skiprows = []
        for i in sorted(skiprows, reverse=True):
            del df.data[i]
        # Get the data type for each column.
        for i in range(df.shape[1]):
            l = []
            for fields in df.get_data():
                l.append(fields[i])
            if all([is_numeric(x) for x in l]):
                df.dtypes.append('numeric')
            else:
                df.dtypes.append('categorical')
        return df
