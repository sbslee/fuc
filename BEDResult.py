import copy

class BEDResult():
    def __init__(self):
        self.head = []
        self.data = []

    @property
    def shape(self):
        """Return the size of the VCFResult."""
        return len(self.data)

    def get_data(self):
        """Return a copy of the data which can be modified safely."""
        return copy.deepcopy(self.data)

    def write(self, bed_path):
        """Write the BEDResult to a file."""
        with open(bed_path, 'w') as f:
            f.write(''.join(self.head))
            for fields in self.get_data():
                f.write('\t'.join(fields) + '\n')

    def intersect(self, other):
        """Return a BEDResult consisting of intersections.

        This method will compute intersections between multiple BED files.
        Header lines and optional BED fields in the other BEDResult will be
        ignored.

        Parameters
        ----------
        other : BEDResult
            Other BEDResult.

        Returns
        -------
        bed_result: BEDResult
            Updated BEDResult.
        """
        def overlap(a, b):
            if b[0] <= a[0] <= b[1]:
                start = a[0]
            elif a[0] <= b[0] <= a[1]:
                start = b[0]
            else:
                return None
            if b[0] <= a[1] <= b[1]:
                end = a[1]
            elif a[0] <= b[1] <= a[1]:
                end = b[1]
            else:
                return None
            return (start, end)
        bed_result = self.__class__()
        bed_result.head = copy.deepcopy(self.head)
        for r1 in self.get_data():
            a = (int(r1[1]), int(r1[2]))
            for r2 in other.data:
                b = (int(r2[1]), int(r2[2]))
                overlap_result = overlap(a, b)
                if overlap_result is not None:
                    r1[1] = str(overlap_result[0])
                    r1[2] = str(overlap_result[1])
                    bed_result.data.append(r1)
        return bed_result

    @classmethod
    def read(cls, bed_path):
        """Create a BEDResult from a file."""
        bed_result = cls()
        with open(bed_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    bed_result.head = line
                    continue
                bed_result.data.append(fields)
        return bed_result
