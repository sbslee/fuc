from typing import List
from dataclasses import dataclass, field

@dataclass
class VcfRecord:
    """Class for storing the information of a single VCF record (row)."""
    chrom : str = ''
    pos : int = 0
    id : str = ''
    ref : str = ''
    alt : List[str] = field(default_factory=list)
    qual : str = ''
    filter : str = ''
    info : str = ''
    format : str = ''
    gt : List[str] = field(default_factory=list)

    @classmethod
    def from_list(cls, l):
        l[1] = int(l[1])
        l[4] = l[4].split(',')
        vf = cls(*l[:9], l[9:])
        return vf

@dataclass
class VcfFrame:
    meta : List[str] = field(default_factory=list)
    head : List[str] = field(default_factory=list)
    data : List[VcfRecord] = field(default_factory=list)

    @classmethod
    def from_file(cls, file_path):
        """Create a VcfFrame from a file."""
        meta = []
        head = []
        data = []
        if file_path.endswith('.gz'):
            f = gzip.open(file_path, 'rb')
        else:
            f = open(file_path)
        for line in f:
            if line.startswith('##'):
                meta.append(line)
            elif line.startswith('#CHROM'):
                head = line.strip().split('\t')
            else:
                fields = line.strip().split('\t')
                record = VcfRecord.from_list(fields)
                data.append(record)
        vf = cls(meta, head, data)
        f.close()
        return vf

    def to_file(self, file_path):
        """Write the VcfFrame to a file."""
        with open(file_path, 'w') as f:
            f.write(''.join(self.meta))
            f.write('\t'.join(self.head) + '\n')
            for record, fields in self.get_data().items():
                f.write('\t'.join(fields) + '\n')
