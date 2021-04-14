from dataclasses import dataclass

class VcfFrame:
    def __init__(self):
        self.meta = []
        self.head = []
        self.data = []

    @classmethod
    def read(cls, vcf_path):
        """Create a VCFResult from a file."""
        vcf_result = cls()
        if vcf_path.endswith('.gz'):
            f = gzip.open(vcf_path, 'rb')
        else:
            f = open(vcf_path)
        for line in f:
            if line.startswith('##'):
                vcf_result.meta.append(line)
            elif line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                vcf_result.head = fields
            else:
                fields = line.strip().split('\t')
                record = VcfRecord.from_list(fields)
                vcf_result.data.append(record)
        f.close()
        return vcf_result

@dataclass
class VcfRecord:
    chrom : str
    pos : int
    id : str
    ref : str
    alt : str
    qual : str
    filter : str
    info : str
    format : str
    fields : list

    @classmethod
    def from_list(cls, l):
        vf = cls(*l[:9])
        return vf
