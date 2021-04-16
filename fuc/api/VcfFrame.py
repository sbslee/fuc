"""
The VcfFrame module is designed for working with VCF files (both zipped
and unzipped).
"""

from typing import List
from dataclasses import dataclass
from copy import deepcopy

from .BedFrame import BedFrame

def has_var(x):
    """Return if the GT field has a variant (e.g. 0/1)."""
    return x.split(':')[0].replace('/', '').replace('.', '').replace('0', '')

@dataclass(eq=False)
class VcfRecord:
    """Class for storing the information of single VCF record."""
    chrom  : str       # CHROM
    pos    : int       # POS
    id     : str       # ID
    ref    : str       # REF
    alt    : List[str] # ALT
    qual   : str       # QUAL
    filter : List[str] # FILTER
    info   : List[str] # INFO
    format : List[str] # FORMAT
    gt     : List[str]

    def __eq__(self, other):
        """Test whether two VcfRecords are equal."""
        return (self.chrom, self.pos, self.ref, self.alt) == (
            other.chrom, other.pos, other.ref, other.alt)

    def to_list(self):
        """Convert the VcfRecord to a list of strings."""
        l = [self.chrom, str(self.pos), self.id, self.ref,
            ','.join(self.alt), self.qual, ';'.join(self.filter),
            ';'.join(self.info), ':'.join(self.format)]
        return l + self.gt

    @classmethod
    def from_list(cls, l):
        """Create a VcfRecord from a list of strings."""
        l[1] = int(l[1])       # POS
        l[4] = l[4].split(',') # ALT
        l[6] = l[6].split(';') # FILTER
        l[7] = l[7].split(';') # INFO
        l[8] = l[8].split(':') # FORMAT
        r = cls(*l[:9], l[9:])
        return r

@dataclass
class VcfFrame:
    HEADERS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT',
               'QUAL', 'FILTER', 'INFO', 'FORMAT']

    """Class for storing the information of single VCF file."""
    meta : List[str]
    head : List[str]
    data : List[VcfRecord]

    @property
    def vmeta(self):
        """Return a view (copy) of the metadata."""
        return deepcopy(self.meta)

    @property
    def vhead(self):
        """Return a view (copy) of the headers."""
        return deepcopy(self.head)

    @property
    def vdata(self):
        """Return a view (copy) of the data."""
        return deepcopy(self.data)

    @property
    def samples(self):
        """Return a list of the sample IDs."""
        return self.head[9:]

    @property
    def shape(self):
        """Return a tuple representing the dimensionality of the VcfFrame."""
        return (len(self.data), len(self.samples))

    def index(self, name):
        """Return the sample index."""
        return self.head.index(name) - 9

    def describe(self):
        """Generate descriptive statistics."""
        print('Samples:')
        print('Name', 'VariantCount')
        for name in self.samples:
            i = self.index(name)
            n = 0
            for r in self.data:
                if has_var(r.gt[i]):
                    n += 1
            print(name, n, sep='\t')

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
                meta.append(line.strip())
            elif line.startswith('#CHROM'):
                head = line.strip().split('\t')
            else:
                data.append(VcfRecord.from_list(line.strip().split('\t')))
        vf = cls(meta, head, data)
        f.close()
        return vf

    def to_file(self, file_path):
        """Write the VcfFrame to a file."""
        with open(file_path, 'w') as f:
            if self.meta:
                f.write('\n'.join(self.meta) + '\n')
            f.write('\t'.join(self.head) + '\n')
            for r in self.data:
                f.write('\t'.join(r.to_list()) + '\n')

    def to_string(self):
        """Render the VcfFrame to a console-friendly tabular output."""
        s = ''
        s += '\n'.join(self.meta) + '\n'
        s += '\t'.join(self.head) + '\n'
        for r in self.data:
            s += '\t'.join(r.to_list()) + '\n'
        return s

    def strip(self, format_subfields=None):
        """Remove unnecessary data from the VcfFrame.

        Parameters
        ----------
        format_subfields : list, optional
            Additional FORMAT subfields (e.g. DP and AD) to be retained
            other than GT, which is included as default.

        Returns
        -------
        vf : VcfFrame
            Stripped VcfFrame.
        """
        if format_subfields is None:
            _format_subfields = ['GT']
        else:
            if 'GT' in format_subfields:
                format_subfields.remove('GT')
            _format_subfields = ['GT'] + format_subfields
        meta = self.vmeta
        head = self.vhead
        data = []
        for r in self.vdata:
            r.id = '.'
            r.into = ['.']
            r.filter = ['.']
            indicies = [r.format.index(x) for x in _format_subfields]
            r.gt = [':'.join([x.split(':')[i] for i in indicies])
                for x in r.gt]
            r.format = _format_subfields
            data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def merge(self, other, format_subfields=None):
        """Merge with the other VcfFrame.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        format_subfields : list, optional
            Additional FORMAT subfields (e.g. DP and AD) to be retained
            other than GT, which is included as default.

        Returns
        -------
        vf : VcfFrame
            Stripped VcfFrame.
        """
        vf1 = self.strip(format_subfields=format_subfields)
        vf2 = other.strip(format_subfields=format_subfields)
        meta = []
        head = vf1.head + vf2.head[9:]
        data = []
        m = './.' # missing value
        for i in range(1, len(vf1.data[0].format)):
            m += ':.'
        for r in vf1.vdata:
            if r not in vf2.data:
                r.gt = r.gt + [m] * vf2.shape[1]
            data.append(r)
        for r in vf2.vdata:
            if r in data:
                data[data.index(r)].gt += r.gt
            else:
                r.gt = [m] * vf1.shape[1] + r.gt
                data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def compare(self, n1, n2):
        """Compare two samples within the VcfFrame.

        Parameters
        ----------
        n1 : string or int
            Test sample or its index in the header row.
        n2 : string or int
            Truth sample or its index in the header row.

        Returns
        -------
        result : tuple
            Comparison result (tp, fp, fn, and tn).
        """
        i1 = self.index(n1) if isinstance(n1, str) else n1
        i2 = self.index(n2) if isinstance(n2, str) else n2
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        for r in self.data:
            a = has_var(r.gt[i1])
            b = has_var(r.gt[i2])
            if a and b:
                tp += 1
            elif a and not b:
                fp += 1
            elif not a and b:
                fn += 1
            else:
                tn += 1
        result = (tp, fp, fn, tn)
        return result

    def multiallelic_sites(self):
        """Return the indicies of multiallelic sites."""
        indicies = []
        n = 0
        for r in self.data:
            if len(r.alt) > 1:
                indicies.append(n)
            n += 1
        return indicies

    def reset_samples(self, samples):
        """Reset the sample list."""
        meta = self.vmeta
        head = deepcopy(self.HEADERS)
        data = []
        indicies = []
        for sample in samples:
            indicies.append(self.index(sample))
        for i in indicies:
            head.append(self.samples[i])
        for r1 in self.vdata:
            r2 = deepcopy(r1)
            r2.gt = []
            for i in indicies:
                r2.gt.append(r1.gt[i])
            data.append(r2)
        vf = self.__class__(meta, head, data)
        return vf

    def add_dp(self):
        """Compute and add the DP subfield of the FORMAT field."""
        meta = self.vmeta
        head = self.vhead
        data = []
        def func(x, i):
            ad = x.split(':')[i].split(',')
            dp = 0
            for depth in ad:
                if depth == '.':
                    return f'{x}:.'
                dp += int(depth)
            return f'{x}:{dp}'
        for r in self.vdata:
            i = r.format.index('AD')
            r.gt = [func(x, i) for x in r.gt]
            r.format += ['DP']
            data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def filter_dp(self, threshold=200):
        """Filter based on the DP subfield of the FORMAT field."""
        meta = self.vmeta
        head = self.vhead
        data = []
        def func(x, i):
            l = x.split(':')
            dp = l[i]
            if dp == '.' or int(dp) < threshold:
                return './.' + ':.' * (len(l)-1)
            return x
        for r in self.vdata:
            i = r.format.index('DP')
            r.gt = [func(x, i) for x in r.gt]
            data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def filter_af(self, threshold=0.1):
        """Filter based on the AF subfield of the FORMAT field."""
        meta = self.vmeta
        head = self.vhead
        data = []
        def func(x, i):
            l = x.split(':')
            af = l[i]
            if af == '.' or float(af) < threshold:
                return './.' + ':.' * (len(l)-1)
            return x
        for r in self.vdata:
            i = r.format.index('AF')
            r.gt = [func(x, i) for x in r.gt]
            data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def filter_empty(self):
        """Filter out VcfRecords that are empty."""
        meta = self.vmeta
        head = self.vhead
        data = []
        for r in self.vdata:
            if all(['.' in x for x in r.gt]):
                continue
            data.append(r)
        vf = self.__class__(meta, head, data)
        return vf

    def filter_bed(self, bed):
        """Filter VcfRecords in the VcfFrame using BED data.

        Parameters
        ----------
        bed : BedFrame or string
            BedFrame or path to a BED file.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        if isinstance(bed, BedFrame):
            bf = bed
        else:
            bf = BedFrame.from_file(bed)
        meta = self.vmeta
        head = self.vhead
        data = []
        for r1 in self.vdata:
            for r2 in bf.data:
                if r1.chrom == r2.chrom and r2.start <= r1.pos <= r2.end:
                    data.append(r1)
                    break
        vf = self.__class__(meta, head, data)
        return vf

    def update(self, other, query_fields, missing_only=False):
        """Copy data from another VcfFrame.

        This method will copy requested data from another VcfFrame for
        overlapping records. You can only request data from the following
        VCF headers: ID, QUAL, FILTER, INFO, and FORMAT. Any other
        requested VCF headers will be ignored.

        Parameters
        ----------
        other : VcfFrame
            Target VcfFrame.
        names : list
            List of VCF headers.
        missing_only : boolean, optional
            If True, only fields with the missing value will be updated.

        Returns
        -------
        vcf_result : VcfFrame
            Updated VcfFrame.
        """
        target_fields = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        for query_field in query_fields:
            if query_field not in target_fields:
                raise ValueError(f'{query_field} is not found '
                    f'in the target fields: {target_fields}')
        meta = self.vmeta
        head = self.vhead
        data = []
        for r1 in self.vdata:
            if r1 in other.vdata:
                r2 = other.data[other.data.index(r1)]
                for query_field in query_fields:
                    v1 = getattr(r1, query_field)
                    v2 = getattr(r2, query_field)
                    if missing_only:
                        if v1 == '.' or v1 == ['.']:
                            setattr(r1, query_field, v2)
                    else:
                        setattr(r1, query_field, v2)
            data.append(r1)
        vf = self.__class__(meta, head, data)
        return vf
        return vcf_result
