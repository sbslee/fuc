"""
The VcfFrame module is designed for working with VCF files (both zipped
and unzipped).
"""

import pandas as pd
import gzip
from copy import deepcopy

class VcfFrame:
    """Class for storing VCF data.

    This class strictly sticks to the standard Variant Call Format
    specification (https://samtools.github.io/hts-specs/VCFv4.3.pdf).

    VCF lines have nine required fields for storing variant data and
    variable-length fields for storing sample genotype data. In all cases,
    missing values are specified with a dot ('.'). The required fields are:
        1. CHROM - An identifier from the reference genome.
        2. POS - The 1-based reference position.
        3. ID - Semicolon-separated list of unique identifiers.
        4. REF - Reference base(s).
        5. ALT - Comma-separated list of alternate non-reference alleles.
        6. QUAL - Phred-scaled quality score for the assertion made in ALT.
        7. FILTER - PASS or a semicolon-separated list of filters that fail.
        8. INFO - Semicolon-separated series of additional information fields.
        9. FORMAT - Colon-separated series of genotype fields.
    """
    def __init__(self, meta, data):
        self.meta = meta
        self.data = data

    @property
    def samples(self):
        """Return a list of the sample IDs."""
        return self.data.columns[9:].to_list()

    @classmethod
    def from_file(cls, file_path):
        """Create a VcfFrame from a VCF file."""
        meta = []
        skip_rows = 0
        if file_path.endswith('.gz'):
            f = gzip.open(file_path, 'rt')
        else:
            f = open(file_path)
        for line in f:
            if line.startswith('##'):
                meta.append(line.strip())
                skip_rows += 1
            else:
                break
        data = pd.read_table(file_path, skiprows=skip_rows)
        vf = cls(meta, data)
        f.close()
        return vf

    def to_file(self, file_path):
        """Write the VcfFrame to a VCF file."""
        with open(file_path, 'w') as f:
            if self.meta:
                f.write('\n'.join(self.meta) + '\n')
            self.data.to_csv(f, sep='\t', index=False)

    def strip(self, format='GT'):
        """Remove unnecessary data from the VcfFrame.

        Parameters
        ----------
        format : str, default: 'GT'
            FORMAT subfields to be retained (e.g. 'GT:AD:DP').

        Returns
        -------
        vf : VcfFrame
            Stripped VcfFrame.
        """
        def outfunc(r):
            idx = [r['FORMAT'].split(':').index(x) for x in format.split(':')]
            infunc = lambda x: ':'.join([x.split(':')[i] for i in idx])
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        data = self.data.copy()
        data[['ID', 'QUAL', 'FILTER', 'INFO']] = '.'
        data = data.apply(outfunc, axis=1)
        data['FORMAT'] = format
        vf = self.__class__([], data)
        return vf

    def merge(self, other, how='inner', format='GT'):
        """Merge with the other VcfFrame.

        This method essentially wraps the `pandas.DataFrame.merge` method.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        how : str, default: 'inner'
            Type of merge to be performed. ['left', 'right', 'outer',
            'inner', 'cross']
        format : str, default: 'GT'
            FORMAT subfields to be retained (e.g. 'GT:AD:DP').

        Returns
        -------
        vf : VcfFrame
            Merged VcfFrame.
        """
        vf1 = self.strip(format=format)
        vf2 = other.strip(format=format)
        dropped = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        shared = ['#CHROM', 'POS', 'REF', 'ALT']
        data = vf1.data.merge(vf2.data.drop(columns=dropped),
            on=shared, how=how)
        data[dropped] = data[dropped].fillna('.')
        data['FORMAT'] = format
        def func(r):
            n = len(r['FORMAT'].split(':'))
            x = '/'.join(['.' for x in r.iloc[9:].dropna()[0].split('/')])
            for i in range(1, n):
                x += ':.'
            r = r.fillna(x)
            return r
        data = data.apply(func, axis=1)
        vf3 = self.__class__([], data)
        return vf3

    def add_dp(self):
        """Compute and add the DP subfield of the FORMAT field."""
        def outfunc(r):
            i = r['FORMAT'].split(':').index('AD')
            def infunc(x):
                ad = x.split(':')[i].split(',')
                dp = 0
                for depth in ad:
                    if depth == '.':
                        return f'{x}:.'
                    dp += int(depth)
                return f'{x}:{dp}'
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            r['FORMAT'] += ':DP'
            return r
        data = self.data.apply(outfunc, axis=1)
        vf = self.__class__(deepcopy(self.meta), data)
        return vf

    def filter_dp(self, threshold=200):
        """Filter rows based on the DP subfield of the FORMAT field."""
        def outfunc(r):
            i = r['FORMAT'].split(':').index('DP')
            def infunc(x):
                l = x.split(':')
                dp = l[i]
                if dp == '.' or int(dp) < threshold:
                    return './.' + ':.' * (len(l)-1)
                return x
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        data = self.data.apply(outfunc, axis=1)
        vf = self.__class__(deepcopy(self.meta), data)
        return vf

    def filter_af(self, threshold=0.1):
        """Filter rows based on the AF subfield of the FORMAT field."""
        def outfunc(r):
            i = r['FORMAT'].split(':').index('AF')
            def infunc(x):
                l = x.split(':')
                af = l[i]
                if af == '.' or float(dp) < threshold:
                    return './.' + ':.' * (len(l)-1)
                return x
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        data = self.data.apply(outfunc, axis=1)
        vf = self.__class__(deepcopy(self.meta), data)
        return vf

    def filter_empty(self):
        """Filter out rows that have no genotypes."""
        def func(r):
            return not all(r.iloc[9:].apply(lambda x: '.' in x.split(':')[0]))
        i = self.data.apply(func, axis=1)
        data = self.data[i]
        vf = self.__class__(deepcopy(self.meta), data)
        return vf
