"""
The VcfFrame module is designed for working with VCF files (both zipped
and unzipped).
"""

import pandas as pd

class VcfFrame:
    """Class for storing VCF data."""
    def __init__(self, meta, data):
        self.meta = meta
        self.data = data

    @classmethod
    def from_file(cls, file_path):
        """Create a VcfFrame from a VCF file."""
        meta = []
        skip_rows = 0
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    meta.append(line.strip())
                    skip_rows += 1
                else:
                    break
        data = pd.read_table(file_path, skiprows=skip_rows)
        vf = cls(meta, data)
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
        def outfunc(r):
            n = len(r['FORMAT'].split(':'))
            x = '/'.join(['.' for x in r.iloc[9:].dropna()[0].split('/')])
            for i in range(1, n):
                x += ':.'
            r = r.fillna(x)
            return r
        data = data.apply(outfunc, axis=1)
        vf3 = self.__class__([], data)
        return vf3
