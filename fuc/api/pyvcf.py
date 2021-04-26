"""
The ``pyvcf`` submodule is designed for working with VCF files (both zipped
and unzipped). It implements ``pyvcf.VcfFrame`` which stores VCF data as a
``pandas.DataFrame`` to allow fast computation and easy manipulation.
"""

import pandas as pd
import gzip
from copy import deepcopy
from . import pybed

CONTIGS = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
    '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M',
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
    'chrX', 'chrY', 'chrM'
]

def read_file(fn):
    """Create a VcfFrame from a VCF file (both zipped and unzipped)."""
    meta = []
    skip_rows = 0
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rt')
    else:
        f = open(fn)
    for line in f:
        if line.startswith('##'):
            meta.append(line.strip())
            skip_rows += 1
        else:
            break
    df = pd.read_table(fn, skiprows=skip_rows)
    df = df.rename(columns={'#CHROM': 'CHROM'})
    vf = VcfFrame(meta, df)
    f.close()
    return vf

def merge(vfs, how='inner', format='GT', sort=True, collapse=False):
    """Return the merged VcfFrame.

    Parameters
    ----------
    vfs : list
        List of VcfFrames to be merged.
    how : str, default: 'inner'
        Type of merge as defined in `pandas.DataFrame.merge`.
    format : str, default: 'GT'
        FORMAT subfields to be retained (e.g. 'GT:AD:DP').
    sort : bool, default: True
        If True, sort the VcfFrame before returning.
    collapse : bool, default: False
        If True, collapse duplicate records.

    Returns
    -------
    merged_vf : VcfFrame
        Merged VcfFrame.
    """
    merged_vf = vfs[0]
    for vf in vfs[1:]:
        merged_vf = merged_vf.merge(vf, how=how, format=format, sort=sort,
            collapse=collapse)
    return merged_vf

def gt_hasvar(x):
    """Return True if the genotype has a variant call (e.g. ``0/1``)."""
    if x.split(':')[0].replace('/', '').replace('.', '').replace('0', ''):
        return True
    else:
        return False

def gt_missing(x):
    """Return True if the genotype has a missing value."""
    return '.' in x.split(':')[0]

def gt_unphase(x):
    """Unphase the genotype (e.g. if ``GT:DP``, ``1|0:20`` to ``0/1:20``)."""
    l = x.split(':')
    gt = l[0]
    if '|' not in gt:
        return x
    l[0] = '/'.join([str(b) for b in sorted([int(a) for a in gt.split('|')])])
    return ':'.join(l)

def row_missing_value(r):
    """Return the missing value for the row (e.g. ``./.:.`` if ``GT:DP``)."""
    if 'X' in r.CHROM or 'Y' in r.CHROM:
        m = '.'
    else:
        m = './.'
    for i in range(1, len(r.FORMAT.split(':'))):
        m += ':.'
    return m

class VcfFrame:
    """Class for storing VCF data.

    This class strictly sticks to the standard Variant Call Format
    specification (https://samtools.github.io/hts-specs/VCFv4.3.pdf).

    VCF lines have nine required fields for storing variant data and
    variable-length fields for storing sample genotype data. In all cases,
    missing values are specified with a dot (``.``). The required fields are:

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
    def __init__(self, meta, df):
        self.meta = meta
        self.df = df

    @property
    def samples(self):
        """Return a list of the sample IDs."""
        return self.df.columns[9:].to_list()

    @property
    def shape(self):
        """Return a tuple representing the dimensionality of the VcfFrame."""
        return (self.df.shape[0], len(self.samples))

    def copy_meta(self):
        """Return a copy of the metadata."""
        return deepcopy(self.meta)

    def copy_df(self):
        """Return a copy of the dataframe."""
        return self.df.copy()

    def to_file(self, file_path):
        """Write the VcfFrame to a VCF file."""
        with open(file_path, 'w') as f:
            if self.meta:
                f.write('\n'.join(self.meta) + '\n')
            self.df.rename(columns={'CHROM': '#CHROM'}).to_csv(
                f, sep='\t', index=False)

    def to_string(self):
        """Render the VcfFrame to a console-friendly tabular output."""
        s = ''
        if self.meta:
            s += '\n'.join(self.meta) + '\n'
        s += self.df.to_csv(index=False, sep='\t')
        return s

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
            idx = [r.FORMAT.split(':').index(x) for x in format.split(':')]
            infunc = lambda x: ':'.join([x.split(':')[i] for i in idx])
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        df = self.df.copy()
        df[['ID', 'QUAL', 'FILTER', 'INFO']] = '.'
        df = df.apply(outfunc, axis=1)
        df.FORMAT = format
        vf = self.__class__([], df)
        return vf

    def merge(self, other, how='inner', format='GT', sort=True,
              collapse=False):
        """Merge with the other VcfFrame.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        how : str, default: 'inner'
            Type of merge as defined in `pandas.DataFrame.merge`.
        format : str, default: 'GT'
            FORMAT subfields to be retained (e.g. 'GT:AD:DP').
        sort : bool, default: True
            If True, sort the VcfFrame before returning.
        collapse : bool, default: False
            If True, collapse duplicate records.

        Returns
        -------
        vf : VcfFrame
            Merged VcfFrame.
        """
        vf1 = self.strip(format=format)
        vf2 = other.strip(format=format)
        dropped = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        shared = ['CHROM', 'POS', 'REF', 'ALT']
        df = vf1.df.merge(vf2.df.drop(columns=dropped), on=shared, how=how)

        # This ensures that the column order is intact when either of the
        # dataframes is empty.
        cols = vf1.df.columns.to_list() + vf2.df.columns[9:].to_list()
        df = df[cols]

        df[dropped] = df[dropped].fillna('.')
        df.FORMAT = format
        def func(r):
            n = len(r.FORMAT.split(':'))
            x = './.'
            for i in range(1, n):
                x += ':.'
            r = r.fillna(x)
            return r
        df = df.apply(func, axis=1)
        vf3 = self.__class__([], df)
        if collapse:
            vf3 = vf3.collapse()
        if sort:
            vf3 = vf3.sort()
        return vf3

    def collapse(self):
        """Collapse duplicate records.

        Duplicate records have the identical values for ``CHROM``, ``POS``,
        and ``REF``. They can result from merging two VCF files.

        .. note::
           The method will sort the order of ALT alleles.

        Returns
        -------
        vf : VcfFrame
            Collapsed VcfFrame.

        Examples
        --------
        Consider the following:

        +-------+-----+-----+-----+--------+---------------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven        | Sara        |
        +=======+=====+=====+=====+========+===============+=============+
        | chr1  | 100 | A   | C   | GT:AD  | 0/1:12,15     | ./.:.       |
        +-------+-----+-----+-----+--------+---------------+-------------+
        | chr1  | 100 | A   | T   | GT:AD  | ./.:.         | 0/1:14,15   |
        +-------+-----+-----+-----+--------+---------------+-------------+

        After collapsing, above will look like this:

        +-------+-----+-----+-----+--------+---------------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven        | Sara        |
        +=======+=====+=====+=====+========+===============+=============+
        | chr1  | 100 | A   | C,T | GT:AD  | 0/1:12,15,0   | 0/2:14,0,15 |
        +-------+-----+-----+-----+--------+---------------+-------------+

        Now here is a slightly more complex scenario:

        +-------+-----+-----+-----+--------+---------------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven        | James       |
        +=======+=====+=====+=====+========+===============+=============+
        | chr1  | 100 | A   | T   | GT:AD  | 0/1:12,15     | ./.:.       |
        +-------+-----+-----+-----+--------+---------------+-------------+
        | chr1  | 100 | A   | T,C | GT:AD  | ./.:.         | 1/2:0,11,17 |
        +-------+-----+-----+-----+--------+---------------+-------------+

        After collapsing, above will look like this:

        +-------+-----+-----+-----+--------+---------------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven        | James       |
        +=======+=====+=====+=====+========+===============+=============+
        | chr1  | 100 | A   | C,T | GT:AD  | 0/2:12,0,15   | 1/2:0,17,11 |
        +-------+-----+-----+-----+--------+---------------+-------------+
        """
        df = self.df.copy()
        dup_idx = df.duplicated(['CHROM', 'POS', 'REF'], keep=False)
        dups = {}
        for i, r in df[dup_idx].iterrows():
            name = f'{r.CHROM}:{r.POS}:{r.REF}'
            if name not in dups:
                dups[name] = []
            dups[name].append(i)

        def collapse_one(df):
            ref_allele = df.REF.unique()[0]
            alt_alleles = []
            for i, r in df.iterrows():
                alt_alleles += r.ALT.split(',')
            alt_alleles = sorted(list(set(alt_alleles)))
            all_alleles = [ref_allele] + alt_alleles

            def infunc(x, r_all_alleles, index_map):
                if gt_missing(x):
                    return ''
                old_fields = x.split(':')
                old_gt = old_fields[0]
                new_gt = '/'.join([str(x) for x in
                    sorted([index_map[int(i)] for i in old_gt.split('/')])
                ])
                new_fields = [new_gt]
                for old_field in old_fields[1:]:
                    old_subfields = old_field.split(',')
                    new_subfields = ['0' for x in all_alleles]
                    if len(old_subfields) == len(r_all_alleles):
                        for i, old_subfield in enumerate(old_subfields):
                            new_subfields[index_map[i]] = old_subfield
                        new_fields.append(','.join(new_subfields))
                    else:
                        new_fields.append(old_field)
                return ':'.join(new_fields)

            def outfunc(r):
                r_alt_alleles = r.ALT.split(',')
                r_all_alleles = [r.REF] + r_alt_alleles
                old_indicies = [i for i in range(len(r_all_alleles))]
                new_indicies = [all_alleles.index(x) for x in r_all_alleles]
                index_map = dict(zip(old_indicies, new_indicies))
                r[9:] = r[9:].apply(infunc, args=(r_all_alleles, index_map))
                return r

            df2 = df.apply(outfunc, axis=1)
            df2 = df2.groupby(['CHROM', 'POS', 'REF']).agg(''.join)
            df2 = df2.reset_index()
            cols = list(df2)
            cols[2], cols[3] = cols[3], cols[2]
            df2 = df2[cols]
            df2.ID = df.ID.unique()[0]
            df2.ALT = ','.join(alt_alleles)
            df2.QUAL = df.QUAL.unique()[0]
            df2.FILTER = df.FILTER.unique()[0]
            df2.INFO = df.INFO.unique()[0]
            df2.FORMAT = df.FORMAT.unique()[0]
            s = df2.squeeze()
            s = s.replace('', row_missing_value(s))
            return s

        for name, i in dups.items():
            df.iloc[i] = collapse_one(df.iloc[i])
        df.drop_duplicates(subset=['CHROM', 'POS', 'REF'], inplace=True)

        vf = self.__class__(self.copy_meta(), df)
        return vf

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
        df = self.df.apply(outfunc, axis=1)
        vf = self.__class__(self.copy_meta(), df)
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
        df = self.df.apply(outfunc, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_af(self, threshold=0.1):
        """Filter rows based on the AF subfield of the FORMAT field."""
        def outfunc(r):
            i = r['FORMAT'].split(':').index('AF')
            def infunc(x):
                l = x.split(':')
                af = l[i]
                if af == '.' or float(af) < threshold:
                    return './.' + ':.' * (len(l)-1)
                return x
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        df = self.df.apply(outfunc, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_empty(self, include=False):
        """Filter out rows that have no genotype calls.

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            empty = all(r.iloc[9:].apply(lambda x: '.' in x.split(':')[0]))
            if include:
                return empty
            else:
                return not empty
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_indels(self, include=False):
        """Filter out rows that have an indel.

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            ref = len(r['REF']) > 1
            alt = max([len(x) for x in r['ALT'].split(',')]) > 1
            has_indel = any([ref, alt])
            if include:
                return has_indel
            else:
                return not has_indel
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_multiallelic(self, include=False):
        """Filter out rows that have multiple ALT alleles.

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            is_multiallelic = ',' in r['ALT']
            if include:
                return is_multiallelic
            else:
                return not is_multiallelic
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_bed(self, bed, include=False):
        """Filter out rows that are not specified in the BED data.

        Parameters
        ----------
        bed : pybed.BedFrame or str
            pybed.BedFrame or path to a BED file.
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        if isinstance(bed, pybed.BedFrame):
            bf = bed
        else:
            bf = pybed.read_file(bed)
        def func(r):
            not_in_bed = bf.gr[r['CHROM'], r['POS']:r['POS']+1].empty
            if include:
                return not_in_bed
            else:
                return not not_in_bed
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_polyploid(self, include=False):
        """Filter out rows that have polyploid genotype (e.g. 0/0/1).

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            has_polyploid = any(
                [x.split(':')[0].count('/') > 1 for x in r[9:]])
            if include:
                return has_polyploid
            else:
                return not has_polyploid
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_sample_counts(self, threshold, include=False):
        """Filter out rows based on the number of samples with the variant.

        Parameters
        ----------
        threshold : int
            Number of samples above which will signify the row to be removed.
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            is_filtered = r[9:].apply(gt_hasvar).sum() > threshold
            if include:
                return is_filtered
            else:
                return not is_filtered
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_by_sample(self, name, include=False):
        """Filter out rows that have a variant call in the target sample.

        Parameters
        ----------
        name : str or int
            Name or index of the sample.
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        name = name if isinstance(name, str) else self.samples[name]
        def func(r):
            is_filtered = gt_hasvar(r[name])
            if include:
                return is_filtered
            else:
                return not is_filtered
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def filter_nonpass(self, include=False):
        """Filter out rows that don't have PASS in the FILTER field.

        Parameters
        ----------
        include : bool, default: False
            If True, include only such rows instead of excluding them.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        def func(r):
            is_filtered = r.FILTER != 'PASS'
            if include:
                return is_filtered
            else:
                return not is_filtered
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def compare(self, n1, n2):
        """Compare two samples within the VcfFrame.

        Parameters
        ----------
        n1 : str or int
            Name or index of the test sample.
        n2 : str or int
            Name or index of the truth sample.

        Returns
        -------
        result : tuple
            Comparison result (tp, fp, fn, tn).
        """
        n1 = n1 if isinstance(n1, str) else self.samples[n1]
        n2 = n2 if isinstance(n2, str) else self.samples[n2]
        def func(r):
            a = gt_hasvar(r[n1])
            b = gt_hasvar(r[n2])
            if a and b:
                return 'tp'
            elif a and not b:
                return 'fp'
            elif not a and b:
                return 'fn'
            else:
                return 'tn'
        d = self.df.apply(func, axis=1).value_counts().to_dict()
        tp = d['tp'] if 'tp' in d else 0
        fp = d['fp'] if 'fp' in d else 0
        fn = d['fn'] if 'fn' in d else 0
        tn = d['tn'] if 'tn' in d else 0
        result = (tp, fp, fn, tn)
        return result

    def combine(self, n1, n2):
        """Return a new column after combining data from the two samples.

        This method is useful when, for example, you are trying to
        consolidate data from multiple replicate samples. When the same
        variant is found in both samples, the method will use the genotype
        data of the first sample.

        Parameters
        ----------
        n1 : str or int
            Name or index of the first (or original) sample.
        n2 : str or int
            Name or index of the second (or replicate) sample.

        Returns
        -------
        s : pandas.Series
            VCF column representing the combined data.
        """
        n1 = n1 if isinstance(n1, str) else self.samples[n1]
        n2 = n2 if isinstance(n2, str) else self.samples[n2]
        def func(r):
            a = gt_hasvar(r[n1])
            b = gt_hasvar(r[n2])
            if a and b:
                return r[n1]
            elif a and not b:
                return r[n1]
            elif not a and b:
                return r[n2]
            else:
                return r[n1]
        s = self.df.apply(func, axis=1)
        return s

    def reset_samples(self, names):
        """Reset the sample list."""
        df = self.df[self.df.columns[:9].to_list() + names]
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def update(self, other, headers=None, missing=True):
        """Copy data from the other VcfFrame.

        This method will copy and paste data from the other VcfFrame for
        overlapping records. By default, the following VCF headers are
        used: ID, QUAL, FILTER, and, INFO.

        Parameters
        ----------
        other : VcfFrame
            Other VcfFrame.
        headers : list, optional
            List of VCF headers to exclude.
        missing : bool, default: True
            If True, only fields with the missing value ('.') will be updated.

        Returns
        -------
        vf : VcfFrame
            Updated VcfFrame.
        """
        targets = ['ID', 'QUAL', 'FILTER', 'INFO']
        if headers is not None:
            for header in headers:
                targets.remove(header)
        def func(r1):
            r2 = other.df[(other.df['CHROM'] == r1['CHROM']) &
                          (other.df['POS'] == r1['POS']) &
                          (other.df['REF'] == r1['REF']) &
                          (other.df['ALT'] == r1['ALT'])]
            if r2.empty:
                return r1
            for target in targets:
                if missing:
                    if r1[target] == '.':
                        r1[target] = r2.iloc[0][target]
                    else:
                        pass
                else:
                    r1[target] = r2.iloc[0][target]
            return r1
        df = self.df.apply(func, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def sort(self):
        """Return the sorted VcfFrame."""
        df = self.df.sort_values(by=['CHROM', 'POS'], ignore_index=True,
            key=lambda col: [CONTIGS.index(x) if isinstance(x, str)
                             else x for x in col])
        vf = self.__class__(self.copy_meta(), df)
        return vf

    def unphase(self):
        """Unphase the genotypes (e.g. 1|0 to 0/1)."""
        def func(r):
            r[9:] = r[9:].apply(gt_unphase)
            return r
        df = self.df.apply(func, axis=1)
        vf = self.__class__(self.copy_meta(), df)
        return vf
