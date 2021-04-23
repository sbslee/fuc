"""
The ``pyvcf`` submodule is designed for working with VCF files (both zipped
and unzipped). It implements `pyvcf.VcfFrame` which stores VCF data as
a ``pandas.DataFrame`` to allow fast computation and easy minipulation.
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

def merge(vfs, how='inner', format='GT', sort=True):
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

    Returns
    -------
    merged_vf : VcfFrame
        Merged VcfFrame.
    """
    merged_vf = vfs[0]
    for vf in vfs[1:]:
        merged_vf = merged_vf.merge(vf, how=how, format=format, sort=sort)
    return merged_vf

def hasvar(x):
    """Return True if the GT field has a variant (e.g. 0/1)."""
    return x.split(':')[0].replace('/', '').replace('.', '').replace('0', '')

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
            idx = [r['FORMAT'].split(':').index(x) for x in format.split(':')]
            infunc = lambda x: ':'.join([x.split(':')[i] for i in idx])
            r.iloc[9:] = r.iloc[9:].apply(infunc)
            return r
        df = self.df.copy()
        df[['ID', 'QUAL', 'FILTER', 'INFO']] = '.'
        df = df.apply(outfunc, axis=1)
        df['FORMAT'] = format
        vf = self.__class__([], df)
        return vf

    def merge(self, other, how='inner', format='GT', sort=True):
        """Merge with the other VcfFrame.

        Notably, this method is smart enough to automatically find and
        collapse duplicate records -- i.e. rows that have the identical
        values for CHROM, POS, and REF. For example, consider these two
        VCF records:

        +-------+-----+-----+-----+--------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven      |
        +=======+=====+=====+=====+========+=============+
        | chr1  | 100 | A   | T   | GT:AD  | 0/1:12,15   |
        +-------+-----+-----+-----+--------+-------------+

        +-------+-----+-----+-----+--------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Sara        |
        +=======+=====+=====+=====+========+=============+
        | chr1  | 100 | A   | T,C | GT:AD  | 1/2:0,14,15 |
        +-------+-----+-----+-----+--------+-------------+

        After merging, above records will be collapsed like this:

        +-------+-----+-----+-----+--------+---------------+-------------+
        | CHROM | POS | REF | ALT | FORMAT | Steven        | Sara        |
        +=======+=====+=====+=====+========+===============+=============+
        | chr1  | 100 | A   | C,T | GT:AD  | 0/2:12,0,15   | 1/2:0,15,14 |
        +-------+-----+-----+-----+--------+---------------+-------------+

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

        Returns
        -------
        vf : VcfFrame
            Merged VcfFrame.
        """
        vf1 = self.strip(format=format)
        vf2 = other.strip(format=format)
        dropped = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        shared = ['CHROM', 'POS', 'REF', 'ALT']
        df = vf1.df.merge(vf2.df.drop(columns=dropped),
            on=shared, how=how)
        df[dropped] = df[dropped].fillna('.')
        df.FORMAT = format

        # Collapse duplicate records.
        def collapse_one(df):
            ref_allele = df.REF.unique()[0]
            alt_alleles = []
            for i, r in df.iterrows():
                alt_alleles += r.ALT.split(',')
            alt_alleles = sorted(list(set(alt_alleles)))
            all_alleles = [ref_allele] + alt_alleles
            id_field = df.ID.unique()[0]
            qual_field = df.QUAL.unique()[0]
            filter_field = df.FILTER.unique()[0]
            info_field = df.INFO.unique()[0]
            format_field = df.FORMAT.unique()[0]

            def infunc(x, r_all_alleles, index_map):
                if pd.isna(x):
                    return x
                if not hasvar(x):
                    return x
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

            collapsed_df = df2 = df.apply(outfunc, axis=1)
            collapsed_df = collapsed_df.fillna('')
            collapsed_df = collapsed_df.groupby(['CHROM', 'POS', 'REF']).agg(''.join)
            collapsed_df = collapsed_df.reset_index()
            collapsed_df.ID = id_field
            collapsed_df.ALT = ','.join(alt_alleles)
            collapsed_df.QUAL = qual_field
            collapsed_df.FILTER = filter_field
            collapsed_df.INFO = info_field
            collapsed_df.FORMAT = format_field
            return collapsed_df.squeeze()

        i = df.duplicated(['CHROM', 'POS', 'REF'], keep=False)
        duplicates = {}
        for j, r in df[i].iterrows():
            name = f'{r.CHROM}:{r.POS}:{r.REF}'
            if name not in duplicates:
                duplicates[name] = []
            duplicates[name].append(j)
        for name in duplicates:
            df.iloc[duplicates[name]] = collapse_one(df.iloc[duplicates[name]])
        df.drop_duplicates(subset=['CHROM', 'POS', 'REF'], inplace=True)

        def func(r):
            n = len(r['FORMAT'].split(':'))
            x = './.'
            for i in range(1, n):
                x += ':.'
            r = r.fillna(x)
            return r
        df = df.apply(func, axis=1)
        vf3 = self.__class__([], df)
        if sort:
            vf3 = vf3.sort()
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
        df = self.df.apply(outfunc, axis=1)
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
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
            a = hasvar(r[n1])
            b = hasvar(r[n2])
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
            a = hasvar(r[n1])
            b = hasvar(r[n2])
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

    def subtract(self, name):
        """Remove rows that have a variant call in the sample.

        Parameters
        ----------
        name : str or int
            Name or index of the sample.

        Returns
        -------
        vf : VcfFrame
            Filtered VcfFrame.
        """
        name = name if isinstance(name, str) else self.samples[name]
        def func(r):
            return not hasvar(r[name])
        i = self.df.apply(func, axis=1)
        df = self.df[i].reset_index(drop=True)
        vf = self.__class__(deepcopy(self.meta), df)
        return vf

    def reset_samples(self, names):
        """Reset the sample list."""
        df = self.df[self.df.columns[:9].to_list() + names]
        vf = self.__class__(deepcopy(self.meta), df)
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
        vf = self.__class__(deepcopy(self.meta), df)
        return vf

    def parse_snpeff(self, idx, sep=' | '):
        """Parse SnpEff annotations.

        SnpEff provides the following functional annotations:

        1. Allele
        2. Annotation
        3. Annotation_Impact
        4. Gene_Name
        5. Gene_ID
        6. Feature_Type
        7. Feature_ID
        8. Transcript_BioType
        9. Rank
        10. HGVS.c
        11. HGVS.p
        12. cDNA.pos / cDNA.length
        13. CDS.pos / CDS.length
        14. AA.pos / AA.length
        15. Distance
        16. ERRORS / WARNINGS
        17. INFO

        Parameters
        ----------
        i : list
            List of annotation indicies.
        sep : str, default: ' | '
            Separator for joining requested annotations.

        Returns
        -------
        s : pandas.Series
            Parsed annotations.
        """
        def func(r):
            ann = [x for x in r['INFO'].split(';') if 'ANN=' in x]
            if not ann:
                return '.'
            ann = ann[0].replace('ANN=', '').split('|')
            ann = sep.join([ann[i] for i in idx])
            return ann
        s = self.df.apply(func, axis=1)
        return s

    def sort(self):
        """Return the sorted VcfFrame."""
        df = self.df.sort_values(by=['CHROM', 'POS'], ignore_index=True,
            key=lambda col: [CONTIGS.index(x) if isinstance(x, str)
                             else x for x in col])
        vf = self.__class__(deepcopy(self.meta), df)
        return vf
