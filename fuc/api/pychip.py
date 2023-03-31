"""
The pychip submodule is designed for working with annotation or manifest
files from the Axiom (Thermo Fisher Scientific) and Infinium (Illumina)
array platforms.
"""

import re
import pandas as pd

from . import common

class AxiomFrame:
    """
    Class for storing Axiom annotation data.

    Parameters
    ----------
    meta : list
        List of metadata lines.
    df : pandas.DataFrame
        DataFrame containing annotation data.
    """
    def __init__(self, meta, df):
        self._meta = meta
        self._df = df.reset_index(drop=True)

    @property
    def meta(self):
        """list : List of metadata lines."""
        return self._meta

    @meta.setter
    def meta(self, value):
        self._meta = value

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing annotation data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """
        Construct AxiomFrame from a CSV file.

        Parameters
        ----------
        fn : str
            CSV file (compressed or uncompressed).

        Returns
        -------
        AxiomFrame
            AxiomFrame object.
        """
        if fn.startswith('~'):
            fn = os.path.expanduser(fn)

        if fn.endswith('.gz'):
            f = gzip.open(fn, 'rt')
        else:
            f = open(fn)

        meta = []
        n = 0
        for line in f:
            if line.startswith('#'):
                meta.append(line)
                n += 1
        f.close()

        df = pd.read_csv(fn, skiprows=n)

        return cls(meta, df)

    def to_vep(self):
        """
        Convert AxiomFrame to the Ensembl VEP format.

        Returns
        -------
        pandas.DataFrame
            Variants in Ensembl VEP format.
        """
        df = self.df[self.df.Chromosome != '---']
        df['Physical Position'] = df['Physical Position'].astype(int)
        df['Position End'] = df['Position End'].astype(int)
        def one_row(r):
            result = []
            nucleotides = ['A', 'C', 'G', 'T']
            chrom = r['Chromosome']
            ref = r['Ref Allele']
            strand = r['Strand']
            start = r['Physical Position']
            end = r['Position End']
            for alt in r['Alt Allele'].split(' // '):
                if ref in nucleotides and alt in nucleotides: # SNV
                    pass
                elif alt == '-': # DEL I
                    pass
                elif len(alt) == len(ref): # MNV
                    pass
                elif len(alt) < len(ref) and ref.startswith(alt): # DEL II
                    start += len(alt)
                    ref = ref[len(alt):]
                    alt = '-'
                elif ref == '-': # INS I
                    start += 1
                    end = start - 1
                elif len(alt) > len(ref) and alt.startswith(ref): # INS II
                    diff = len(alt) - len(ref)
                    start += diff
                    end = start - 1
                    ref = '-'
                    alt = alt[diff:]
                else:
                    pass
                line = [chrom, start, end, f'{ref}/{alt}', strand]
                result.append('|'.join([str(x) for x in line]))
            return ','.join(result)
        s = df.apply(one_row, axis=1)
        s = ','.join(s)
        data = [x.split('|') for x in s.split(',')]
        df = pd.DataFrame(data).drop_duplicates()
        df.iloc[:, 1] = df.iloc[:, 1].astype(int)
        df.iloc[:, 2] = df.iloc[:, 2].astype(int)
        df = df.sort_values(by=[0, 1])
        return df

class InfiniumFrame:
    """
    Class for storing Infinium manifest data.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing manifest data.
    """
    def __init__(self, df):
        self._df = df.reset_index(drop=True)

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing manifest data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """
        Construct InfiniumFrame from a CSV file.

        Parameters
        ----------
        fn : str
            CSV file (compressed or uncompressed).

        Returns
        -------
        InfiniumFrame
            InfiniumFrame object.
        """
        if fn.startswith('~'):
            fn = os.path.expanduser(fn)

        if fn.endswith('.gz'):
            f = gzip.open(fn, 'rt')
        else:
            f = open(fn)

        lines = f.readlines()
        f.close()

        for i, line in enumerate(lines):
            if line.startswith('[Assay]'):
                start = i
                headers = lines[i+1].strip().split(',')
            elif line.startswith('[Controls]'):
                end = i

        lines = lines[start+2:end]
        lines = [x.strip().split(',') for x in lines]

        df = pd.DataFrame(lines, columns=headers)

        return cls(df)

    def to_vep(self, fasta):
        """
        Convert InfiniumFrame to the Ensembl VEP format.

        Returns
        -------
        pandas.DataFrame
            Variants in Ensembl VEP format.

        Parameters
        ----------
        region : str
            Region ('chrom:start-end').
        """
        nucleotides = ['A', 'C', 'G', 'T']
        n = 15
        df = self.df[(self.df.Chr != 'XY') & (self.df.Chr != '0')]
        df.MapInfo = df.MapInfo.astype(int)
        df = df.head(1000)
        def one_row(r):
            if r.Chr == 'MT':
                chrom = 'chrM'
            else:
                chrom = f'chr{r.Chr}'
            pos = r.MapInfo
            matches = re.findall(r'\[([^\]]+)\]', r.SourceSeq)
            if not matches:
                raise ValueError(f'Something went wrong: {r}')
            a1, a2 = matches[0].split('/')
            left_seq = r.SourceSeq.split('[')[0]
            right_seq = r.SourceSeq.split(']')[1]
            locus  = f'{chrom}:{pos}-{pos}'
            if a1 in nucleotides and a2 in nucleotides: # SNV
                ref = common.extract_sequence(fasta, locus)
                if ref == a1:
                    pass
                elif ref == a2:
                    a1, a2 = a2, a1
                else:
                    raise ValueError(f'Reference allele not found: {locus}')
            elif a1 == '-' or a2 == '-':
                affected = a1 if a2 == '-' else a2
                c = len(affected)
                if c == 1:
                    left_query = common.extract_sequence(fasta, f'{chrom}:{pos-n}-{pos-1}')
                    right_query = common.extract_sequence(fasta, f'{chrom}:{pos+1}-{pos+n}')
                    print(locus, left_seq[-n:], left_query, left_seq[-n:] == left_query, 'DEL I')
                    print(locus, right_seq[:n], right_query, right_seq[:n] == right_query, 'DEL I')
                    print()
                else:
                    left_query = common.extract_sequence(fasta, f'{chrom}:{pos-n}-{pos}') # -1
                    right_query = common.extract_sequence(fasta, f'{chrom}:{pos+c}-{pos+c+n-1}')
                    if left_seq[-n:] == left_query[:-1] and right_seq[:n] == right_query: # DEL I
                        print(locus, left_seq[-n:], left_query[:-1], left_seq[-n:] == left_query[:-1], 'DEL II')
                        print(locus, right_seq[:n], right_query, right_seq[:n] == right_query, 'DEL II')
                        print()
                    else:
                        #left_query = common.extract_sequence(fasta, f'{chrom}:{pos-n}-{pos}')
                        #left_query = left_query[len(affected):] + affected
                        print(locus, left_seq[-n:], left_query, left_seq[-n:] == left_query, 'Other')
                        print(locus, right_seq[:n], right_query, right_seq[:n] == right_query, 'Other')
                        print()
            else:
                raise ValueError(f'Unsupported format: {locus}')
            data = pd.Series([r.Chr, r.MapInfo, r.MapInfo, f'{a1}/{a2}', '+'])
            return data
        df = df.apply(one_row, axis=1)
        df = df.sort_values(by=[0, 1])
        df = df.drop_duplicates()
        return df
