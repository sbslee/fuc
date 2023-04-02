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
        n = 20
        k = 100
        df = self.df[(self.df.Chr != 'XY') & (self.df.Chr != '0')]
        df.MapInfo = df.MapInfo.astype(int)
        df = df.head(5000)

        def compare_seq(seq1, seq2):
            if len(seq1) != len(seq2):
                return False
            for i in range(len(seq1)):
                if seq1[i] != seq2[i]:
                    if seq1[i] == 'N' or seq2[i] == 'N':
                        continue
                    else:
                        return False
            return True

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
            left_seq = r.SourceSeq.split('[')[0][-n:].upper()
            right_seq = r.SourceSeq.split(']')[1][:n].upper()
            locus  = f'{chrom}:{pos}-{pos}'
            if a1 in nucleotides and a2 in nucleotides: # SNV
                ref = common.extract_sequence(fasta, locus)
                if ref == a1:
                    pass
                elif ref == a2:
                    a1, a2 = a2, a1
                elif ref == common.reverse_complement(a1):
                    a1, a2 = ref, common.reverse_complement(a2)
                elif ref == common.reverse_complement(a2):
                    a1, a2 = ref, common.reverse_complement(a1)
                else:
                    raise ValueError(f'Reference allele not found: {[locus, a1, a2, ref]}')
            elif a1 == '-' or a2 == '-':
                affected = a1 if a2 == '-' else a2
                c = len(affected)
                ref_seq = common.extract_sequence(fasta, f'{chrom}:{pos-k}-{pos+k}')
                left_query = ref_seq[:k][-n:]
                right_query = ref_seq[k+c:][:n]
                if compare_seq(left_seq, left_query) and compare_seq(right_seq, right_query):
                    print(locus, left_seq, 'DEL I')
                else:
                    left_query = ref_seq[:k+c+1][-n:]
                    right_query = ref_seq[k+c+1:][:n]
                    if compare_seq(left_seq, left_query) and compare_seq(right_seq, right_query):
                        print(locus, left_seq, 'INS')
                    else:
                        left_query = ref_seq[:k][-n:]
                        right_query = ref_seq[k+c+1:][:n]
                        if left_query in r.SourceSeq.split('[')[0]:
                            repeats = r.SourceSeq.split('[')[0].split(left_query)[1]
                            left_query += repeats
                            left_query = left_query[len(repeats):]
                            right_query = ref_seq[k+c+len(repeats):][:n]
                            if compare_seq(left_seq, left_query) and compare_seq(right_seq, right_query):
                                print(locus, left_seq, 'DEL II')
                        else:
                            print(locus, left_seq, left_query, left_seq == left_query)
                            print(locus, right_seq, right_query, right_seq == right_query)
                            print()
            else:
                raise ValueError(f'Unsupported format: {locus}')
            data = pd.Series([r.Chr, r.MapInfo, r.MapInfo, f'{a1}/{a2}', '+'])
            return data
        df = df.apply(one_row, axis=1)
        df = df.sort_values(by=[0, 1])
        df = df.drop_duplicates()
        return df
