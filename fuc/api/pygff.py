"""
The pygff submodule is designed for working with GFF (and GTF) files. It
implements ``pygff.GffFrame`` which stores GFF data as ``pandas.DataFrame``
to allow fast computation and easy manipulation. The submodule strictly
adheres to the standard `GFF specification
<https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`_.

A GFF file contains nine columns as follows:

+-----+------------+------------------------------------+------------------------+
| No. | Name       | Description                        | Examples               |
+=====+============+====================================+========================+
| 1   | Seqid      | Chromosome or contig identifier    | 'NC_000001.10'         |
+-----+------------+------------------------------------+------------------------+
| 2   | Source     | 1-based reference position         | 10041, 23042           |
+-----+------------+------------------------------------+------------------------+
| 3   | Type       | ';'-separated variant identifiers  | '.', 'rs35', 'rs9;rs53'|
+-----+------------+------------------------------------+------------------------+
| 4   | Start      | Reference allele                   | 'A', 'GT'              |
+-----+------------+------------------------------------+------------------------+
| 5   | End        | ','-separated alternate alleles    | 'T', 'ACT', 'C,T'      |
+-----+------------+------------------------------------+------------------------+
| 6   | Score      | Phred-scaled quality score for ALT | '.', 67, 12            |
+-----+------------+------------------------------------+------------------------+
| 7   | Strand     | ';'-separated filters that failed  | '.', 'PASS', 'q10;s50' |
+-----+------------+------------------------------------+------------------------+
| 8   | Phase      | ';'-separated information fields   | '.', 'DP=14;AF=0.5;DB' |
+-----+------------+------------------------------------+------------------------+
| 9   | Attributes | ':'-separated genotype fields      | 'GT', 'GT:AD:DP'       |
+-----+------------+------------------------------------+------------------------+
"""

import gzip

import pandas as pd

class GffFrame:
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
        """pandas.DataFrame : DataFrame containing GFF data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @classmethod
    def from_file(cls, fn):
        """
        Construct a GffFrame from a GFF file.

        Parameters
        ----------
        fn : str
            GFF file path (zipped or unzipped).

        Returns
        -------
        GffFrame
            GffFrame object.
        """
        meta = []
        skip_rows = 0

        if fn.startswith('~'):
            fn = os.path.expanduser(fn)

        if fn.endswith('.gz'):
            f = gzip.open(fn, 'rt')
        else:
            f = open(fn)

        for line in f:
            if line.startswith('#'):
                meta.append(line.strip())
                skip_rows += 1
            else:
                break

        f.close()

        names = [
            'Seqid', 'Source', 'Type', 'Start', 'End',
            'Score', 'Strand', 'Phase', 'Attributes'
        ]

        df = pd.read_table(fn, skiprows=skip_rows, names=names, header=None)

        return cls(meta, df)

    def protein_length(self, gene, name=None):
        temp = self.df[self.df.Type == 'CDS'].query(f"Attributes.str.contains('{gene}')")
        def one_row(r):
            l = r.Attributes.split(';')
            l = [x for x in l if x.startswith('Name=')]
            name = l[0].replace('Name=', '')
            return name
        temp['Protein_ID'] = temp.apply(one_row, axis=1)
        ids = temp.Protein_ID.unique()
        if name is None and len(ids) != 1:
            raise ValueError(f'Multiple protein sequences available: {ids}')
        elif name is None and len(ids) == 1:
            name = ids[0]
        temp = temp[temp.Protein_ID == name]
        length = (temp.End - temp.Start).sum() + temp.shape[0] - 3
        return int(length / 3)
