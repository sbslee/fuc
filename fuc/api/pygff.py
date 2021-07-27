"""
The pygff submodule is designed for working with GFF/GTF files. It implements
``pygff.GffFrame`` which stores GFF/GTF data as ``pandas.DataFrame`` to allow
fast computation and easy manipulation. The submodule strictly adheres to the
standard `GFF specification
<https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`_.

A GFF/GTF file contains nine columns as follows:

+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| No. | Name       | Description              | Examples                                                                                                              |
+=====+============+==========================+=======================================================================================================================+
| 1   | Seqid      | Landmark ID              | 'NC_000001.10', 'NC_012920.1'                                                                                         |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 2   | Source     | Feature source           | 'RefSeq', 'BestRefSeq', 'Genescan', 'Genebank'                                                                        |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 3   | Type       | Feature type             | 'transcript', 'exon', 'gene'                                                                                          |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 4   | Start      | Start coordinate         | 11874, 14409                                                                                                          |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 5   | End        | End coordinate           | 11874, 14409                                                                                                          |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 6   | Score      | Feature score            | '.', '1730.55', '1070'                                                                                                |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 7   | Strand     | Feature strand           | '.', '-', '+', '?'                                                                                                    |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 8   | Phase      | CDS phase                | '.', '0', '1', '2'                                                                                                    |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
| 9   | Attributes | ';'-separated attributes | 'ID=NC_000001.10:1..249250621;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA' |
+-----+------------+--------------------------+-----------------------------------------------------------------------------------------------------------------------+
"""

import gzip

import pandas as pd

class GffFrame:
    """
    Class for storing GFF/GTF data.

    Parameters
    ----------
    meta : list
        List of metadata lines.
    df : pandas.DataFrame
        DataFrame containing GFF/GTF data.
    fasta : str
        FASTA sequence lines.
    """
    def __init__(self, meta, df, fasta):
        self._meta = meta
        self._df = df.reset_index(drop=True)
        self._fasta = fasta

    @property
    def meta(self):
        """list : List of metadata lines."""
        return self._meta

    @meta.setter
    def meta(self, value):
        self._meta = value

    @property
    def df(self):
        """pandas.DataFrame : DataFrame containing GFF/GTF data."""
        return self._df

    @df.setter
    def df(self, value):
        self._df = value.reset_index(drop=True)

    @property
    def fasta(self):
        """dict : FASTA sequence lines."""
        return self._fasta

    @fasta.setter
    def fasta(self, value):
        self._fasta = value

    @classmethod
    def from_file(cls, fn):
        """
        Construct GffFrame from a GFF/GTF file.

        Parameters
        ----------
        fn : str
            GFF/GTF file (zipped or unzipped).

        Returns
        -------
        GffFrame
            GffFrame object.
        """
        if fn.startswith('~'):
            fn = os.path.expanduser(fn)

        if fn.endswith('.gz'):
            f = gzip.open(fn, 'rt')
        else:
            f = open(fn)

        meta = []
        data = []
        fasta = {}

        fasta_mode = False
        current_fasta = None

        for line in f:
            if '##FASTA' in line:
                fasta_mode = True
            elif fasta_mode and line.startswith('>'):
                name = line.strip().replace('>', '')
                fasta[name] = []
                current_fasta = name
            elif fasta_mode:
                fasta[current_fasta].append(line.strip())
            elif line.startswith('#'):
                meta.append(line.strip())
            else:
                fields = line.strip().split('\t')
                data.append(fields)

        f.close()

        columns = [
            'Seqid', 'Source', 'Type', 'Start', 'End',
            'Score', 'Strand', 'Phase', 'Attributes'
        ]

        df = pd.DataFrame(data, columns=columns)
        df.Start = df.Start.astype(int)
        df.End = df.End.astype(int)

        return cls(meta, df, fasta)

    def protein_length(self, gene, name=None):
        """
        Return the protein length of a gene.

        Parameters
        ----------
        gene : str
            Name of the gene.
        name : str, optional
            Protein sequence ID (e.g. 'NP_005219.2'). Required when the gene
            has multiple protein sequences available.

        Returns
        -------
        int
            Protein length.
        """
        df = self.df[self.df.Type == 'CDS'].query(
            f"Attributes.str.contains('{gene}')")
        def one_row(r):
            l = r.Attributes.split(';')
            l = [x for x in l if x.startswith('Name=')]
            name = l[0].replace('Name=', '')
            return name
        df['Protein_ID'] = df.apply(one_row, axis=1)
        ids = df.Protein_ID.unique()
        if name is None and len(ids) != 1:
            raise ValueError(f'Multiple protein sequences available: {ids}')
        elif name is None and len(ids) == 1:
            name = ids[0]
        df = df[df.Protein_ID == name]
        length = (df.End - df.Start).sum() + df.shape[0] - 3
        return int(length / 3)
