"""
The pymaf submodule is designed for working with MAF files. It implements
``pymaf.MafFrame`` which stores MAF data as ``pandas.DataFrame`` to allow
fast computation and easy manipulation. The ``pymaf.MafFrame`` class also
contains many useful plotting methods such as ``MafFrame.plot_oncoplot`` and
``MafFrame.plot_summary``. The submodule strictly adheres to the
standard `MAF specification
<https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_.

A typical MAF file contains many columns ranging from gene symbol to
protein change. However, most of the analysis in pymaf uses the
following columns:

+-----+------------------------+----------------------+-------------------------------+
| No. | Name                   | Description          | Examples                      |
+=====+========================+======================+===============================+
| 1   | Hugo_Symbol            | HUGO gene symbol     | 'TP53', 'Unknown'             |
+-----+------------------------+----------------------+-------------------------------+
| 2   | Chromosome             | Chromosome name      | 'chr1', '1', 'X'              |
+-----+------------------------+----------------------+-------------------------------+
| 3   | Start_Position         | Start coordinate     | 119031351                     |
+-----+------------------------+----------------------+-------------------------------+
| 4   | End_Position           | End coordinate       | 44079555                      |
+-----+------------------------+----------------------+-------------------------------+
| 5   | Variant_Classification | Translational effect | 'Missense_Mutation', 'Silent' |
+-----+------------------------+----------------------+-------------------------------+
| 6   | Variant_Type           | Mutation type        | 'SNP', 'INS', 'DEL'           |
+-----+------------------------+----------------------+-------------------------------+
| 7   | Reference_Allele       | Reference allele     | 'T', '-', 'ACAA'              |
+-----+------------------------+----------------------+-------------------------------+
| 8   | Tumor_Seq_Allele1      | First tumor allele   | 'A', '-', 'TCA'               |
+-----+------------------------+----------------------+-------------------------------+
| 9   | Tumor_Seq_Allele2      | Second tumor allele  | 'A', '-', 'TCA'               |
+-----+------------------------+----------------------+-------------------------------+
| 10  | Tumor_Sample_Barcode   | Sample ID            | 'TCGA-AB-3002'                |
+-----+------------------------+----------------------+-------------------------------+
| 11  | Protein_Change         | Protein change       | 'p.L558Q'                     |
+-----+------------------------+----------------------+-------------------------------+

It is also recommended to include additional custom columns such as variant
allele frequecy (VAF) and transcript name.

If sample annotation data are available for a given MAF file, use
the :class:`fuc.api.common.AnnFrame` class to import the data.

There are nine nonsynonymous variant classifcations that pymaf primarily
uses: Missense_Mutation, Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del,
In_Frame_Ins, Nonsense_Mutation, Nonstop_Mutation, Splice_Site, and
Translation_Start_Site.
"""

import re
import copy
import warnings
import itertools

from . import pyvcf, common

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec


CHROM_LENGTHS = {
    'hg18': [
        247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
        158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
        114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
        63811651, 62435964, 46944323, 49691432, 154913754, 57772954
    ],
    'hg19': [
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
        115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566, 155270560, 59373566
    ],
    'hg38': [
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
        159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
        58617616, 64444167, 46709983, 50818468, 156040895, 57227415
    ],
}

COMMON_COLUMNS = [
    'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome',
    'Start_Position', 'End_Position', 'Strand', 'Variant_Classification',
    'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1',
    'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'Protein_Change'
]

# Below is the list of calculated variant consequences from Ensembl VEP:
# https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
# (accessed on 2021-05-31)
#
# Note that both frameshift_variant and protein_altering_variant require
# additional information to find their correct Variant_Classification.

VEP_CONSEQUENCES = {
    'transcript_ablation':                'Splice_Site',
    'splice_acceptor_variant':            'Splice_Site',
    'splice_donor_variant':               'Splice_Site',
    'stop_gained':                        'Nonsense_Mutation',
    'frameshift_variant':                 'AMBIGUOUS',
    'stop_lost':                          'Nonstop_Mutation',
    'start_lost':                         'Translation_Start_Site',
    'transcript_amplification':           'Intron',
    'inframe_insertion':                  'In_Frame_Ins',
    'inframe_deletion':                   'In_Frame_Del',
    'missense_variant':                   'Missense_Mutation',
    'protein_altering_variant':           'AMBIGUOUS',
    'splice_region_variant':              'Splice_Region',
    'incomplete_terminal_codon_variant':  'Silent',
    'start_retained_variant':             'Silent',
    'stop_retained_variant':              'Silent',
    'synonymous_variant':                 'Silent',
    'coding_sequence_variant':            'Missense_Mutation',
    'mature_miRNA_variant':               'RNA',
    '5_prime_UTR_variant':                "5'UTR",
    '3_prime_UTR_variant':                "3'UTR",
    'non_coding_transcript_exon_variant': 'RNA',
    'intron_variant':                     'Intron',
    'NMD_transcript_variant':             'Silent',
    'non_coding_transcript_variant':      'RNA',
    'upstream_gene_variant':              "5'Flank",
    'downstream_gene_variant':            "3'Flank",
    'TFBS_ablation':                      'Targeted_Region',
    'TFBS_amplification':                 'Targeted_Region',
    'TF_binding_site_variant':            'IGR',
    'regulatory_region_ablation':         'Targeted_Region',
    'regulatory_region_amplification':    'Targeted_Region',
    'feature_elongation':                 'Targeted_Region',
    'regulatory_region_variant':          'IGR',
    'feature_truncation':                 'Targeted_Region',
    'intergenic_variant':                 'IGR',
}

VARCLS_LIST = [
    "3'Flank",
    "3'UTR",
    "5'Flank",
    "5'UTR",
    'De_novo_Start_InFrame',
    'De_novo_Start_OutOfFrame',
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'IGR',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Intron',
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Nonstop_Mutation',
    'RNA',
    'Silent',
    'Splice_Region',
    'Splice_Site',
    'Start_Codon_Ins',
    'Start_Codon_SNP',
    'Stop_Codon_Del',
    'Targeted_Region',
    'Translation_Start_Site',
    'lincRNA',
]

NONSYN_NAMES = [
    'Missense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
    'In_Frame_Del', 'In_Frame_Ins', 'Nonsense_Mutation',
    'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site'
]

NONSYN_COLORS = [
    'tab:green', 'tab:blue', 'tab:purple', 'tab:olive', 'tab:red',
    'tab:cyan', 'tab:pink', 'tab:orange', 'tab:brown'
]

SNV_CLASSES = {
    'A>C': {'class': 'T>G', 'type': 'Tv'},
    'A>G': {'class': 'T>C', 'type': 'Ti'},
    'A>T': {'class': 'T>A', 'type': 'Tv'},
    'C>A': {'class': 'C>A', 'type': 'Tv'},
    'C>G': {'class': 'C>G', 'type': 'Tv'},
    'C>T': {'class': 'C>T', 'type': 'Ti'},
    'G>A': {'class': 'C>T', 'type': 'Ti'},
    'G>C': {'class': 'C>G', 'type': 'Tv'},
    'G>T': {'class': 'C>A', 'type': 'Tv'},
    'T>A': {'class': 'T>A', 'type': 'Tv'},
    'T>C': {'class': 'T>C', 'type': 'Ti'},
    'T>G': {'class': 'T>G', 'type': 'Tv'},
}

SNV_CLASS_ORDER = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

class MafFrame:
    """Class for storing MAF data.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing MAF data.

    See Also
    --------
    MafFrame.from_file
        Construct MafFrame from a MAF file.
    """
    def __init__(self, df):
        self.df = df.reset_index(drop=True)

    @property
    def shape(self):
        """tuple : Dimensionality of MafFrame (variants, samples)."""
        return (self.df.shape[0], len(self.samples))

    @property
    def samples(self):
        """list : List of the sample names."""
        return list(self.df.Tumor_Sample_Barcode.unique())

    @property
    def genes(self):
        """list : List of the genes."""
        return list(self.df.Hugo_Symbol.unique())

    def compute_clonality(self, vaf_col, threshold=0.25):
        """
        Compute the clonality of variants based on
        :ref:`VAF <glossary:Variant allele frequency (VAF)>`.

        A mutation will be defined as "Subclonal" if the VAF is less than the
        threshold percentage (e.g. 25%) of the highest VAF in the sample and
        is defined as "Clonal" if it is equal to or above this threshold.

        Parameters
        ----------
        vaf_col : str
            MafFrame column containing VAF data.
        threshold : float
            Minimum VAF to be considered as "Clonal".

        Returns
        -------
        panda.Series
            Clonality for each variant.

        Examples
        --------

        >>> import matplotlib.pyplot as plt
        >>> from fuc import common, pymaf
        >>> common.load_dataset('tcga-laml')
        >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
        >>> mf = pymaf.MafFrame.from_file(maf_file)
        >>> mf.df['Clonality'] = mf.compute_clonality('i_TumorVAF_WU')
        >>> mf.df['Clonality'][:10]
        0    Clonal
        1    Clonal
        2    Clonal
        3    Clonal
        4    Clonal
        5    Clonal
        6    Clonal
        7    Clonal
        8    Clonal
        9    Clonal
        Name: Clonality, dtype: object
        """
        d = self.df.groupby('Tumor_Sample_Barcode')[vaf_col].max().to_dict()
        def one_row(r):
            m = d[r.Tumor_Sample_Barcode]
            if r[vaf_col] < m * threshold:
                result = 'Subclonal'
            else:
                result = 'Clonal'
            return result
        s = self.df.copy().apply(one_row, axis=1)
        return s

    @classmethod
    def from_file(cls, fn):
        """
        Construct MafFrame from a MAF file.

        Parameters
        ----------
        fn : str
            MAF file (zipped or unzipped).

        Returns
        -------
        MafFrame
            MafFrame object.

        See Also
        --------
        MafFrame
            MafFrame object creation using constructor.

        Examples
        --------

        >>> from fuc import common, pymaf
        >>> common.load_dataset('tcga-laml')
        >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
        >>> mf = pymaf.MafFrame.from_file(maf_file)
        """
        # Read the input MAF file.
        df = pd.read_table(fn)

        # Check the letter case of column names. This will help distinguish
        # missing columns from columns with incorrect letter case (e.g.
        # 'End_Position' vs. 'End_position').
        lower_names = [x.lower() for x in COMMON_COLUMNS]
        for col in df.columns:
            if col.lower() in lower_names:
                i = lower_names.index(col.lower())
                if col != COMMON_COLUMNS[i]:
                    message = (
                        f"Input column '{col}' will be renamed "
                        f"as '{COMMON_COLUMNS[i]}'."
                    )
                    warnings.warn(message)
                    df = df.rename(columns={col: COMMON_COLUMNS[i]})

        # Set the data type of chromosomes as string (e.g. 'chr1' vs. '1').
        if 'Chromosome' in df.columns:
            df.Chromosome = df.Chromosome.astype(str)

        return cls(df)

    @classmethod
    def from_vcf(cls, vcf, keys=None, names=None):
        """
        Construct MafFrame from a VCF file or VcfFrame.

        It is recommended that the input VCF data be functionally annotated
        by an annotation tool such as Ensembl VEP, SnpEff, and ANNOVAR;
        however, the method can handle unannotated VCF data as well.

        The preferred tool for functional annotation is Ensembl VEP with
        "RefSeq transcripts" as the transcript database and the filtering
        option "Show one selected consequence per variant".

        Parameters
        ----------
        vcf : str or VcfFrame
            VCF file or VcfFrame.
        keys : str or list
            Genotype key (e.g. 'AD', 'AF') or list of genotype keys to be
            added to the MafFrame.
        names : str or list
            Column name or list of column names for ``keys`` (must be the
            same length). By default, the genotype keys will be used as
            column names.

        Examples
        --------
        Below is a simple example:

        >>> from fuc import pyvcf, pymaf
        >>> data = {
        ...     'CHROM': ['chr1', 'chr2'],
        ...     'POS': [100, 101],
        ...     'ID': ['.', '.'],
        ...     'REF': ['G', 'T'],
        ...     'ALT': ['A', 'C'],
        ...     'QUAL': ['.', '.'],
        ...     'FILTER': ['.', '.'],
        ...     'INFO': ['CSQ=T|missense_variant|MODERATE|MTOR|2475|Transcript|NM_001386500.1|protein_coding|47/58||||6792|6644|2215|S/Y|tCt/tAt|rs587777894&COSV63868278&COSV63868313||-1||EntrezGene||||||||G|G||deleterious(0)|possibly_damaging(0.876)||||||||||||||||||likely_pathogenic&pathogenic|0&1&1|1&1&1|26619011&27159400&24631838&26018084&27830187|||||', 'CSQ=C|splice_donor_variant|HIGH|MTOR|2475|Transcript|NM_001386500.1|protein_coding||46/57||||||||||-1||EntrezGene||||||||A|A|||||||||||||||||||||||||||||'],
        ...     'FORMAT': ['GT:AD:DP:AF', 'GT:AD:DP:AF'],
        ...     'A': ['0/1:176,37:213:0.174', '0/1:966,98:1064:0.092']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID REF ALT QUAL FILTER                                               INFO       FORMAT                      A
        0  chr1  100  .   G   A    .      .  CSQ=T|missense_variant|MODERATE|MTOR|2475|Tran...  GT:AD:DP:AF   0/1:176,37:213:0.174
        1  chr2  101  .   T   C    .      .  CSQ=C|splice_donor_variant|HIGH|MTOR|2475|Tran...  GT:AD:DP:AF  0/1:966,98:1064:0.092
        >>> mf = pymaf.MafFrame.from_vcf(vf)
        >>> mf.df
          Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome  Start_Position  End_Position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Protein_Change Tumor_Sample_Barcode
        0        MTOR           2475      .          .       chr1             100           100      -      Missense_Mutation          SNP                G                 A                 A       p.S2215Y                    A
        1        MTOR           2475      .          .       chr2             101           101      -            Splice_Site          SNP                T                 C                 C              .                    A

        We can add genotype keys such as AD and AF:

        >>> mf = pymaf.MafFrame.from_vcf(vf, keys=['AD', 'AF'])
        >>> mf.df
          Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome  Start_Position  End_Position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Protein_Change Tumor_Sample_Barcode      AD     AF
        0        MTOR           2475      .          .       chr1             100           100      -      Missense_Mutation          SNP                G                 A                 A       p.S2215Y                    A  176,37  0.174
        1        MTOR           2475      .          .       chr2             101           101      -            Splice_Site          SNP                T                 C                 C              .                    A  966,98  0.092

        The method can accept a VCF file as input instead of VcfFrame:

        >>> mf = pymaf.MafFrame.from_vcf('annotated.vcf')

        The method can also handle unannotated VCF data:

        >>> data = {
        ...     'CHROM': ['chr1', 'chr1', 'chr1'],
        ...     'POS': [100, 200, 300],
        ...     'ID': ['.', '.', '.'],
        ...     'REF': ['G', 'C', 'TTC'],
        ...     'ALT': ['A', 'CAG', 'T'],
        ...     'QUAL': ['.', '.', '.'],
        ...     'FILTER': ['.', '.', '.'],
        ...     'INFO': ['.', '.', '.'],
        ...     'FORMAT': ['GT', 'GT', 'GT'],
        ...     'A': ['0/1', '0/1', '0/1']
        ... }
        >>> vf = pyvcf.VcfFrame.from_dict([], data)
        >>> vf.df
          CHROM  POS ID  REF  ALT QUAL FILTER INFO FORMAT    A
        0  chr1  100  .    G    A    .      .    .     GT  0/1
        1  chr1  200  .    C  CAG    .      .    .     GT  0/1
        2  chr1  300  .  TTC    T    .      .    .     GT  0/1
        >>> mf = pymaf.MafFrame.from_vcf(vf)
        >>> mf.df
          Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome  Start_Position  End_Position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Protein_Change Tumor_Sample_Barcode
        0           .              .      .          .       chr1             100           100      .                      .          SNP                G                 A                 A              .                    A
        1           .              .      .          .       chr1             200           201      .                      .          INS                -                AG                AG              .                    A
        2           .              .      .          .       chr1             301           302      .                      .          DEL               TC                 -                 -              .                    A
        """
        # Parse the input VCF.
        if isinstance(vcf, str):
            vf = pyvcf.VcfFrame.from_file(vcf)
        else:
            vf = vcf

        # Set some default values in case the VCF is not annotated.
        ncbi_build = '.'

        # Get the NCBI_Build data.
        for line in vf.meta:
            if line.startswith('##VEP'):
                ncbi_build = re.search(r'assembly="(.*?)"', line).group(1)
                break

        # Define the conversion algorithm.
        def one_row(r):
            has_annot = 'CSQ=' in r.INFO

            # Set some default values in case the VCF is not annotated.
            strand = '.'
            variant_classification = '.'
            protein_change = '.'
            hugo_symbol = '.'
            entrez_gene_id = '.'

            # Get the sequence data.
            inframe = abs(len(r.REF) - len(r.ALT)) / 3 == 0
            if len(r.REF) == len(r.ALT) == 1:
                variant_type = 'SNP'
                start_position = r.POS
                end_position = r.POS
                reference_allele = r.REF
                tumor_seq_allele1 = r.ALT
                tumor_seq_allele2 = r.ALT
            elif len(r.REF) > len(r.ALT):
                variant_type = 'DEL'
                start_position = r.POS + 1
                end_position = r.POS + len(r.REF) - len(r.ALT)
                reference_allele = r.REF[1:]
                tumor_seq_allele1 = '-'
                tumor_seq_allele2 = '-'
            else:
                variant_type = 'INS'
                start_position = r.POS
                end_position = r.POS + 1
                reference_allele = '-'
                tumor_seq_allele1 = r.ALT[1:]
                tumor_seq_allele2 = r.ALT[1:]

            fields = r.INFO.replace('CSQ=', '').split(',')[0].split('|')

            # Get the Strand data.
            if has_annot:
                strand = '+' if fields[19] == '1' else '-'

            # Get the Variant_Classification data.
            if has_annot:
                consequence = fields[1].split('&')[0]
                if consequence == 'frameshift_variant':
                    if variant_type == 'DEL':
                        variant_classification = 'Frame_Shift_Del'
                    else:
                        variant_classification = 'Frame_Shift_Ins'
                elif consequence == 'protein_altering_variant':
                    if inframe:
                        if variant_type == 'DEL':
                            variant_classification = 'In_Frame_Del'
                        else:
                            variant_classification = 'In_Frame_Ins'
                    else:
                        if variant_type == 'DEL':
                            variant_classification = 'Frame_Shift_Del'
                        else:
                            variant_classification = 'Frame_Shift_Ins'
                elif consequence in VEP_CONSEQUENCES:
                    variant_classification = VEP_CONSEQUENCES[consequence]
                else:
                    m = f'Found unknown Ensembl VEP consequence: {consequence}'
                    raise ValueError(m)

            # Get the Tumor_Sample_Barcode data.
            s = r[9:].apply(pyvcf.gt_hasvar)
            tumor_sample_barcode = ','.join(s[s].index.to_list())

            # Get the Protein_Change data.
            if has_annot:
                pos = fields[14]
                aa = fields[15].split('/')
                if len(aa) > 1:
                    protein_change = f'p.{aa[0]}{pos}{aa[1]}'

            # Get other data.
            if has_annot:
                hugo_symbol = fields[3]
                entrez_gene_id = fields[4]

            d = dict(
                Hugo_Symbol = hugo_symbol,
                Entrez_Gene_Id = entrez_gene_id,
                Center = '.',
                NCBI_Build = ncbi_build,
                Chromosome = r.CHROM,
                Start_Position = start_position,
                End_Position = end_position,
                Strand = strand,
                Variant_Classification = variant_classification,
                Variant_Type = variant_type,
                Reference_Allele = reference_allele,
                Tumor_Seq_Allele1 = tumor_seq_allele1,
                Tumor_Seq_Allele2 = tumor_seq_allele2,
                Tumor_Sample_Barcode = tumor_sample_barcode,
                Protein_Change = protein_change,
                CHROM = r.CHROM, # will be dropped
                POS = r.POS,     # will be dropped
                REF = r.REF,     # will be dropped
                ALT = r.ALT,     # will be dropped
            )

            return pd.Series(d)

        # Apply the conversion algorithm.
        df = vf.df.apply(one_row, axis=1)

        # Expand the Tumor_Sample_Barcode column to multiple rows.
        s = df['Tumor_Sample_Barcode'].str.split(',').apply(
            pd.Series, 1).stack()
        s.index = s.index.droplevel(-1)
        s.name = 'Tumor_Sample_Barcode'
        del df['Tumor_Sample_Barcode']
        df = df.join(s)

        # Append extra genotype keys, if necessary.
        if keys is not None:
            if names is None:
                names = keys
            if isinstance(keys, str):
                keys = [keys]
            if isinstance(names, str):
                names = [names]
            for i, key in enumerate(keys):
                temp_df = vf.extract(key)
                temp_df = pd.concat([vf.df.iloc[:, :9], temp_df], axis=1)
                temp_df = temp_df.drop(
                    columns=['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
                temp_df = pd.melt(
                    temp_df,
                    id_vars=['CHROM', 'POS', 'REF', 'ALT'],
                    var_name='Tumor_Sample_Barcode',
                )
                temp_df = temp_df[temp_df.value != '.']
                df = df.merge(temp_df,
                    on=['CHROM', 'POS', 'REF', 'ALT', 'Tumor_Sample_Barcode'])
                df = df.rename(columns={'value': names[i]})

        # Drop the extra columns.
        df = df.drop(columns=['CHROM', 'POS', 'REF', 'ALT'])

        return cls(df)

    def matrix_prevalence(self):
        """
        Compute a matrix of variant counts with a shape of (genes, samples).

        Returns
        -------
        pandas.DataFrame
            The said matrix.
        """
        s = self.df.groupby(
            'Hugo_Symbol')['Tumor_Sample_Barcode'].value_counts()
        s.name = 'Count'
        df = s.to_frame().reset_index()
        df = df.pivot(index='Hugo_Symbol',
            columns='Tumor_Sample_Barcode', values='Count')
        df.columns.name = ''
        df = df.fillna(0)
        return df

    def matrix_genes(self, mode='variants', count=10):
        """
        Compute a matrix of counts with a shape of (genes, variant
        classifications).

        This method only considers the nine nonsynonymous variant
        classifications.

        Parameters
        ----------
        mode : {'variants', 'samples'}, default: 'variants'
            Determines how to identify top mutated genes:

            * 'variants': Count the number of observed variants.
            * 'samples': Count the number of affected samples. Using this
              option will create an additional variant classification called
              'Multi_Hit'.

        count : int, default: 10
            Number of top mutated genes to include.

        Returns
        -------
        pandas.DataFrame
            The said matrix.
        """
        if mode == 'variants':
            df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]
            df = df.groupby('Hugo_Symbol')[
                'Variant_Classification'].value_counts().to_frame()
            df.columns = ['Count']
            df = df.reset_index()
            df = df.pivot(index='Hugo_Symbol', columns='Variant_Classification',
                values='Count')
            df = df.fillna(0)
            for varcls in NONSYN_NAMES:
                if varcls not in df.columns:
                    df[varcls] = 0
            i = df.sum(axis=1).sort_values(ascending=False).index
            df = df.reindex(index=i)
            df = df[NONSYN_NAMES]
            df = df[:count]
            df = df.rename_axis(None, axis=1)
        elif mode == 'samples':
            df = self.matrix_waterfall(count)
            df = df.apply(lambda r: r.value_counts(), axis=1)
            for varcls in NONSYN_NAMES + ['Multi_Hit']:
                if varcls not in df.columns:
                    df[varcls] = np.nan
            df = df[NONSYN_NAMES + ['Multi_Hit']]
            df = df.fillna(0)
        else:
            raise ValueError(f'Found incorrect mode: {mode}')
        return df

    def matrix_tmb(self):
        """
        Compute a matrix of variant counts with a shape of (samples, variant
        classifications).

        Returns
        -------
        pandas.DataFrame
            The said matrix.
        """
        df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]
        df = df.groupby('Tumor_Sample_Barcode')[
            'Variant_Classification'].value_counts().to_frame()
        df.columns = ['Count']
        df = df.reset_index()
        df = df.pivot(index='Tumor_Sample_Barcode',
            columns='Variant_Classification', values='Count')
        df = df.fillna(0)
        for varcls in NONSYN_NAMES:
            if varcls not in df.columns:
                df[varcls] = 0
        i = df.sum(axis=1).sort_values(ascending=False).index
        df = df.reindex(index=i)
        df = df[NONSYN_NAMES]
        df = df.rename_axis(None, axis=1)
        return df

    def matrix_waterfall(self, count=10, keep_empty=False):
        """
        Compute a matrix of variant classifications with a shape of
        (genes, samples).

        If there are multiple variant classifications available for a given
        cell, they will be replaced as 'Multi_Hit'.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to include.
        keep_empty : bool, default: False
            If True, keep samples with all ``NaN``'s.

        Returns
        -------
        pandas.DataFrame
            The said matrix.
        """
        df = self.df[self.df.Variant_Classification.isin(NONSYN_NAMES)]

        f = lambda x: ''.join(x) if len(x) == 1 else 'Multi_Hit'
        df = df.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode'])[
            'Variant_Classification'].apply(f).to_frame()
        df = df.reset_index()
        df = df.pivot(index='Hugo_Symbol', columns='Tumor_Sample_Barcode',
            values='Variant_Classification')

        # Sort the rows (genes).
        i = df.isnull().sum(axis=1).sort_values(ascending=True).index
        df = df.reindex(index=i)

        # Select the top mutated genes.
        df = df[:count]

        # Drop samples with all NaN's.
        if not keep_empty:
            df = df.dropna(axis=1, how='all')

        # Sort the columns (samples).
        c = df.applymap(lambda x: 0 if pd.isnull(x) else 1).sort_values(
            df.index.to_list(), axis=1, ascending=False).columns
        df = df[c]
        df = df.fillna('None')
        df = df.rename_axis(None, axis=1)

        return df

    def plot_genes(
        self, mode='variants', count=10, flip=False, ax=None, figsize=None,
        **kwargs
    ):
        """
        Create a bar plot showing variant distirbution for top mutated genes.

        Parameters
        ----------
        mode : {'variants', 'samples'}, default: 'variants'
            Determines how to identify top mutated genes:

            * 'variants': Count the number of observed variants.
            * 'samples': Count the number of affected samples. Using this
              option will create an additional variant classification called
              'Multi_Hit'.
        count : int, default: 10
            Number of top mutated genes to display.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.bar` or
            :meth:`pandas.DataFrame.plot.barh`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        By default (``mode='variants'``), the method identifies top mutated
        genes by counting the number of observed variants:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_genes()
            >>> plt.tight_layout()

        We can also identify top mutated genes by counting the number of
        affected samples:

        .. plot::
            :context: close-figs

            >>> mf.plot_genes(mode='samples')
            >>> plt.tight_layout()
        """
        if mode == 'variants':
            colors = NONSYN_COLORS
        elif mode == 'samples':
            colors = NONSYN_COLORS + ['k']
        else:
            raise ValueError(f'Found incorrect mode: {mode}')
        df = self.matrix_genes(count=count, mode=mode)
        df = df.iloc[::-1]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            df = df.iloc[::-1]
            kind = 'bar'
            xlabel, ylabel = '', 'Count'
        else:
            kind = 'barh'
            xlabel, ylabel = 'Count', ''

        df.plot(
            kind=kind, ax=ax, stacked=True, legend=False,
            color=colors, **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_oncoplot(
        self, count=10, keep_empty=False, figsize=(15, 10), label_fontsize=15,
        ticklabels_fontsize=15, legend_fontsize=15
    ):
        """
        Create an oncoplot.

        See this :ref:`tutorial <tutorials:Create customized oncoplots>` to
        learn how to create customized oncoplots.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.
        keep_empty : bool, default: False
            If True, display samples that do not have any mutations.
        figsize : tuple, default: (15, 10)
            Width, height in inches. Format: (float, float).
        label_fontsize : float, default: 15
            Font size of labels.
        ticklabels_fontsize : float, default: 15
            Font size of tick labels.
        legend_fontsize : float, default: 15
            Font size of legend texts.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_oncoplot()
        """
        g = {'height_ratios': [1, 10, 1], 'width_ratios': [10, 1]}
        fig, axes = plt.subplots(3, 2, figsize=figsize, gridspec_kw=g)
        [[ax1, ax2], [ax3, ax4], [ax5, ax6]] = axes

        # Create the TMB plot.
        samples = list(self.matrix_waterfall(count=count,
            keep_empty=keep_empty).columns)
        self.plot_tmb(ax=ax1, samples=samples, width=0.95)
        ax1.set_xlabel('')
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.set_xlim(-0.5, len(samples)-0.5)
        ax1.set_ylabel('TMB', fontsize=label_fontsize)
        ax1.set_yticks([0, self.matrix_tmb().sum(axis=1).max()])
        ax1.tick_params(axis='y', which='major',
                        labelsize=ticklabels_fontsize)

        # Remove the top right plot.
        ax2.remove()

        # Create the waterfall plot.
        self.plot_waterfall(count=count, ax=ax3, linewidths=1, keep_empty=keep_empty)
        ax3.set_xlabel('')
        ax3.tick_params(axis='y', which='major', labelrotation=0,
                        labelsize=ticklabels_fontsize)

        # Create the genes plot.
        self.plot_genes(count=count, ax=ax4, mode='samples', width=0.95)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax4.set_yticks([])
        ax4.set_xlabel('Samples', fontsize=label_fontsize)
        ax4.set_xticks([0, self.matrix_genes(
            count=10, mode='samples').sum(axis=1).max()])
        ax4.set_ylim(-0.5, count-0.5)
        ax4.tick_params(axis='x', which='major',
                        labelsize=ticklabels_fontsize)

        # Create the legend.
        handles = common.legend_handles(NONSYN_NAMES+['Multi_Hit'],
            colors=NONSYN_COLORS+['k'])
        ax5.legend(
            handles=handles,
            title='Variant_Classification',
            loc='upper center',
            ncol=4,
            fontsize=legend_fontsize,
            title_fontsize=legend_fontsize
        )
        ax5.axis('off')

        # Remove the bottom right plot.
        ax6.remove()

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.01, hspace=0.01)

    def plot_oncoplot_matched(
        self, af, patient_col, group_col, group_order, colors='Set2',
        figsize=(15, 10), label_fontsize=12, ticklabels_fontsize=12,
        legend_fontsize=12
    ):
        """
        Create an oncoplot for mached samples.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        colors : str
            Colormap name for the sample groups.
        figsize : tuple, default: (15, 10)
            Width, height in inches. Format: (float, float).
        label_fontsize : float, default: 12
            Font size of labels.
        ticklabels_fontsize : float, default: 12
            Font size of tick labels.
        legend_fontsize : float, default: 12
            Font size of legend texts.
        """
        fig, axes = plt.subplots(3, 2, figsize=figsize,
            gridspec_kw={'height_ratios': [1, 10, 1.5], 'width_ratios': [10, 1]}
        )

        [[ax1, ax2], [ax3, ax4], [ax5, ax6]] = axes

        patients = self.matrix_waterfall_matched(af, patient_col, group_col, group_order).columns

        self.plot_tmb_matched(
            af, patient_col, group_col, group_order=group_order, ax=ax1,
            legend=False, patients=patients, width=0.90,
            color=sns.color_palette(colors)[:3]
        )
        ax1.set_xticks([])
        ax1.set_xlim(-0.5, 53-0.5)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.set_ylabel('TMB', fontsize=label_fontsize)
        ax1.tick_params(axis='y', which='major',
                        labelsize=ticklabels_fontsize)

        ax2.remove()

        self.plot_waterfall_matched(af, patient_col, group_col, group_order=group_order, ax=ax3)
        ax3.set_xticks([])
        ax3.tick_params(axis='y', which='major', labelrotation=0,
                        labelsize=ticklabels_fontsize)

        self.plot_mutated_matched(
            af, patient_col, group_col, group_order=group_order, ax=ax4, palette=colors
        )
        ax4.set_yticks([])
        ax4.legend().remove()
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax4.tick_params(axis='x', which='major',
                        labelsize=ticklabels_fontsize)
        ax4.set_xlabel('Patients', fontsize=label_fontsize)

        # Create the legends.
        handles1 = common.legend_handles(NONSYN_NAMES+['Multi_Hit'],
            colors=NONSYN_COLORS+['k'])
        handles2 = common.legend_handles(group_order, colors=colors)
        leg1 = ax5.legend(handles=handles1, loc=(0, 0), title='Variant_Classification', ncol=4, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
        leg2 = ax5.legend(handles=handles2, loc=(0.8, 0), title=group_col, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
        ax5.add_artist(leg1)
        ax5.add_artist(leg2)
        ax5.axis('off')

        # Remove the bottom right plot.
        ax6.remove()

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.01, hspace=0.01)

    def plot_clonality(
        self, vaf_col, af=None, group_col=None, group_order=None, count=10,
        threshold=0.25, subclonal=False, ax=None, figsize=None
    ):
        """
        Create a bar plot summarizing the clonality of variants in top
        mutated genes.

        Clonality will be calculated based on VAF using
        :meth:`MafFrame.compute_clonality`.

        Parameters
        ----------
        vaf_col : str
            MafFrame column containing VAF data.
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        count : int, defualt: 10
            Number of top mutated genes to display.
        threshold : float, default: 0.25
            VAF threshold percentage.
        subclonal : bool, default: False
            If True, display subclonality (1 - clonality).
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        MafFrame.compute_clonality
            Compute the clonality of variants based on VAF.

        Examples
        --------

        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_clonality('i_TumorVAF_WU')
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_clonality('i_TumorVAF_WU',
            ...                   af=af,
            ...                   group_col='FAB_classification',
            ...                   group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        df = self.df.copy()
        df['Clonality'] = self.compute_clonality(vaf_col, threshold=threshold)

        if group_col is None:
            s = df.groupby('Hugo_Symbol')['Clonality'].value_counts()
            s.name = 'Count'
            df = s.to_frame().reset_index()
            df = df.pivot(index='Hugo_Symbol', columns='Clonality', values='Count')
        else:
            df = df.merge(af.df[group_col], left_on='Tumor_Sample_Barcode', right_index=True)
            s = df.groupby(['Hugo_Symbol', group_col])['Clonality'].value_counts()
            s.name = 'Count'
            df = s.to_frame().reset_index()
            df = df.pivot(index=['Hugo_Symbol', group_col], columns='Clonality', values='Count')

        df = df.reset_index()
        df = df.fillna(0)
        l = ['Clonal', 'Subclonal']
        df[l] = df[l].div(df[l].sum(axis=1), axis=0)
        genes = self.matrix_genes(count=count).index

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if subclonal:
            y = 'Subclonal'
        else:
            y = 'Clonal'

        sns.barplot(
            x='Hugo_Symbol', y=y, data=df, order=genes, hue=group_col,
            hue_order=group_order, ax=ax
        )

        ax.set_xlabel('')

        return ax

    def plot_evolution(
        self, samples, vaf_col, anchor=None, normalize=True, count=5,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a line plot visualizing changes in VAF between specified
        samples.

        Parameters
        ----------
        samples : list
            List of samples to display.
        vaf_col : str
            MafFrame column containing VAF data.
        anchor : str, optional
            Sample to use as the anchor. If absent, use the first sample in
            the list.
        normalize : bool, default: True
            If False, do not normalize VAF by the maximum value.
        count : int, default: 5
            Number of top variants to display.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.lineplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        df = self.df[self.df.Tumor_Sample_Barcode.isin(samples)]

        if df.empty:
            message = f'No variants to display for the samples: {samples}.'
            raise ValueError(message)

        df = df[df.Variant_Classification.isin(NONSYN_NAMES)]

        def one_row(r):
            if r.Protein_Change == '.':
                variant_name = f'{r.Hugo_Symbol} ({r.Variant_Classification})'
            else:
                variant_name = f'{r.Hugo_Symbol} ({r.Protein_Change})'
            return variant_name

        df['Variant_Name'] = df.apply(one_row, axis=1)
        df = df.pivot(index=['Variant_Name'],
            columns=['Tumor_Sample_Barcode'], values=[vaf_col])
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''
        df = df.fillna(0)

        for sample in samples:
            if sample not in df.columns:
                df[sample] = 0

        df = df[samples]

        if anchor is None:
            anchor = samples[0]

        df = df.sort_values(by=anchor, ascending=False)
        if normalize:
            df = df / df.max()
        df = df.fillna(0)
        df = df.iloc[:count, :].T
        df = df.loc[samples]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.lineplot(data=df, ax=ax, **kwargs)

        ax.set_ylabel('Fraction')

        return ax

    def plot_genepair(
        self, x, y, vaf_col, af=None, group_col=None, group_order=None,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a scatter plot of VAF between Gene X and Gene Y.

        Parameters
        ----------
        x, y : str
            Gene names.
        vaf_col : str
            MafFrame column containing VAF data.
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.scatterplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_genepair('DNMT3A', 'FLT3', 'i_TumorVAF_WU')
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_genepair('DNMT3A', 'FLT3', 'i_TumorVAF_WU',
            ...                  af=af,
            ...                  group_col='FAB_classification')
            >>> plt.tight_layout()
        """
        df = self.df[self.df.Hugo_Symbol.isin([x, y])]
        df = df[['Tumor_Sample_Barcode', 'Hugo_Symbol', vaf_col]]
        df = df.sort_values(vaf_col, ascending=False)
        df = df.drop_duplicates(subset=['Tumor_Sample_Barcode', 'Hugo_Symbol'])
        df = df.pivot(index='Tumor_Sample_Barcode',
            columns='Hugo_Symbol', values=vaf_col)
        df = df.fillna(0)

        if group_col is not None:
            df = df.merge(af.df[group_col], left_index=True, right_index=True)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.scatterplot(
            x=x, y=y, data=df, ax=ax, hue=group_col, hue_order=group_order,
            **kwargs
        )

        # Print summary statistics including R-squared and p-value.
        results = smf.ols(f'{y} ~ {x}', data=df).fit()
        print(f'Results for {y} ~ {x}:')
        print(f'R^2 = {results.rsquared:.2f}')
        print(f'  P = {results.f_pvalue:.2e}')

        return ax

    def plot_regplot(
        self, af, group_col, a, b, a_size=None, b_size=None, genes=None,
        count=10, to_csv=None, ax=None, figsize=None, **kwargs
    ):
        """
        Create a scatter plot with a linear regression model fit visualizing
        correlation between gene mutation frequencies in two sample groups
        A and B.

        Each point in the plot represents a gene.

        The method will automatically calculate and print summary statistics
        including R-squared and p-value.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        group_col : str
            AnnFrame column containing sample group information.
        a, b : str
            Sample group names.
        a_size, b_size : int, optional
            Sample group sizes to use as denominator. By default, these are
            inferred from the MafFrame and AnnFrame objects.
        genes : list, optional
            Genes to display. When absent, top mutated genes (``count``) will
            be used.
        count : int, defualt: 10
            Number of top mutated genes to display. Ignored if ``genes`` is
            specified.
        to_csv : str, optional
            Write the plot's data to a CSV file.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.regplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_regplot(af, 'FAB_classification', 'M1', 'M2')
            Results for M2 ~ M1:
            R^2 = 0.43
              P = 3.96e-02
            >>> plt.tight_layout()
        """
        df1 = self.matrix_prevalence()
        df2 = af.df[af.df.index.isin(df1.columns)]
        i_a = df2[df2[group_col] == a].index
        i_b = df2[df2[group_col] == b].index

        # Determine which genes to display.
        if genes is None:
            genes = self.matrix_genes(count=count).index.to_list()

        # Determine each group's sample size.
        if a_size is None:
            a_size = len(i_a)
        if b_size is None:
            b_size = len(i_b)

        f = lambda x: 0 if x == 0 else 1
        s_a = df1.T.loc[i_a].applymap(f).sum().loc[genes] / a_size
        s_b = df1.T.loc[i_b].applymap(f).sum().loc[genes] / b_size
        df3 = pd.concat([s_a, s_b], axis=1)
        df3.columns = [a, b]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        # Draw the main plot.
        sns.regplot(x=a, y=b, data=df3, ax=ax, **kwargs)

        # Write the DataFrame to a CSV file.
        if to_csv is not None:
            df3.to_csv(to_csv)

        # Print summary statistics including R-squared and p-value.
        results = smf.ols(f'{b} ~ {a}', data=df3).fit()
        print(f'Results for {b} ~ {a}:')
        print(f'R^2 = {results.rsquared:.2f}')
        print(f'  P = {results.f_pvalue:.2e}')

        return ax

    def plot_interactions(
        self, count=10, cmap=None, ax=None, figsize=None, **kwargs
    ):
        """
        Create a heatmap representing mutually exclusive or co-occurring set
        of genes.

        This method performs pair-wise Fisher’s Exact test to detect such
        significant pair of genes.

        Parameters
        ----------
        count : int, defualt: 10
            Number of top mutated genes to display.
        cmap : str, optional
            Color map.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.heatmap`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_interactions(count=25, cmap='BrBG')
            >>> plt.tight_layout()
        """
        df = self.matrix_prevalence()
        genes = self.matrix_genes(count=count, mode='samples').index.to_list()
        df = df.loc[genes]
        df = df.applymap(lambda x: True if x else False)
        df = df.T
        pairs = list(itertools.combinations(genes, 2))
        data = []

        def one_pair(a, b):
            s_a = df[a].to_list()
            s_b = df[b].to_list()
            ab = 0
            AB = 0
            aB = 0
            Ab = 0
            for i in range(len(s_a)):
                if s_a[i] and s_b[i]:
                    AB += 1
                elif s_a[i] and not s_b[i]:
                    Ab += 1
                elif not s_a[i] and s_b[i]:
                    aB += 1
                else:
                    ab += 1
            return (ab, AB, aB, Ab)

        for pair in pairs:
            a = pair[0]
            b = pair[1]
            ab, AB, aB, Ab = one_pair(a, b)
            event = 'Co_Occurence' if AB else 'Mutually_Exclusive'
            data.append([a, b, ab, AB, aB, Ab, event])

        df = pd.DataFrame(data,
            columns=['A', 'B', 'ab', 'AB', 'aB', 'Ab', 'Event'])

        def one_row(r):
            oddsr, p = fisher_exact([[r.AB, r.aB], [r.Ab, r.ab]],
                alternative='two-sided')
            return pd.Series([oddsr, p], index=['Odds_Ratio', 'P_Value'])

        df = pd.concat([df.apply(one_row, axis=1), df], axis=1)
        df = df.sort_values('P_Value')

        def one_row(r):
            r['Log_P_Value'] = -np.log10(r.P_Value)
            if r.P_Value < 0.05:
                r['Label'] = '*'
            elif r.P_Value < 0.1:
                r['Label'] = '.'
            else:
                r['Label'] = ''
            if r.Event == 'Mutually_Exclusive':
                r.Log_P_Value *= -1
            return r

        df = df.apply(one_row, axis=1)

        annot = df.pivot(index='A', columns='B', values='Label')
        annot = annot.fillna('')

        df = df.pivot(index='A', columns='B', values='Log_P_Value')
        df = df.fillna(0)

        for gene in genes:
            if gene not in df.columns:
                df[gene] = 0
            if gene not in annot.columns:
                annot[gene] = ''

        df = df.T
        annot = annot.T

        for gene in genes:
            if gene not in df.columns:
                df[gene] = 0
            if gene not in annot.columns:
                annot[gene] = ''

        annot = annot[genes]
        annot = annot.loc[genes]

        df = df[genes]
        df = df.loc[genes]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        # Create a mask for the heatmap.
        corr = np.corrcoef(np.random.randn(count, 200))
        mask = np.zeros_like(corr)
        mask[np.triu_indices_from(mask)] = True

        sns.heatmap(
            df, annot=annot, fmt='', cmap=cmap, mask=mask, vmax=3, vmin=-3,
            center=0, ax=ax, **kwargs
        )

        ax.set_xlabel('')
        ax.set_ylabel('')

        return ax

    def plot_lollipop(
        self, gene, alpha=0.7, ax=None, figsize=None, legend=True
    ):
        """
        Create a lollipop or stem plot showing amino acid changes of a gene.

        Parameters
        ----------
        gene : str
            Name of the gene.
        alpha : float, default: 0.7
            Set the color transparency. Must be within the 0-1 range,
            inclusive.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_lollipop('DNMT3A')
            >>> plt.tight_layout()
        """
        # Only select variants from the gene.
        df1 = self.df[self.df.Hugo_Symbol == gene]

        # Raise an error if there are no SNVs to plot.
        if df1.empty:
            raise ValueError(f"No variants to plot for the gene: '{gene}'.")

        # Count each amino acid change.
        df2 = df1.Protein_Change.value_counts().to_frame().reset_index()
        df2.columns = ['Protein_Change', 'Count']

        # Identify variant classification for each amino acid change.
        df3 = df1[['Protein_Change', 'Variant_Classification']
            ].drop_duplicates(subset=['Protein_Change'])
        df4 = pd.merge(df2, df3, on='Protein_Change')

        # Extract amino acid positions. Sort the counts by position.
        def one_row(r):
            digits = [x for x in r.Protein_Change if x.isdigit()]
            if not digits:
                return np.nan
            return int(''.join(digits))
        df4['Protein_Position'] = df4.apply(one_row, axis=1)
        df4 = df4.dropna(subset=['Protein_Position'])
        df4 = df4.sort_values(['Protein_Position'])

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        for i, nonsyn_name in enumerate(NONSYN_NAMES):
            temp = df4[df4.Variant_Classification == nonsyn_name]
            color = NONSYN_COLORS[i]
            ax.vlines(temp.Protein_Position, ymin=0, ymax=temp.Count,
                alpha=alpha, color=color)
            ax.plot(temp.Protein_Position, temp.Count, 'o', alpha=alpha,
                color=color, label=nonsyn_name)

        ax.set_xlabel('Position')
        ax.set_ylabel('Count')

        if legend:
            ax.legend()

        return ax

    def plot_mutated(
        self, af=None, group_col=None, group_order=None, genes=None,
        count=10, ax=None, figsize=None
    ):
        """
        Create a bar plot visualizing the mutation prevalence of top
        mutated genes.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        genes : list, optional
            Genes to display. When absent, top mutated genes (``count``) will
            be used.
        count : int, defualt: 10
            Number of top mutated genes to display. Ignored if ``genes`` is
            specified.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_mutated()
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_mutated(af=af,
            ...                 group_col='FAB_classification',
            ...                 group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        df = self.matrix_prevalence()

        # Determine which genes to display.
        if genes is None:
            genes = self.matrix_genes(count=count).index.to_list()

        df = df.loc[genes]
        df = df.applymap(lambda x: True if x else False)
        if group_col is None:
            df = (df.sum(axis=1) / df.shape[1]).to_frame().reset_index()
            df.columns.values[1] = 'Prevalence'
        else:
            df = df.T
            df = pd.merge(df, af.df[group_col], left_index=True, right_index=True)
            df = df.groupby([group_col]).mean().reset_index()
            df = df.melt(id_vars=[group_col])
            df.columns = [group_col, 'Hugo_Symbol', 'Prevalence']

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.barplot(
            x='Hugo_Symbol', y='Prevalence', data=df, hue=group_col,
            hue_order=group_order, ax=ax
        )

        ax.set_xlabel('')

        return ax

    def plot_mutated_matched(
        self, af, patient_col, group_col, group_order, ax=None, figsize=None,
        **kwargs
    ):
        """
        Create a bar plot visualizing the mutation prevalence of top
        mutated genes.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list
            List of sample group names.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        df = self.matrix_waterfall_matched(af, patient_col, group_col, group_order)
        df = df.applymap(lambda x: 0 if x == 'None' else 1)
        s = df.sum(axis=1) / len(df.columns) * 100
        s.name = 'Count'
        df = s.to_frame().reset_index()

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.barplot(
            x='Count', y='Gene', hue='Group', data=df, hue_order=group_order,
            orient='h', ax=ax, **kwargs
        )

        ax.set_xlabel('Patients (%)')
        ax.set_ylabel('')

        return ax

    def plot_rainfall(
        self, sample, palette=None, legend='auto', ax=None, figsize=None,
        **kwargs
    ):
        """
        Create a rainfall plot visualizing inter-variant distance on a linear
        genomic scale for single sample.

        Parameters
        ----------
        sample : str
            Name of the sample.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        legend : {'auto', 'brief', 'full', False}, default: 'auto'
            Display setting of the legend according to
            :meth:`seaborn.scatterplot`.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.scatterplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pymaf
            >>> common.load_dataset('brca')
            >>> maf_file = '~/fuc-data/brca/brca.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_rainfall('TCGA-A8-A08B',
            ...                  figsize=(14, 7),
            ...                  palette=sns.color_palette('Set2')[:6])
            >>> plt.tight_layout()
        """
        # Select variants from the sample.
        df = self.df[self.df.Tumor_Sample_Barcode == sample]

        # Remove indels.
        df = df[df.Variant_Type == 'SNP']

        # Raise an error if there are no SNVs to plot.
        if df.empty:
            message = (
                'There are no SNVs to be drawn '
                f"for the sample: '{sample}'."
            )
            raise ValueError(message)

        # Get SNV class for each variant.
        def one_row(r):
            change = r.Reference_Allele + '>' + r.Tumor_Seq_Allele2
            return SNV_CLASSES[change]['class']
        df['SNV_Class'] = df.apply(one_row, axis=1)

        # Convert string chromosomes to integers for ordering.
        def one_row(r):
            r.Chromosome = int(r.Chromosome.replace(
                'chr', '').replace('X', '23').replace('Y', '24'))
            return r
        df = df.apply(one_row, axis=1)
        df = df[['Chromosome', 'Start_Position', 'SNV_Class']]
        df = df.sort_values(['Chromosome', 'Start_Position'])

        # Update positions as if all chromosomes are one long molecule.
        def one_row(r):
            if r.Chromosome == 1:
                return r
            r.Start_Position += sum(CHROM_LENGTHS['hg19'][:r.Chromosome-1])
            return r
        df = df.apply(one_row, axis=1)
        s = np.diff(df.Start_Position)
        s = np.insert(s, 0, 0)
        s = np.log10(s + 1)
        df['Interevent_Distance'] = s
        df = df.reset_index(drop=True)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        bounds = [0] + df.drop_duplicates(subset=['Chromosome'],
            keep='last').index.to_list()

        xticks = []
        for i, bound in enumerate(bounds):
            if i == 0:
                continue
            elif i == 1:
                xticks.append(bound / 2)
            else:
                xticks.append(bounds[i-1] + (bound - bounds[i-1]) / 2)

        for bound in bounds:
            ax.axvline(x=bound, color='lightgray', zorder=1)

        sns.scatterplot(
            x=df.index, y='Interevent_Distance', data=df, hue='SNV_Class',
            hue_order=SNV_CLASS_ORDER, palette=palette, ax=ax, legend=legend,
            zorder=2, **kwargs
        )

        ax.set_xlabel('Chromosomes')
        ax.set_ylabel('Interevent distance')
        ax.set_xticks(xticks)
        ax.set_xticklabels(['X' if x == 23 else 'Y' if x == 24 else x
            for x in df.Chromosome.unique()])

        return ax

    def plot_snvclsc(
        self, af=None, group_col=None, group_order=None, palette=None,
        flip=False, ax=None, figsize=None, **kwargs
    ):
        """
        Create a bar plot summarizing the count distrubtions of the six
        :ref:`glossary:SNV classes` for all samples.

        A grouped bar plot can be created with ``group_col`` (requires an AnnFrame).

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        MafFrame.plot_snvclsp
            Create a box plot summarizing the proportion distrubtions of
            the six :ref:`glossary:SNV classes` for all sample.
        MafFrame.plot_snvclss
            Create a bar plot showing the proportions of the six
            :ref:`glossary:SNV classes` for individual samples.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_snvclsc(palette=sns.color_palette('Dark2'))
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_snvclsc(af=af,
            ...                 group_col='FAB_classification',
            ...                 group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        # Add the SNV_Class column.
        df = self.df[self.df.Variant_Type == 'SNP']
        def one_row(r):
            change = r.Reference_Allele + '>' + r.Tumor_Seq_Allele2
            return SNV_CLASSES[change]['class']
        s = df.apply(one_row, axis=1)
        s.name = 'SNV_Class'
        df = pd.concat([df, s], axis=1)

        # Count the occurance of each SNV class.
        if group_col is not None:
            df = pd.merge(df, af.df[group_col], left_on='Tumor_Sample_Barcode',
                right_index=True)
            s = df.groupby([group_col]).SNV_Class.value_counts()
            df = s.to_frame().rename(columns={'SNV_Class': 'Count'}
                ).reset_index()
        else:
            s = df.SNV_Class.value_counts()
            df = s.to_frame().reset_index()
            df.columns = ['SNV_Class', 'Count']

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = 'Count', 'SNV_Class'
            xlabel, ylabel = 'Count', ''
        else:
            x, y = 'SNV_Class', 'Count'
            xlabel, ylabel = '', 'Count'

        sns.barplot(
            x=x, y=y, data=df, ax=ax, hue=group_col, hue_order=group_order,
            palette=palette, order=SNV_CLASS_ORDER, **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_snvclsp(
        self, af=None, group_col=None, group_order=None, palette=None, flip=False,
        ax=None, figsize=None, **kwargs
    ):
        """
        Create a box plot summarizing the proportion distrubtions of the six
        :ref:`glossary:SNV classes` for all sample.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.boxplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        MafFrame.plot_snvclsc
            Create a bar plot summarizing the count distrubtions of the six
            :ref:`glossary:SNV classes` for all samples.
        MafFrame.plot_snvclss
            Create a bar plot showing the proportions of the six
            :ref:`glossary:SNV classes` for individual samples.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import seaborn as sns
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_snvclsp(palette=sns.color_palette('Set2'))
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_snvclsp(af=af,
            ...                 group_col='FAB_classification',
            ...                 group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        # Add the SNV_Class column.
        df = self.df[self.df.Variant_Type == 'SNP']
        def one_row(r):
            change = r.Reference_Allele + '>' + r.Tumor_Seq_Allele2
            return SNV_CLASSES[change]['class']
        s = df.apply(one_row, axis=1)
        s.name = 'SNV_Class'
        df = pd.concat([df, s], axis=1)

        # Compute the proportions of SNV classes in each sample.
        s = df.groupby('Tumor_Sample_Barcode')['SNV_Class'].value_counts()
        s.name = 'Count'
        df = s.to_frame().reset_index()
        df = df.pivot(index='Tumor_Sample_Barcode', columns='SNV_Class')
        df = df.fillna(0)
        df = df.apply(lambda r: r/r.sum(), axis=1)
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''

        if group_col is None:
            df = pd.melt(df, var_name='SNV_Class', value_name='Proportion')
        else:
            df = pd.merge(df, af.df[group_col], left_index=True, right_index=True)
            df = pd.melt(df, id_vars=[group_col], var_name='SNV_Class',
                value_name='Proportion')

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = 'Proportion', 'SNV_Class'
            xlabel, ylabel = 'Proportion', ''
        else:
            x, y = 'SNV_Class', 'Proportion'
            xlabel, ylabel = '', 'Proportion'

        sns.boxplot(
            x=x, y=y, data=df, hue=group_col, hue_order=group_order,
            palette=palette, ax=ax, **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_snvclss(
        self, samples=None, color=None, colormap=None, width=0.8,
        legend=True, flip=False, to_csv=None, ax=None, figsize=None, **kwargs
    ):
        """
        Create a bar plot showing the proportions of the six
        :ref:`glossary:SNV classes` for individual samples.

        Parameters
        ----------
        samples : list, optional
            List of samples to display (in that order too). If samples that
            are absent in the MafFrame are provided, the method will give a
            warning but still draw an empty bar for those samples.
        color : list, optional
            List of color tuples. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        colormap : str or matplotlib colormap object, optional
            Colormap to select colors from. See the :ref:`tutorials:Control
            plot colors` tutorial for details.
        width : float, default: 0.8
            The width of the bars.
        legend : bool, default: True
            Place legend on axis subplots.
        flip : bool, default: False
            If True, flip the x and y axes.
        to_csv : str, optional
            Write the plot's data to a CSV file.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.bar` or
            :meth:`pandas.DataFrame.plot.barh`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        MafFrame.plot_snvclsc
            Create a bar plot summarizing the count distrubtions of the six
            :ref:`glossary:SNV classes` for all samples.
        MafFrame.plot_snvclsp
            Create a box plot summarizing the proportion distrubtions of
            the six :ref:`glossary:SNV classes` for all sample.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> ax = mf.plot_snvclss(width=1, color=plt.get_cmap('Set2').colors)
            >>> ax.legend(loc='upper right')
            >>> plt.tight_layout()
        """
        # Add the SNV_Class column.
        df = self.df[self.df.Variant_Type == 'SNP']
        def one_row(r):
            change = r.Reference_Allele + '>' + r.Tumor_Seq_Allele2
            return SNV_CLASSES[change]['class']
        s = df.apply(one_row, axis=1)
        s.name = 'SNV_Class'
        df = pd.concat([df, s], axis=1)

        # Compute the proportions of SNV classes in each sample.
        s = df.groupby('Tumor_Sample_Barcode')['SNV_Class'].value_counts()
        s.name = 'Count'
        df = s.to_frame().reset_index()
        df = df.pivot(index='Tumor_Sample_Barcode', columns='SNV_Class')
        df = df.fillna(0)
        df = df.apply(lambda r: r/r.sum(), axis=1)
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''
        for x in SNV_CLASS_ORDER:
            if x not in df.columns:
                df[x] = 0
        df = df[SNV_CLASS_ORDER]

        # Determine which samples should be displayed.
        if samples is not None:
            missing_samples = []
            missing_data = []
            for sample in samples:
                if sample not in df.index:
                    missing_samples.append(sample)
                    missing_data.append([0] * 6)
            if missing_samples:
                message = (
                    'Although the following samples are absent in the '
                    'MafFrame, they will still be displayed as empty bar: '
                    f'{missing_samples}.'
                )
                warnings.warn(message)
                temp = pd.DataFrame(missing_data)
                temp.index = missing_samples
                temp.columns = SNV_CLASS_ORDER
                df = pd.concat([df, temp]).loc[samples]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            kind = 'barh'
            xlabel, ylabel = 'Proportion', 'Samples'
        else:
            kind = 'bar'
            xlabel, ylabel = 'Samples', 'Proportion'

        df.plot(
            kind=kind, ax=ax, stacked=True, legend=legend, width=width,
            color=color, colormap=colormap, **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if flip:
            ax.set_yticks([])
        else:
            ax.set_xticks([])

        # Write the DataFrame to a CSV file.
        if to_csv is not None:
            df.to_csv(to_csv)

        return ax

    def plot_titv(
        self, af=None, group_col=None, group_order=None, flip=False, ax=None,
        figsize=None, **kwargs
    ):
        """
        Create a box plot showing the :ref:`Ti/Tv <glossary:Transitions (Ti)
        and transversions (Tv)>` proportions of samples.

        Parameters
        ----------
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.boxplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        See Also
        --------
        fuc.api.pyvcf.VcfFrame.plot_titv
            Similar method for the :class:`fuc.api.pyvcf.VcfFrame` class.

        Examples
        --------
        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_titv()
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_titv(af=af,
            ...              group_col='FAB_classification',
            ...              group_order=['M0', 'M1', 'M2'])
            >>> plt.tight_layout()
        """
        df = self.df[self.df.Variant_Type == 'SNP']
        def one_row(r):
            change = r.Reference_Allele + '>' + r.Tumor_Seq_Allele2
            return SNV_CLASSES[change]['type']
        s = df.apply(one_row, axis=1)
        s.name = 'SNV_Type'
        df = pd.concat([df, s], axis=1)
        s = df.groupby('Tumor_Sample_Barcode')['SNV_Type'].value_counts()
        s.name = 'Count'
        df = s.to_frame().reset_index()
        df = df.pivot(index='Tumor_Sample_Barcode', columns='SNV_Type')

        df = df.fillna(0)
        df = df.apply(lambda r: r/r.sum(), axis=1)
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''

        if group_col is not None:
            df = pd.merge(df, af.df[group_col], left_index=True, right_index=True)
            df = df.reset_index(drop=True)
            df = df.set_index(group_col)
            df = pd.melt(df, var_name='SNV_Type', value_name='Proportion',
                ignore_index=False)
            df = df.reset_index()
        else:
            df = pd.melt(df, var_name='SNV_Type', value_name='Proportion')

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = 'Proportion', 'SNV_Type'
            xlabel, ylabel = 'Proportion', ''
        else:
            x, y = 'SNV_Type', 'Proportion'
            xlabel, ylabel = '', 'Proportion'

        sns.boxplot(
            x=x, y=y, data=df, hue=group_col, hue_order=group_order, ax=ax,
            **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_summary(
        self, figsize=(15, 10), title_fontsize=16, ticklabels_fontsize=12,
        legend_fontsize=12

    ):
        """Create a summary figure for MafFrame.

        Parameters
        ----------
        figsize : tuple, default: (15, 10)
            Width, height in inches. Format: (float, float).
        title_fontsize : float, default: 16
            Font size of subplot titles.
        ticklabels_fontsize : float, default: 12
            Font size of tick labels.
        legend_fontsize : float, default: 12
            Font size of legend texts.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_summary()
        """
        g = {'height_ratios': [10, 10, 1]}
        fig, axes = plt.subplots(3, 3, figsize=figsize, gridspec_kw=g)
        [[ax1, ax2, ax3], [ax4, ax5, ax6], [ax7, ax8, ax9]] = axes
        gs = axes[2, 0].get_gridspec()
        for ax in axes[2, :]:
            ax.remove()
        axbig = fig.add_subplot(gs[2, :])

        # Create the 'Variant classification (variants)' figure.
        self.plot_varcls(ax=ax1)
        ax1.set_yticks([])
        ax1.set_title('Variant classification (variants)',
                      fontsize=title_fontsize)
        ax1.set_xlabel('')
        ax1.tick_params(axis='x', which='major',
                        labelsize=ticklabels_fontsize)

        # Create the 'Variant type' figure.
        self.plot_vartype(ax=ax2, palette='Pastel1', flip=True)
        ax2.set_title('Variant type', fontsize=title_fontsize)
        ax2.set_xlabel('')
        ax2.tick_params(axis='both', which='major',
                        labelsize=ticklabels_fontsize)

        # Create the 'SNV class' figure.
        self.plot_snvclsc(ax=ax3, flip=True,
            palette=sns.color_palette('Set2'))
        ax3.set_title('SNV class', fontsize=title_fontsize)
        ax3.set_xlabel('')
        ax3.tick_params(axis='both', which='major',
                        labelsize=ticklabels_fontsize)

        # Create the 'Variants per sample' figure.
        median = self.matrix_tmb().sum(axis=1).median()
        self.plot_tmb(ax=ax4, width=1)
        ax4.set_title(f'Variants per sample (median={median:.1f})',
                      fontsize=title_fontsize)
        ax4.set_xlabel('')
        ax4.set_ylabel('')
        ax4.tick_params(axis='y', which='major',
                        labelsize=ticklabels_fontsize)

        ax4.axhline(y=median, color='red', linestyle='dashed')

        # Create the 'Variant classification (samples)' figure.
        self.plot_varsum(ax=ax5)
        ax5.set_title('Variant classification (samples)',
                      fontsize=title_fontsize)
        ax5.set_yticks([])
        ax5.set_xlabel('')
        ax5.tick_params(axis='x', which='major',
                        labelsize=ticklabels_fontsize)

        # Create the 'Top 10 mutated genes' figure.
        self.plot_genes(ax=ax6)
        ax6.set_title('Top 10 mutated genes', fontsize=title_fontsize)
        ax6.set_xlabel('')
        ax6.tick_params(axis='both', which='major',
                        labelsize=ticklabels_fontsize)

        # Add the legend.
        axbig.legend(
            handles=common.legend_handles(NONSYN_NAMES, colors=NONSYN_COLORS),
            title='Variant_Classification',
            loc='upper center',
            ncol=3,
            fontsize=legend_fontsize,
            title_fontsize=legend_fontsize
        )
        axbig.axis('off')

        plt.tight_layout()

    def plot_tmb(
        self, samples=None, width=0.8, ax=None, figsize=None, **kwargs
    ):
        """
        Create a bar plot showing the :ref:`TMB <glossary:Tumor mutational
        burden (TMB)>` distributions of samples.

        Parameters
        ----------
        samples : list, optional
            List of samples to display (in that order too). If samples that
            are absent in the MafFrame are provided, the method will give a
            warning but still draw an empty bar for those samples.
        width : float, default: 0.8
            The width of the bars.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.bar`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------
        Below is a simple example:

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_tmb(width=1)
            >>> plt.tight_layout()
        """
        df = self.matrix_tmb()

        if samples is not None:
            df = df.T
            missing_samples = []
            for sample in samples:
                if sample not in df.columns:
                    missing_samples.append(sample)
            if missing_samples:
                message = (
                    'Although the following samples are absent in the '
                    'MafFrame, they will still be displayed as empty bar: '
                    f'{missing_samples}.'
                )
                warnings.warn(message)
                for missing_sample in missing_samples:
                    df[missing_sample] = 0
            df = df[samples]
            df = df.T

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        df.plot.bar(stacked=True, ax=ax, width=width, legend=False,
            color=NONSYN_COLORS, **kwargs)

        ax.set_xlabel('Samples')
        ax.set_ylabel('Count')
        ax.set_xticks([])

        return ax

    def plot_tmb_matched(
        self, af, patient_col, group_col, group_order=None, patients=None,
        legend=True, ax=None, figsize=None, **kwargs
    ):
        """
        Create a grouped bar plot showing TMB distributions for different
        group levels in each patient.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        patients : list, optional
            List of patient names.
        legend : bool, default: True
            Place legend on axis subplots.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`pandas.DataFrame.plot.bar`

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        df = self.matrix_tmb().T

        for sample in af.samples:
            if sample not in df.columns:
                df[sample] = 0

        s = df.sum()
        s.name = 'TMB'
        df = pd.concat([s, af.df[[patient_col, group_col]]], axis=1)

        df = df.pivot(index=patient_col, columns=group_col, values='TMB')

        if group_order is not None:
            df = df[group_order]

        i = df.sum(axis=1).sort_values(ascending=False).index
        df = df.loc[i]

        if patients is not None:
            df = df.loc[patients]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        df.plot(ax=ax, kind='bar', stacked=True, legend=legend, **kwargs)

        ax.set_xlabel('')
        ax.set_ylabel('TMB')

        return ax

    def plot_vaf(
        self, vaf_col, count=10, af=None, group_col=None, group_order=None,
        flip=False, sort=True, ax=None, figsize=None, **kwargs
    ):
        """
        Create a box plot showing the :ref:`VAF <glossary:Variant allele
        frequency (VAF)>` distributions of top mutated genes.

        A grouped box plot can be created with ``group_col`` (requires an
        AnnFrame).

        Parameters
        ----------
        vaf_col : str
            MafFrame column containing VAF data.
        count : int, default: 10
            Number of top mutated genes to display.
        af : AnnFrame, optional
            AnnFrame containing sample annotation data.
        group_col : str, optional
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        flip : bool, default: False
            If True, flip the x and y axes.
        sort : bool, default: True
            If False, do not sort the genes by median value.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.boxplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        Below is a simple example:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_vaf('i_TumorVAF_WU')
            >>> plt.tight_layout()

        We can create a grouped bar plot based on FAB classification:

        .. plot::
            :context: close-figs

            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_vaf('i_TumorVAF_WU',
            ...             af=af,
            ...             group_col='FAB_classification',
            ...             group_order=['M1', 'M2', 'M3'],
            ...             count=5)
            >>> plt.tight_layout()
        """
        genes = self.matrix_genes(count=count).index.to_list()

        if sort:
            medians = self.df.groupby('Hugo_Symbol')[vaf_col].median()
            genes = medians[genes].sort_values(
                ascending=False).index.to_list()

        df = self.df[self.df.Hugo_Symbol.isin(genes)]

        if group_col is not None:
            df = pd.merge(df, af.df, left_on='Tumor_Sample_Barcode',
                right_index=True)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = vaf_col, 'Hugo_Symbol'
            xlabel, ylabel = 'VAF', ''
        else:
            x, y = 'Hugo_Symbol', vaf_col
            xlabel, ylabel = '', 'VAF'

        sns.boxplot(
            x=x, y=y, data=df, ax=ax, order=genes, hue=group_col,
            hue_order=group_order, **kwargs
        )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_varcls(self, ax=None, figsize=None, **kwargs):
        """Create a bar plot for the nonsynonymous variant classes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`matplotlib.axes.Axes.bar` and :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_varcls()
            >>> plt.tight_layout()
        """
        d = self.df.Variant_Classification.value_counts().to_dict()
        counts = {}
        for varcls in NONSYN_NAMES:
            if varcls in d:
                counts[varcls] = d[varcls]
            else:
                counts[varcls] = 0
        s = pd.Series(counts).reindex(index=NONSYN_NAMES)
        df = s.to_frame().reset_index()
        df.columns = ['Variant_Classification', 'Count']

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.barplot(x='Count', y='Variant_Classification', data=df,
                    ax=ax, palette=NONSYN_COLORS, **kwargs)

        ax.set_ylabel('')

        return ax

    def plot_matrixg(
        self, gene, af, group_col, group_order=None, cbar=True, ax=None,
        figsize=None, **kwargs
    ):
        """
        Create a heatmap of count matrix with a shape of (sample groups,
        protein changes).

        Parameters
        ----------
        gene : str
            Name of the gene.
        af : AnnFrame
            AnnFrame containing sample annotation data.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list, optional
            List of sample group names.
        cbar : bool, default: True
            Whether to draw a colorbar.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.heatmap`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> annot_file = '~/fuc-data/tcga-laml/tcga_laml_annot.tsv'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> af = common.AnnFrame.from_file(annot_file, sample_col=0)
            >>> mf.plot_matrixg('IDH1', af, 'FAB_classification', linewidth=0.5, square=True, annot=True)
            >>> plt.tight_layout()
        """
        df = self.df[self.df.Hugo_Symbol == gene]

        if df.empty:
            raise ValueError(f'No protein changes were found: {gene}')

        df = df[['Tumor_Sample_Barcode', 'Protein_Change']]
        df = df[df.Protein_Change != '.']
        df = df.merge(af.df[group_col], left_on='Tumor_Sample_Barcode', right_index=True)
        s = df.groupby(group_col)['Protein_Change'].value_counts()
        s.name = 'Count'
        df = s.to_frame().reset_index()
        df = df.pivot(index=group_col,
                      columns='Protein_Change',
                      values='Count')
        df = df.fillna(0)

        if group_order is not None:
            missing_groups = []
            for group in group_order:
                if group not in df.index:
                    missing_groups.append(group)
            if missing_groups:
                message = (
                    'Although the following sample groups are absent in the '
                    'MafFrame, they will still be displayed as empty bar: '
                    f'{missing_groups}.'
                )
                warnings.warn(message)
                df = df.T
                for missing_group in missing_groups:
                    df[missing_group] = 0
                df = df[groups]
                df = df.T

        # Sort protein changes by position.
        f = lambda s: int(''.join([x for x in list(s) if x.isdigit()]))
        l = df.columns
        s = pd.Series([f(x) for x in l], index=l).sort_values()
        df = df[s.index]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(
            df, ax=ax, cbar=cbar, **kwargs
        )

        ax.set_xlabel('')
        ax.set_ylabel('')

        return ax

    def plot_matrixs(
        self, gene, samples=None, c0='lightgray', c1='red', l0='0', l1='1',
        cbar=True, square=False, ax=None, figsize=None, **kwargs
    ):
        """
        Create a heatmap of presence/absence matrix with a shape of (samples,
        protein changes).

        Parameters
        ----------
        gene : str
            Name of the gene.
        samples : list, optional
            List of samples to display (in that order too). If samples that
            are absent in the MafFrame are provided, the method will give a
            warning but still draw an empty bar for those samples.
        c0 : str, default: 'lightgray'
            Color for absence.
        c1 : str, default: 'red'
            Color for presence.
        l0 : str, default: '0'
            Label for absence.
        l1 : str, default: '1'
            Label for presence.
        cbar : bool, default: True
            Whether to draw a colorbar.
        square : bool, default: False
            If True, set the Axes aspect to "equal" so each cell will be
            square-shaped.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.heatmap`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_matrixs('KRAS', linewidth=0.5, square=True)
            >>> plt.tight_layout()
        """
        df = self.df[self.df.Hugo_Symbol == gene]

        if df.empty:
            raise ValueError(f'No protein changes were found: {gene}')

        df = df[['Tumor_Sample_Barcode', 'Protein_Change']]
        df = df[df.Protein_Change != '.']
        df['Presence'] = 1
        df = df.pivot(index='Tumor_Sample_Barcode',
                      columns='Protein_Change',
                      values='Presence')
        df = df.fillna(0)

        if samples is not None:
            missing_samples = []
            for sample in samples:
                if sample not in df.index:
                    missing_samples.append(sample)
            if missing_samples:
                message = (
                    'Although the following samples are absent in the '
                    'MafFrame, they will still be displayed as empty bar: '
                    f'{missing_samples}.'
                )
                warnings.warn(message)
                df = df.T
                for missing_sample in missing_samples:
                    df[missing_sample] = 0
                df = df[samples]
                df = df.T

        # Sort protein changes by position.
        f = lambda s: int(''.join([x for x in list(s) if x.isdigit()]))
        l = df.columns
        s = pd.Series([f(x) for x in l], index=l).sort_values()
        df = df[s.index]

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if len(np.unique(df.values)) == 1:
            cmap = [c1]
            cbar_ticklabels = [l1]
        else:
            cmap = [c0, c1]
            cbar_ticklabels = [l0, l1]

        sns.heatmap(
            df, cbar=cbar, cmap=cmap, square=square, ax=ax, **kwargs
        )

        if cbar:
            colorbar = ax.collections[0].colorbar
            n=len(cmap)
            r = colorbar.vmax - colorbar.vmin
            colorbar.set_ticks([colorbar.vmin + r/n * (0.5+i) for i in range(n)])
            colorbar.set_ticklabels(cbar_ticklabels)

        ax.set_xlabel('')
        ax.set_ylabel('')

        return ax

    def plot_varsum(self, flip=False, ax=None, figsize=None):
        """
        Create a summary box plot for variant classifications.

        Parameters
        ----------
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_varsum()
            >>> plt.tight_layout()
        """
        df = self.matrix_tmb()
        df = pd.melt(df, value_vars=df.columns)

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = 'variable', 'value'
            xlabel, ylabel = '', 'Samples'
        else:
            x, y = 'value', 'variable'
            xlabel, ylabel = 'Samples', ''

        sns.boxplot(x=x, y=y, data=df, ax=ax, showfliers=False,
            palette=NONSYN_COLORS)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return ax

    def plot_vartype(self, palette=None, flip=False, ax=None, figsize=None, **kwargs):
        """
        Create a bar plot summarizing the count distrubtions of viaration
        types for all samples.

        Parameters
        ----------
        palette : str, optional
            Name of the seaborn palette. See the :ref:`tutorials:Control plot
            colors` tutorial for details.
        flip : bool, default: False
            If True, flip the x and y axes.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.barplot`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_vartype()
            >>> plt.tight_layout()
        """
        s = self.df.Variant_Type.value_counts()
        df = s.to_frame().reset_index()

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        if flip:
            x, y = 'Variant_Type', 'index'
            xlabel, ylabel = 'Count', ''
        else:
            x, y = 'index', 'Variant_Type'
            xlabel, ylabel = '', 'Count'

        sns.barplot(x=x, y=y, data=df, ax=ax, palette=palette, **kwargs)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        return ax

    def plot_waterfall(
        self, count=10, keep_empty=False, samples=None, ax=None,
        figsize=None, **kwargs
    ):
        """
        Create a waterfall plot (oncoplot).

        See this :ref:`tutorial <tutorials:Create customized oncoplots>` to
        learn how to create customized oncoplots.

        Parameters
        ----------
        count : int, default: 10
            Number of top mutated genes to display.
        keep_empty : bool, default: False
            If True, display samples that do not have any mutations.
        samples : list, optional
            List of samples to display (in that order too). If samples that
            are absent in the MafFrame are provided, the method will give a
            warning but still draw an empty bar for those samples.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).
        kwargs
            Other keyword arguments will be passed down to
            :meth:`seaborn.heatmap`.

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.

        Examples
        --------

        .. plot::

            >>> import matplotlib.pyplot as plt
            >>> from fuc import common, pymaf
            >>> common.load_dataset('tcga-laml')
            >>> maf_file = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
            >>> mf = pymaf.MafFrame.from_file(maf_file)
            >>> mf.plot_waterfall(linewidths=0.5)
            >>> plt.tight_layout()
        """
        df = self.matrix_waterfall(count=count, keep_empty=keep_empty)

        if samples is not None:
            missing_samples = []
            for sample in samples:
                if sample not in df.columns:
                    missing_samples.append(sample)
            if missing_samples:
                message = (
                    'Although the following samples are absent in the '
                    'MafFrame, they will still be displayed as empty bar: '
                    f'{missing_samples}.'
                )
                warnings.warn(message)
                for missing_sample in missing_samples:
                    df[missing_sample] = 'None'
            df = df[samples]

        # Apply the mapping between items and integers.
        l = reversed(NONSYN_NAMES + ['Multi_Hit', 'None'])
        d = {k: v for v, k in enumerate(l)}
        df = df.applymap(lambda x: d[x])

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        colors = list(reversed(NONSYN_COLORS + ['k', 'lightgray']))

        sns.heatmap(df, cmap=colors, ax=ax, xticklabels=False,
            cbar=False, **kwargs)

        ax.set_xlabel('Samples')
        ax.set_ylabel('')

        return ax

    def plot_waterfall_matched(
        self, af, patient_col, group_col, group_order, count=10, ax=None,
        figsize=None
    ):
        """
        Create a waterfall plot using matched samples from each patient.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list
            List of sample group names.
        count : int, default: 10
            Number of top mutated genes to include.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        df = self.matrix_waterfall_matched(af, patient_col,
            group_col, group_order, count=count)
        genes = df.index.get_level_values(0).unique().to_list()

        l = reversed(NONSYN_NAMES + ['Multi_Hit', 'None'])
        d = {k: v for v, k in enumerate(l)}
        df = df.applymap(lambda x: d[x])

        colors = list(reversed(NONSYN_COLORS + ['k', 'lightgray']))

        # Determine which matplotlib axes to plot on.
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(df, cmap=colors, xticklabels=True, cbar=False, ax=ax)

        n = len(group_order)
        i = n / 2
        yticks = []
        for gene in genes:
            yticks.append(i)
            i += n

        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_yticks(yticks)
        ax.set_yticklabels(genes)

        # Add horizontal lines.
        for i, gene in enumerate(genes, start=1):
            ax.axhline(i*n, color='white')

        # Add vertical lines.
        for i, sample in enumerate(af.samples, start=1):
            ax.axvline(i, color='white')

        return ax

    def matrix_waterfall_matched(
        self, af, patient_col, group_col, group_order, count=10
    ):
        """
        Compute a matrix of variant classifications with a shape of
        (gene-group pairs, patients).

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        patient_col : str
            AnnFrame column containing patient information.
        group_col : str
            AnnFrame column containing sample group information.
        group_order : list
            List of sample group names.
        count : int, default: 10
            Number of top mutated genes to include.
        ax : matplotlib.axes.Axes, optional
            Pre-existing axes for the plot. Otherwise, crete a new one.
        figsize : tuple, optional
            Width, height in inches. Format: (float, float).

        Returns
        -------
        matplotlib.axes.Axes
            The matplotlib axes containing the plot.
        """
        df = self.matrix_waterfall(count=count)

        genes = df.index

        missing_samples = {}

        for sample in af.samples:
            if sample not in df.columns:
                missing_samples[sample] = ['None'] * len(genes)

        df = pd.concat(
            [df, pd.DataFrame(missing_samples, index=genes)], axis=1)

        df = df[af.samples].T
        df = df.merge(af.df[[patient_col, group_col]],
            left_index=True, right_index=True)
        df = df.reset_index(drop=True)

        temps = []

        for group in group_order:
            temp = df[df[group_col] == group].set_index(patient_col)[genes]
            tuples = [(x, group) for x in genes]
            mi = pd.MultiIndex.from_tuples(tuples, names=['Gene', 'Group'])
            temp.columns = mi
            temps.append(temp)

        df = pd.concat(temps, axis=1)

        tuples = []

        for gene in genes:
            for group in group_order:
                tuples.append((gene, group))

        df = df[tuples]
        df = df.T

        c = df.applymap(lambda x: 0 if x == 'None' else 1).sort_values(
            df.index.to_list(), axis=1, ascending=False).columns
        df = df[c]

        return df


    def to_vcf(
        self, fasta=None, ignore_indels=False, cols=None, names=None
    ):
        """
        Write the MafFrame to a sorted VcfFrame.

        Converting from MAF to VCF is pretty straightforward for SNVs, but it
        can be challenging for INDELs and complex events involving multiple
        nucleotides (e.g. 'AAGG' → 'CCCG'). This is because, for the latter
        case we need to identify the "anchor" nucleotide for each event,
        which is crucial for constructing a properly formatted VCF. For
        example, a deletion event 'AGT' → '-' in MAF would have to be
        converted to 'CAGT' → 'C' in the VCF where 'C' is our anchor
        nucleotide. The position should be shifted by one as well.

        In order to tackle this issue, the method makes use of a reference
        assembly (i.e. FASTA file). If SNVs are your only concern, then you
        do not need a FASTA file and can just set ``ignore_indels`` as True.
        If you are going to provide a FASTA file, please make sure to select
        the appropriate one (e.g. one that matches the genome assembly). For
        example, if your MAF is in hg19/GRCh37, use the 'hs37d5.fa' file
        which can be freely downloaded from the 1000 Genomes Project.

        Parameters
        ----------
        fasta : str, optional
            FASTA file. Required if ``ignore_indels`` is False.
        ignore_indels : bool, default: False
            If True, do not include INDELs in the VcfFrame. Useful when
            a FASTA file is not available.
        cols : str or list, optional
            Column(s) in the MafFrame which contain additional genotype
            data of interest. If provided, these data will be added to
            individual sample genotypes (e.g. '0/1:0.23').
        names : str or list, optional
            Name(s) to be displayed in the FORMAT field (e.g. AD, AF, DP).
            If not provided, the original column name(s) will be displayed.

        Returns
        -------
        VcfFrame
            The VcfFrame object.

        Examples
        --------

        >>> from fuc import pymaf
        >>> mf = pymaf.MafFrame.from_file('in.maf')
        >>> vf = mf.to_vcf(fasta='hs37d5.fa')
        >>> vf = mf.to_vcf(ignore_indels=True)
        >>> vf = mf.to_vcf(fasta='hs37d5.fa', cols='i_TumorVAF_WU', names='AF')
        """
        if not ignore_indels and fasta is None:
            raise ValueError("A FASTA file is required when 'ignore_indels' "
                             "argument is False.")

        if cols is None:
            cols = []
        if names is None:
            names = []

        if isinstance(cols, str):
            cols = [cols]
        if isinstance(names, str):
            names = [names]

        if cols and not names:
            names = cols
        if len(cols) != len(names):
            raise ValueError("Arguments 'cols' and 'names' "
                             "have different lengths.")

        # Create the minimal VCF.
        index_cols = ['Chromosome', 'Start_Position',
                      'Reference_Allele', 'Tumor_Seq_Allele2']
        df = self.df.pivot(index=index_cols,
                           columns='Tumor_Sample_Barcode',
                           values='Tumor_Seq_Allele2')
        f = lambda x: '0/0' if pd.isnull(x) else '0/1'
        df = df.applymap(f)
        df.columns.name = None
        df = df.reset_index()
        df = df.rename(columns={'Chromosome': 'CHROM',
                                'Start_Position': 'POS',
                                'Reference_Allele': 'REF',
                                'Tumor_Seq_Allele2': 'ALT'})
        df['ID'] = '.'
        df['QUAL'] = '.'
        df['FILTER'] = '.'
        df['INFO'] = '.'
        df['FORMAT'] = 'GT'
        df = df[pyvcf.HEADERS + self.samples]

        # Add requested genotype information.
        f = lambda x: '.' if pd.isnull(x) else str(x)
        for i, col in enumerate(cols):
            _ = self.df.pivot(index=index_cols,
                              columns='Tumor_Sample_Barcode',
                              values='i_TumorVAF_WU')
            _ = _.reset_index()
            _ = _.drop(index_cols, axis=1)
            _ = _[self.samples]
            _ = _.applymap(f)
            df.iloc[:, 9:] = df.iloc[:, 9:] + ':' + _
            df.FORMAT = df.FORMAT + ':' + names[i]

        # Handle INDELs.
        l = ['A', 'C', 'G', 'T']
        if ignore_indels:
            i = (df.REF.isin(l)) & (df.ALT.isin(l))
            df = df[i]
        else:
            def one_row(r):
                if r.REF in l and r.ALT in l:
                    return r
                region = f'{r.CHROM}:{r.POS-1}-{r.POS-1}'
                anchor = common.extract_sequence(fasta, region)
                if not anchor:
                    return r
                r.POS = r.POS - 1
                if r.ALT == '-':
                    r.REF = anchor + r.REF
                    r.ALT = anchor
                elif r.REF == '-':
                    r.REF = anchor
                    r.ALT = anchor + r.ALT
                else:
                    r.REF = anchor + r.REF
                    r.ALT = anchor + r.ALT
                return r
            df = df.apply(one_row, axis=1)

        # Create the metadata.
        meta = [
            '##fileformat=VCFv4.3',
            '##source=fuc.api.pymaf.MafFrame.to_vcf',
        ]

        # Create the VcfFrame.
        vf = pyvcf.VcfFrame(meta, df)
        vf = vf.sort()

        return vf

    def to_file(self, fn):
        """Write MafFrame to a MAF file.

        Parameters
        ----------
        fn : str
            MAF file path.
        """
        with open(fn, 'w') as f:
            f.write(self.to_string())

    def to_string(self):
        """Render MafFrame to a console-friendly tabular output.

        Returns
        -------
        str
            String representation of MafFrame.
        """
        return self.df.to_csv(index=False, sep='\t')

    def filter_annot(self, af, expr):
        """
        Filter the MafFrame using sample annotation data.

        Samples are selected by querying the columns of an AnnFrame with a
        boolean expression. Samples not present in the MafFrame will be
        excluded automatically.

        Parameters
        ----------
        af : AnnFrame
            AnnFrame containing sample annotation data.
        expr : str
            Query expression to evaluate.

        Returns
        -------
        MafFrame
            Filtered MafFrame.

        Examples
        --------

        >>> from fuc import common, pymaf
        >>> common.load_dataset('tcga-laml')
        >>> mf = pymaf.MafFrame.from_file('~/fuc-data/tcga-laml/tcga_laml.maf.gz')
        >>> af = common.AnnFrame.from_file('~/fuc-data/tcga-laml/tcga_laml_annot.tsv', sample_col=0)
        >>> filtered_mf = mf.filter_annot(af, "FAB_classification == 'M4'")
        """
        samples = af.df.query(expr).index
        i = self.df.Tumor_Sample_Barcode.isin(samples)
        df = self.df[i]
        mf = self.__class__(df)
        return mf
