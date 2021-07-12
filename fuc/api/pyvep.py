"""
The pyvep submodule is designed for parsing VCF annotation data from the
`Ensembl VEP <https://asia.ensembl.org/info/docs/tools/vep/index.html>`_
program. It should be used with ``pyvcf.VcfFrame``.

The input VCF should already contain functional annotation data from
Ensembl VEP. The recommended method is Ensembl VEP's
`web interface <https://asia.ensembl.org/Tools/VEP>`_ with
“RefSeq transcripts” as the transcript database and the filtering option
“Show one selected consequence per variant”.

A typical VEP-annotated VCF file contains many fields ranging from gene
symbol to protein change. However, most of the analysis in pyvep uses the
following fields:

+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| No. | Name               | Description         | Examples                                                                                                 |
+=====+====================+=====================+==========================================================================================================+
| 1   | Allele             | Variant allele      | 'C', 'T', 'CG'                                                                                           |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 2   | Consequence        | Consequence type    | 'missense_variant', 'stop_gained'                                                                        |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 3   | IMPACT             | Consequence impact  | 'MODERATE', 'HIGH'                                                                                       |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 4   | SYMBOL             | Gene symbol         | 'TP53', 'APC'                                                                                            |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 5   | Gene               | Ensembl ID          | 7157, 55294                                                                                              |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 6   | BIOTYPE            | Transcript biotype  | 'protein_coding', 'pseudogene', 'miRNA'                                                                  |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 7   | EXON               | Exon number         | '8/17', '17/17'                                                                                          |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 8   | Protein_position   | Amino acid position | 234, 1510                                                                                                |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 9   | Amino_acids        | Amino acid changes  | 'R/\*', 'A/X', 'V/A'                                                                                     |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 10  | Existing_variation | Variant identifier  | 'rs17851045', 'rs786201856&COSV57325157'                                                                 |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 11  | STRAND             | Genomic strand      | 1, -1                                                                                                    |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 12  | SIFT               | SIFT prediction     | 'tolerated(0.6)', 'tolerated_low_confidence(0.37)', 'deleterious_low_confidence(0.03)', 'deleterious(0)' |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 13  | PolyPhen           | PolyPhen prediction | 'benign(0.121)', 'possibly_damaging(0.459)', 'probably_damaging(0.999)'                                  |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
| 14  | CLIN_SIG           | ClinVar prediction  | 'pathogenic', 'likely_pathogenic&pathogenic'                                                             |
+-----+--------------------+---------------------+----------------------------------------------------------------------------------------------------------+
"""

import re
import pandas as pd
import numpy as np
from . import pyvcf

# https://m.ensembl.org/info/docs/tools/vep/vep_formats.html

DATA_TYPES = {
    'Allele': 'str',
    'Consequence': 'str',
    'IMPACT': 'str',
    'SYMBOL': 'str',
    'Gene': 'str',
    'Feature_type': 'str',
    'Feature': 'str',
    'BIOTYPE': 'str',
    'EXON': 'str',
    'INTRON': 'str',
    'HGVSc': 'str',
    'HGVSp': 'str',
    'cDNA_position': 'str',
    'CDS_position': 'str',
    'Protein_position': 'str',
    'Amino_acids': 'str',
    'Codons': 'str',
    'Existing_variation': 'str',
    'DISTANCE': 'str',
    'STRAND': 'str',
    'FLAGS': 'str',
    'SYMBOL_SOURCE': 'str',
    'HGNC_ID': 'str',
    'MANE_SELECT': 'str',
    'MANE_PLUS_CLINICAL': 'str',
    'TSL': 'str',
    'APPRIS': 'str',
    'REFSEQ_MATCH': 'str',
    'REFSEQ_OFFSET': 'str',
    'GIVEN_REF': 'str',
    'USED_REF': 'str',
    'BAM_EDIT': 'str',
    'SIFT': 'str',
    'PolyPhen': 'str',
    'AF': 'float64',
    'AFR_AF': 'float64',
    'AMR_AF': 'float64',
    'EAS_AF': 'float64',
    'EUR_AF': 'float64',
    'SAS_AF': 'float64',
    'AA_AF': 'float64',
    'EA_AF': 'float64',
    'gnomAD_AF': 'float64',
    'gnomAD_AFR_AF': 'float64',
    'gnomAD_AMR_AF': 'float64',
    'gnomAD_ASJ_AF': 'float64',
    'gnomAD_EAS_AF': 'float64',
    'gnomAD_FIN_AF': 'float64',
    'gnomAD_NFE_AF': 'float64',
    'gnomAD_OTH_AF': 'float64',
    'gnomAD_SAS_AF': 'float64',
    'CLIN_SIG': 'str',
    'SOMATIC': 'str',
    'PHENO': 'str',
    'PUBMED': 'str',
    'MOTIF_NAME': 'str',
    'MOTIF_POS': 'str',
    'HIGH_INF_POS': 'str',
    'MOTIF_SCORE_CHANGE': 'str',
    'TRANSCRIPTION_FACTORS': 'str',
}

SEVERITIY = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
]

def row_firstann(r):
    """Return the first result in the row.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    str
        Annotation result. Empty string if annotation is not found.

    See Also
    --------
    row_mostsevere
        Similar method that returns result with the most severe
        consequence in the row.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
    ...     'POS': [100, 101, 200, 201],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'C', 'G', 'CCCA'],
    ...     'ALT': ['C', 'T', 'A', 'C'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': [
    ...         'CSQ=C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375733.2|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375743.4|protein_coding|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|missense_variant|MODERATE|SPEN|ENSG00000065526|Transcript|ENST00000375759.3|protein_coding|12/15||||10322|10118|3373|Q/P|cAg/cCg|||1||HGNC|17575|||||tolerated_low_confidence(0.08)|possibly_damaging(0.718)||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000440560.1|protein_coding|||||||||||4615|-1|cds_start_NF|HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000462525.1|processed_transcript|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|upstream_gene_variant|MODIFIER|SPEN|ENSG00000065526|Transcript|ENST00000487496.1|processed_transcript|||||||||||1329|1||HGNC|17575||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000537142.1|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000922256|promoter_flanking_region|||||||||||||||||||||||||||||||',
    ...         'CSQ=T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000261443.5|protein_coding||||||||||COSV99795232|3453|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000339438.6|protein_coding||||||||||COSV99795232|3071|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000358528.4|protein_coding||||||||||COSV99795232|3075|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000369530.1|protein_coding||||||||||COSV99795232|3470|-1||HGNC|29905|||||||||1|1||||||,T|missense_variant|MODERATE|NRAS|ENSG00000213281|Transcript|ENST00000369535.4|protein_coding|3/7||||502|248|83|A/D|gCc/gAc|COSV99795232||-1||HGNC|7989|||||deleterious(0.02)|probably_damaging(0.946)|||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000438362.2|protein_coding||||||||||COSV99795232|3074|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000483407.1|processed_transcript||||||||||COSV99795232|3960|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000530886.1|protein_coding||||||||||COSV99795232|3463|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000534699.1|protein_coding||||||||||COSV99795232|4115|-1||HGNC|29905|||||||||1|1||||||',
    ...         'CSQ=A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000358660.3|protein_coding|10/16||||1296|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.999)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000368196.3|protein_coding|10/16||||1375|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|possibly_damaging(0.894)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000392302.2|protein_coding|11/17||||1339|1165|389|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.983)|||1|1||||||,A|downstream_gene_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000489021.2|processed_transcript||||||||||COSV62328771|1037|1||HGNC|8031|||||||||1|1||||||,A|stop_gained&NMD_transcript_variant|HIGH|NTRK1|ENSG00000198400|Transcript|ENST00000497019.2|nonsense_mediated_decay|10/16||||1183|1032|344|W/*|tgG/tgA|COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000524377.1|protein_coding|11/17||||1314|1273|425|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.996)|||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000530298.1|retained_intron|12/17||||1313|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000534682.1|retained_intron|2/4||||496|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001787924|open_chromatin_region||||||||||COSV62328771||||||||||||||1|1||||||',
    ...         'CSQ=-|upstream_gene_variant|MODIFIER|LRRC56|ENSG00000161328|Transcript|ENST00000270115.7|protein_coding||||||||||rs1164486792|3230|1||HGNC|25430||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000311189.7|protein_coding|2/6||||198-200|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397594.1|protein_coding|1/5||||79-81|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397596.2|protein_coding|2/5||||162-164|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000417302.1|protein_coding|2/6||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000451590.1|protein_coding|2/5||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000462734.1|processed_transcript||||||||||rs1164486792|700|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000468682.2|processed_transcript|2/3||||514-516|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000478324.1|processed_transcript||||||||||rs1164486792|683|-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000479482.1|processed_transcript||||||||||rs1164486792|319|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000482021.1|processed_transcript|2/2||||149-151|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion&NMD_transcript_variant|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000493230.1|nonsense_mediated_decay|2/7||||202-204|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|downstream_gene_variant|MODIFIER|RP13-46H24.1|ENSG00000254739|Transcript|ENST00000526431.1|antisense||||||||||rs1164486792|4636|1||Clone_based_vega_gene|||||||||||1||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000035647|promoter||||||||||rs1164486792|||||||||||||||1||||||'
    ...     ],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict(meta, data)
    >>> vf.df.apply(pyvep.row_firstann, axis=1)
    0    C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG...
    1    T|downstream_gene_variant|MODIFIER|CSDE1|ENSG0...
    2    A|missense_variant|MODERATE|NTRK1|ENSG00000198...
    3    -|upstream_gene_variant|MODIFIER|LRRC56|ENSG00...
    dtype: object
    """
    results = pyvcf.row_parseinfo(r, 'CSQ')
    if not results:
        return ''
    return results.split(',')[0]

def row_mostsevere(r):
    """Return result with the most severe consequence in the row.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    str
        Annotation result. Empty string if annotation is not found.

    See Also
    --------
    row_firstann
        Similar method that returns the first result in the row.

    Examples
    --------

    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
    ...     'POS': [100, 101, 200, 201],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'C', 'G', 'CCCA'],
    ...     'ALT': ['C', 'T', 'A', 'C'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': [
    ...         'CSQ=C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375733.2|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375743.4|protein_coding|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|missense_variant|MODERATE|SPEN|ENSG00000065526|Transcript|ENST00000375759.3|protein_coding|12/15||||10322|10118|3373|Q/P|cAg/cCg|||1||HGNC|17575|||||tolerated_low_confidence(0.08)|possibly_damaging(0.718)||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000440560.1|protein_coding|||||||||||4615|-1|cds_start_NF|HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000462525.1|processed_transcript|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|upstream_gene_variant|MODIFIER|SPEN|ENSG00000065526|Transcript|ENST00000487496.1|processed_transcript|||||||||||1329|1||HGNC|17575||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000537142.1|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000922256|promoter_flanking_region|||||||||||||||||||||||||||||||',
    ...         'CSQ=T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000261443.5|protein_coding||||||||||COSV99795232|3453|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000339438.6|protein_coding||||||||||COSV99795232|3071|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000358528.4|protein_coding||||||||||COSV99795232|3075|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000369530.1|protein_coding||||||||||COSV99795232|3470|-1||HGNC|29905|||||||||1|1||||||,T|missense_variant|MODERATE|NRAS|ENSG00000213281|Transcript|ENST00000369535.4|protein_coding|3/7||||502|248|83|A/D|gCc/gAc|COSV99795232||-1||HGNC|7989|||||deleterious(0.02)|probably_damaging(0.946)|||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000438362.2|protein_coding||||||||||COSV99795232|3074|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000483407.1|processed_transcript||||||||||COSV99795232|3960|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000530886.1|protein_coding||||||||||COSV99795232|3463|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000534699.1|protein_coding||||||||||COSV99795232|4115|-1||HGNC|29905|||||||||1|1||||||',
    ...         'CSQ=A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000358660.3|protein_coding|10/16||||1296|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.999)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000368196.3|protein_coding|10/16||||1375|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|possibly_damaging(0.894)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000392302.2|protein_coding|11/17||||1339|1165|389|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.983)|||1|1||||||,A|downstream_gene_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000489021.2|processed_transcript||||||||||COSV62328771|1037|1||HGNC|8031|||||||||1|1||||||,A|stop_gained&NMD_transcript_variant|HIGH|NTRK1|ENSG00000198400|Transcript|ENST00000497019.2|nonsense_mediated_decay|10/16||||1183|1032|344|W/*|tgG/tgA|COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000524377.1|protein_coding|11/17||||1314|1273|425|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.996)|||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000530298.1|retained_intron|12/17||||1313|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000534682.1|retained_intron|2/4||||496|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001787924|open_chromatin_region||||||||||COSV62328771||||||||||||||1|1||||||',
    ...         'CSQ=-|upstream_gene_variant|MODIFIER|LRRC56|ENSG00000161328|Transcript|ENST00000270115.7|protein_coding||||||||||rs1164486792|3230|1||HGNC|25430||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000311189.7|protein_coding|2/6||||198-200|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397594.1|protein_coding|1/5||||79-81|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397596.2|protein_coding|2/5||||162-164|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000417302.1|protein_coding|2/6||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000451590.1|protein_coding|2/5||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000462734.1|processed_transcript||||||||||rs1164486792|700|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000468682.2|processed_transcript|2/3||||514-516|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000478324.1|processed_transcript||||||||||rs1164486792|683|-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000479482.1|processed_transcript||||||||||rs1164486792|319|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000482021.1|processed_transcript|2/2||||149-151|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion&NMD_transcript_variant|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000493230.1|nonsense_mediated_decay|2/7||||202-204|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|downstream_gene_variant|MODIFIER|RP13-46H24.1|ENSG00000254739|Transcript|ENST00000526431.1|antisense||||||||||rs1164486792|4636|1||Clone_based_vega_gene|||||||||||1||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000035647|promoter||||||||||rs1164486792|||||||||||||||1||||||'
    ...     ],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict(meta, data)
    >>> vf.df.apply(pyvep.row_mostsevere, axis=1)
    0    C|missense_variant|MODERATE|SPEN|ENSG000000655...
    1    T|missense_variant|MODERATE|NRAS|ENSG000002132...
    2    A|stop_gained&NMD_transcript_variant|HIGH|NTRK...
    3    -|inframe_deletion|MODERATE|HRAS|ENSG000001747...
    dtype: object
    """
    results = pyvcf.row_parseinfo(r, 'CSQ')
    if not results:
        return ''
    f = lambda x: SEVERITIY.index(x.split('|')[1].split('&')[0])
    return sorted(results.split(','), key=f)[0]

def parseann(vf, targets, sep=' | ', as_series=False):
    """Parse variant annotations in VcfFrame.

    If there are multiple results per row, the method will use the first
    one for extracting requested data.

    Parameters
    ----------
    vf : VcfFrame
        VcfFrame.
    targets : list
        List of subfield IDs or indicies or both.
    sep : str, default: ' | '
        Separator for joining requested annotations.
    as_series : bool, default: False
        If True, return 1D array instead of VcfFrame.

    Returns
    -------
    VcfFrame or pandas.Series
        Parsed annotations in VcfFrame or 1D array.

    Examples
    --------
    Assume we have the following data:

    >>> meta = [
    ...     '##fileformat=VCFv4.1',
    ...     '##VEP="v104" time="2021-05-20 10:50:12" cache="/net/isilonP/public/ro/ensweb-data/latest/tools/grch37/e104/vep/cache/homo_sapiens/104_GRCh37" db="homo_sapiens_core_104_37@hh-mysql-ens-grch37-web" 1000genomes="phase3" COSMIC="92" ClinVar="202012" HGMD-PUBLIC="20204" assembly="GRCh37.p13" dbSNP="154" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"',
    ...     '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|SIFT|PolyPhen|AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">'
    ... ]
    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
    ...     'POS': [100, 101, 200, 201],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'C', 'CAG', 'C'],
    ...     'ALT': ['T', 'T', 'C', 'T'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': [
    ...         'CSQ=T|missense_variant|MODERATE|MTOR|ENSG00000198793|Transcript|ENST00000361445.4|protein_coding|47/58||||6721|6644|2215|S/Y|tCt/tAt|rs587777894&COSV63868278&COSV63868313||-1||HGNC|3942|||||deleterious(0)|possibly_damaging(0.876)||likely_pathogenic&pathogenic|0&1&1|1&1&1|26619011&27159400&24631838&26018084&27830187|||||',
    ...         'CSQ=T|synonymous_variant|LOW|MTOR|ENSG00000198793|Transcript|ENST00000361445.4|protein_coding|49/58||||6986|6909|2303|L|ctG/ctA|rs11121691&COSV63870864||-1||HGNC|3942|||||||0.2206|benign|0&1|1&1|24996771|||||',
    ...         'CSQ=-|frameshift_variant|HIGH|BRCA2|ENSG00000139618|Transcript|ENST00000380152.3|protein_coding|18/27||||8479-8480|8246-8247|2749|Q/X|cAG/c|rs80359701||1||HGNC|1101||||||||pathogenic||1|26467025&26295337&15340362|||||',
    ...         'CSQ=T|missense_variant|MODERATE|MTOR|ENSG00000198793|Transcript|ENST00000361445.4|protein_coding|30/58||||4516|4439|1480|R/H|cGc/cAc|rs780930764&COSV63868373||-1||HGNC|3942|||||tolerated(0.13)|benign(0)||likely_benign|0&1|1&1||||||'
    ...     ],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict(meta, data)

    We can extract gene symbols:

    >>> pyvep.parseann(vf, ['SYMBOL']).df
      CHROM  POS ID  REF ALT QUAL FILTER   INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .   MTOR     GT    0/1
    1  chr1  101  .    C   T    .      .   MTOR     GT    0/1
    2  chr2  200  .  CAG   C    .      .  BRCA2     GT    0/1
    3  chr2  201  .    C   T    .      .   MTOR     GT    0/1

    If we want to add variant impact:

    >>> pyvep.parseann(vf, ['SYMBOL', 'IMPACT']).df
      CHROM  POS ID  REF ALT QUAL FILTER             INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .  MTOR | MODERATE     GT    0/1
    1  chr1  101  .    C   T    .      .       MTOR | LOW     GT    0/1
    2  chr2  200  .  CAG   C    .      .     BRCA2 | HIGH     GT    0/1
    3  chr2  201  .    C   T    .      .  MTOR | MODERATE     GT    0/1

    If we want to change the delimter:

    >>> pyvep.parseann(vf, ['SYMBOL', 'IMPACT'], sep=',').df
      CHROM  POS ID  REF ALT QUAL FILTER           INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .  MTOR,MODERATE     GT    0/1
    1  chr1  101  .    C   T    .      .       MTOR,LOW     GT    0/1
    2  chr2  200  .  CAG   C    .      .     BRCA2,HIGH     GT    0/1
    3  chr2  201  .    C   T    .      .  MTOR,MODERATE     GT    0/1

    Finally, we can return a 1D array instead of VcfFrame:

    >>> pyvep.parseann(vf, ['SYMBOL'], as_series=True)
    0     MTOR
    1     MTOR
    2    BRCA2
    3     MTOR
    dtype: object
    """
    _targets = [x if isinstance(x, int) else annot_names(vf).index(x) for x in targets]
    def func(r):
        ann = row_firstann(r)
        if not ann:
            return '.'
        ann = ann.split('|')
        ann = sep.join([ann[i] for i in _targets])
        return ann
    s = vf.df.apply(func, axis=1)
    if as_series:
        return s
    new_vf = vf.copy()
    new_vf.df.INFO = s
    return new_vf

def to_frame(vf, as_zero=False):
    """
    Create a DataFrame containing analysis-ready VEP annotation data.

    Parameters
    ----------
    vf : VcfFrame
        VcfFrame object.
    as_zero : bool, default: False
        If True, missing values will be converted to zeros instead of ``NaN``.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing VEP annotation data.
    """
    f = lambda r: pd.Series(row_firstann(r).split('|'))
    df = vf.df.apply(f, axis=1)
    if as_zero:
        df = df.replace('', 0)
    else:
        df = df.replace('', np.nan)
    df.columns = annot_names(vf)
    d = {}
    for col in df.columns:
        if col in DATA_TYPES:
            d[col] = DATA_TYPES[col]
    df = df.astype(d)
    return df

def pick_result(vf, mode='mostsevere'):
    """Return a new VcfFrame after picking one result per row.

    Parameters
    ----------
    vf : VcfFrame
        VcfFrame object.
    mode : {'mostsevere', 'firstann'}, default: 'mostsevere'
        Selection mode:

            - mostsevere: use ``row_mostsevere`` method
            - firstann: use ``row_firstann`` method

    Returns
    -------
    VcfFrame
        New VcfFrame.

    Examples
    --------

    >>> meta = [
    ...     '##fileformat=VCFv4.1',
    ...     '##VEP="v104" time="2021-05-20 10:50:12" cache="/net/isilonP/public/ro/ensweb-data/latest/tools/grch37/e104/vep/cache/homo_sapiens/104_GRCh37" db="homo_sapiens_core_104_37@hh-mysql-ens-grch37-web" 1000genomes="phase3" COSMIC="92" ClinVar="202012" HGMD-PUBLIC="20204" assembly="GRCh37.p13" dbSNP="154" gencode="GENCODE 19" genebuild="2011-04" gnomAD="r2.1" polyphen="2.2.2" regbuild="1.0" sift="sift5.2.2"',
    ...     '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|SIFT|PolyPhen|AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">'
    ... ]
    >>> data = {
    ...     'CHROM': ['chr1', 'chr1', 'chr2', 'chr2'],
    ...     'POS': [100, 101, 200, 201],
    ...     'ID': ['.', '.', '.', '.'],
    ...     'REF': ['G', 'C', 'G', 'CCCA'],
    ...     'ALT': ['C', 'T', 'A', 'C'],
    ...     'QUAL': ['.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.'],
    ...     'INFO': [
    ...         'CSQ=C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375733.2|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000375743.4|protein_coding|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|missense_variant|MODERATE|SPEN|ENSG00000065526|Transcript|ENST00000375759.3|protein_coding|12/15||||10322|10118|3373|Q/P|cAg/cCg|||1||HGNC|17575|||||tolerated_low_confidence(0.08)|possibly_damaging(0.718)||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000440560.1|protein_coding|||||||||||4615|-1|cds_start_NF|HGNC|12936||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000462525.1|processed_transcript|||||||||||4615|-1||HGNC|12936||||||||||||||||,C|upstream_gene_variant|MODIFIER|SPEN|ENSG00000065526|Transcript|ENST00000487496.1|processed_transcript|||||||||||1329|1||HGNC|17575||||||||||||||||,C|downstream_gene_variant|MODIFIER|ZBTB17|ENSG00000116809|Transcript|ENST00000537142.1|protein_coding|||||||||||4617|-1||HGNC|12936||||||||||||||||,C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000922256|promoter_flanking_region|||||||||||||||||||||||||||||||',
    ...         'CSQ=T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000261443.5|protein_coding||||||||||COSV99795232|3453|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000339438.6|protein_coding||||||||||COSV99795232|3071|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000358528.4|protein_coding||||||||||COSV99795232|3075|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000369530.1|protein_coding||||||||||COSV99795232|3470|-1||HGNC|29905|||||||||1|1||||||,T|missense_variant|MODERATE|NRAS|ENSG00000213281|Transcript|ENST00000369535.4|protein_coding|3/7||||502|248|83|A/D|gCc/gAc|COSV99795232||-1||HGNC|7989|||||deleterious(0.02)|probably_damaging(0.946)|||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000438362.2|protein_coding||||||||||COSV99795232|3074|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000483407.1|processed_transcript||||||||||COSV99795232|3960|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000530886.1|protein_coding||||||||||COSV99795232|3463|-1||HGNC|29905|||||||||1|1||||||,T|downstream_gene_variant|MODIFIER|CSDE1|ENSG00000009307|Transcript|ENST00000534699.1|protein_coding||||||||||COSV99795232|4115|-1||HGNC|29905|||||||||1|1||||||',
    ...         'CSQ=A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000358660.3|protein_coding|10/16||||1296|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.999)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000368196.3|protein_coding|10/16||||1375|1255|419|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|possibly_damaging(0.894)|||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000392302.2|protein_coding|11/17||||1339|1165|389|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.983)|||1|1||||||,A|downstream_gene_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000489021.2|processed_transcript||||||||||COSV62328771|1037|1||HGNC|8031|||||||||1|1||||||,A|stop_gained&NMD_transcript_variant|HIGH|NTRK1|ENSG00000198400|Transcript|ENST00000497019.2|nonsense_mediated_decay|10/16||||1183|1032|344|W/*|tgG/tgA|COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|missense_variant|MODERATE|NTRK1|ENSG00000198400|Transcript|ENST00000524377.1|protein_coding|11/17||||1314|1273|425|A/T|Gcc/Acc|COSV62328771||1||HGNC|8031|||||deleterious(0.01)|probably_damaging(0.996)|||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000530298.1|retained_intron|12/17||||1313|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|non_coding_transcript_exon_variant|MODIFIER|NTRK1|ENSG00000198400|Transcript|ENST00000534682.1|retained_intron|2/4||||496|||||COSV62328771||1||HGNC|8031|||||||||1|1||||||,A|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00001787924|open_chromatin_region||||||||||COSV62328771||||||||||||||1|1||||||',
    ...         'CSQ=-|upstream_gene_variant|MODIFIER|LRRC56|ENSG00000161328|Transcript|ENST00000270115.7|protein_coding||||||||||rs1164486792|3230|1||HGNC|25430||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000311189.7|protein_coding|2/6||||198-200|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397594.1|protein_coding|1/5||||79-81|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000397596.2|protein_coding|2/5||||162-164|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000417302.1|protein_coding|2/6||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000451590.1|protein_coding|2/5||||214-216|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000462734.1|processed_transcript||||||||||rs1164486792|700|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000468682.2|processed_transcript|2/3||||514-516|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000478324.1|processed_transcript||||||||||rs1164486792|683|-1||HGNC|5173||||||||||1||||||,-|upstream_gene_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000479482.1|processed_transcript||||||||||rs1164486792|319|-1||HGNC|5173||||||||||1||||||,-|non_coding_transcript_exon_variant|MODIFIER|HRAS|ENSG00000174775|Transcript|ENST00000482021.1|processed_transcript|2/2||||149-151|||||rs1164486792||-1||HGNC|5173||||||||||1||||||,-|inframe_deletion&NMD_transcript_variant|MODERATE|HRAS|ENSG00000174775|Transcript|ENST00000493230.1|nonsense_mediated_decay|2/7||||202-204|26-28|9-10|VG/G|gTGGgc/ggc|rs1164486792||-1||HGNC|5173||||||||||1||||||,-|downstream_gene_variant|MODIFIER|RP13-46H24.1|ENSG00000254739|Transcript|ENST00000526431.1|antisense||||||||||rs1164486792|4636|1||Clone_based_vega_gene|||||||||||1||||||,-|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000035647|promoter||||||||||rs1164486792|||||||||||||||1||||||'
    ...     ],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT'],
    ...     'Steven': ['0/1', '0/1', '0/1', '0/1'],
    ... }
    >>> vf = pyvcf.VcfFrame.from_dict(meta, data)
    >>> vf1 = pyvep.pick_result(vf, mode='mostsevere')
    >>> vf2 = pyvep.pick_result(vf, mode='firstann')
    >>> vf1.df.INFO
    0    CSQ=C|missense_variant|MODERATE|SPEN|ENSG00000...
    1    CSQ=T|missense_variant|MODERATE|NRAS|ENSG00000...
    2    CSQ=A|stop_gained&NMD_transcript_variant|HIGH|...
    3    CSQ=-|inframe_deletion|MODERATE|HRAS|ENSG00000...
    Name: INFO, dtype: object
    >>> vf2.df.INFO
    0    CSQ=C|downstream_gene_variant|MODIFIER|ZBTB17|...
    1    CSQ=T|downstream_gene_variant|MODIFIER|CSDE1|E...
    2    CSQ=A|missense_variant|MODERATE|NTRK1|ENSG0000...
    3    CSQ=-|upstream_gene_variant|MODIFIER|LRRC56|EN...
    Name: INFO, dtype: object
    """
    funcs = {'mostsevere': row_mostsevere, 'firstann': row_firstann}
    new_vf = vf.copy()
    one_row = lambda r: pyvcf.row_updateinfo(r, 'CSQ', funcs[mode](r))
    new_vf.df.INFO = vf.df.apply(one_row, axis=1)
    return new_vf

def annot_names(vf):
    """Return the list of avaialble consequence annotations in the VcfFrame.

    Parameters
    ----------
    vf : VcfFrame
        VcfFrame.

    Returns
    -------
    list
        List of consequence annotations.
    """
    l = []
    for i, line in enumerate(vf.meta):
        if 'ID=CSQ' in line:
            l = re.search(r'Format: (.*?)">', vf.meta[i]).group(1).split('|')
    return l

def filter_query(vf, expr, opposite=None, as_index=False, as_zero=False):
    """
    Select rows that satisfy the query expression.

    This method essentially wraps the :meth:`pandas.DataFrame.query` method.

    Parameters
    ----------
    expr : str
        The query string to evaluate.
    opposite : bool, default: False
        If True, return rows that don't meet the said criteria.
    as_index : bool, default: False
        If True, return boolean index array instead of VcfFrame.
    as_zero : bool, default: False
        If True, missing values will be converted to zeros instead of ``NaN``.

    Returns
    -------
    VcfFrame or pandas.Series
        Filtered VcfFrame or boolean index array.

    Examples
    --------

    >>> from fuc import common, pyvcf, pyvep
    >>> common.load_dataset('tcga-laml')
    >>> fn = '~/fuc-data/tcga-laml/tcga_laml_vep.vcf'
    >>> vf = pyvcf.VcfFrame.from_file(fn)
    >>> filtered_vf = pyvep.filter_query(vf, "SYMBOL == 'TP53'")
    >>> filtered_vf = pyvep.filter_query(vf, "SYMBOL != 'TP53'")
    >>> filtered_vf = pyvep.filter_query(vf, "SYMBOL != 'TP53'", opposite=True)
    >>> filtered_vf = pyvep.filter_query(vf, "Consequence in ['splice_donor_variant', 'stop_gained']")
    >>> filtered_vf = pyvep.filter_query(vf, "(SYMBOL == 'TP53') and (Consequence.str.contains('stop_gained'))")
    >>> filtered_vf = pyvep.filter_query(vf, "gnomAD_AF < 0.001")
    >>> filtered_vf = pyvep.filter_query(vf, "gnomAD_AF < 0.001", as_zero=True)
    >>> filtered_vf = pyvep.filter_query(vf, "(IMPACT == 'HIGH') or (Consequence.str.contains('missense_variant') and (PolyPhen.str.contains('damaging') or SIFT.str.contains('deleterious')))")
    """
    df = to_frame(vf, as_zero=as_zero)
    i = vf.df.index.isin(df.query(expr).index)
    if opposite:
        i = ~i
    if as_index:
        return i
    return vf.__class__(vf.copy_meta(), vf.df[i])
