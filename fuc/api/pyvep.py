"""
The pyvep submodule is designed for parsing VCF annotation data from
the `Ensembl Variant Effect Predictor (VEP)
<https://asia.ensembl.org/info/docs/tools/vep/index.html>`_. It is
designed to be used with ``pyvcf.VcfFrame``.
"""

import re
import pandas as pd
from . import pyvcf

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

def _get_keys(vf):
    """Return existing annotation keys (e.g. Allele, IMPACT)."""
    l = []
    for i, line in enumerate(vf.meta):
        if 'ID=CSQ' in line:
            l = re.search(r'Format: (.*?)">', vf.meta[i]).group(1).split('|')
    return l

def row_firstann(r):
    """Return the first result in the row.

    Parameters
    ----------
    r : pandas.Series
        VCF row.

    Returns
    -------
    str
        VEP result.
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
        VEP result.
    """
    results = pyvcf.row_parseinfo(r, 'CSQ')
    if not results:
        return ''
    f = lambda x: SEVERITIY.index(x.split('|')[1].split('&')[0])
    return sorted(results.split(','), key=f)[0]

def filter_nothas(vf, s, include=False):
    """Filter out rows based on the VEP annotations.

    Parameters
    ----------
    vf : fuc.api.pyvcf.VcfFrame
        Input VcfFrame.
    targets : list
        List of annotations (e.g. ['missense_variant', 'stop_gained']).
    include : bool, default: False
        If True, include only such rows instead of excluding them.

    Returns
    -------
    vf : VcfFrame
        Filtered VcfFrame.
    """
    def func(r):
        ann = row_firstann(r)
        return s not in ann
    i = vf.df.apply(func, axis=1)
    if include:
        i = ~i
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def filter_clinsig(vf, whitelist=None, blacklist=None, opposite=False,
                   index=False):
    """Select rows based on the given CLIN_SIG values.

    List of CLIN_SIG values:

        - benign
        - likely_benign
        - pathogenic
        - likely_pathogenic
        - drug_response
        - risk_factor
        - uncertain_significance
        - conflicting_interpretations_of_pathogenicity
        - not_provided
        - other
        - benign/likely_benign
        - pathogenic/likely_pathogenic
        - ...

    Parameters
    ----------
    whilelist : list, default: None
        If these CLIN_SIG values are present, select the row.
    blacklist : list, default: None
        If these CLIN_SIG values are present, do not select the row.
    opposite : bool, default: False
        If True, return rows that don't meet the said criteria.
    index : bool, default: False
        If True, return boolean index array instead of VcfFrame.

    Returns
    -------
    VcfFrame or pandas.Series
        Filtered VcfFrame or boolean index array.

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
    >>> pyvep.parseann(vf, ['CLIN_SIG']).df
      CHROM  POS ID  REF ALT QUAL FILTER                          INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .  likely_pathogenic&pathogenic     GT    0/1
    1  chr1  101  .    C   T    .      .                        benign     GT    0/1
    2  chr2  200  .  CAG   C    .      .                    pathogenic     GT    0/1
    3  chr2  201  .    C   T    .      .                 likely_benign     GT    0/1

    We can select rows with pathogenic or likely_pathogenic:

    >>> whitelist=['pathogenic', 'likely_pathogenic']
    >>> temp_vf = pyvep.filter_clinsig(vf, whitelist=whitelist)
    >>> pyvep.parseann(temp_vf, ['CLIN_SIG']).df
      CHROM  POS ID  REF ALT QUAL FILTER                          INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .  likely_pathogenic&pathogenic     GT    0/1
    1  chr2  200  .  CAG   C    .      .                    pathogenic     GT    0/1

    We can also remove those rows:

    >>> temp_vf = pyvep.filter_clinsig(vf, whitelist=whitelist, opposite=True)
    >>> pyvep.parseann(temp_vf, ['CLIN_SIG']).df
      CHROM  POS ID REF ALT QUAL FILTER           INFO FORMAT Steven
    0  chr1  101  .   C   T    .      .         benign     GT    0/1
    1  chr2  201  .   C   T    .      .  likely_benign     GT    0/1

    Finally, we can return boolean index array from the filtering:

    >>> pyvep.filter_clinsig(vf, whitelist=whitelist, index=True)
    0     True
    1    False
    2     True
    3    False
    dtype: bool
    """
    if whitelist is None and blacklist is None:
        raise ValueError('must provide either whitelist or blcklist')
    if whitelist is None:
        whitelist = []
    if blacklist is None:
        blacklist = []
    def func(r):
        ann = row_firstann(r)
        values = ann.split('|')[get_index(vf, 'CLIN_SIG')].split('&')
        if not list(set(values) & set(whitelist)):
            return False
        if list(set(values) & set(blacklist)):
            return False
        return True
    i = vf.df.apply(func, axis=1)
    if opposite:
        i = ~i
    if index:
        return i
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def filter_impact(vf, values, opposite=False, index=False):
    """Select rows based on the given IMPACT values.

    List of IMPACT values:

        - LOW
        - MODERATE
        - HIGH
        - MODIFIER

    Parameters
    ----------
    values : list, default: None
        If any one of the IMPACT values is present, select the row.
    opposite : bool, default: False
        If True, return rows that don't meet the said criteria.
    index : bool, default: False
        If True, return boolean index array instead of VcfFrame.

    Returns
    -------
    VcfFrame or pandas.Series
        Filtered VcfFrame or boolean index array.

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
    >>> pyvep.parseann(vf, ['IMPACT']).df
      CHROM  POS ID  REF ALT QUAL FILTER      INFO FORMAT Steven
    0  chr1  100  .    G   T    .      .  MODERATE     GT    0/1
    1  chr1  101  .    C   T    .      .       LOW     GT    0/1
    2  chr2  200  .  CAG   C    .      .      HIGH     GT    0/1
    3  chr2  201  .    C   T    .      .  MODERATE     GT    0/1

    We can select rows with either LOW or HIGH:

    >>> temp_vf = pyvep.filter_impact(vf, ['LOW', 'HIGH'])
    >>> pyvep.parseann(temp_vf, ['IMPACT']).df
      CHROM  POS ID REF ALT QUAL FILTER      INFO FORMAT Steven
    0  chr1  100  .   G   T    .      .  MODERATE     GT    0/1
    1  chr2  201  .   C   T    .      .  MODERATE     GT    0/1

    We can also remove those rows:

    >>> temp_vf = pyvep.filter_impact(vf, ['LOW', 'HIGH'], opposite=True)
    >>> pyvep.parseann(temp_vf, ['IMPACT']).df
      CHROM  POS ID REF ALT QUAL FILTER      INFO FORMAT Steven
    0  chr1  100  .   G   T    .      .  MODERATE     GT    0/1
    1  chr2  201  .   C   T    .      .  MODERATE     GT    0/1

    Finally, we can return boolean index array from the filtering:

    >>> pyvep.filter_impact(vf, ['LOW', 'HIGH'], index=True)
    0     True
    1    False
    2    False
    3     True
    dtype: bool

    """
    def func(r):
        ann = row_firstann(r)
        impact = ann.split('|')[get_index(vf, 'IMPACT')]
        return impact in values
    i = vf.df.apply(func, axis=1)
    if opposite:
        i = ~i
    if index:
        return i
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def parseann(vf, targets, sep=' | ', as_series=False):
    """Parse VEP annotations in VcfFrame.

    Parameters
    ----------
    vf : fuc.pyvcf.VcfFrame
        VcfFrame containing VEP annotations.
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
    """
    _targets = [x if isinstance(x, int) else _get_keys(vf).index(x) for x in targets]
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

def get_index(vf, target):
    """Return the index of the target field (e.g. CLIN_SIG)."""
    headers = _get_keys(vf)
    return headers.index(target)

def get_table(vf):
    """Write the VcfFrame as a tab-delimited text file."""
    df = vf.df.copy()
    headers = _get_keys(vf)
    def func(r):
        ann = row_firstann(r)
        if ann:
            s = ['.' if x == '' else x for x in ann.split('|')]
        else:
            s = ['.' for x in headers]
        s = pd.Series(s, index=headers)
        s = pd.concat([r[:9], s, r[9:]])
        return s
    df = df.apply(func, axis=1)
    return df

def write_table(vf, fn):
    """Write the VcfFrame as a tab-delimited text file."""
    df = get_table(vf)
    df.to_csv(fn, sep='\t', index=False)
