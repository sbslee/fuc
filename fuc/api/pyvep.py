"""
The ``vep`` submodule is designed for parsing VCF annotation data from
the Ensembl Variant Effect Predictor (VEP). It should be used with
``pyvcf.VcfFrame``.

1. Allele
2. Consequence
3. IMPACT
4. SYMBOL
5. Gene
6. Feature_type
7. Feature
8. BIOTYPE
9. EXON
10. INTRON
11. HGVSc
12. HGVSp
13. cDNA_position
14. CDS_position
15. Protein_position
16. Amino_acids
17. Codons
18. Existing_variation
19. DISTANCE
20. STRAND
21. FLAGS
22. SYMBOL_SOURCE
23. HGNC_ID
24. MANE_SELECT
25. MANE_PLUS_CLINICAL
26. TSL
27. APPRIS
28. SIFT
29. PolyPhen
30. AF
31. CLIN_SIG
32. SOMATIC
33. PHENO
34. PUBMED
35. MOTIF_NAME
36. MOTIF_POS
37. HIGH_INF_POS
38. MOTIF_SCORE_CHANGE
39. TRANSCRIPTION_FACTORS
"""

import re
import pandas as pd

def row_first_ann(r):
    """Return the first VEP annotation for the row."""
    ann = [x for x in r.INFO.split(';') if 'CSQ=' in x]
    if not ann:
        return ''
    ann = ann[0].replace('CSQ=', '').split(',')[0]
    return ann

def filter_ann(vf, targets, include=True):
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
        ann = row_first_ann(r)
        if not ann:
            return False
        ann = ann.split('|')[1]
        has_ann = ann in targets
        if include:
            return has_ann
        else:
            return not has_ann
    i = vf.df.apply(func, axis=1)
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def filter_clinsig(vf, whitelist=None, blacklist=None, include=False,
                   index=False):
    """Filter out rows whose variant is deemed as not clincally significant.

    Parameters
    ----------
    whilelist : list, default: None
        CLIN_SIG values that signifiy a variant is clincally significant.
        By default, it includes ``pathogenic`` and ``likely_pathogenic``
    blacklist : list, default: None
        CLIN_SIG values that signifiy a variant is not clincally significant.
        By default, it includes ``benign`` and ``likely_benign``.
    include : bool, default: False
        If True, include only such rows instead of excluding them.
    index : bool, default: False
        If True, return the boolean index instead of a new VcfFrame.

    Returns
    -------
    vf : VcfFrame or pandas.Series
        Filtered VcfFrame or boolean index.
    """
    if whitelist is None:
        whitelist = ['pathogenic', 'likely_pathogenic']
    if blacklist is None:
        blacklist = ['benign', 'likely_benign']
    def func(r):
        ann = row_first_ann(r)
        values = ann.split('|')[get_index(vf, 'CLIN_SIG')].split('&')
        if not list(set(values) & set(whitelist)):
            return False
        if list(set(values) & set(blacklist)):
            return False
        return True
    i = vf.df.apply(func, axis=1)
    if include:
        i = ~i
    if index:
        return i
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def filter_impact(vf, values, include=False, index=False):
    """Filter out rows based on the IMPACT field."""
    def func(r):
        ann = row_first_ann(r)
        impact = ann.split('|')[get_index(vf, 'IMPACT')]
        return impact not in values
    i = vf.df.apply(func, axis=1)
    if include:
        i = ~i
    if index:
        return i
    df = vf.df[i].reset_index(drop=True)
    vf = vf.__class__(vf.copy_meta(), df)
    return vf

def parse_ann(vf, idx, sep=' | '):
    """Parse VEP annotations.

    Parameters
    ----------
    vf : fuc.api.pyvcf.VcfFrame
        Input VcfFrame.
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
        ann = row_first_ann(r)
        if not ann:
            return '.'
        ann = ann.split('|')
        ann = sep.join([ann[i] for i in idx])
        return ann
    s = vf.df.apply(func, axis=1)
    return s

def get_index(vf, target):
    """Return the index of the target field (e.g. CLIN_SIG)."""
    headers = get_headers(vf)
    return headers.index(target)

def get_headers(vf):
    """Return the list of field IDs."""
    headers = []
    for i, line in enumerate(vf.meta):
        if 'ID=CSQ' in line:
            headers = re.search(r'Format: (.*?)">',
                                vf.meta[i]).group(1).split('|')
    return headers

def get_table(vf):
    """Write the VcfFrame as a tab-delimited text file."""
    df = vf.df.copy()
    headers = get_headers(vf)
    def func(r):
        ann = row_first_ann(r)
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
