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
    i = 0
    for j, line in enumerate(vf.meta):
        if 'ID=CSQ' in line:
            i = j
    if not i:
        return 0
    s = re.search(r'Format: (.*?)">', vf.meta[i]).group(1)
    return s.split('|').index(target)
