"""
The pysnpeff submodule is designed for parsing VCF annotation data from
the `SnpEff <https://pcingola.github.io/SnpEff/>`_ program. It should be
used with ``pyvcf.VcfFrame``.

One VCF record can have several SnpEff annotations if, for example,
the record is a multiallelic site or the variant is shared by
multiple genes. When more than one annotations are reported, SnpEff
will sort them by their importance. For more details, visit the `official
website <https://pcingola.github.io/SnpEff/>`_.

For each annotation, SnpEff provides the following data:

1. Allele - ALT allele.
2. Annotation - Sequence Ontology terms concatenated using '&'.
3. Annotation_Impact - HIGH, MODERATE, LOW, or MODIFIER.
4. Gene_Name - Common gene name (HGNC).
5. Gene_ID - Gene ID.
6. Feature_Type - Which type of feature is in the next field.
7. Feature_ID - Transcript ID, Motif ID, miRNA, ChipSeq peak, etc.
8. Transcript_BioType - Coding or noncoding.
9. Rank - Exon or Intron rank / total number of exons or introns.
10. HGVS.c - Variant using HGVS notation (DNA level).
11. HGVS.p - Variant using HGVS notation (Protein level).
12. cDNA.pos / cDNA.length - Position in cDNA and trancript's cDNA length.
13. CDS.pos / CDS.length - Position and number of coding bases.
14. AA.pos / AA.length - Position and number of AA.
15. Distance - All items in this field are options.
16. ERRORS / WARNINGS - Messages that can affect annotation accuracy.
17. INFO - Additional information.
"""

def row_firstann(r):
    """Return the first SnpEff annotation for the row."""
    ann = [x for x in r.INFO.split(';') if 'ANN=' in x]
    if not ann:
        return ''
    ann = ann[0].replace('ANN=', '').split(',')[0]
    return ann

def filter_ann(vf, targets, include=True):
    """Filter out rows based on the SnpEff annotations.

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

def parseann(vf, idx, sep=' | '):
    """Parse SnpEff annotations.

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
        ann = row_firstann(r)
        if not ann:
            return '.'
        ann = ann.split('|')
        ann = sep.join([ann[i] for i in idx])
        return ann
    s = vf.df.apply(func, axis=1)
    return s
