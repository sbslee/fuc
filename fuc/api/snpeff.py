"""
One VCF record can have several SnpEff annotations if, for example,
the record is a multiallelic site or the variant is shared by
multiple genes. When more than one annotations are reported, SnpEff
will sort them by their importance. For more details, visit the official
website (https://pcingola.github.io/SnpEff/).

For each annotation, SnpEff provides the following data:

1. Allele - ALT allele.
2. Annotation - Annotated using Sequence Ontology terms. Multiple effects can be concatenated using ``&`` (e.g. ``intron_variant&nc_transcript_variant``).
3. Annotation_Impact - A simple estimation of putative impact or deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}.
4. Gene_Name - Common gene name (HGNC). Optional: use closest gene when the variant is 'intergenic'.
5. Gene_ID - Gene ID.
6. Feature_Type - Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.). It is preferred to use Sequence Ontology (SO) terms, but 'custom' (user defined) are allowed.
7. Feature_ID - Depending on the annotation, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID).
8. Transcript_BioType - The bare minimum is at least a description on whether the transcript is {"Coding", "Noncoding"}. Whenever possible, use ENSEMBL biotypes.
9. Rank - Exon or Intron rank / total number of exons or introns.
10. HGVS.c - Variant using HGVS notation (DNA level).
11. HGVS.p - If variant is coding, this field describes the variant using HGVS notation (Protein level). Since transcript ID is already mentioned in 'feature ID', it may be omitted here.
12. cDNA.pos / cDNA.length - Position in cDNA and trancript's cDNA length (one based).
13. CDS.pos / CDS.length - Position and number of coding bases (one based includes START and STOP codons).
14. AA.pos / AA.length - Position and number of AA (one based, including START, but not STOP).
15. Distance - All items in this field are options, so the field could be empty.
16. ERRORS / WARNINGS - Add errors, warnings or informative message that can affect annotation accuracy.
17. INFO - Additional information.
"""

def row_first_ann(r):
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

def parse_snpeff(vf, idx, sep=' | '):
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
        ann = row_first_ann(r)
        if not ann:
            return '.'
        ann = ann.split('|')
        ann = sep.join([ann[i] for i in idx])
        return ann
    s = vf.df.apply(func, axis=1)
    return s
