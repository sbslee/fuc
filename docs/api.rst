API
***

Introduction
============

This section describes application programming interface (API) for the ``fuc`` package.

Below is the list of submodules available in API:

- **common** : The ``common`` submodule is used by other ``fuc`` submodules such as ``pyvcf`` and ``pybed``. It also provides many useful methods.
- **pybed** : The ``pybed`` submodule is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **pyfq** : The ``pyfq`` submodule is designed for working with FASTQ files (both zipped and unzipped).
- **pyvcf** : The ``pyvcf`` submodule is designed for working with VCF files (both zipped and unzipped). It implements `pyvcf.VcfFrame` which stores VCF data as a ``pandas.DataFrame`` to allow fast computation and easy minipulation.
- **snpeff** : One VCF record can have several SnpEff annotations if, for example, the record is a multiallelic site or the variant is shared by multiple genes. When more than one annotations are reported, SnpEff will sort them by their importance. For more details, visit the official website (https://pcingola.github.io/SnpEff/).  For each annotation, SnpEff provides the following data:  1. Allele - ALT allele. 2. Annotation - Annotated using Sequence Ontology terms. Multiple effects can be concatenated using ``&`` (e.g. ``intron_variant&nc_transcript_variant``). 3. Annotation_Impact - A simple estimation of putative impact or deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}. 4. Gene_Name - Common gene name (HGNC). Optional: use closest gene when the variant is 'intergenic'. 5. Gene_ID - Gene ID. 6. Feature_Type - Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.). It is preferred to use Sequence Ontology (SO) terms, but 'custom' (user defined) are allowed. 7. Feature_ID - Depending on the annotation, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID). 8. Transcript_BioType - The bare minimum is at least a description on whether the transcript is {"Coding", "Noncoding"}. Whenever possible, use ENSEMBL biotypes. 9. Rank - Exon or Intron rank / total number of exons or introns. 10. HGVS.c - Variant using HGVS notation (DNA level). 11. HGVS.p - If variant is coding, this field describes the variant using HGVS notation (Protein level). Since transcript ID is already mentioned in 'feature ID', it may be omitted here. 12. cDNA.pos / cDNA.length - Position in cDNA and trancript's cDNA length (one based). 13. CDS.pos / CDS.length - Position and number of coding bases (one based includes START and STOP codons). 14. AA.pos / AA.length - Position and number of AA (one based, including START, but not STOP). 15. Distance - All items in this field are options, so the field could be empty. 16. ERRORS / WARNINGS - Add errors, warnings or informative message that can affect annotation accuracy. 17. INFO - Additional information.

For getting help on a specific module (e.g. ``pyvcf``):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

fuc.api.common
==============

.. automodule:: fuc.api.common
   :members:

fuc.api.pybed
=============

.. automodule:: fuc.api.pybed
   :members:

fuc.api.pyfq
============

.. automodule:: fuc.api.pyfq
   :members:

fuc.api.pyvcf
=============

.. automodule:: fuc.api.pyvcf
   :members:

fuc.api.snpeff
==============

.. automodule:: fuc.api.snpeff
   :members:

