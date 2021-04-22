API
***

Introduction
============

This section describes application programming interface (API) for the ``fuc`` package.

Below is the list of modules available in API:

- **BedFrame** : The BedFrame module is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **FastqFrame** : The FastqFrame module is designed for working with FASTQ files (both zipped and unzipped).
- **VcfFrame** : The VcfFrame module is designed for working with VCF files (both zipped and unzipped).
- **common** : The common module is used by other ``fuc`` modules such as `VcfFrame` and `BedFrame`. It also provides many useful methods.

For getting help on a specific module (e.g. `VcfFrame`):

.. code:: python3

   from fuc.api import VcfFrame
   help(VcfFrame)

To give:

.. parsed-literal::

   Python Library Documentation: module fuc.api.VcfFrame in fuc.api
   
   NAME
       fuc.api.VcfFrame
   
   DESCRIPTION
       The VcfFrame module is designed for working with VCF files (both zipped
       and unzipped).
   
   CLASSES
       builtins.object
           VcfFrame
       
       class VcfFrame(builtins.object)
        |  VcfFrame(meta, df)
        |  
        |  Class for storing VCF data.
        |  
        |  This class strictly sticks to the standard Variant Call Format
        |  specification (https://samtools.github.io/hts-specs/VCFv4.3.pdf).
        |  
        |  VCF lines have nine required fields for storing variant data and
        |  variable-length fields for storing sample genotype data. In all cases,
        |  missing values are specified with a dot (``.``). The required fields are:
        |  
        |  1. CHROM - An identifier from the reference genome.
        |  2. POS - The 1-based reference position.
        |  3. ID - Semicolon-separated list of unique identifiers.
        |  4. REF - Reference base(s).
        |  5. ALT - Comma-separated list of alternate non-reference alleles.
        |  6. QUAL - Phred-scaled quality score for the assertion made in ALT.
        |  7. FILTER - PASS or a semicolon-separated list of filters that fail.
        |  8. INFO - Semicolon-separated series of additional information fields.
        |  9. FORMAT - Colon-separated series of genotype fields.
        |  
        |  Methods defined here:
        |  
        |  __init__(self, meta, df)
        |      Initialize self.  See help(type(self)) for accurate signature.
        |  
        |  add_dp(self)
        |      Compute and add the DP subfield of the FORMAT field.
        |  
        |  combine(self, n1, n2)
        |      Return a new column after combining data from the two samples.
        |      
        |      This method is useful when, for example, you are trying to
        |      consolidate data from multiple replicate samples. When the same
        |      variant is found in both samples, the method will use the genotype
        |      data of the first sample.
        |      
        |      Parameters
        |      ----------
        |      n1 : str or int
        |          Name or index of the first (or original) sample.
        |      n2 : str or int
        |          Name or index of the second (or replicate) sample.
        |      
        |      Returns
        |      -------
        |      s : pandas.Series
        |          VCF column representing the combined data.
        |  
        |  compare(self, n1, n2)
        |      Compare two samples within the VcfFrame.
        |      
        |      Parameters
        |      ----------
        |      n1 : str or int
        |          Name or index of the test sample.
        |      n2 : str or int
        |          Name or index of the truth sample.
        |      
        |      Returns
        |      -------
        |      result : tuple
        |          Comparison result (tp, fp, fn, tn).
        |  
        |  filter_af(self, threshold=0.1)
        |      Filter rows based on the AF subfield of the FORMAT field.
        |  
        |  filter_bed(self, bed)
        |      Filter rows based on BED data.
        |      
        |      Parameters
        |      ----------
        |      bed : BedFrame or str
        |          BedFrame or path to a BED file.
        |      
        |      Returns
        |      -------
        |      vf : VcfFrame
        |          Filtered VcfFrame.
        |  
        |  filter_dp(self, threshold=200)
        |      Filter rows based on the DP subfield of the FORMAT field.
        |  
        |  filter_empty(self)
        |      Filter out rows that have no genotype calls.
        |  
        |  filter_multiallelic(self)
        |      Filter out rows that have multiple alternative alleles.
        |  
        |  merge(self, other, how='inner', format='GT')
        |      Merge with the other VcfFrame.
        |      
        |      This method essentially wraps the `pandas.DataFrame.merge` method.
        |      
        |      Parameters
        |      ----------
        |      other : VcfFrame
        |          Other VcfFrame.
        |      how : str, default: 'inner'
        |          Type of merge to be performed. ['left', 'right', 'outer',
        |          'inner', 'cross']
        |      format : str, default: 'GT'
        |          FORMAT subfields to be retained (e.g. 'GT:AD:DP').
        |      
        |      Returns
        |      -------
        |      vf : VcfFrame
        |          Merged VcfFrame.
        |  
        |  parse_snpeff(self, idx, sep=' | ')
        |      Parse SnpEff annotations.
        |      
        |      SnpEff provides the following functional annotations:
        |      
        |      1. Allele
        |      2. Annotation
        |      3. Annotation_Impact
        |      4. Gene_Name
        |      5. Gene_ID
        |      6. Feature_Type
        |      7. Feature_ID
        |      8. Transcript_BioType
        |      9. Rank
        |      10. HGVS.c
        |      11. HGVS.p
        |      12. cDNA.pos / cDNA.length
        |      13. CDS.pos / CDS.length
        |      14. AA.pos / AA.length
        |      15. Distance
        |      16. ERRORS / WARNINGS
        |      17. INFO
        |      
        |      Parameters
        |      ----------
        |      i : list
        |          List of annotation indicies.
        |      sep : str, default: ' | '
        |          Separator for joining requested annotations.
        |      
        |      Returns
        |      -------
        |      s : pandas.Series
        |          Parsed annotations.
        |  
        |  reset_samples(self, names)
        |      Reset the sample list.
        |  
        |  strip(self, format='GT')
        |      Remove unnecessary data from the VcfFrame.
        |      
        |      Parameters
        |      ----------
        |      format : str, default: 'GT'
        |          FORMAT subfields to be retained (e.g. 'GT:AD:DP').
        |      
        |      Returns
        |      -------
        |      vf : VcfFrame
        |          Stripped VcfFrame.
        |  
        |  subtract(self, name)
        |      Remove rows that have a variant call in the sample.
        |      
        |      Parameters
        |      ----------
        |      name : str or int
        |          Name or index of the sample.
        |      
        |      Returns
        |      -------
        |      vf : VcfFrame
        |          Filtered VcfFrame.
        |  
        |  to_file(self, file_path)
        |      Write the VcfFrame to a VCF file.
        |  
        |  to_string(self)
        |      Render the VcfFrame to a console-friendly tabular output.
        |  
        |  update(self, other, headers=None, missing=True)
        |      Copy data from the other VcfFrame.
        |      
        |      This method will copy and paste data from the other VcfFrame for
        |      overlapping records. By default, the following VCF headers are
        |      used: ID, QUAL, FILTER, and, INFO.
        |      
        |      Parameters
        |      ----------
        |      other : VcfFrame
        |          Other VcfFrame.
        |      headers : list, optional
        |          List of VCF headers to exclude.
        |      missing : bool, default: True
        |          If True, only fields with the missing value ('.') will be updated.
        |      
        |      Returns
        |      -------
        |      vf : VcfFrame
        |          Updated VcfFrame.
        |  
        |  ----------------------------------------------------------------------
        |  Class methods defined here:
        |  
        |  from_file(file_path) from builtins.type
        |      Create a VcfFrame from a VCF file.
        |  
        |  ----------------------------------------------------------------------
        |  Readonly properties defined here:
        |  
        |  samples
        |      Return a list of the sample IDs.
        |  
        |  ----------------------------------------------------------------------
        |  Data descriptors defined here:
        |  
        |  __dict__
        |      dictionary for instance variables (if defined)
        |  
        |  __weakref__
        |      list of weak references to the object (if defined)
   
   FILE
       /Users/sbslee/Desktop/fuc/fuc/api/VcfFrame.py
   

fuc.api.BedFrame
================

.. automodule:: fuc.api.BedFrame
   :members:

fuc.api.FastqFrame
==================

.. automodule:: fuc.api.FastqFrame
   :members:

fuc.api.VcfFrame
================

.. automodule:: fuc.api.VcfFrame
   :members:

fuc.api.common
==============

.. automodule:: fuc.api.common
   :members:

