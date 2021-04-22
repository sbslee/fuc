README
******

.. image:: https://badge.fury.io/py/fuc.svg
    :target: https://badge.fury.io/py/fuc

.. image:: https://readthedocs.org/projects/sbslee-fuc/badge/?version=latest
   :target: https://sbslee-fuc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Introduction
============

The main goal of the ``fuc`` package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use ``fuc`` for both command line interface (CLI) and application programming interface (API) whose documentations are available at `Read the Docs <https://sbslee-fuc.readthedocs.io/en/latest/>`_.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

| Author: Seung-been "Steven" Lee
| Email: sbstevenlee@gmail.com
| License: MIT License

Examples
========

To merge VCF files with CLI:

.. code-block:: console

   $ fuc vfmerge 1.vcf 2.vcf 3.vcf > merged.vcf

To filter a VCF file based on a BED file using API:

.. code:: python3

   from fuc import pyvcf
   vf = pyvcf.read_file('original.vcf')
   filtered_vf = vf.filter_bed('targets.bed')
   filtered_vf.to_file('filtered.vcf')

Required Packages
=================

The following packages are required to run ``fuc``:

.. parsed-literal::

   numpy
   pandas
   pyranges

Getting Started
===============

There are various ways you can install ``fuc``. The easiest one would be to use ``pip``:

.. code-block:: console

   $ pip install fuc

Above will automatically download and install all the dependencies as well.

Alternatively, you can clone the GitHub repository and then install ``fuc`` this way:

.. code-block:: console

   $ git clone https://github.com/sbslee/fuc
   $ cd fuc
   $ pip install .

Above will also allow you to install a development version that's not available in PyPI.

For getting help on CLI:

.. code-block:: console

   $ fuc -h
   usage: fuc [-h] [-v] COMMAND ...
   
   positional arguments:
     COMMAND        name of the command
       bfintxn      [BED] find intersection of two or more BED files
       bfsum        [BED] summarize a BED file
       dfmerge      [TABLE] merge two text files
       dfsum        [TABLE] summarize a text file
       fuccompf     [FUC] compare two files
       fucexist     [FUC] check whether files/dirs exist
       qfcount      [FASTQ] count sequence reads in FASTQ files
       qfsum        [FASTQ] summarize a FASTQ file
       vfmerge      [VCF] merge two or more VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     -v, --version  show the version number and exit

For getting help on a specific command (e.g. `vfmerge`):

.. code-block:: console

   $ fuc vfmerge -h
   usage: fuc vfmerge [-h] [--how TEXT] [--format TEXT] vcf_files [vcf_files ...]
   
   This command will merge multiple VCF files (both zipped and unzipped). By
   default, only the GT subfield of the FORMAT field will be included in the
   merged VCF. Use '--format' to include additional FORMAT subfields such as AD
   and DP.
   
   positional arguments:
     vcf_files      VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     --how TEXT     type of merge as defined in `pandas.DataFrame.merge`
                    (default: 'inner')
     --format TEXT  FORMAT subfields to be retained (e.g. 'GT:AD:DP') (default:
                    'GT')

Below is the list of modules available in API:

- **common** : The common module is used by other ``fuc`` modules such as `VcfFrame` and `BedFrame`. It also provides many useful methods.
- **pybed** : The ``pybed`` module is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **pyfq** : The ``pyfq`` module is designed for working with FASTQ files (both zipped and unzipped).
- **pyvcf** : The ``pyvcf`` module is designed for working with VCF files (both zipped and unzipped).

For getting help on a specific module (e.g. `pyvcf`):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

To give:

.. parsed-literal::

   Python Library Documentation: module fuc.api.pyvcf in fuc.api
   
   NAME
       fuc.api.pyvcf
   
   DESCRIPTION
       The ``pyvcf`` module is designed for working with VCF files (both zipped
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
        |      bed : pybed.BedFrame or str
        |          pybed.BedFrame or path to a BED file.
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
        |      Parameters
        |      ----------
        |      other : VcfFrame
        |          Other VcfFrame.
        |      how : str, default: 'inner'
        |          Type of merge as defined in `pandas.DataFrame.merge`.
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
   
   FUNCTIONS
       merge(vfs, how='inner', format='GT')
           Return the merged VcfFrame.
           
           Parameters
           ----------
           vfs : list
               List of VcfFrames to be merged.
           how : str, default: 'inner'
               Type of merge as defined in `pandas.DataFrame.merge`.
           format : str, default: 'GT'
               FORMAT subfields to be retained (e.g. 'GT:AD:DP').
           
           Returns
           -------
           merged_vf : VcfFrame
               Merged VcfFrame.
       
       read_file(fn)
           Create a VcfFrame from a VCF file (both zipped and unzipped).
   
   FILE
       /Users/sbslee/Desktop/fuc/fuc/api/pyvcf.py
   
