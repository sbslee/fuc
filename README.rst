README
******

.. image:: https://badge.fury.io/py/fuc.svg
    :target: https://badge.fury.io/py/fuc

.. image:: https://readthedocs.org/projects/sbslee-fuc/badge/?version=latest
   :target: https://sbslee-fuc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Introduction
============

The main goal of the fuc package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use fuc for both command line interface (CLI) and application programming interface (API) whose documentations are available at `Read the Docs <https://sbslee-fuc.readthedocs.io/en/latest/>`_.

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

The following packages are required to run fuc:

.. parsed-literal::

   numpy
   pandas
   pyranges

Getting Started
===============

There are various ways you can install fuc. The easiest one would be to use pip:

.. code-block:: console

   $ pip install fuc

Above will automatically download and install all the dependencies as well.

Alternatively, you can clone the GitHub repository and then install fuc this way:

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
       fuccompf     [FUC] compare contents of two files
       fucexist     [FUC] check whether files/dirs exist
       qfcount      [FASTQ] count sequence reads in FASTQ files
       qfsum        [FASTQ] summarize a FASTQ file
       vfmerge      [VCF] merge two or more VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     -v, --version  show the version number and exit

For getting help on a specific command (e.g. vfmerge):

.. code-block:: console

   $ fuc vfmerge -h

Below is the list of submodules available in API:

- **common** : The common submodule is used by other fuc submodules such as pyvcf and pybed. It also provides many day-to-day actions used in the field of bioinformatics.
- **pybed** : The pybed submodule is designed for working with BED files. It implements ``pybed.BedFrame`` which stores BED data as ``pyranges.PyRanges`` to allow fast computation and easy manipulation.
- **pyfq** : The pyfq submodule is designed for working with FASTQ files (both zipped and unzipped). It implements ``pyfq.FqFrame`` which stores FASTQ data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
- **pysnpeff** : The pysnpeff submodule is designed for parsing VCF annotation data from the SnpEff program. It should be used with ``pyvcf.VcfFrame``.
- **pyvcf** : The pyvcf submodule is designed for working with VCF files (both zipped and unzipped). It implements ``pyvcf.VcfFrame`` which stores VCF data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
- **pyvep** : The pyvep submodule is designed for parsing VCF annotation data from the Ensembl Variant Effect Predictor (VEP). It should be used with ``pyvcf.VcfFrame``.

For getting help on a specific module (e.g. pyvcf):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

