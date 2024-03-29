..
   This file was automatically generated by docs/create.py.

API
***

Introduction
============

This section describes application programming interface (API) for the fuc package.

Below is the list of submodules available in the fuc API:

- **common** : The common submodule is used by other fuc submodules such as pyvcf and pybed. It also provides many day-to-day actions used in the field of bioinformatics.
- **pybam** : The pybam submodule is designed for working with sequence alignment files (SAM/BAM/CRAM). It essentially wraps the `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast computation and easy manipulation. If you are mainly interested in working with depth of coverage data, please check out the pycov submodule which is specifically designed for the task.
- **pybed** : The pybed submodule is designed for working with BED files. It implements ``pybed.BedFrame`` which stores BED data as ``pandas.DataFrame`` via the `pyranges <https://github.com/biocore-ntnu/pyranges>`_ package to allow fast computation and easy manipulation. The submodule strictly adheres to the standard `BED specification <https://genome.ucsc.edu/FAQ/FAQformat.html>`_.
- **pychip** : The pychip submodule is designed for working with annotation or manifest files from the Axiom (Thermo Fisher Scientific) and Infinium (Illumina) array platforms.
- **pycov** : The pycov submodule is designed for working with depth of coverage data from sequence alingment files (SAM/BAM/CRAM). It implements ``pycov.CovFrame`` which stores read depth data as ``pandas.DataFrame`` via the `pysam <https://pysam.readthedocs.io/en/latest/api.html>`_ package to allow fast computation and easy manipulation. The ``pycov.CovFrame`` class also contains many useful plotting methods such as ``CovFrame.plot_region`` and ``CovFrame.plot_uniformity``.
- **pyfq** : The pyfq submodule is designed for working with FASTQ files. It implements ``pyfq.FqFrame`` which stores FASTQ data as ``pandas.DataFrame`` to allow fast computation and easy manipulation.
- **pygff** : The pygff submodule is designed for working with GFF/GTF files. It implements ``pygff.GffFrame`` which stores GFF/GTF data as ``pandas.DataFrame`` to allow fast computation and easy manipulation. The submodule strictly adheres to the standard `GFF specification <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>`_.
- **pykallisto** : The pykallisto submodule is designed for working with RNAseq quantification data from Kallisto. It implements ``pykallisto.KallistoFrame`` which stores Kallisto's output data as ``pandas.DataFrame`` to allow fast computation and easy manipulation. The ``pykallisto.KallistoFrame`` class also contains many useful plotting methods such as ``KallistoFrame.plot_differential_abundance``.
- **pymaf** : The pymaf submodule is designed for working with MAF files. It implements ``pymaf.MafFrame`` which stores MAF data as ``pandas.DataFrame`` to allow fast computation and easy manipulation. The ``pymaf.MafFrame`` class also contains many useful plotting methods such as ``MafFrame.plot_oncoplot`` and ``MafFrame.plot_summary``. The submodule strictly adheres to the standard `MAF specification <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_.
- **pysnpeff** : The pysnpeff submodule is designed for parsing VCF annotation data from the `SnpEff <https://pcingola.github.io/SnpEff/>`_ program. It should be used with ``pyvcf.VcfFrame``.
- **pyvcf** : The pyvcf submodule is designed for working with VCF files. It implements ``pyvcf.VcfFrame`` which stores VCF data as ``pandas.DataFrame`` to allow fast computation and easy manipulation. The ``pyvcf.VcfFrame`` class also contains many useful plotting methods such as ``VcfFrame.plot_comparison`` and ``VcfFrame.plot_tmb``. The submodule strictly adheres to the standard `VCF specification <https://samtools.github.io/hts-specs/VCFv4.3.pdf>`_.
- **pyvep** : The pyvep submodule is designed for parsing VCF annotation data from the `Ensembl VEP <https://asia.ensembl.org/info/docs/tools/vep/index.html>`_ program. It should be used with ``pyvcf.VcfFrame``.

For getting help on a specific submodule (e.g. pyvcf):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

fuc.common
==========

.. automodule:: fuc.api.common
   :members:

fuc.pybam
=========

.. automodule:: fuc.api.pybam
   :members:

fuc.pybed
=========

.. automodule:: fuc.api.pybed
   :members:

fuc.pychip
==========

.. automodule:: fuc.api.pychip
   :members:

fuc.pycov
=========

.. automodule:: fuc.api.pycov
   :members:

fuc.pyfq
========

.. automodule:: fuc.api.pyfq
   :members:

fuc.pygff
=========

.. automodule:: fuc.api.pygff
   :members:

fuc.pykallisto
==============

.. automodule:: fuc.api.pykallisto
   :members:

fuc.pymaf
=========

.. automodule:: fuc.api.pymaf
   :members:

fuc.pysnpeff
============

.. automodule:: fuc.api.pysnpeff
   :members:

fuc.pyvcf
=========

.. automodule:: fuc.api.pyvcf
   :members:

fuc.pyvep
=========

.. automodule:: fuc.api.pyvep
   :members:

