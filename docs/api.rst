API
***

Introduction
============

This section describes application programming interface (API) for the ``fuc`` package.

Below is the list of submodules available in API:

- **common** : The ``common`` submodule is used by other ``fuc`` submodules such as
- **pybed** : The ``pybed`` submodule is designed for working with BED files. For example,
- **pyfq** : The ``pyfq`` submodule is designed for working with FASTQ files (both zipped
- **pyvcf** : The ``pyvcf`` submodule is designed for working with VCF files (both zipped
- **snpeff** : One VCF record can have several SnpEff annotations if, for example,

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

