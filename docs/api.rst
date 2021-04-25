API
***

Introduction
============

This section describes application programming interface (API) for the ``fuc`` package.

Below is the list of submodules available in API:

- **common** : The ``common`` submodule is used by other ``fuc`` submodules such as ``pyvcf`` and ``pybed``. It also provides many useful methods.
- **pybed** : The ``pybed`` submodule is designed for working with BED files. For example, it can be used to find the intersection between multiple BED files.
- **pyfq** : The ``pyfq`` submodule is designed for working with FASTQ files (both zipped and unzipped).
- **pyvcf** : The ``pyvcf`` submodule is designed for working with VCF files (both zipped and unzipped). It implements ``pyvcf.VcfFrame`` which stores VCF data as a ``pandas.DataFrame`` to allow fast computation and easy minipulation.
- **snpeff** : The ``snpeff`` submodule is designed for parsing VCF annotation data from the SnpEff program.

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

