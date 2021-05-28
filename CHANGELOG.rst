Changelog
*********

0.9.0 (in development)
----------------------

* Add new submodule ``pymaf``.
* Deprecate ``fuc.api.pyvcf.read_file`` method and add ``fuc.api.pyvcf.VcfFrame.from_file`` method.
* Deprecate ``fuc.api.pybed.read_file`` method and add ``fuc.api.pybed.BedFrame.from_file`` method.
* Deprecate ``fuc.api.pyfq.read_file`` method and add ``fuc.api.pyfq.FqFrame.from_file`` method.

0.8.0 (2021-05-27)
------------------

* Add ``pysam`` package as dependency for working with SAM/BAM/CRAM files.
* Add new submodules ``pybam`` and ``pycov``.
* Rename the commands (e.g. ``vfmerge`` to ``vcf_merge``).
* Add new commands ``bam_head/index/rename/slice``.

0.7.0 (2021-05-23)
------------------

* Add ``lxml`` package as dependency for parsing HTML files.
* Add ``matplotlib`` and ``seaborn`` packages as dependency for creating figures.
* Add ``fucdemux`` command.
* Add ``fuc.api.pyvcf.VcfFrame.filter_phased`` method.
* Add ``fuc.api.pyvcf.VcfFrame.meta_keys`` method.
* Update ``fuc.api.pyvep.filter_clinsig`` method.
* Update ``fuc.api.pyvep.filter_impact`` method.
* Add ``as_nan`` argument to ``fuc.api.pyvcf.VcfFrame.markmiss_ad/af/dp`` methods.
* Deprecate ``fuc.api.pyvcf.update`` method.
* Add ``fuc.api.pyvcf.row_updateinfo/row_parseinfo`` methods.
* The ``fuc`` package is now available on `Bioconda <https://anaconda.org/bioconda/fuc>`__.

0.6.0 (2021-05-16)
------------------

* Update Read the Docs.
* Add ``fuc.api.pyvcf.VcfFrame.markmiss_ad`` method.
* Add ``full`` argument to ``fuc.api.pyvcf.VcfFrame.markmiss_ad/af/dp`` methods.
* Add ``fucfind`` command.
* Update ``dfsum`` command.

0.5.0 (2021-05-06)
------------------

* Add ``bioconda`` package as dependency for working with BGZF compressed files.
* Update ``fuc.api.pyvcf.read_file`` and ``fuc.api.pyvcf.VcfFrame.to_file`` methods to support BGZF compressed files.
* Update Read the Docs.
* Add ``fuc.api.pyvcf.VcfFrame.slice`` method.
* Add ``vfslice`` command.

0.4.1 (2021-05-03)
------------------

* Update Read the Docs.
* Add new methods to ``fuc.api.pyvcf.VcfFrame``.
* :issue:`6`: Add ``sphinx.ext.linkcode`` extension to Read the Docs.

0.3.2 (2021-04-30)
------------------

* Rename ``snpeff`` submodule to ``pysnpeff``.
* Add new submodule ``pyvep``.
* Update ``fuc.api.pyvcf.VcfFrame`` class.
* Add ``autodocsumm`` extension to Read the Docs.
* Add contents to Read the Docs.

0.2.0 (2021-04-26)
------------------

* :issue:`2`: Fix Read the Docs automodule not working properly.
* :issue:`3`: Add ``sphinx-issues`` extension to Read the Docs.
* Rename submodules ``fuc.api.BedFrame``, ``fuc.api.FastqFrame``, and ``fuc.api.VcfFrame`` to ``fuc.api.pybed``, ``fuc.api.pyfq``, and ``fuc.api.pyvcf``, respectively.
* Add new methods to ``fuc.api.pyvcf`` submodule.
* Add new methods to ``fuc.api.pyvcf.VcfFrame`` class.
* Add new submodule ``fuc.api.snpeff``.

0.1.4 (2021-04-21)
------------------

* Initial release.
