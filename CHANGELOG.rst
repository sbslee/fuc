Changelog
*********

0.9.0 (2021-06-01)
------------------

* Add new submodule ``pymaf``.
* Deprecate :meth:`pyvcf.read_file` method and add :meth:`pyvcf.VcfFrame.from_file` method.
* Deprecate :meth:`pybed.read_file` method and add :meth:`pybed.BedFrame.from_file` method.
* Deprecate :meth:`pyfq.read_file` method and add :meth:`pyfq.FqFrame.from_file` method.
* Deprecate :meth:`pycov.read_file` method and add :meth:`pycov.CovFrame.from_file` method.
* Add new method :meth:`common.parse_region`.

0.8.0 (2021-05-27)
------------------

* Add ``pysam`` package as dependency for working with SAM/BAM/CRAM files.
* Add new submodules ``pybam`` and ``pycov``.
* Rename the commands (e.g. ``vfmerge`` to ``vcf_merge``).
* :issue:`11`: Add new command ``bam_slice``.
* Add new commands ``bam_head/index/rename``.

0.7.0 (2021-05-23)
------------------

* Add ``lxml`` package as dependency for parsing HTML files.
* Add ``matplotlib`` and ``seaborn`` packages as dependency for creating figures.
* Add ``fucdemux`` command.
* Add :meth:`pyvcf.VcfFrame.filter_phased` method.
* Add :meth:`pyvcf.VcfFrame.meta_keys` method.
* Update :meth:`pyvep.filter_clinsig` method.
* Update :meth:`pyvep.filter_impact` method.
* Add ``as_nan`` argument to :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` methods.
* Deprecate :meth:`pyvcf.update` method.
* Add :meth:`pyvcf.row_updateinfo/parseinfo` methods.
* The ``fuc`` package is now available on `Bioconda <https://anaconda.org/bioconda/fuc>`__.

0.6.0 (2021-05-16)
------------------

* Update Read the Docs.
* Add :meth:`pyvcf.VcfFrame.markmiss_ad` method.
* Add ``full`` argument to :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` methods.
* Add ``fucfind`` command.
* Update ``dfsum`` command.

0.5.0 (2021-05-06)
------------------

* Add ``biopython`` package as dependency for working with BGZF compressed files.
* Update :meth:`pyvcf.read_file` method and :meth:`pyvcf.VcfFrame.to_file` method to support BGZF compressed files.
* Update Read the Docs.
* Add :meth:`pyvcf.VcfFrame.slice` method.
* Add ``vfslice`` command.

0.4.1 (2021-05-03)
------------------

* Update Read the Docs.
* Add new methods to :class:`pyvcf.VcfFrame` class.
* :issue:`6`: Add ``sphinx.ext.linkcode`` extension to Read the Docs.

0.3.2 (2021-04-30)
------------------

* Rename ``snpeff`` submodule to ``pysnpeff``.
* Add new submodule ``pyvep``.
* Update :class:`pyvcf.VcfFrame` class.
* Add ``autodocsumm`` extension to Read the Docs.
* Add contents to Read the Docs.

0.2.0 (2021-04-26)
------------------

* :issue:`2`: Fix Read the Docs automodule not working properly.
* :issue:`3`: Add ``sphinx-issues`` extension to Read the Docs.
* Rename submodules ``BedFrame``, ``FastqFrame``, and ``VcfFrame`` to ``pybed``, ``pyfq``, and ``pyvcf``, respectively.
* Add new methods to ``pyvcf`` submodule.
* Add new methods to :class:`pyvcf.VcfFrame` class.
* Add new submodule ``snpeff``.

0.1.4 (2021-04-21)
------------------

* Initial release.
