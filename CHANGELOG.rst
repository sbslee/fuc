Changelog
*********

0.14.0 (in development)
-----------------------

* Deprecate methods :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` and add new method :meth:`pyvcf.VcfFrame.markmiss`.

0.13.0 (2021-06-16)
-------------------

* Add new method :meth:`pymaf.MafFrame.to_vcf`.
* :issue:`21`: Add new command :command:`maf_maf2vcf`.
* Add new method :meth:`pyvcf.VcfFrame.rename`.
* Add new command :command:`vcf_rename`.
* Add new method :meth:`pymaf.MafFrame.plot_vaf`.
* Update :meth:`pyvcf.VcfFrame.slice` method.
* Update :command:`vcf_slice` command.

0.12.0 (2021-06-12)
-------------------

* Add new method :meth:`pyvcf.VcfFrame.add_af`.
* Add new method :meth:`pyvcf.VcfFrame.extract`.
* Deprecate methods :meth:`pyvep.filter_af/biotype/nothas/impact`.
* Add new method :meth:`pyvep.filter_query`.
* :issue:`19`: Add new command :command:`vcf_vep`.
* Rename :meth:`pyvcf.VcfFrame.plot_histplot` to :meth:`pyvcf.VcfFrame.plot_tmb`.
* Add ``scipy`` package as dependency for performing statistical analysis.
* Add new method :meth:`pyvcf.VcfFrame.plot_hist`.

0.11.0 (2021-06-10)
-------------------

* :issue:`16`: Add new method :meth:`pyvcf.VcfFrame.cfilter_empty`.
* Add new methods :meth:`pyvep.filter_af/lof`.
* Add ``matplotlib-venn`` package as dependency for plotting Venn diagrams.
* Add new methods :meth:`pyvcf.plot_comparison/regplot/histplot`.
* :issue:`17`: Add new method :meth:`pyvep.filter_biotype`.
* Add new class :class:`pyvcf.AnnFrame`.

0.10.0 (2021-06-03)
-------------------

* Add new methods :meth:`pymaf.plot_summary/varsum`.
* Add new command :command:`maf_sumplt`.
* Add new method :meth:`pymaf.MafFrame.to_string`.
* Update :command:`maf_oncoplt` command.
* Add new method :meth:`pyvcf.VcfFrame.filter_qual`.
* Deprecate :meth:`pymaf.plot_legend` method and add :meth:`pymaf.legend_handles` method.
* Add new methods :meth:`pymaf.AnnFrame.legend_handles/plot_annot`.
* Add new method :meth:`pyvcf.VcfFrame.expand`.
* Rename methods :meth:`pyvcf.gt_missing/haspolyp` to :meth:`pyvcf.gt_miss/polyp`.
* Add new method :meth:`pybed.BedFrame.from_frame`.
* :issue:`14`: Add new method :meth:`pyvcf.VcfFrame.to_bed` and new command :command:`vcf_vcf2bed`.

0.9.0 (2021-06-01)
------------------

* Add new submodule ``pymaf``.
* Deprecate :meth:`pyvcf.read_file` method and add :meth:`pyvcf.VcfFrame.from_file` method.
* Deprecate :meth:`pybed.read_file` method and add :meth:`pybed.BedFrame.from_file` method.
* Deprecate :meth:`pyfq.read_file` method and add :meth:`pyfq.FqFrame.from_file` method.
* Deprecate :meth:`pycov.read_file` method and add :meth:`pycov.CovFrame.from_file` method.
* Add new method :meth:`common.parse_region`.
* Add new commands :command:`maf_oncoplt/vcf2maf`.

0.8.0 (2021-05-27)
------------------

* Add ``pysam`` package as dependency for working with SAM/BAM/CRAM files.
* Add new submodules ``pybam`` and ``pycov``.
* Rename the commands (e.g. :command:`vfmerge` to :command:`vcf_merge`).
* :issue:`11`: Add new command :command:`bam_slice`.
* Add new commands :command:`bam_head/index/rename`.

0.7.0 (2021-05-23)
------------------

* Add ``lxml`` package as dependency for parsing HTML files.
* Add ``matplotlib`` and ``seaborn`` packages as dependency for creating figures.
* Add new command :command:`fucdemux`.
* Add new method :meth:`pyvcf.VcfFrame.filter_phased`.
* Add new method :meth:`pyvcf.VcfFrame.meta_keys`.
* Update :meth:`pyvep.filter_clinsig` method.
* Update :meth:`pyvep.filter_impact` method.
* Add ``as_nan`` argument to :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` methods.
* Deprecate :meth:`pyvcf.update` method.
* Add new methods :meth:`pyvcf.row_updateinfo/parseinfo`.
* The ``fuc`` package is now available on `Bioconda <https://anaconda.org/bioconda/fuc>`__.

0.6.0 (2021-05-16)
------------------

* Update Read the Docs.
* Add new method :meth:`pyvcf.VcfFrame.markmiss_ad`.
* Add ``full`` argument to :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` methods.
* Add new command :command:`fucfind`.
* Update :command:`dfsum` command.

0.5.0 (2021-05-06)
------------------

* Add ``biopython`` package as dependency for working with BGZF compressed files.
* Update :meth:`pyvcf.read_file` method and :meth:`pyvcf.VcfFrame.to_file` method to support BGZF compressed files.
* Update Read the Docs.
* Add new method :meth:`pyvcf.VcfFrame.slice`.
* Add new command :command:`vfslice`.

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
