Changelog
*********

0.18.0 (in development)
-----------------------

* Update :command:`fq-count` command to run significantly faster.
* Update :command:`fuc-find` command to support pattern matching that is more robust than just file extension.
* Update :meth:`pyvcf.VcfFrame.subset` method to take string as input in addition to list.
* Fix a minor bug in :meth:`pymaf.MafFrame.plot_snvclss` method when certain SNV classes are entirely missing.
* Add ``statsmodels`` package as dependency for performing statistical analysis.
* Update :meth:`pymaf.MafFrame.plot_regplot` method to calculate and print summary statistics as well.
* Update :meth:`pyvcf.VcfFrame.plot_regplot` method to calculate and print summary statistics as well.
* Add :meth:`pyvcf.VcfFrame.miss2ref` method.
* Update :meth:`pymaf.MafFrame.plot_tmb` method to draw empty bars with warning when specified samples do not exist.
* Update :meth:`pymaf.MafFrame.plot_waterfall` method to draw empty bars with warning when specified samples do not exist.
* Add ``flip`` argument to :meth:`pymaf.MafFrame.plot_genes` method.

0.17.0 (2021-07-08)
-------------------

* Add :meth:`pymaf.MafFrame.plot_lollipop` method.
* :issue:`30`: Add :meth:`pymaf.MafFrame.plot_rainfall` method.
* :issue:`30`: Add :meth:`pyvcf.VcfFrame.plot_rainfall` method.
* Update :meth:`pymaf.MafFrame.to_vcf` method to output sorted VCF.
* Add :meth:`pymaf.MafFrame.matrix_prevalence` method.
* Add :meth:`pymaf.MafFrame.plot_regplot` method.
* Add ``samples`` argument to :meth:`pymaf.MafFrame.plot_snvclss` method.
* Add :meth:`pymaf.MafFrame.plot_evolution` method.
* Add new submodule ``pygff``.

0.16.0 (2021-07-02)
-------------------

* Rename the commands (e.g. :command:`vcf_merge` to :command:`vcf-merge`).
* Add ``flip`` argument to :meth:`pymaf.MafFrame.plot_vaf` method.
* Update :meth:`pymaf.MafFrame.plot_vaf` method to support creation of a grouped bar plot.
* Factor out ``count`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclsc`.
* Factor out ``proportion`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclsp`.
* Factor out ``samples`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclss`.
* Factor out ``titv`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_titv`.
* Deprecate :meth:`pymaf.MafFrame.plot_snvcls` method.
* Add ``hue_order`` argument to :meth:`pyvcf.VcfFrame.plot_hist` method.
* Update aesthetic aspect of :meth:`pymaf.MafFrame.plot_oncoplot` method.
* Add ``width`` argument to :meth:`pymaf.MafFrame.plot_tmb` method.
* Add ``palette`` and ``flip`` arguments to :meth:`pymaf.MafFrame.plot_vartype` method.
* Update :meth:`pymaf.MafFrame.plot_snvclsc` method to support creation of a grouped bar plot.
* Update :meth:`pymaf.MafFrame.plot_snvclsp` method to support creation of a grouped box plot.
* Add :meth:`pyvcf.VcfFrame.plot_snvclsc` method (simply wraps :meth:`pymaf.MafFrame.plot_snvclsc` method).
* Add :meth:`pyvcf.VcfFrame.plot_snvclsp` method (simply wraps :meth:`pymaf.MafFrame.plot_snvclsp` method).
* Add :meth:`pyvcf.VcfFrame.plot_snvclss` method (simply wraps :meth:`pymaf.MafFrame.plot_snvclss` method).
* Add :meth:`pyvcf.VcfFrame.plot_titv` method (simply wraps :meth:`pymaf.MafFrame.plot_titv` method).
* :issue:`28`: Update :meth:`pymaf.MafFrame.from_vcf` method to handle unannotated VCF data.

0.15.0 (2021-06-24)
-------------------

* Update :command:`vcf_filter` command.
* Update :command:`tbl_sum` command.
* Add ``samples`` and ``shape`` attributes to :class:`pymaf.AnnFrame` class.
* Rename :meth:`pymaf.MafFrame.compute_genes/tmb/waterfall` methods to :meth:`pymaf.MafFrame.matrix_genes/tmb/waterfall`.
* Add ``keep_empty`` argument to :meth:`pymaf.MafFrame.matrix_waterfall/plot_oncoplot/plot_waterfall` methods.
* Add :meth:`pymaf.MafFrame.filter_annot` method.
* Add :meth:`pymaf.AnnFrame.sorted_samples` method.
* Fix minor bug in :meth:`pymaf.MafFrame.to_frame` method.
* Deprecate :meth:`pyvep.filter_lof/clinsig` methods.
* Update :meth:`pymaf.MafFrame.from_vcf` method to extract genotype keys (e.g. DP, AD, AF).
* Update :command:`bam_slice` and :command:`bam_rename` commands.
* Deprecate :meth:`pybam.rename` method.

0.14.0 (2021-06-20)
-------------------

* :issue:`23`: Deprecate methods :meth:`pyvcf.VcfFrame.markmiss_ad/af/dp` and add new method :meth:`pyvcf.VcfFrame.markmiss`.
* Add new command :command:`vcf_filter`.
* Update methods :meth:`pycov.CovFrame.slice/plot_region`.
* :issue:`24`: Add new method :meth:`pyvcf.VcfFrame.drop_duplicates`.
* Update :meth:`pymaf.MafFrame.plot_snvcls` method to support various plotting modes.
* Rename ``horizontal`` argument of :meth:`pymaf.MafFrame.plot_varsum` method to ``flip``.

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
