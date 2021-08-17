Changelog
*********

0.22.0 (in development)
-----------------------

* Update :meth:`pyvcf.VcfFrame.from_file` method to be more memory efficient by pre-specifying data type for each VCF column.

0.21.0 (2021-08-16)
-------------------

* Rename :meth:`pyvcf.VcfFrame.compare` method to :meth:`pyvcf.VcfFrame.calculate_concordance`.
* Add new method :meth:`pyvcf.VcfFrame.compare`.
* Add new property ``contigs`` to :class:`pyvcf.VcfFrame`.
* Add new method :meth:`pyvcf.VcfFrame.plot_region`.
* Add special genotype keys to :meth:`pyvcf.VcfFrame.extract` method.
* :issue:`39`: Update :meth:`pyvcf.VcfFrame.extract` method to ignore rows where the genotype key of interest is not present.
* Rename :meth:`pyvcf.VcfFrame.extract` method to :meth:`pyvcf.VcfFrame.extract_format`.
* Rename :meth:`pyvcf.VcfFrame.plot_hist` method to :meth:`pyvcf.VcfFrame.plot_hist_format`.
* Add new method :meth:`pyvcf.VcfFrame.extract_info`.
* Add new method :meth:`pyvcf.VcfFrame.plot_hist_info`.
* Add new method :meth:`common.plot_exons`.
* Update :meth:`pycov.CovFrame.plot_uniformity` method to accept a list of x positions.
* Add new command :command:`ngs-fq2bam`.
* Add new command :command:`ngs-hc`.
* Add new command :command:`ngs-pon`.
* Add new command :command:`ngs-m2`.
* Add new method :meth:`common.conda_env`.
* Add new argument ``meta_only`` to :meth:`pyvcf.VcfFrame.from_file` method.
* Add new argument ``nrows`` to :meth:`pyvcf.VcfFrame.from_file` method.
* Deprecate :meth:`pybam.header` method.
* Add new method :meth:`pybam.has_chr`.

0.20.0 (2021-08-07)
-------------------

* Add new method :meth:`pymaf.MafFrame.matrix_waterfall_matched`.
* Add new method :meth:`pymaf.MafFrame.plot_waterfall_matched`.
* Add new method :meth:`pymaf.MafFrame.plot_tmb_matched`.
* Add new method :meth:`pymaf.MafFrame.plot_mutated_matched`.
* Add new method :meth:`pymaf.MafFrame.plot_oncoplot_matched`.
* Deprecate method :meth:`pymaf.MafFrame.legend_handles`.
* Add new method :meth:`common.legend_handles`.
* Deprecate classes :class:`pyvcf.AnnFrame` and :class:`pymaf.AnnFrame`. Add new class :class:`common.AnnFrame`.
* Rename :meth:`common.file2list` method to :meth:`convert_file2list`.
* Add new method :meth:`common.convert_num2cat`.
* Add new command :command:`fuc-undetm`.
* Add new method :meth:`common.plot_annot_matched`.
* Add new argument ``sheet`` to :command:`fuc-demux` command.
* Add new class :class:`common.Variant`.
* Add new method :meth:`pyvcf.rescue_filtered_variants`.
* Add new arguments ``a_size`` and ``b_size`` to :meth:`pymaf.MafFrame.plot_regplot` method.
* Rename ``hue`` and ``hue_order`` arguments in plotting methods to ``group_col`` and ``group_order``, respectively.

0.19.0 (2021-07-31)
-------------------

* Fix bug in :meth:`pymaf.MafFrame.plot_mutated` when using the ``hue`` option.
* Add new argument ``sort`` to :meth:`pymaf.MafFrame.plot_vaf` method.
* Add new method :meth:`pymaf.MafFrame.plot_matrixs`.
* Add new method :meth:`pymaf.MafFrame.plot_matrixg`.
* Add new method :meth:`pymaf.MafFrame.compute_clonality`.
* Add new method :meth:`pymaf.MafFrame.plot_clonality`.
* Fix bug in :meth:`pymaf.MafFrame.plot_evolution` when there are no variants to display for the specified samples.
* :issue:`34`: Add new method :meth:`pymaf.MafFrame.plot_genepair`.
* :issue:`34`: Add new method :meth:`pymaf.MafFrame.plot_interactions`.
* Update the :command:`fuc-demux` command to output a better figure.
* Add new method :meth:`common.plot_cytobands`.
* Add new method :meth:`pycov.CovFrame.plot_uniformity`.
* Add new method :meth:`pycov.CovFrame.plot_distribution`.
* Rename :meth:`pycov.CovFrame.from_file` method to :meth:`pycov.CovFrame.from_bam`.
* Add new method :meth:`pycov.CovFrame.from_file`.
* Add new command :command:`fuc-depth`.
* Add new method :meth:`common.file2list`.
* Add new method :meth:`pyvcf.VcfFrame.chr_prefix`.
* Fix bug in :meth:`pyvcf.gt_unphase` when '.|.' is provided.
* Update :meth:`pyvcf.VcfFrame.compare` method to only consider biallelic sites.
* Update :meth:`pyvcf.VcfFrame.compare` method to support comparison of SNVs only and INDELs only.
* Update :meth:`pymaf.MafFrame.from_vcf` method so that ``names`` argument is no longer required when ``keys`` argument is used.

0.18.0 (2021-07-20)
-------------------

* Update :command:`fq-count` command to run significantly faster.
* Update :command:`fuc-find` command to support pattern matching that is more robust than just file extension.
* Update :meth:`pyvcf.VcfFrame.subset` method to take string as input in addition to list.
* Fix bug in :meth:`pymaf.MafFrame.plot_snvclss` method when certain SNV classes are entirely missing.
* Add new package ``statsmodels`` as dependency for performing statistical analysis.
* Update :meth:`pymaf.MafFrame.plot_regplot` method to calculate and print summary statistics as well.
* Update :meth:`pyvcf.VcfFrame.plot_regplot` method to calculate and print summary statistics as well.
* :issue:`32`: Add :meth:`pyvcf.VcfFrame.miss2ref` method.
* Update :meth:`pymaf.MafFrame.plot_tmb` method to draw empty bars with warning when specified samples do not exist.
* Update :meth:`pymaf.MafFrame.plot_waterfall` method to draw empty bars with warning when specified samples do not exist.
* Add ``flip`` argument to :meth:`pymaf.MafFrame.plot_genes` method.
* Add new method :meth:`pymaf.MafFrame.plot_mutated`.

0.17.0 (2021-07-08)
-------------------

* Add new method :meth:`pymaf.MafFrame.plot_lollipop`.
* :issue:`30`: Add :meth:`pymaf.MafFrame.plot_rainfall` method.
* :issue:`30`: Add :meth:`pyvcf.VcfFrame.plot_rainfall` method.
* Update :meth:`pymaf.MafFrame.to_vcf` method to output sorted VCF.
* Add new method :meth:`pymaf.MafFrame.matrix_prevalence`.
* Add new method :meth:`pymaf.MafFrame.plot_regplot`.
* Add new argument ``samples`` to :meth:`pymaf.MafFrame.plot_snvclss` method.
* Add new method :meth:`pymaf.MafFrame.plot_evolution`.
* Add new submodule ``pygff``.

0.16.0 (2021-07-02)
-------------------

* Rename the commands (e.g. :command:`vcf_merge` to :command:`vcf-merge`).
* Add new argument ``flip`` to :meth:`pymaf.MafFrame.plot_vaf` method.
* Update :meth:`pymaf.MafFrame.plot_vaf` method to support creation of a grouped bar plot.
* Factor out ``count`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclsc`.
* Factor out ``proportion`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclsp`.
* Factor out ``samples`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_snvclss`.
* Factor out ``titv`` mode of :meth:`pymaf.MafFrame.plot_snvcls` method to new method :meth:`pymaf.MafFrame.plot_titv`.
* Deprecate method :meth:`pymaf.MafFrame.plot_snvcls`.
* Add new argument ``hue_order`` to :meth:`pyvcf.VcfFrame.plot_hist` method.
* Update aesthetic aspect of :meth:`pymaf.MafFrame.plot_oncoplot` method.
* Add new argument ``width`` to :meth:`pymaf.MafFrame.plot_tmb` method.
* Add new arguments ``palette`` and ``flip`` to :meth:`pymaf.MafFrame.plot_vartype` method.
* Update :meth:`pymaf.MafFrame.plot_snvclsc` method to support creation of a grouped bar plot.
* Update :meth:`pymaf.MafFrame.plot_snvclsp` method to support creation of a grouped box plot.
* Add new method :meth:`pyvcf.VcfFrame.plot_snvclsc` (simply wraps :meth:`pymaf.MafFrame.plot_snvclsc` method).
* Add new method :meth:`pyvcf.VcfFrame.plot_snvclsp` (simply wraps :meth:`pymaf.MafFrame.plot_snvclsp` method).
* Add new method :meth:`pyvcf.VcfFrame.plot_snvclss` (simply wraps :meth:`pymaf.MafFrame.plot_snvclss` method).
* Add new method :meth:`pyvcf.VcfFrame.plot_titv` (simply wraps :meth:`pymaf.MafFrame.plot_titv` method).
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
* Fix bug in :meth:`pymaf.MafFrame.to_frame` method.
* Deprecate methods :meth:`pyvep.filter_lof/clinsig`.
* Update :meth:`pymaf.MafFrame.from_vcf` method to extract genotype keys (e.g. DP, AD, AF).
* Update :command:`bam_slice` and :command:`bam_rename` commands.
* Deprecate method :meth:`pybam.rename`.

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
* Deprecate method :meth:`pymaf.plot_legend` and add :meth:`pymaf.legend_handles` method.
* Add new methods :meth:`pymaf.AnnFrame.legend_handles/plot_annot`.
* Add new method :meth:`pyvcf.VcfFrame.expand`.
* Rename methods :meth:`pyvcf.gt_missing/haspolyp` to :meth:`pyvcf.gt_miss/polyp`.
* Add new method :meth:`pybed.BedFrame.from_frame`.
* :issue:`14`: Add new method :meth:`pyvcf.VcfFrame.to_bed` and new command :command:`vcf_vcf2bed`.

0.9.0 (2021-06-01)
------------------

* Add new submodule ``pymaf``.
* Deprecate method :meth:`pyvcf.read_file` and add :meth:`pyvcf.VcfFrame.from_file` method.
* Deprecate method :meth:`pybed.read_file` and add :meth:`pybed.BedFrame.from_file` method.
* Deprecate method :meth:`pyfq.read_file` and add :meth:`pyfq.FqFrame.from_file` method.
* Deprecate method :meth:`pycov.read_file` and add :meth:`pycov.CovFrame.from_file` method.
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
* Deprecate method :meth:`pyvcf.update`.
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
* :issue:`6`: Add new extension ``sphinx.ext.linkcode`` to Read the Docs.

0.3.2 (2021-04-30)
------------------

* Rename ``snpeff`` submodule to ``pysnpeff``.
* Add new submodule ``pyvep``.
* Update :class:`pyvcf.VcfFrame` class.
* Add new extension ``autodocsumm`` to Read the Docs.
* Add contents to Read the Docs.

0.2.0 (2021-04-26)
------------------

* :issue:`2`: Fix Read the Docs automodule not working properly.
* :issue:`3`: Add new extension ``sphinx-issues`` to Read the Docs.
* Rename submodules ``BedFrame``, ``FastqFrame``, and ``VcfFrame`` to ``pybed``, ``pyfq``, and ``pyvcf``, respectively.
* Add new methods to ``pyvcf`` submodule.
* Add new methods to :class:`pyvcf.VcfFrame` class.
* Add new submodule ``snpeff``.

0.1.4 (2021-04-21)
------------------

* Initial release.
