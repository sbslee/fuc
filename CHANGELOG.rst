Changelog
*********

0.38.0 (in development)
-----------------------

* Update :meth:`pyvcf.has_chr_prefix` method to ignore the HLA contigs for GRCh38.

0.37.0 (2023-09-09)
-------------------

* :issue:`67`: Fix bug in :meth:`pymaf.MafFrame.plot_waterfall` method where ``count=1`` was causing color mismatch.
* Add new submodule ``pychip``.
* Add new method :meth:`common.reverse_complement`.
* Fix bug in :meth:`common.extract_sequence` method where a long DNA sequence output was truncated.
* :issue:`68`: Refresh the variant consequences database from Ensembl VEP. The database's latest update was on May 31, 2021.

0.36.0 (2022-08-12)
-------------------

* ``fuc`` now has a citation! Please refer to the publication “`ClinPharmSeq: A targeted sequencing panel for clinical pharmacogenetics implementation <https://doi.org/10.1371/journal.pone.0272129>`__” by Lee et al., 2022 (Steven is the first author). Fore more details, see the Citation section in README.
* Update ``pyvcf`` submodule to accept "sites-only" VCF.
* Add new method :meth:`pyvcf.VcfFrame.filter_gsa`.
* Add new method :meth:`pyvcf.VcfFrame.duplicated`.
* Add new optional argument ``to_csv`` to :meth:`pymaf.MafFrame.plot_regplot_tmb` method.
* Add new optional argument ``count`` to :meth:`pymaf.MafFrame.plot_mutated_matched` method.

0.35.0 (2022-07-12)
-------------------

* Fix bug in :meth:`pyvcf.VcfFrame.pseudophase` method.
* Add new methods :meth:`pyvcf.VcfFrame.diploidize` and :meth:`pyvcf.gt_diploidize`.
* Update :meth:`pyvcf.VcfFrame.get_af` method to handle situations where there are multiple records with the same ``REF`` allele.
* Add new method :meth:`pymaf.MafFrame.plot_regplot_tmb`.
* Rename :meth:`pyvcf.VcfFrame.plot_regplot` method to :meth:`pyvcf.VcfFrame.plot_regplot_tmb` and :meth:`pymaf.MafFrame.plot_regplot` method to :meth:`pymaf.MafFrame.plot_regplot_gene`.

0.34.0 (2022-06-08)
-------------------

* Add new optional argument ``--stranded`` to :command:`ngs-quant` command.
* Add new method :meth:`pycov.CovFrame.merge`.
* Add new method :meth:`pycov.merge`.
* :issue:`61`: Update :meth:`pymaf.MafFrame.from_vcf` method to automatically detect CSQ field in INFO column (thanks `@lbeltrame <https://github.com/lbeltrame>`__).
* :issue:`63`: Update :meth:`pyvcf.VcfFrame.sort` method to handle contigs that are not pre-defined.

0.33.1 (2022-05-03)
-------------------

* Add new method :meth:`pybam.index` which simply wraps :meth:`pysam.index` method.
* Update :command:`bam-index` command to use :meth:`pybam.index` method.
* Add new method :meth:`pybam.slice`.
* Update :command:`bam-slice` command to use :meth:`pybam.slice` method.
* Update :command:`ngs-bam2fq` and :command:`ngs-fq2bam` commands to allow users to run in local environment.
* Update :command:`ngs-fq2bam` command to handle cases where input FASTQ does not have information on flowcell and barcode.
* Update :meth:`pyvcf.call` method to run more efficiently.

0.32.0 (2022-04-02)
-------------------

* Add new optional argument ``filter_off`` to :class:`pykallisto.KallistoFrame` constructor, which is useful for generating a simple count or tpm matrix.
* Add new optional argument ``--dir-path`` to :command:`vcf-call` command for storing intermediate files.
* Add new optional argument ``--gap_frac`` to :command:`vcf-call` command so that users can control indel calling sensitivity.
* Add new optional argument ``--group-samples`` to :command:`vcf-call` command so that users can group samples into populations and apply the HWE assumption within but not across the populations.
* Fix minor bug in :meth:`pyvcf.call` method when ``pybed.BedFrame`` object is given as ``regions``.

0.31.0 (2022-03-01)
-------------------

* Fix bug in :meth:`pykallisto.KallistoFrame.compute_fold_change` method.
* Add new method :meth:`pyvcf.call` and new command :command:`vcf-call`.
* Combine optional arguments ``bam`` and ``fn`` into single positional argument ``bams`` for :meth:`pycov.CovFrame.from_bam` method. The same goes for :command:`bam-depth` command (combine ``--bam`` and ``--fn`` into ``bams``).
* Combine optional arguments ``bed`` and ``region`` into single optional argument ``regions`` for :meth:`pycov.CovFrame.from_bam` method. The same goes for :command:`bam-depth` command (combine ``--bed`` and ``--region`` into ``--regions``).
* Update :meth:`pycov.CovFrame.from_bam` method and :command:`bam-depth` command to automatically handle the 'chr' string.
* Rename :meth:`pyvcf.VcfFrame.variants` method to :meth:`pyvcf.VcfFrame.to_variants`.
* Add new optional arguments ``force`` and ``missing`` to :meth:`pyvcf.row_updateinfo` method.
* Add new method :meth:`pyvcf.gt_ploidy`.
* Update :meth:`pyvcf.gt_polyp` method to use :meth:`pyvcf.gt_ploidy` method internally.
* :issue:`53`: Add new methods to compute AC/AN/AF in the INFO column: :meth:`pyvcf.row_computeinfo` and :meth:`pyvcf.VcfFrame.compute_info`.
* :issue:`54`: Update :meth:`pyvcf.VcfFrame.cfilter_empty` method so that users can control missingness threshold for filtering samples.
* Rename :meth:`pyvcf.VcfFrame.cfilter_empty` method to :meth:`pyvcf.VcfFrame.empty_samples`.
* Update :meth:`common.sort_regions` method to support regions containing an ALT contig (e.g. chr16_KI270854v1_alt).

0.30.0 (2022-02-05)
-------------------

* Update :command:`fuc-find` command to allow users to control whether to use recursive retrieving.
* Add new command :command:`ngs-trim`.
* Add new command :command:`ngs-quant`.
* Add new submodule ``pykallisto``.
* Update :meth:`pycov.CovFrame.from_bam` method to use filename as sample name when the SM tag is missing.
* Add new method :meth:`pyvcf.row_phased`. From now on, it's used to get the ``pyvcf.VcfFrame.phased`` property.
* Add new method :meth:`pyvcf.split` and :command:`vcf-split` command for splitting VCF by individual.
* Update :meth:`pyvcf.merge` method, :meth:`pyvcf.VcfFrame.merge` method, and :command:`vcf-merge` command to automatically handle the 'chr' string.

0.29.0 (2021-12-19)
-------------------

* Add new property ``pyvcf.VcfFrame.phased``.
* Update :meth:`pyvcf.VcfFrame.slice` method to automatically handle the 'chr' string.
* Add new argument ``--thread`` to :command:`ngs-hc` command. This argument will be used to set ``--native-pair-hmm-threads`` for GATK's :command:`HaplotypeCaller` command, ``--reader-threads`` for GATK's :command:`GenomicsDBImport` command, and ``-XX:ParallelGCThreads`` and ``-XX:ConcGCThreads`` for Java.
* Add new argument ``--batch`` to :command:`ngs-hc` command. This argument will be used to set ``--batch-size`` for GATK's :command:`GenomicsDBImport` command.
* Update :command:`ngs-bam2fq` command to fix the SGE issue that outputs an error like ``Unable to run job: denied: "XXXXX" is not a valid object name (cannot start with a digit)``.
* Update :command:`ngs-hc` command so that when ``--posix`` is set, it will use ``--genomicsdb-shared-posixfs-optimizations`` argument from GATK's :command:`GenomicsDBImport` command in addition to exporting relevant shell variable (i.e. ``export TILEDB_DISABLE_FILE_LOCKING=1``).
* Add new argument ``--job`` to :command:`ngs-fq2bam` command.
* Update :command:`ngs-fq2bam` command so that BAM creation step and BAM processing step are now in one step.
* Update :command:`ngs-fq2bam` command so that ``--thread`` is now also used to set ``-XX:ParallelGCThreads`` and ``-XX:ConcGCThreads`` for Java.
* Add new method :meth:`common.parse_list_or_file`.

0.28.0 (2021-12-05)
-------------------

* Update :meth:`pyvcf.VcfFrame.filter_empty` method so that users can choose a varying number of missing genotypes as threshold.
* Add new method :meth:`pyvcf.plot_af_correlation`.
* Update :command:`bam-slice` command to support BED file as input for specifying regions. Additionally, from now on, the command will automatically handle the annoying 'chr' prefix.
* Add new method :meth:`pycov.CovFrame.matrix_uniformity`.
* Fix bug in :meth:`pyvcf.slice` method when the input region is missing start or end.
* Add new command :command:`ngs-bam2fq`.
* Add new command :command:`fa-filter`.
* Update :meth:`pycov.CovFrame.plot_region` and :meth:`pyvcf.VcfFrame.plot_region` methods to raise an error if the CovFrame/VcfFrame is empty.
* Update :meth:`pyvcf.VcfFrame.filter_*` methods so that they don't raise an error when the VcfFrame is empty (i.e. will return the empty VcfFrame).
* Update :meth:`common.plot_exons` method to not italicize text by default (use ``name='$text$'`` to italicize).
* Add new argument ``--posix`` to :command:`ngs-hc` command.
* Add new method :meth:`common.AnnFrame.subset`.
* Update :meth:`common.AnnFrame.plot_annot` method to raise an error if user provides an invalid group in ``group_order``.
* Add new method :meth:`pymaf.MafFrame.get_gene_concordance`.

0.27.0 (2021-11-20)
-------------------

* Rename ``file`` argument to ``vcf`` for :command:`vcf-slice` command.
* Add new command :command:`vcf-index`.
* Add new method :meth:`pyvcf.has_chr_prefix`.
* Add new command :meth:`common.update_chr_prefix`.
* Update :meth:`pyvcf.slice` method to automatically handle the 'chr' prefix.
* Fix bug caused by a typo in :meth:`pyvcf.VcfFrame.filter_sampany` method.

0.26.0 (2021-10-24)
-------------------

* Add new method :meth:`pybam.count_allelic_depth`.
* Update :meth:`common.parse_variant` method to handle position-only strings as input (e.g. '22-42127941-G-A' vs. '22-42127941').
* Add new command :command:`bam-aldepth`.
* Rename :meth:`pybam.has_chr` method to :meth:`pybam.has_chr_prefix`.
* Rename :meth:`pybed.BedFrame.chr_prefix`, :meth:`pycov.CovFrame.chr_prefix`, :meth:`pyvcf.VcfFrame.chr_prefix` methods to :meth:`pybed.BedFrame.update_chr_prefix`, :meth:`pycov.CovFrame.update_chr_prefix`, :meth:`pyvcf.VcfFrame.update_chr_prefix`.
* Add new properties ``pybed.BedFrame.has_chr_prefix``, ``pycov.CovFrame.has_chr_prefix``, ``pyvcf.VcfFrame.has_chr_prefix``.
* Add new method :meth:`pyvcf.slice`.
* Add new method :meth:`pyvcf.VcfFrame.from_string`.
* Remove ``nrows`` argument from :meth:`pyvcf.VcfFrame.from_file` method.
* Add new argument ``regions`` to :meth:`pyvcf.VcfFrame.from_file` method.
* Add new property ``pybed.BedFrame.shape``.
* Add new method :meth:`pybed.BedFrame.to_regions`.
* Add new method :meth:`pybed.BedFrame.from_regions`.
* Update :meth:`pyvcf.VcfFrame.from_file` method to accept BED data to specify regions of interest.
* Update :command:`vcf-slice` command to run significantly faster by allowing random access.
* Add new method :meth:`common.sort_regions`.
* Fix minor bug in :meth:`pyvcf.VcfFrame.get_af` method when the variant of interest does not exist in VcfFrame.

0.25.0 (2021-10-09)
-------------------

* Add new method :meth:`common.sort_variants`.
* Add new method :meth:`pyvcf.VcfFrame.variants`.
* Add new method :meth:`pymaf.MafFrame.variants`.
* Add new method :meth:`pymaf.MafFrame.subset`.
* Add new method :meth:`pymaf.MafFrame.calculate_concordance`.
* Add new method :meth:`pymaf.MafFrame.copy`.
* Add new method :meth:`pymaf.MafFrame.filter_indel`.
* Add new method :meth:`pymaf.MafFrame.plot_comparison`.

0.24.0 (2021-10-02)
-------------------

* Add new command :command:`fuc-bgzip`.
* Add new command :command:`tabix-index`.
* Fix bug in :meth:`pyvcf.VcfFrame.from_file` method when ``meta_only`` is ``True``.
* Update :meth:`pyvcf.VcfFrame.from_file` method to extract VCF headers as well when ``meta_only`` is ``True``.
* Add new command :command:`tabix-slice`.
* Update :meth:`pyvcf.VcfFrame.chr_prefix`, :meth:`pybed.BedFrame.chr_prefix`, and :meth:`pycov.CovFrame.chr_prefix` methods to skip lines that already have ``chr`` string when ``mode='add'``.
* Add new methods :meth:`common.rename` and :meth:`pycov.CovFrame.rename`.
* Add new command :command:`cov-rename`.
* Add new method :meth:`pyvcf.gt_het`.
* Add new method :meth:`pyvcf.gt_pseudophase`.

0.23.0 (2021-09-21)
-------------------

* Update :class:`pycov.CovFrame` class to ensure that the ``Chromosome`` column is always string.
* Update :meth:`pycov.CovFrame.from_file` method to accept file-like object as input as well.
* Add new argument ``metadata`` to :meth:`pyvcf.VcfFrame.strip` method.
* Update :meth:`pyvcf.VcfFrame.from_file` method to accept file-like object as input as well.
* Add new method :meth:`pycov.CovFrame.mask_bed`.
* Add new method :meth:`pycov.CovFrame.chr_prefix`.
* Add new property ``contigs`` to :class:`pybed.BedFrame` class.
* Add new method :meth:`pybed.BedFrame.chr_prefix`.
* Add new methods :meth:`pybed.BedFrame.copy_meta` and :meth:`pybed.BedFrame.sort`.
* Add new method :meth:`pybed.BedFrame.merge`.
* Add new property ``empty`` to :class:`pyvcf.VcfFrame` class.
* Fix minor bug in :meth:`pyvcf.VcfFrame.strip` method when sample genotypes don't have the same number of fields as FORMAT.
* Add new method :meth:`pycov.CovFrame.subset` method.
* Add new method :meth:`common.color_print`.
* Add new method :meth:`pycov.concat`.
* Add new command :command:`cov-concat`.
* Update :class:`pyvcf.VcfFrame` to enforce the dtypes.
* Update :meth:`pyvcf.VcfFrame.add_af` method to output allele fraction for each ALT allele.
* Fix bug in :meth:`pyvcf.VcfFrame.add_af` method when the sum of allelic depths is 0.
* Add new method :meth:`pyvcf.VcfFrame.get_af`.

0.22.0 (2021-09-04)
-------------------

* Update :meth:`pyvcf.VcfFrame.from_file` method to be more memory efficient by pre-specifying data type for each VCF column.
* Update :meth:`pyvcf.VcfFrame.from_file` method to raise error if one or more VCF columns are missing, except for the FORMAT column (i.e. "sites-only" VCFs).
* Add new property ``sites_only`` to :class:`pyvcf.VcfFrame`.
* Update :meth:`pyvcf.VcfFrame.merge` method to handle sites-only VCFs.
* Add new method :meth:`pyvcf.VcfFrame.filter_vcf`.
* Add new arguments ``--bed`` and ``--vcf`` to :command:`vcf-slice` command.
* Update :meth:`common.parse_region` method to output ``NaN`` instead of 0.
* Add new method :meth:`common.parse_variant`.
* Update :meth:`pycov.CovFrame.from_file` method to be more memory efficient by pre-specifying data type for each of the columns in the input text file.
* Update :meth:`pycov.CovFrame.from_file` method to raise error if 'Chromosome' or 'Position' column is missing.
* Add new method :meth:`pyvcf.VcfFrame.fetch`.
* Update :meth:`pyvcf.VcfFrame.strip` method to handle cases where one or more specified FORMAT keys are missing in a row.
* Add new method :meth:`pyvcf.VcfFrame.pseudophase`.
* Update :meth:`pyvcf.VcfFrame.filter_vcf` method to also use REF and ALT (previously it only used CHROM and POS).
* Add new argument ``--zero`` to :command:`bam-depth` command.
* Update :meth:`pycov.CovFrame.plot_region` method: 1) New argument ``label`` has been added. 2) Argument ``names`` has been deprecated. 3) New argument ``sample`` has been added. 4) From now on, by default the method will plot profile for single sample specified by ``sample`` as opposed to all samples at once. 5) From now on, argument ``region`` can be omitted if there is only one contig.
* Add new property ``contigs`` to :class:`pyvcf.CovFrame`.
* Add new methods :meth:`pyvcf.CovFrame.copy` and :meth:`pyvcf.CovFrame.copy_df`.
* Update :meth:`pyvcf.CovFrame.from_file` method to accept GZIP compressed files. Also add new argument ``compression``.
* Add new methods :meth:`pyvcf.CovFrame.to_string` and :meth:`pyvcf.CovFrame.to_file`.

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
