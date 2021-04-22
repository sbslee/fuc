CLI
***

Introduction
============

This section describes command line interface (CLI) for the ``fuc`` package.

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
       fuccompf     [FUC] compare two files
       fucexist     [FUC] check whether files/dirs exist
       qfcount      [FASTQ] count sequence reads in FASTQ files
       qfsum        [FASTQ] summarize a FASTQ file
       vfmerge      [VCF] merge two or more VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     -v, --version  show the version number and exit

For getting help on a specific command (e.g. `vfmerge`):

.. code-block:: console

   $ fuc vfmerge -h
   usage: fuc vfmerge [-h] [--how TEXT] [--format TEXT] vcf_files [vcf_files ...]
   
   This command will merge multiple VCF files (both zipped and unzipped). By
   default, only the GT subfield of the FORMAT field will be included in the
   merged VCF. Use '--format' to include additional FORMAT subfields such as AD
   and DP.
   
   positional arguments:
     vcf_files      VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     --how TEXT     type of merge as defined in `pandas.DataFrame.merge`
                    (default: 'inner')
     --format TEXT  FORMAT subfields to be retained (e.g. 'GT:AD:DP') (default:
                    'GT')

bfintxn
=======

.. code-block:: console

   $ fuc bfintxn -h
   usage: fuc bfintxn [-h] bed_files [bed_files ...]
   
   This command will compute intersections beween multiple BED files. It
   essentially wraps the `pyranges.PyRanges.intersect` method.
   
   positional arguments:
     bed_files   BED files
   
   optional arguments:
     -h, --help  show this help message and exit

bfsum
=====

.. code-block:: console

   $ fuc bfsum -h
   usage: fuc bfsum [-h] [--bases INTEGER] [--decimals INTEGER] bed_file
   
   This command will compute summary statstics of the BED file. This includes the
   total numbers of probes and covered base pairs for each chromosome. By
   default, covered base paris are displayed in bp, but if you prefer you can,
   for example, use '--bases 1000' to display base pairs in kb.
   
   positional arguments:
     bed_file            input BED file
   
   optional arguments:
     -h, --help          show this help message and exit
     --bases INTEGER     number used to divide the bases (default: 1)
     --decimals INTEGER  maximum number of decimals (default: 0)

dfmerge
=======

.. code-block:: console

   $ fuc dfmerge -h
   usage: fuc dfmerge [-h] [--how TEXT] [--on TEXT [TEXT ...]]
                      [--left_delimiter TEXT] [--right_delimiter TEXT]
                      [--output_delimiter TEXT]
                      left_file right_file output_file
   
   This command will merge two text files using one or more shared columns. This
   essentially wraps the `pandas.DataFrame.merge` method.
   
   positional arguments:
     left_file             left file
     right_file            right file
     output_file           output file
   
   optional arguments:
     -h, --help            show this help message and exit
     --how TEXT            type of merge to be performed ['left', 'right',
                           'outer', 'inner', 'cross'] (default: 'inner')
     --on TEXT [TEXT ...]  column names to join on
     --left_delimiter TEXT
                           left delimiter (default: '\t')
     --right_delimiter TEXT
                           right delimiter (default: '\t')
     --output_delimiter TEXT
                           output delimiter (default: '\t')

dfsum
=====

.. code-block:: console

   $ fuc dfsum -h
   usage: fuc dfsum [-h] [--delimiter TEXT] text_file
   
   This command will summarize a text file. It essentially wraps the
   `pandas.DataFrame.describe` method.
   
   positional arguments:
     text_file         text file
   
   optional arguments:
     -h, --help        show this help message and exit
     --delimiter TEXT  delimiter (default: '\t')

fuccompf
========

.. code-block:: console

   $ fuc fuccompf -h
   usage: fuc fuccompf [-h] file1 file2
   
   This command will compare two files.
   
   positional arguments:
     file1       first file
     file2       second file
   
   optional arguments:
     -h, --help  show this help message and exit

fucexist
========

.. code-block:: console

   $ fuc fucexist -h
   usage: fuc fucexist [-h] [paths ...]
   
   This command will check whether files/dirs exist. It will look for stdin if
   there are no arguments (e.g. $ cat files.list | fuc fucexist).
   
   positional arguments:
     paths       file/dir paths (default: stdin)
   
   optional arguments:
     -h, --help  show this help message and exit

qfcount
=======

.. code-block:: console

   $ fuc qfcount -h
   usage: fuc qfcount [-h] [paths ...]
   
   This command will count sequence reads in FASTQ files (both zipped and
   unzipped). It will look for stdin if there are no arguments (e.g. $ cat
   files.list | fuc qfcount).
   
   positional arguments:
     paths       FASTQ file paths (default: stdin)
   
   optional arguments:
     -h, --help  show this help message and exit

qfsum
=====

.. code-block:: console

   $ fuc qfsum -h
   usage: fuc qfsum [-h] fastq_file
   
   This command will output a summary of the input FASTQ file (both zipped and
   unqzipped). The summary includes the total number of sequence reads, the
   distribution of read lengths, and the numbers of unique and duplicate
   sequences.
   
   positional arguments:
     fastq_file  input FASTQ file
   
   optional arguments:
     -h, --help  show this help message and exit

vfmerge
=======

.. code-block:: console

   $ fuc vfmerge -h
   usage: fuc vfmerge [-h] [--how TEXT] [--format TEXT] vcf_files [vcf_files ...]
   
   This command will merge multiple VCF files (both zipped and unzipped). By
   default, only the GT subfield of the FORMAT field will be included in the
   merged VCF. Use '--format' to include additional FORMAT subfields such as AD
   and DP.
   
   positional arguments:
     vcf_files      VCF files
   
   optional arguments:
     -h, --help     show this help message and exit
     --how TEXT     type of merge as defined in `pandas.DataFrame.merge`
                    (default: 'inner')
     --format TEXT  FORMAT subfields to be retained (e.g. 'GT:AD:DP') (default:
                    'GT')

