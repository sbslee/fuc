# CLI

## Table of contents

* [Introduction](#Introduction)
* [Commands](#Commands)
	* [bfintxn](#bfintxn) 
	* [bfsum](#bfsum) 
	* [dfmerge](#dfmerge) 
	* [dfsum](#dfsum) 
	* [dfsumcol](#dfsumcol) 
	* [fuccompf](#fuccompf) 
	* [fucexist](#fucexist) 
	* [qfcount](#qfcount) 
	* [qfsum](#qfsum) 
	* [vfmerge](#vfmerge) 

## Introduction <a name="Introduction"></a>

For getting help on CLI:

```
$ fuc -h
usage: fuc [-h] [-v] COMMAND ...

positional arguments:
  COMMAND        name of the command
    bfintxn      [BED] find intersection of two or more BED files
    bfsum        [BED] summarize a BED file
    dfmerge      [TABLE] merge two text files
    dfsum        [TABLE] summarize a text file
    dfsumcol     [TABLE] summarize a column in a text file
    fuccompf     [FUC] compare two files
    fucexist     [FUC] check whether files/dirs exist
    qfcount      [FASTQ] count sequence reads in a FASTQ file
    qfsum        [FASTQ] summarize a FASTQ file
    vfmerge      [VCF] merge two or more VCF files

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show the version number and exit
```

For getting help on a specific command (e.g. `qfcount`):

```
$ fuc qfcount -h
usage: fuc qfcount [-h] fastq_file

This command will count sequence reads in a FASTQ file (both zipped and
unzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

## Commands <a name="Commands"></a>

### bfintxn <a name="bfintxn"></a>

```
usage: fuc bfintxn [-h] input_bed [input_bed ...] output_bed

This command will compute intersections beween multiple BED files.

positional arguments:
  input_bed   input BED files
  output_bed  output BED file

optional arguments:
  -h, --help  show this help message and exit
```

### bfsum <a name="bfsum"></a>

```
usage: fuc bfsum [-h] [--bases BASES] [--decimals DECIMALS] bed_file

This command computes summary statstics of the given BED file. This includes
the total numbers of probes and covered base pairs for each chromosome. By
default, covered base paris are displayed in bp, but if you prefer you can,
for example, use '--bases 1000' to display base pairs in kb.

positional arguments:
  bed_file             input BED file

optional arguments:
  -h, --help           show this help message and exit
  --bases BASES        number used to divide the bases (default: 1)
  --decimals DECIMALS  maximum number of decimals (default: 0)
```

### dfmerge <a name="dfmerge"></a>

```
usage: fuc dfmerge [-h] [--left_delimiter LEFT_DELIMITER]
                   [--right_delimiter RIGHT_DELIMITER]
                   [--output_delimiter OUTPUT_DELIMITER]
                   left_file right_file output_file on [on ...]

This command will merge two text files using one or more shared columns.

positional arguments:
  left_file             left file
  right_file            right file
  output_file           merged file
  on                    column names to join on

optional arguments:
  -h, --help            show this help message and exit
  --left_delimiter LEFT_DELIMITER
                        delimiter for the left file (default: '\t')
  --right_delimiter RIGHT_DELIMITER
                        delimiter for the right file (default: '\t')
  --output_delimiter OUTPUT_DELIMITER
                        delimiter for the output file (default: '\t')
```

### dfsum <a name="dfsum"></a>

```
usage: fuc dfsum [-h] [--delimiter DELIMITER]
                 [--columns COLUMNS [COLUMNS ...]] [--exclude_columns]
                 [--rows ROWS] [--exclude_rows]
                 table_file

This command will output a summary of the input text file. For each column, it
will return the counts of unique records for categorical data and the summary
statistics (minimum, maximum, mean, and median) for numeric data. You can use
'--columns' to specify which columns should be displayed. For filtering, you
can use '--rows' to express SQLite WHERE clause which will select rows that
meet certain criteria.

positional arguments:
  table_file            input table file

optional arguments:
  -h, --help            show this help message and exit
  --delimiter DELIMITER
                        delimiter for the table (default: '\t')
  --columns COLUMNS [COLUMNS ...]
                        specify which columns to summarize
  --exclude_columns     use this tag to exclude specified columns
  --rows ROWS           SQLite WHERE clause specifying which rows to summarize
  --exclude_rows        use this tag to exclude specified rows
```

### dfsumcol <a name="dfsumcol"></a>

```
usage: fuc dfsumcol [-h] [--group_col GROUP_COL] [--delimiter DELIMITER]
                    [--skiprows SKIPROWS [SKIPROWS ...]] [--rows ROWS]
                    [--exclude_rows]
                    table_file target_col

This command will output a summary table for the target column in the input
text file. The target column must be categorical. You can also use '--
group_col' to group the observations by another categorical column. For
filtering, you can use '--rows' to express SQLite WHERE clause which will
select rows that meet certain criteria.

positional arguments:
  table_file            input table file
  target_col            target column

optional arguments:
  -h, --help            show this help message and exit
  --group_col GROUP_COL
                        column to group by
  --delimiter DELIMITER
                        delimiter for the table (default: '\t')
  --skiprows SKIPROWS [SKIPROWS ...]
                        line numbers to skip
  --rows ROWS           SQLite WHERE clause specifying which rows to summarize
  --exclude_rows        use this tag to exclude specified rows
```

### fuccompf <a name="fuccompf"></a>

```
usage: fuc fuccompf [-h] file1 file2

This command will compare two files.

positional arguments:
  file1       first file
  file2       second file

optional arguments:
  -h, --help  show this help message and exit
```

### fucexist <a name="fucexist"></a>

```
usage: fuc fucexist [-h] [paths [paths ...]]

This command will check whether files/dirs exist. It will look for stdin if
there are no arguments (e.g. $ cat files.list | fuc fucexist).

positional arguments:
  paths       file/dir paths (default: stdin)

optional arguments:
  -h, --help  show this help message and exit
```

### qfcount <a name="qfcount"></a>

```
usage: fuc qfcount [-h] fastq_file

This command will count sequence reads in a FASTQ file (both zipped and
unzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

### qfsum <a name="qfsum"></a>

```
usage: fuc qfsum [-h] fastq_file

This command will output a summary of the input FASTQ file (both zipped and
unqzipped). The summary includes the total number of sequence reads, the
distribution of read lengths, and the numbers of unique and duplicate
sequences.

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

### vfmerge <a name="vfmerge"></a>

```
usage: fuc vfmerge [-h]
                   [--format_subfields FORMAT_SUBFIELDS [FORMAT_SUBFIELDS ...]]
                   input_vcf [input_vcf ...]

This command will merge multiple VCF files (both zipped and unzipped). By
default, only the GT subfield of the FORMAT field will be included in the
merged VCF. Use '--format_subfields' to include additional FORMAT subfields
such as AD and DP.

positional arguments:
  input_vcf             input VCF files

optional arguments:
  -h, --help            show this help message and exit
  --format_subfields FORMAT_SUBFIELDS [FORMAT_SUBFIELDS ...]
                        FORMAT subfields
```

