# CLI

## Table of contents

* [bed_intersect_files](#bed_intersect_files) 
* [bed_summarize_file](#bed_summarize_file) 
* [fastq_compute_read_lengths](#fastq_compute_read_lengths) 
* [fastq_count_reads](#fastq_count_reads) 
* [table_check_files](#table_check_files) 
* [table_compare_files](#table_compare_files) 
* [table_merge_files](#table_merge_files) 
* [table_summarize_column](#table_summarize_column) 
* [table_summarize_file](#table_summarize_file) 
* [vcf_merge_files](#vcf_merge_files) 

## bed_intersect_files <a name="bed_intersect_files"></a>

```
usage: fuc bed_intersect_files [-h] input_bed [input_bed ...] output_bed

This command will compute intersections beween multiple BED files.

positional arguments:
  input_bed   input BED files
  output_bed  output BED file

optional arguments:
  -h, --help  show this help message and exit
```

## bed_summarize_file <a name="bed_summarize_file"></a>

```
usage: fuc bed_summarize_file [-h] [--bases BASES] [--decimals DECIMALS]
                              bed_file

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

## fastq_compute_read_lengths <a name="fastq_compute_read_lengths"></a>

```
usage: fuc fastq_compute_read_lengths [-h] fastq_file

This command will compute the distribution of sequence read lengths for a
FASTQ file (both zipped and unqzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

## fastq_count_reads <a name="fastq_count_reads"></a>

```
usage: fuc fastq_count_reads [-h] fastq_file

This command will count sequence reads from a FASTQ file (both zipped and
unzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

## table_check_files <a name="table_check_files"></a>

```
usage: fuc table_check_files [-h] [--delimiter DELIMITER]
                             input_file column_header

This command will checks whether or not files in the given list exist inthe
operating system.

positional arguments:
  input_file            input file containing the list of file paths
  column_header         column header

optional arguments:
  -h, --help            show this help message and exit
  --delimiter DELIMITER
                        column delimiter (default: ',')
```

## table_compare_files <a name="table_compare_files"></a>

```
usage: fuc table_compare_files [-h] file1 file2

This command will compare two files.

positional arguments:
  file1       first file
  file2       second file

optional arguments:
  -h, --help  show this help message and exit
```

## table_merge_files <a name="table_merge_files"></a>

```
usage: fuc table_merge_files [-h] [--left_delimiter LEFT_DELIMITER]
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

## table_summarize_column <a name="table_summarize_column"></a>

```
usage: fuc table_summarize_column [-h] [--group_col GROUP_COL]
                                  [--delimiter DELIMITER]
                                  [--skiprows SKIPROWS [SKIPROWS ...]]
                                  [--rows ROWS] [--exclude_rows]
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

## table_summarize_file <a name="table_summarize_file"></a>

```
usage: fuc table_summarize_file [-h] [--delimiter DELIMITER]
                                [--columns COLUMNS [COLUMNS ...]]
                                [--exclude_columns] [--rows ROWS]
                                [--exclude_rows]
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

## vcf_merge_files <a name="vcf_merge_files"></a>

```
usage: fuc vcf_merge_files [-h] [--subfield SUBFIELD [SUBFIELD ...]]
                           input_vcf [input_vcf ...] output_vcf

This command merges multiple VCF files (both zipped and unzipped). By default,
only GT subfield of FORMAT field is included. Use '--subfield' to include
additional subfields such as AD and DP.

positional arguments:
  input_vcf             input VCF files
  output_vcf            output VCF file

optional arguments:
  -h, --help            show this help message and exit
  --subfield SUBFIELD [SUBFIELD ...]
                        FORMAT subfields
```

