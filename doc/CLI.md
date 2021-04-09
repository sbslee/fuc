# API

## Table of contents

* [check_files.py](#check_files.py) 
* [merge_files.py](#merge_files.py) 
* [summarize_file.py](#summarize_file.py) 
* [merge_vcfs.py](#merge_vcfs.py) 
* [summarize_bed.py](#summarize_bed.py) 
* [intersect_beds.py](#intersect_beds.py) 
* [count_reads.py](#count_reads.py) 

## check_files.py

```
usage: check_files.py [-h] [--delimiter DELIMITER] input_file column_header

This command checks whether or not files in the given list exist in the
operating system.

positional arguments:
  input_file            input file containing the list of file paths
  column_header         column header

optional arguments:
  -h, --help            show this help message and exit
  --delimiter DELIMITER
                        column delimiter (default: ',')
```

## merge_files.py

```
usage: merge_files.py [-h] [--left_delimiter LEFT_DELIMITER]
                      [--right_delimiter RIGHT_DELIMITER]
                      [--output_delimiter OUTPUT_DELIMITER]
                      left_file right_file output_file on [on ...]

This command will merge two text files on a shared column.

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

## summarize_file.py

```
usage: summarize_file.py [-h] [--delimiter DELIMITER]
                         [--columns COLUMNS [COLUMNS ...]] [--exclude_columns]
                         [--rows ROWS] [--exclude_rows]
                         table_file

This command will output a summary of the input text file. For each column, it
will return the counts of unique records for categorical data and the summary
statistics (minimum, maximum, mean, and median) for numeric data. You can use
'--columns' to specify which columns should be displayed. Additionally, you
can use '--rows' to filter for rows that meet certain criteria.

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

## merge_vcfs.py

```
usage: merge_vcfs.py [-h] [--subfield SUBFIELD [SUBFIELD ...]]
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

## summarize_bed.py

```
usage: summarize_bed.py [-h] [--bases BASES] [--decimals DECIMALS] bed_file

This command computes summary statstics of the given BED file. This includes
the total numbers of probes and covered base pairs for each chromosome. By
default, covered base paris are displayed in bp, but if you prefer you can,
for example, use '--bases 1000' to display base pairs in kb.

positional arguments:
  bed_file             input BED file

optional arguments:
  -h, --help           show this help message and exit
  --bases BASES        number used to divide the bases (default: 1)
  --decimals DECIMALS  maximum number of decimals (default: 10)
```

## intersect_beds.py

```
usage: intersect_beds.py [-h] input_bed [input_bed ...] output_bed

This command computes intersections between multiple BED files.

positional arguments:
  input_bed   input BED files
  output_bed  output BED file

optional arguments:
  -h, --help  show this help message and exit
```

## count_reads.py

```
usage: count_reads.py [-h] fastq_file

This command will count sequence reads from a FASTQ file (both zipped and
unzipped).

positional arguments:
  fastq_file  input FASTQ file

optional arguments:
  -h, --help  show this help message and exit
```

