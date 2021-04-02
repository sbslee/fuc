# fuc
This repository is my attempt to wrap some of the most frequently used commands in the field of bioinformatics with pure Python 3 code, only using standard libraries. That means no installation for external packages, not even the popular ones like `numpy` and `pandas`. The motivation for not relying on external packages is quite simple: I just got tired of managing different Python environments for doing simple things.

You can use Python scripts in this repository for both command line interface (CLI) and application programming interface (API). Note that I'm not calling this repository as a Python package as it's merely a collection of Python scripts. In other words, there is no installation (i.e. no `setup.py`). Just clone this repository and you are good to go. This repository not being a package gives me the flexibility to be as creative as I want. However, this may change in the future.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

# Usage

To clone this repository:

```
$ git clone https://github.com/sbslee/fuc
```

To run a command from the `fuc` repository (e.g. `check_files.py`):

```
$ python3 path/to/fuc/check_files.py input_file column_header
```

To access a module from the `fuc` repository (e.g. `VCFResult.py`):

```
import sys
sys.path.append('/path/to/fuc')
from VCFResult import VCFResult
```

# CLI

## check_files.py

```
$ python3 check_files.py -h
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
$ python3 merge_files.py -h
usage: merge_files.py [-h] left_file right_file output_file on

This command will merge two text files on a shared column.

positional arguments:
  left_file    left file
  right_file   right file
  output_file  merged file
  on           column name to join on

optional arguments:
  -h, --help   show this help message and exit
```

## merge_vcfs.py

```
$ python3 merge_vcfs.py -h
usage: merge_vcfs.py [-h] [--subfield SUBFIELD [SUBFIELD ...]]
                     input_vcf [input_vcf ...] output_vcf

This command merges multiple VCF files. By default, only GT subfield of FORMAT
field is included. Use '--subfield' to include additional subfields such as AD
and DP.

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
$ python3 summarize_bed.py -h
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
$ python3 intersect_beds.py -h
usage: intersect_beds.py [-h] input_bed [input_bed ...] output_bed

This command computes intersections between multiple BED files.

positional arguments:
  input_bed   input BED files
  output_bed  output BED file

optional arguments:
  -h, --help  show this help message and exit
```

# API

## DataFrame.py

The `DataFrame` module provides a suite of tools to work with text files.

## VCFResult.py

The `VCFResult` module provides a suite of tools to work with VCF files.

## BEDResult.py

The `BEDResult` module provides a suite of tools to work with BED files.
