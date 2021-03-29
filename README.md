# fuc
Frequently used commands

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
