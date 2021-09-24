import sys

from .. import api

from Bio import bgzf

description = f"""
#################################
# Write a BGZF compressed file. #
#################################

BGZF (Blocked GNU Zip Format) is a modified form of gzip compression which can be applied to any file format to provide compression with efficient random access.

In addition to being required for random access to and writing of BAM files, the BGZF format can also be used for most of the sequence data formats (e.g. FASTA, FASTQ, GenBank, VCF, MAF).

Usage examples:
  $ fuc {api.common._script_name()} in.vcf > in.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Write a BGZF compressed file.',
        description=description,
    )
    parser.add_argument(
        'file',
        help='File to be compressed.'
    )

def main(args):
    w = bgzf.BgzfWriter(fileobj=sys.stdout.buffer)

    with open(args.file) as f:
        data = f.read()

    w.write(data)
    w.close()
