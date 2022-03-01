import sys

from .. import api

from Bio import bgzf

description = """
Write a BGZF compressed file.

BGZF (Blocked GNU Zip Format) is a modified form of gzip compression which
can be applied to any file format to provide compression with efficient
random access. In addition to being required for random access to and writing
of BAM files, the BGZF format can also be used for most of the sequence data
formats (e.g. FASTA, FASTQ, GenBank, VCF, MAF).
"""

epilog = f"""
[Example] When the input is a VCF file:
  $ fuc {api.common._script_name()} in.vcf > out.vcf.gz

[Example] When the input is stdin:
  $ cat in.vcf | fuc {api.common._script_name()} > out.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Write a BGZF compressed file."""
    )
    parser.add_argument(
        'file',
        nargs='*',
        help=
"""Input file to be compressed (default: stdin)."""
    )

def main(args):
    if args.file:
        with open(args.file[0]) as f:
            data = f.read()
    elif not sys.stdin.isatty():
        data = sys.stdin.read()
    else:
        raise ValueError('No input data detected')

    w = bgzf.BgzfWriter(fileobj=sys.stdout.buffer)
    w.write(data)
    w.close()
