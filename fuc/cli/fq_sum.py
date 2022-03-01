from .. import api

description = """
Summarize a FASTQ file.

This command will output a summary of the input FASTQ file. The summary
includes the total number of sequence reads, the distribution of read
lengths, and the numbers of unique and duplicate sequences.
"""

epilog = f"""
[Example] Summarize a FASTQ file:
  $ fuc {api.common._script_name()} in.fastq
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Summarize a FASTQ file."""
    )
    parser.add_argument(
        'fastq',
        help=
"""Input FASTQ file (compressed or uncompressed)."""
    )

def main(args):
    qf = api.pyfq.FqFrame.from_file(args.fastq)
    print(f'# Total: {qf.shape[0]:,}')
    unique = qf.df.SEQ.nunique()
    duplicate = qf.shape[0] - unique
    print(f'# Unique: {unique:,}')
    print(f'# Duplicate: {duplicate:,}')
    lengths = qf.readlen()
    print('# Read length and count:')
    for length in sorted(lengths):
        print(length, f'{lengths[length]:,}', sep='\t')
