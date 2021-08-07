from .. import api

description = f"""
This command will output a summary of the input FASTQ file (both zipped and unqzipped).

The summary includes the total number of sequence reads, the distribution of read lengths, and the numbers of unique and duplicate sequences.

Usage examples:
  $ fuc {api.common._script_name()} in.fastq
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Summarize a FASTQ file.',
        description=description,
    )
    parser.add_argument('fastq', help='FASTQ file.')

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
