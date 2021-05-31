from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[FASTQ] summarize a FASTQ file',
        description=
            'This command will output a summary of the input FASTQ '
            'file (both zipped and unqzipped). The summary includes '
            'the total number of sequence reads, the distribution of read '
            'lengths, and the numbers of unique and duplicate sequences.'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')

def main(args):
    qf = api.pyfq.FqFrame.from_file(args.fastq_file)
    print(f'# Total: {qf.shape[0]:,}')
    unique = qf.df.SEQ.nunique()
    duplicate = qf.shape[0] - unique
    print(f'# Unique: {unique:,}')
    print(f'# Duplicate: {duplicate:,}')
    lengths = qf.readlen()
    print('# Read length and count:')
    for length in sorted(lengths):
        print(length, f'{lengths[length]:,}', sep='\t')
