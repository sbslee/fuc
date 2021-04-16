from fuc.api.common import get_script_name
from fuc.api.FastqFrame import FastqFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FASTQ] summarize a FASTQ file',
        description='This command will output a summary of the input FASTQ '
            'file (both zipped and unqzipped). The summary includes '
            'the total number of sequence reads, the distribution of read '
            'lengths, and the numbers of unique and duplicate sequences.'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')

def main(args):
    qf = FastqFrame.from_file(args.fastq_file)
    print(f'# Total: {qf.shape:,}')
    unique = len(set(qf.data))
    duplicate = qf.shape - unique
    print(f'# Unique: {unique:,}')
    print(f'# Duplicate: {duplicate:,}')
    lengths = qf.readlen()
    print('# Read length and count:')
    for length in sorted(lengths):
        print(length, f'{lengths[length]:,}', sep='\t')
