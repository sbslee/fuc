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
    print(f'# Total: {qf.shape}')
    unique = []
    duplicate = []
    for r in qf.data:
        if r in unique:
            duplicate.append(r)
        else:
            unique.append(r)
    print(f'# Unique: {len(unique)}')
    print(f'# Duplicate: {len(duplicate)}')
    lengths = qf.readlen()
    print('# Read length and count:')
    for length in sorted(lengths):
        print(length, lengths[length], sep='\t')
