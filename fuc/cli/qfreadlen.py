from fuc.api.common import get_script_name
from fuc.api.FastqFrame import FastqFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FASTQ] compute read lengths for a FASTQ file',
        description='This command will compute the distribution of sequence '
                    'read lengths for a FASTQ file (both zipped and '
                    'unqzipped).'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')

def main(args):
    fastq_result = FastqFrame.from_file(args.fastq_file)
    lengths = fastq_result.get_read_length()
    for length in sorted(lengths):
        print(length, lengths[length], sep='\t')
