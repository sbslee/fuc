from fuc.api.common import get_script_name
from fuc.api.FastqFrame import FastqFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FASTQ] count sequence reads in a FASTQ file',
        description='This command will count sequence reads in a FASTQ '
                    'file (both zipped and unzipped).'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')
    return parser

def main(args):
    fastq_result = FastqFrame.from_file(args.fastq_file)
    print(fastq_result.shape)
