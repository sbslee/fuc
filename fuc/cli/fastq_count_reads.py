from fuc.api.common import get_script_name
from fuc.api.FASTQResult import FASTQResult

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='count sequence reads in FASTQ',
        description='This command will count sequence reads from a FASTQ '
                    'file (both zipped and unzipped).'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')
    return parser

def main(args):
    fastq_result = FASTQResult.read(args.fastq_file)
    print(fastq_result.shape)
