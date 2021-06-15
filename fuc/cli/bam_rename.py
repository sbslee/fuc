from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script(),
        help='[BAM] rename the samples in a SAM/BAM/CRAM file',
        description='This command will rename the sample in a SAM/BAM/CRAM file.'
    )
    parser.add_argument('input_bam', help='SAM/BAM/CRAM file')
    parser.add_argument('name', help='sample name')
    parser.add_argument('output_bam', help='output BAM file')

def main(args):
    api.pybam.rename(args.input_bam, args.name, args.output_bam)
