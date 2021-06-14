from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[MAF] convert a MAF file to a VCF file',
        description='This command will convert a MAF file to a VCF file.'
    )
    parser.add_argument('maf', help='MAF file')
    parser.add_argument('--fasta', metavar='PATH', help='FASTA file')
    parser.add_argument('--ignore_indels', action='store_true', help='use this tag to exclude INDELs from the output')

def main(args):
    mf = api.pymaf.MafFrame.from_file(args.maf)
    vf = mf.to_vcf(fasta=args.fasta, ignore_indels=args.ignore_indels)
    print(vf.to_string())
