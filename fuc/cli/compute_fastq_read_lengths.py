from fuc.api.FASTQResult import FASTQResult

def create_parser(subparsers):
    parser = subparsers.add_parser(
        'compute_fastq_read_lengths',
        help='compute read lengths in FASTQ',
        description='This command will compute the distribution of sequence '
                    'read lengths for a FASTQ file (both zipped and '
                    'unqzipped).'
    )
    parser.add_argument('fastq_file', help='input FASTQ file')

def main(args):
    fastq_result = FASTQResult.read(args.fastq_file)
    lengths = fastq_result.get_read_length()
    for length in sorted(lengths):
        print(length, lengths[length], sep='\t')
