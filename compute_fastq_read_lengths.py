import argparse
from api.FASTQResult import FASTQResult

def main():
    parser = argparse.ArgumentParser(description='This command will compute '
        'the distribution of sequence read lengths for a FASTQ file (both '
        'zipped and unzipped).')
    parser.add_argument('fastq_file', help='input FASTQ file')
    args = parser.parse_args()
    fastq_result = FASTQResult.read(args.fastq_file)
    lengths = fastq_result.get_read_length()
    for length in sorted(lengths):
        print(length, lengths[length], sep='\t')

if __name__ == '__main__':
    main()
