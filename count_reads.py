import argparse
from FASTQResult import FASTQResult

def main():
    parser = argparse.ArgumentParser(description='This command will count '
        'sequence reads from a FASTQ file (both zipped and unzipped).')
    parser.add_argument('fastq_file', help='input FASTQ file')
    args = parser.parse_args()
    fastq_result = FASTQResult.read(args.fastq_file)
    print('Read count:', fastq_result.shape, sep='\t')

if __name__ == '__main__':
    main()
