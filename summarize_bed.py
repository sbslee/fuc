import argparse
from BEDResult import BEDResult

def main():
    parser = argparse.ArgumentParser(description='This command computes '
        'summary statstics for a BED file.')
    parser.add_argument('bed_file', help='input BED file')
    args = parser.parse_args()
    bed_result = BEDResult.read(args.bed_file)
    chrom_dict = {}
    total = [0, 0]
    for fields in bed_result.get_data():
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        bases = end - start
        if chrom not in chrom_dict:
            chrom_dict[chrom] = [0, 0]
        chrom_dict[chrom][0] += 1
        chrom_dict[chrom][1] += bases
        total[0] += 1
        total[1] += bases
    print('Chrom', 'Count', 'Bases', sep='\t')
    for chrom in chrom_dict:
        results = chrom_dict[chrom]
        print(chrom, results[0], results[1], sep='\t')
    print('Total', total[0], total[1], sep='\t')

if __name__ == '__main__':
    main()
