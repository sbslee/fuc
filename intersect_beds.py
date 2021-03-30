import argparse
from BEDResult import BEDResult

def main():
    parser = argparse.ArgumentParser(description='This command computes '
        'intersections between multiple BED files.')
    parser.add_argument('input_bed', help='input BED files', nargs='+')
    parser.add_argument('output_bed', help='output BED file')
    args = parser.parse_args()
    bed_list = [BEDResult.read(x) for x in args.input_bed]
    final_bed = bed_list[0]
    for bed in bed_list[1:]:
        final_bed = final_bed.intersect(bed)
    final_bed.write(args.output_bed)

if __name__ == '__main__':
    main()
