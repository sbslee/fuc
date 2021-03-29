import argparse
from common import VCFResult, read_vcf, merge_vcfs, write_vcf

def main():
    parser = argparse.ArgumentParser(description='This command merges '
        'multiple VCF files. By default, only GT subfield of FORMAT field '
        "is included. Use '--subfield' to include additional subfields "
        'such as AD and DP.')
    parser.add_argument('input_vcf', help='input VCF files', nargs='+')
    parser.add_argument('output_vcf', help='output VCF file')
    parser.add_argument('--subfield', help='FORMAT subfields',
        nargs='+')
    args = parser.parse_args()
    vcf_list = [read_vcf(x) for x in args.input_vcf]
    merged_vcf = VCFResult()
    for vcf in vcf_list:
        merged_vcf = merge_vcfs(merged_vcf, vcf, args.subfield)
    write_vcf(merged_vcf, args.output_vcf)

if __name__ == '__main__':
    main()
