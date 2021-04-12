import argparse
from api.VCFResult import VCFResult

def main():
    parser = argparse.ArgumentParser(description='This command merges '
        'multiple VCF files (both zipped and unzipped). By default, only '
        "GT subfield of FORMAT field is included. Use '--subfield' "
        'to include additional subfields such as AD and DP.')
    parser.add_argument('input_vcf', help='input VCF files', nargs='+')
    parser.add_argument('output_vcf', help='output VCF file')
    parser.add_argument('--subfield', help='FORMAT subfields',
        nargs='+')
    args = parser.parse_args()
    vcf_list = [VCFResult.read(x) for x in args.input_vcf]
    merged_vcf = VCFResult()
    for vcf in vcf_list:
        merged_vcf = merged_vcf.merge(vcf, args.subfield)
    merged_vcf.write(args.output_vcf)

if __name__ == '__main__':
    main()
