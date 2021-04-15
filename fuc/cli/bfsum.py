from fuc.api.common import get_script_name
from fuc.api.BedFrame import BedFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[BED] summarize a BED file',
        description='This command computes summary statstics of the given '
            'BED file. This includes the total numbers of probes and '
            'covered base pairs for each chromosome. By default, covered '
            'base paris are displayed in bp, but if you prefer you can, '
            "for example, use '--bases 1000' to display base pairs in kb."
    )
    parser.add_argument('bed_file', help='input BED file')
    parser.add_argument('--bases', type=int, default=1, help='number used '
        'to divide the bases (default: 1)')
    parser.add_argument('--decimals', type=int, default=0, help='maximum '
        'number of decimals (default: 0)')

def main(args):
    bed_result = BedFrame.from_file(args.bed_file)
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
    print('Chrom', 'Probes', 'Bases', sep='\t')
    for chrom in chrom_dict:
        results = chrom_dict[chrom]
        probes = results[0]
        bases = f'{results[1]/args.bases:.{args.decimals}f}'
        print(chrom, probes, bases, sep='\t')
    probes = total[0]
    bases = f'{total[1]/args.bases:.{args.decimals}f}'
    print('Total', probes, bases, sep='\t')
