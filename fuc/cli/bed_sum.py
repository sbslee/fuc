from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[BED] summarize a BED file',
        description=
            'This command will compute summary statstics of the '
            'BED file. This includes the total numbers of probes and '
            'covered base pairs for each chromosome. By default, covered '
            'base paris are displayed in bp, but if you prefer you can, '
            "for example, use '--bases 1000' to display base pairs in kb."
    )
    parser.add_argument('bed_file', help='input BED file')
    parser.add_argument('--bases', metavar='INTEGER', type=int, default=1,
        help='number used to divide the bases (default: 1)')
    parser.add_argument('--decimals', metavar='INTEGER', type=int, default=0,
        help='maximum number of decimals (default: 0)')

def main(args):
    bf = api.pybed.BedFrame.from_file(args.bed_file)
    chrom_dict = {}
    total = [0, 0]
    for i, r in bf.gr.df.iterrows():
        chrom = r['Chromosome']
        start = r['Start']
        end = r['End']
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
