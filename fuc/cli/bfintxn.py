from fuc.api.common import get_script_name
from fuc.api.BedFrame import BedFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[BED] find intersection of two or more BED files',
        description='This command will compute intersections beween '
            'multiple BED files.'
    )
    parser.add_argument('input_bed', help='input BED files', nargs='+')
    parser.add_argument('output_bed', help='output BED file')

def main(args):
    bed_list = [BedFrame.from_read(x) for x in args.input_bed]
    final_bed = bed_list[0]
    for bed in bed_list[1:]:
        final_bed = final_bed.intersect(bed)
    final_bed.write(args.output_bed)