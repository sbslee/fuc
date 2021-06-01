from .. import api
import matplotlib.pyplot as plt

DOC_URL = 'https://sbslee-fuc.readthedocs.io/en/latest/api.html#fuc.api.pymaf.MafFrame.plot_oncoplot'

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[MAF] create an oncoplot from a MAF file',
        description=
            'This command will create an oncoplot from a MAF file. '
            'The format of output image (PDF/PNG/JPEG/SVG) will be '
            "automatically determined by the output file's extension. "
            'This command essentially wraps the `pymaf.plot_oncoplot` '
            f"method. Visit the method's documentation ({DOC_URL}) to see "
            'example plots.'
    )
    parser.add_argument(
        'maf_file',
        help='input MAF file'
    )
    parser.add_argument(
        'output_file',
        help='output inage file'
    )
    parser.add_argument(
        '--count',
        metavar='INTEGER',
        type=int,
        default=10,
        help='number of top mutated genes to display (default: 10)'
    )
    parser.add_argument(
        '--figsize',
        metavar='FLOAT',
        type=float,
        default=[15, 10],
        nargs=2,
        help='width, height in inches (default: [15, 10])'
    )
    parser.add_argument(
        '--label_fontsize',
        metavar='FLOAT',
        type=float,
        default=15,
        help='font size of labels (default: 15)'
    )
    parser.add_argument(
        '--ticklabels_fontsize',
        metavar='FLOAT',
        type=float,
        default=15,
        help='font size of tick labels (default: 15)'
    )
    parser.add_argument(
        '--legend_fontsize',
        metavar='FLOAT',
        type=float,
        default=15,
        help='font size of legend texts (default: 15)'
    )

def main(args):
    mf = api.pymaf.MafFrame.from_file(args.maf_file)
    mf.plot_oncoplot(
        count=args.count,
        figsize=args.figsize,
        label_fontsize=args.label_fontsize,
        ticklabels_fontsize=args.ticklabels_fontsize,
        legend_fontsize=args.legend_fontsize
    )
    plt.savefig(args.output_file)
