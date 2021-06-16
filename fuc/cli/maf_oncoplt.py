from .. import api
import matplotlib.pyplot as plt

description = f"""
This command will create an oncoplot with a MAF file. It essentially wraps
the 'pymaf.plot_oncoplot' method from the fuc API. Visit the method's
documentation to see example plots.

The format of output image (PDF/PNG/JPEG/SVG) will be automatically
determined by the output file's extension.

usage examples:
  $ fuc {api.common._script_name()} in.maf out.png
  $ fuc {api.common._script_name()} in.maf out.pdf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[MAF] create an oncoplot with a MAF file',
        description=description,
    )
    parser.add_argument('maf', help='MAF file')
    parser.add_argument('out', help='image file')
    parser.add_argument(
        '--count',
        metavar='INT',
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
    mf = api.pymaf.MafFrame.from_file(args.maf)
    mf.plot_oncoplot(
        count=args.count,
        figsize=args.figsize,
        label_fontsize=args.label_fontsize,
        ticklabels_fontsize=args.ticklabels_fontsize,
        legend_fontsize=args.legend_fontsize
    )
    plt.savefig(args.out)
