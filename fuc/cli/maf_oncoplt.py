from .. import api
import matplotlib.pyplot as plt

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[MAF] create an oncoplot from a MAF file',
        description='This command will create an oncoplot from a MAF file.'
    )
    parser.add_argument('maf_file', help='input MAF file')
    parser.add_argument('plot_file', help='output plot file')

def main(args):
    mf = api.pymaf.MafFrame.from_file(args.maf_file)
    mf.plot_oncoplot(fontsize=14)
    plt.savefig(args.plot_file)
