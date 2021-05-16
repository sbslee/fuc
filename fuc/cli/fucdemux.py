from fuc.api.common import get_script_name
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FUC] parse the Reports directory from bcl2fastq or bcl2fastq2',
        description=
            'This command will parse the Reports directory from the '
            'bcl2fastq or bcl2fastq2 prograrm. In the output directory, the '
            'command will create four files: flowcell_summary.csv, '
            'lane_summary.csv, top_unknown_barcodes.csv, and reports.pdf.'
    )
    parser.add_argument('reports_dir', help='Reports directory')
    parser.add_argument('output_dir', help='output directory')

def main(args):
    for path in Path(args.reports_dir).rglob('all/all/all/laneBarcode.html'):
        html_file = path.absolute()

    dfs = pd.read_html(html_file)

    d = {
        '% of thelane': '% of the lane',
        '% Perfectbarcode': '% Perfect barcode',
        '% One mismatchbarcode': '% One mismatch barcode',
        '% PFClusters': '% PF Clusters',
        '% >= Q30bases': '% >= Q30 bases',
        'Mean QualityScore': 'Mean Quality Score',
    }

    df1 = dfs[1]
    df2 = dfs[2].rename(columns=d)
    df3 = dfs[3].dropna()

    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(15, 15))

    sns.histplot(data=df2[:-1], x='PF Clusters', ax=ax1, kde=True)
    sns.histplot(data=df2[:-1], x='% of the lane', ax=ax2, kde=True)
    sns.histplot(data=df2[:-1], x='% >= Q30 bases', ax=ax3, kde=True)
    sns.histplot(data=df2[:-1], x='Mean Quality Score', ax=ax4, kde=True)

    os.mkdir(args.output_dir)

    plt.savefig(f'{args.output_dir}/reports.pdf')

    df1_string = df1.to_csv(index=False)
    df2_string = df2.to_csv(index=False, na_rep='NaN')
    df3_string = df3.to_csv(index=False)

    with open(f'{args.output_dir}/flowcell_summary.csv', 'w') as f:
        f.write(df1_string)
    with open(f'{args.output_dir}/lane_summary.csv', 'w') as f:
        f.write(df2_string)
    with open(f'{args.output_dir}/top_unknown_barcodes.csv', 'w') as f:
        f.write(df3_string)
