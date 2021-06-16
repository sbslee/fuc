from .. import api
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

description = f"""
This command will parse the Reports directory from the bcl2fastq or
bcl2fastq2 prograrm. In the output directory, the command will create four
files:

- flowcell_summary.csv
- lane_summary.csv
- top_unknown_barcodes.csv
- reports.pdf

usage examples:
  $ fuc {api.common._script_name()} reports_dir output_dir
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] parse Reports directory from bcl2fastq or bcl2fastq2',
        description=description,
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

    fig, axes = plt.subplots(3, 2, figsize=(14, 17))

    kwargs = dict(data=df2[:-1], kde=True)

    sns.histplot(ax=axes[0][0], x='PF Clusters', **kwargs)
    sns.histplot(ax=axes[0][1], x='% of the lane', **kwargs)
    sns.histplot(ax=axes[1][0], x='Yield (Mbases)', **kwargs)
    sns.histplot(ax=axes[1][1], x='% PF Clusters', **kwargs)
    sns.histplot(ax=axes[2][0], x='% >= Q30 bases', **kwargs)
    sns.histplot(ax=axes[2][1], x='Mean Quality Score', **kwargs)

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
