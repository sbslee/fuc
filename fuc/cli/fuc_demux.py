import os
from pathlib import Path

from .. import api

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

description = """
Parse the Reports directory from bcl2fastq.

This command will parse, and extract various statistics from, HTML files in
the Reports directory created by the bcl2fastq or bcl2fastq2 prograrm. After
creating an output directory, the command will write the following files:
  - flowcell-summary.csv
  - lane-summary.csv
  - top-unknown-barcodes.csv
  - reports.pdf

Use --sheet to sort samples in the lane-summary.csv file in the same order
as your SampleSheet.csv file. You can also provide a modified version of your
SampleSheet.csv file to subset samples for the lane-summary.csv and
reports.pdf files.
"""

epilog = f"""
[Example] Parse a Reports directory:
  $ fuc {api.common._script_name()} Reports output

[Example] Sort and/or subset samples:
  $ fuc {api.common._script_name()} Reports output --sheet SampleSheet.csv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        help=
"""Parse the Reports directory from bcl2fastq."""
    )
    parser.add_argument(
        'reports',
        help=
"""Input Reports directory."""
    )
    parser.add_argument(
        'output',
        help=
"""Output directory (will be created)."""
    )
    parser.add_argument(
        '--sheet',
        metavar='PATH',
        help=
"""SampleSheet.csv file. Used for sorting and/or subsetting
samples."""
    )

def main(args):
    for path in Path(args.reports).rglob('all/all/all/laneBarcode.html'):
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
    df2.Sample = df2.Sample.astype(str)
    df3 = dfs[3].dropna()

    if args.sheet is None:
        df4 = pd.DataFrame()
    else:
        with open(args.sheet) as f:
            for i, line in enumerate(f):
                if line.startswith('[Data]'):
                    n = i + 1
        df4 = pd.read_csv(args.sheet, skiprows=n)
        df4.Sample_Name = df4.Sample_Name.astype(str)
        sample_order = df4.Sample_Name.to_list() + ['Undetermined']
        df2.index = df2.Sample
        df2 = df2.loc[sample_order]
        df2 = df2.reset_index(drop=True)

    fig, axes = plt.subplots(3, 2, figsize=(14, 17))

    kwargs = dict(data=df2[:-1], kde=True)

    sns.histplot(ax=axes[0][0], x='PF Clusters', **kwargs)
    sns.histplot(ax=axes[0][1], x='% of the lane', **kwargs)
    sns.histplot(ax=axes[1][0], x='Yield (Mbases)', **kwargs)
    sns.histplot(ax=axes[1][1], x='% PF Clusters', **kwargs)
    sns.histplot(ax=axes[2][0], x='% >= Q30 bases', **kwargs)
    sns.histplot(ax=axes[2][1], x='Mean Quality Score', **kwargs)

    os.mkdir(args.output)

    plt.tight_layout()
    plt.savefig(f'{args.output}/reports.pdf')

    df1_string = df1.to_csv(index=False)
    df2_string = df2.to_csv(index=False, na_rep='NaN')
    df3_string = df3.to_csv(index=False)

    with open(f'{args.output}/flowcell-summary.csv', 'w') as f:
        f.write(df1_string)
    with open(f'{args.output}/lane-summary.csv', 'w') as f:
        f.write(df2_string)
    with open(f'{args.output}/top-unknown-barcodes.csv', 'w') as f:
        f.write(df3_string)
