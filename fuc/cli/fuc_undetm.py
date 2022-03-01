import gzip

from .. import api

description = """
Compute top unknown barcodes using undertermined FASTQ from bcl2fastq.

This command will compute top unknown barcodes using undertermined FASTQ from
the bcl2fastq or bcl2fastq2 prograrm.
"""

epilog = f"""
[Example] Compute top unknown barcodes:
  $ fuc {api.common._script_name()} Undetermined_S0_R1_001.fastq.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compute top unknown barcodes using undertermined FASTQ from
bcl2fastq."""
    )
    parser.add_argument(
        'fastq',
        help=
"""Undertermined FASTQ (compressed or uncompressed)."""
    )
    parser.add_argument(
        '--count',
        metavar='INT',
        type=int,
        default=30,
        help=
"""Number of top unknown barcodes to return (default: 30)."""
    )

def main(args):
    if args.fastq.endswith('.gz'):
        f = gzip.open(args.fastq, 'rt')
    else:
        f = open(args.fastq)

    data = {}

    for line in f:
        if not line.startswith('@'):
            continue
        fields = line.strip().split(':')
        name = fields[-1]
        if name not in data:
            data[name] = 0
        data[name] += 1

    f.close()

    data = dict(sorted(data.items(), key=lambda x: x[1], reverse=True))

    for i, name in enumerate(data):
        if i < 30:
            print(name, data[name])
