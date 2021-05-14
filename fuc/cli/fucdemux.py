from fuc.api.common import get_script_name
import json

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FUC] parse Stats.json from bcl2fastq/bcl2fastq2',
        description='This command will parse the Stats.json file from '
            'the bcl2fastq/bcl2fastq2 program.'
    )
    parser.add_argument('json_file', help='Stats.json file')

def main(args):
    with open(args.json_file) as f:
        data = json.load(f)
    for sample in data['ConversionResults'][0]['DemuxResults']:
        print(sample)
