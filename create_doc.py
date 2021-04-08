import subprocess

from common import fuc_dir

commands = ['check_files.py', 'merge_files.py', 'summarize_file.py',
            'merge_vcfs.py', 'summarize_bed.py', 'intersect_beds.py',
            'count_reads.py']

cli_file = f'{fuc_dir()}/doc/CLI.md'

with open(cli_file, 'w') as f:
    f.write('# API\n')
    f.write('\n')

for command in commands:
    args = ['python3', f'{fuc_dir()}/{command}', '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    with open(cli_file, 'a') as f:
        f.write(f'## {command}\n')
        f.write('\n')
        f.write('```\n')
        f.write(result.stdout)
        f.write('```\n')
        f.write('\n')
