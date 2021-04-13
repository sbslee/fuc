import subprocess

from fuc.api.common import fuc_dir
from fuc.cli import commands

cli_file = f'{fuc_dir()}/doc/CLI.md'

with open(cli_file, 'w') as f:
    f.write('# CLI\n')
    f.write('\n')
    f.write('## Table of contents\n')
    f.write('\n')
    for command in commands:
        f.write(f'* [{command}](#{command}) \n')
    f.write('\n')

for command in commands:
    args = ['fuc', command, '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    with open(cli_file, 'a') as f:
        f.write(f'## {command} <a name="{command}"></a>\n')
        f.write('\n')
        f.write('```\n')
        f.write(result.stdout)
        f.write('```\n')
        f.write('\n')
