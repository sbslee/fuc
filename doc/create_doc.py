import subprocess

from fuc.api.common import fuc_dir
from fuc.cli import commands

readme_file = f'{fuc_dir()}/README.md'

with open(readme_file, 'w') as f:
    f.write('# README\n')
    f.write('\n')
    f.write('The `fuc` package is my attempt to wrap some of the most '
            'frequently used commands in the field of bioinformatics with '
            'pure Python 3 code, only using standard libraries. That means '
            'no installation for external Python packages, not even the '
            'popular ones like `numpy` and `pandas`. This also includes '
            'many famous bioinformatics tools such as `samtools` and '
            '`bedtools`. The motivation for not relying on external '
            'packages or programs is quite simple: I just got tired of '
            'managing different Python environments for doing simple '
            'things.\n')
    f.write('\n')
    f.write('You can use `fuc` for both command line interface (CLI) '
            'and application programming interface (API). Click '
            '[here](doc/CLI.md) to see the CLI documentation and '
            '[here](doc/API.md) to see the API documentation.\n')
    f.write('\n')
    f.write('To install `fuc`, enter the following in your terminal:\n')
    f.write('```\n')
    f.write('$ git clone https://github.com/sbslee/fuc\n')
    f.write('$ cd fuc\n')
    f.write('$ pip install .\n')
    f.write('```\n')
    f.write('\n')
    f.write('Your contributions (e.g. feature ideas, pull requests) '
            'are most welcome.\n')
    f.write('\n')
    f.write('Author: Seung-been "Steven" Lee<br/>\n')
    f.write('Email: sbstevenlee@gmail.com<br/>\n')
    f.write('License: MIT License\n')
    f.write('\n')

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
