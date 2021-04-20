import subprocess
import pydoc

from fuc.api.common import fuc_dir
from fuc.cli import commands
from fuc import VcfFrame
import fuc

modules = [x for x in dir(fuc) if x not in ['api', 'cli'] and '__' not in x]
readme_file = f'{fuc_dir()}/README.md'
cli_file = f'{fuc_dir()}/doc/CLI.md'
api_file = f'{fuc_dir()}/doc/API.md'

# -- README.md ---------------------------------------------------------------

with open(readme_file, 'w') as f:
    f.write('# README\n')
    f.write('\n')
    f.write('## Introduction\n')
    f.write('\n')
    f.write('The main goal of the `fuc` package is to wrap some of the most '
            'frequently used commands in the field of bioinformatics '
            'into one place.\n')
    f.write('\n')
    f.write('You can use `fuc` for both command line interface (CLI) '
            'and application programming interface (API). Click '
            '[here](doc/CLI.md) to see the CLI documentation and '
            '[here](doc/API.md) to see the API documentation.\n')
    f.write('\n')
    f.write('Your contributions (e.g. feature ideas, pull requests) '
            'are most welcome.\n')
    f.write('\n')
    f.write('Author: Seung-been "Steven" Lee<br/>\n')
    f.write('Email: sbstevenlee@gmail.com<br/>\n')
    f.write('License: MIT License\n')
    f.write('\n')
    f.write('## Examples\n')
    f.write('\n')
    f.write('To merge VCF files with CLI:\n')
    f.write('\n')
    f.write('```\n')
    f.write('$ fuc vfmerge 1.vcf 2.vcf 3.vcf > merged.vcf\n')
    f.write('```\n')
    f.write('\n')
    f.write('To filter a VCF file based on a BED file using API:\n')
    f.write('\n')
    f.write('```\n')
    f.write('from fuc.api.VcfFrame import VcfFrame\n')
    f.write("vf = VcfFrame.from_file('original.vcf')\n")
    f.write("filtered_vf = vf.filter_bed('1.bed')\n")
    f.write("filtered_vf.to_file('filtered.vcf')\n")
    f.write('```\n')
    f.write('\n')
    f.write('## Required Packages\n')
    f.write('\n')
    f.write('The following packages are required to run `fuc`:\n')
    f.write('\n')
    f.write('```\n')
    f.write('numpy\n')
    f.write('pandas\n')
    f.write('pyranges\n')
    f.write('```\n')
    f.write('\n')
    f.write('## Getting Started\n')
    f.write('\n')
    f.write('To install `fuc`, enter the following in your terminal:\n')
    f.write('\n')
    f.write('```\n')
    f.write('$ git clone https://github.com/sbslee/fuc\n')
    f.write('$ cd fuc\n')
    f.write('$ pip install .\n')
    f.write('```\n')
    f.write('\n')
    f.write('For getting help on CLI:\n')
    f.write('\n')
    f.write('```\n')
    f.write('$ fuc -h\n')
    args = ['fuc', '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    f.write(result.stdout)
    f.write('```\n')
    f.write('\n')
    example = 'vfmerge'
    f.write(f'For getting help on a specific command (e.g. `{example}`):\n')
    f.write('\n')
    f.write('```\n')
    f.write(f'$ fuc {example} -h\n')
    args = ['fuc', example, '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    f.write(result.stdout)
    f.write('```\n')
    f.write('\n')
    f.write('Below is the list of modules available in API:\n')
    f.write('\n')
    for module in modules:
        description = pydoc.getdoc(getattr(fuc, module)).replace('\n', ' ')
        f.write(f'- **{module}** : {description}\n')
    f.write('\n')
    f.write('For getting help on a specific module (e.g. `VcfFrame`):\n')
    f.write('\n')
    f.write('```\n')
    f.write('from fuc import VcfFrame\n')
    f.write('help(VcfFrame)\n')
    f.write('```\n')
    f.write('\n')
    f.write('To give:\n')
    f.write('\n')
    f.write('```\n')
    result = pydoc.render_doc(VcfFrame, renderer=pydoc.plaintext)
    f.write(result)
    f.write('```\n')

# -- CLI.md ------------------------------------------------------------------

with open(cli_file, 'w') as f:
    f.write('# CLI\n')
    f.write('\n')
    f.write('## Table of contents\n')
    f.write('\n')
    f.write('* [Introduction](#Introduction)\n')
    f.write('* [Commands](#Commands)\n')
    for command in commands:
        f.write(f'\t* [{command}](#{command}) \n')
    f.write('\n')
    f.write('## Introduction <a name="Introduction"></a>\n')
    f.write('\n')
    f.write('For getting help on CLI:\n')
    f.write('\n')
    f.write('```\n')
    f.write('$ fuc -h\n')
    args = ['fuc', '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    f.write(result.stdout)
    f.write('```\n')
    f.write('\n')
    example = 'qfcount'
    f.write(f'For getting help on a specific command (e.g. `{example}`):\n')
    f.write('\n')
    f.write('```\n')
    f.write(f'$ fuc {example} -h\n')
    args = ['fuc', example, '-h']
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    f.write(result.stdout)
    f.write('```\n')
    f.write('\n')
    f.write('## Commands <a name="Commands"></a>\n')
    f.write('\n')
    for command in commands:
        args = ['fuc', command, '-h']
        result = subprocess.run(args, capture_output=True,
            text=True, check=True)
        f.write(f'### {command} <a name="{command}"></a>\n')
        f.write('\n')
        f.write('```\n')
        f.write(result.stdout)
        f.write('```\n')
        f.write('\n')

# -- API.md ------------------------------------------------------------------

with open(api_file, 'w') as f:
    f.write('# API\n')
    f.write('\n')
    f.write('## Table of contents\n')
    f.write('\n')
    for module in modules:
        f.write(f'* [{module}](#{module}) \n')
    f.write('\n')
    for module in modules:
        f.write(f'## {module} <a name="{module}"></a>\n')
        f.write('\n')
        f.write('```\n')
        doc = pydoc.render_doc(getattr(fuc, module), renderer=pydoc.plaintext)
        f.write(f'{doc}')
        f.write('```\n')
        f.write('\n')
