import subprocess
import pydoc

from fuc.api.common import fuc_dir
from fuc.cli import commands
from fuc import VcfFrame
import fuc

modules = [x for x in dir(fuc) if x not in ['api', 'cli'] and '__' not in x]

# -- README.rst ---------------------------------------------------------------

readme_file = f'{fuc_dir()}/README.rst'

fuc_help = subprocess.run(['fuc', '-h'], capture_output=True, text=True, check=True).stdout
fuc_help = '\n'.join(['   ' + x for x in fuc_help.splitlines()])

vfmerge_help = subprocess.run(['fuc', 'vfmerge', '-h'], capture_output=True, text=True, check=True).stdout
vfmerge_help = '\n'.join(['   ' + x for x in fuc_help.splitlines()])

module_help = ''
for module in modules:
    description = pydoc.getdoc(getattr(fuc, module)).replace('\n', ' ')
    module_help += f'- **{module}** : {description}\n'

vcfframe_help = pydoc.render_doc(VcfFrame, renderer=pydoc.plaintext)
vcfframe_help = '\n'.join(['   ' + x for x in vcfframe_help.splitlines()])

d = dict(fuc_help=fuc_help, vfmerge_help=vfmerge_help, module_help=module_help, vcfframe_help=vcfframe_help)

readme = """
README
******

.. image:: https://readthedocs.org/projects/sbslee-fuc/badge/?version=latest
   :target: https://sbslee-fuc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Introduction
============

The main goal of the ``fuc`` package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use ``fuc`` for both command line interface (CLI) and application programming interface (API) whose documentations are available at `Read the Docs <https://sbslee-fuc.readthedocs.io/en/latest/>`_.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

| Author: Seung-been "Steven" Lee
| Email: sbstevenlee@gmail.com
| License: MIT License

Examples
========

To merge VCF files with CLI:

.. code-block:: console

   $ fuc vfmerge 1.vcf 2.vcf 3.vcf > merged.vcf

To filter a VCF file based on a BED file using API:

.. code:: python3

   from fuc.api.VcfFrame import VcfFrame
   vf = VcfFrame.from_file('original.vcf')
   filtered_vf = vf.filter_bed('targets.bed')
   filtered_vf.to_file('filtered.vcf')

Required Packages
=================

The following packages are required to run ``fuc``:

.. parsed-literal::

   numpy
   pandas
   pyranges

Getting Started
===============

To install ``fuc``, enter the following in your terminal:

.. code-block:: console

   $ git clone https://github.com/sbslee/fuc
   $ cd fuc
   $ pip install .

For getting help on CLI:

.. code-block:: console

   $ fuc -h
{fuc_help}

For getting help on a specific command (e.g. `vfmerge`):

.. code-block:: console

   $ fuc vfmerge -h
{vfmerge_help}

Below is the list of modules available in API:

{module_help}
For getting help on a specific module (e.g. `VcfFrame`):

.. code:: python3

   from fuc.api import VcfFrame
   help(VcfFrame)

To give:

.. parsed-literal::

{vcfframe_help}
""".format(**d)

with open(readme_file, 'w') as f:
    f.write(readme.lstrip())

# -- cli.rst -----------------------------------------------------------------

cli_file = f'{fuc_dir()}/docs/cli.rst'

cli = """
CLI
***

Introduction
============

This section describes command line interface (CLI) for the ``fuc`` package.

For getting help on CLI:

.. code-block:: console

   $ fuc -h
{fuc_help}

For getting help on a specific command (e.g. `vfmerge`):

.. code-block:: console

   $ fuc vfmerge -h
{vfmerge_help}

""".format(**d)

for command in commands:
    s = f'{command}\n'
    s += '=' * (len(s)-1) + '\n'
    s += '\n'
    s += '.. code-block:: console\n'
    s += '\n'
    s += f'   $ fuc {command} -h\n'
    command_help = subprocess.run(['fuc', command, '-h'], capture_output=True, text=True, check=True).stdout
    command_help = '\n'.join(['   ' + x for x in command_help.splitlines()])
    s += command_help + '\n'
    s += '\n'
    cli += s

with open(cli_file, 'w') as f:
    f.write(cli.lstrip())

# -- api.rst -----------------------------------------------------------------

api_file = f'{fuc_dir()}/docs/api.rst'

api = """
API
***

Introduction
============

This section describes application programming interface (API) for the ``fuc`` package.

Below is the list of modules available in API:

{module_help}
For getting help on a specific module (e.g. `VcfFrame`):

.. code:: python3

   from fuc.api import VcfFrame
   help(VcfFrame)

To give:

.. parsed-literal::

{vcfframe_help}

""".format(**d)

for module in modules:
    s = f'fuc.api.{module}\n'
    s += '=' * (len(s)-1) + '\n'
    s += '\n'
    s += f'.. automodule:: fuc.api.{module}\n'
    s += '   :members:\n'
    s += '\n'
    api += s

with open(api_file, 'w') as f:
    f.write(api.lstrip())
