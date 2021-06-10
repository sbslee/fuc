import subprocess
import pydoc

from fuc.api.common import fuc_dir
from fuc.cli import commands
from fuc import pyvcf
import fuc

modules = [x for x in dir(fuc) if x not in ['api', 'cli'] and '__' not in x]

# -- README.rst ---------------------------------------------------------------

credit = """
..
   This file was automatically generated by docs/create.py.
"""

readme_file = f'{fuc_dir()}/README.rst'

fuc_help = subprocess.run(['fuc', '-h'], capture_output=True, text=True, check=True).stdout
fuc_help = '\n'.join(['   ' + x for x in fuc_help.splitlines()])

module_help = ''
for module in modules:
    description = pydoc.getdoc(getattr(fuc, module)).split('\n\n')[0].replace('\n', ' ')
    module_help += f'- **{module}** : {description}\n'

d = dict(credit=credit, fuc_help=fuc_help, module_help=module_help)

readme = """
{credit}
README
******

.. image:: https://badge.fury.io/py/fuc.svg
    :target: https://badge.fury.io/py/fuc

.. image:: https://readthedocs.org/projects/sbslee-fuc/badge/?version=latest
   :target: https://sbslee-fuc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://anaconda.org/bioconda/fuc/badges/version.svg
   :target: https://anaconda.org/bioconda/fuc

.. image:: https://anaconda.org/bioconda/fuc/badges/license.svg
   :target: https://github.com/sbslee/fuc/blob/main/LICENSE

.. image:: https://anaconda.org/bioconda/fuc/badges/downloads.svg
   :target: https://anaconda.org/bioconda/fuc/files

.. image:: https://anaconda.org/bioconda/fuc/badges/installer/conda.svg
   :target: https://conda.anaconda.org/bioconda

Introduction
============

The main goal of the fuc package is to wrap some of the most frequently used commands in the field of bioinformatics into one place.

You can use fuc for both command line interface (CLI) and application programming interface (API) whose documentations are available at `Read the Docs <https://sbslee-fuc.readthedocs.io/en/latest/>`_.

Currently, the following file formats are supported by fuc:

- Sequence Alignment/Map (SAM)
- Binary Alignment/Map (BAM)
- CRAM
- Variant Call Format (VCF)
- Mutation Annotation Format (MAF)
- Browser Extensible Data (BED)
- FASTQ
- delimiter-separated values format (e.g. comma-separated values or CSV format)

Additionally, you can use fuc to parse output data from the following programs:

- Ensembl Variant Effect Predictor (VEP)
- SnpEff
- bcl2fastq and bcl2fastq2

Your contributions (e.g. feature ideas, pull requests) are most welcome.

| Author: Seung-been "Steven" Lee
| Email: sbstevenlee@gmail.com
| License: MIT License

CLI Examples
============

SAM/BAM/CRAM
------------

To print the header of a BAM file:

.. code-block:: console

   $ fuc bam_head example.bam

BED
---

To find intersection between BED files:

.. code-block:: console

   $ fuc bed_intxn 1.bed 2.bed 3.bed > intersect.bed

FASTQ
-----

To count sequence reads in a FASTQ file:

.. code-block:: console

   $ fuc fq_count example.fastq

FUC
---

To check whether a file exists in the operating system:

.. code-block:: console

   $ fuc fuc_exist example.txt

To find all VCF files within the current directory recursively:

.. code-block:: console

   $ fuc fuc_find . vcf

TABLE
-----

To merge two tab-delimited files:

.. code-block:: console

   $ fuc tbl_merge left.txt right.txt > merged.txt

VCF
---

To merge VCF files:

.. code-block:: console

   $ fuc vcf_merge 1.vcf 2.vcf 3.vcf > merged.vcf

API Examples
============

VCF
---

To filter a VCF file based on a BED file:

.. code:: python3

   >>> from fuc import pyvcf
   >>> vf = pyvcf.VcfFrame.from_file('original.vcf')
   >>> filtered_vf = vf.filter_bed('targets.bed')
   >>> filtered_vf.to_file('filtered.vcf')

To remove indels from a VCF file:

.. code:: python3

   >>> from fuc import pyvcf
   >>> vf = pyvcf.VcfFrame.from_file('with_indels.vcf')
   >>> filtered_vf = vf.filter_indel()
   >>> filtered_vf.to_file('no_indels.vcf')

To create a Venn diagram showing genotype concordance between groups:

.. code:: python3

    >>> from fuc import pyvcf, common
    >>> common.load_dataset('pyvcf')
    >>> f = '~/fuc-data/pyvcf/plot_comparison.vcf'
    >>> vf = pyvcf.VcfFrame.from_file(f)
    >>> a = ['Steven_A', 'John_A', 'Sara_A']
    >>> b = ['Steven_B', 'John_B', 'Sara_B']
    >>> c = ['Steven_C', 'John_C', 'Sara_C']
    >>> vf.plot_comparison(a, b, c)

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/plot_comparison.png

To create a histogram of tumor mutational burden (TMB) distribution:

.. code:: python3

    >>> from fuc import pyvcf
    >>> vcf_data = {{
    ...     'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
    ...     'POS': [100, 101, 102, 103, 103],
    ...     'ID': ['.', '.', '.', '.', '.'],
    ...     'REF': ['T', 'T', 'T', 'T', 'T'],
    ...     'ALT': ['C', 'C', 'C', 'C', 'C'],
    ...     'QUAL': ['.', '.', '.', '.', '.'],
    ...     'FILTER': ['.', '.', '.', '.', '.'],
    ...     'INFO': ['.', '.', '.', '.', '.'],
    ...     'FORMAT': ['GT', 'GT', 'GT', 'GT', 'GT'],
    ...     'Steven_N': ['0/0', '0/0', '0/1', '0/0', '0/0'],
    ...     'Steven_T': ['0/0', '0/1', '0/1', '0/1', '0/1'],
    ...     'Sara_N': ['0/0', '0/1', '0/0', '0/0', '0/0'],
    ...     'Sara_T': ['0/0', '0/0', '1/1', '1/1', '0/1'],
    ...     'John_N': ['0/0', '0/0', '0/0', '0/0', '0/0'],
    ...     'John_T': ['0/1', '0/0', '1/1', '1/1', '0/1'],
    ...     'Rachel_N': ['0/0', '0/0', '0/0', '0/0', '0/0'],
    ...     'Rachel_T': ['0/1', '0/1', '0/0', '0/1', '0/1'],
    ... }}
    >>> annot_data = {{
    ...     'Sample': ['Steven_N', 'Steven_T', 'Sara_N', 'Sara_T', 'John_N', 'John_T', 'Rachel_N', 'Rachel_T'],
    ...     'Subject': ['Steven', 'Steven', 'Sara', 'Sara', 'John', 'John', 'Rachel', 'Rachel'],
    ...     'Type': ['Normal', 'Tumor', 'Normal', 'Tumor', 'Normal', 'Tumor', 'Normal', 'Tumor'],
    ... }}
    >>> vf = pyvcf.VcfFrame.from_dict([], vcf_data)
    >>> af = pyvcf.AnnFrame.from_dict(annot_data, 'Sample')
    >>> vf.plot_histplot(hue='Type', af=af)

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/plot_histplot.png

MAF
---

To create an oncoplot with a MAF file:

.. code:: python3

    >>> from fuc import common, pymaf
    >>> common.load_dataset('tcga-laml')
    >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
    >>> mf = pymaf.MafFrame.from_file(f)
    >>> mf.plot_oncoplot()

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/oncoplot.png

To create a customized oncoplot with a MAF file, see the 'Create customized oncoplot' tutorial:

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/customized_oncoplot.png

To create a summary figure for a MAF file:

.. code:: python3

    >>> from fuc import common, pymaf
    >>> common.load_dataset('tcga-laml')
    >>> f = '~/fuc-data/tcga-laml/tcga_laml.maf.gz'
    >>> mf = pymaf.MafFrame.from_file(f)
    >>> mf.plot_summary()

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/maf_summary.png

SAM/BAM/CRAM
------------

To create read depth profile of a region from a CRAM file:

.. code:: python3

    >>> from fuc import pycov
    >>> cf = pycov.CovFrame.from_file('HG00525.final.cram', zero=True,
    ...    region='chr12:21161194-21239796', names=['HG00525'])
    >>> cf.plot_region('chr12', start=21161194, end=21239796)

.. image:: https://raw.githubusercontent.com/sbslee/fuc-data/main/images/coverage.png

Installation
============

The following packages are required to run fuc:

.. parsed-literal::

   biopython
   lxml
   matplotlib
   numpy
   pandas
   pyranges
   pysam
   seaborn

There are various ways you can install fuc. The recommended way is via conda:

.. code-block:: console

   $ conda install -c bioconda fuc

Above will automatically download and install all the dependencies as well. Alternatively, you can use pip to install fuc and all of its dependencies:

.. code-block:: console

   $ pip install fuc

Finally, you can clone the GitHub repository and then install fuc this way:

.. code-block:: console

   $ git clone https://github.com/sbslee/fuc
   $ cd fuc
   $ pip install .

The nice thing about this approach is that you will have access to development versions that are not available in Anaconda or PyPI. For example, you can access a development branch with the ``git checkout`` command. When you do this, please make sure your environment already has all the dependencies installed.

Getting Help
============

For detailed documentations on fuc's CLI and API, please refer to the `Read the Docs <https://sbslee-fuc.readthedocs.io/en/latest/>`_.

For getting help on CLI:

.. code-block:: console

   $ fuc -h
{fuc_help}

For getting help on a specific command (e.g. vcf_merge):

.. code-block:: console

   $ fuc vcf_merge -h

Below is the list of submodules available in API:

{module_help}
For getting help on a specific module (e.g. pyvcf):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

""".format(**d)

with open(readme_file, 'w') as f:
    f.write(readme.lstrip())

# -- cli.rst -----------------------------------------------------------------

cli_file = f'{fuc_dir()}/docs/cli.rst'

cli = """
{credit}
CLI
***

Introduction
============

This section describes command line interface (CLI) for the fuc package.

For getting help on CLI:

.. code-block:: console

   $ fuc -h
{fuc_help}

For getting help on a specific command (e.g. vcf_merge):

.. code-block:: console

   $ fuc vcf_merge -h

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
{credit}
API
***

Introduction
============

This section describes application programming interface (API) for the fuc package.

Below is the list of submodules available in API:

{module_help}
For getting help on a specific module (e.g. pyvcf):

.. code:: python3

   from fuc import pyvcf
   help(pyvcf)

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
