# README
This repository is my attempt to wrap some of the most frequently used commands in the field of bioinformatics with pure Python 3 code, only using standard libraries. That means no installation for external Python packages, not even the popular ones like `numpy` and `pandas`. This also includes many famous bioinformatics tools such as `samtools` and `bedtools`. The motivation for not relying on external packages or programs is quite simple: I just got tired of managing different Python environments for doing simple things.

You can use Python scripts in this repository for both command line interface (CLI) and application programming interface (API). Note that I'm not calling this repository as a Python package as it's merely a collection of Python scripts. In other words, there is no installation (i.e. no `setup.py`). Just clone this repository and you are good to go. This repository not being a package gives me the flexibility to be as creative as I want. However, this may change in the future.

Your contributions (e.g. feature ideas, pull requests) are most welcome.

Click [here](doc/CLI.md) to see the CLI documentation and [here](doc/API.md) to see the API documentation.

Author: Seung-been "Steven" Lee<br/>
Email: sbstevenlee@gmail.com<br/>
License: MIT License

# Usage

To clone this repository:

```
$ git clone https://github.com/sbslee/fuc
```

To run a command from the `fuc` repository (e.g. `check_files.py`):

```
$ python3 path/to/fuc/check_files.py input_file column_header
```

To access a module from the `fuc` repository (e.g. `VCFResult.py`):

```
import sys
sys.path.append('/path/to/fuc')
from VCFResult import VCFResult
```
