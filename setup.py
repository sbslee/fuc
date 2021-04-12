from setuptools import setup, find_packages

exec(open('fuc/version.py').read())

setup(
    name='fuc',
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email='sbstevenlee@gmail.com',
    description=('Frequently used commands in bioinformatics with '
                 'pure Python 3 code'),
    url='https://github.com/sbslee/fuc',
    packages=find_packages(),
    entry_points={'console_scripts': ['fuc=fuc.__main__:main']}
)
