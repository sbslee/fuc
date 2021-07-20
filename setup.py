from setuptools import setup, find_packages

exec(open('fuc/version.py').read())

requirements = [
    'biopython', 'lxml', 'matplotlib', 'matplotlib-venn', 'numpy', 'pandas',
    'pyranges', 'pysam', 'scipy', 'seaborn', 'statsmodels'
]

setup(
    name='fuc',
    version=__version__,
    author='Seung-been "Steven" Lee',
    author_email='sbstevenlee@gmail.com',
    description=('Frequently used commands in bioinformatics'),
    long_description=open('README.rst').read(),
    long_description_content_type='text/x-rst',
    url='https://github.com/sbslee/fuc',
    packages=find_packages(),
    license='MIT',
    install_requires=requirements,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    entry_points={'console_scripts': ['fuc=fuc.__main__:main']}
)
