from setuptools import setup, find_packages

setup(
    name='gfetchseq',
    version='0.1',
    description='Fetch DNA gh38 sequences into a formatted fasta file from a gatk_interval file',
    url='http://github.com/lexxxxxxa/gfetchseq',
    author='Lexa',
    author_email='alexper.recovery@gmail.com',
    license='MUNI',
    install_requires=['bioblend', 'json', 'csv', 'time', 'datetime', 'logging',
    'pandas', 'termcolor', 'Bio'],
    packages=find_packages(),
    entry_points=dict(
        console_scripts=['rq=src.main']
    )
)