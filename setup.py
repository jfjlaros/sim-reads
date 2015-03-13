import sys
from setuptools import setup

if sys.version_info < (2, 6):
    raise Exception('fastools requires Python 2.6 or higher.')

# Todo: How does this play with pip freeze requirement files?
requires = ['SOAPpy']

# Python 2.6 does not include the argparse module.
try:
    import argparse
except ImportError:
    requires.append('argparse')

import fastools as distmeta

setup(
    name='sim-reads',
    version=distmeta.__version__,
    description='Generate simulated reads.',
    long_description=distmeta.__doc__,
    author=distmeta.__author__,
    author_email=distmeta.__contact__,
    url=distmeta.__homepage__,
    license='MIT License',
    platforms=['any'],
    packages=['sim-reads'],
    install_requires=requires,
    entry_points = {
        'console_scripts': [
            'sim_reads = sim_reads.make_fastq:main',
            'fetch_reference = sim_reads.fetch_reference:main'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
    ],
    keywords='bioinformatics'
)
