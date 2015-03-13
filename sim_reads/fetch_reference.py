#!/usr/bin/env python

"""
Fetch a reference sequence from the NCBI.


The e-mail parameter is used to connect to the NCBI.
"""

import argparse
import sys

from Bio import Entrez

from . import version


def fetch_ref(accno, output, email):
    """
    Fetch a reference sequence.

    :arg str accno: Accession number of a reference sequence.
    :arg stream output: Open handle to the output file.
    :arg str email: E-mail address used to connect to the NCBI.
    """
    Entrez.email = email

    handle = Entrez.efetch(db='nuccore', rettype='fasta', retmode='text',
        id=accno)

    output.write(handle.read())
    handle.close()


def main():
    """
    Main entry point.
    """
    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])

    parser.add_argument('input', metavar='INPUT', type=str,
        help='accession number of a reference sequence')
    parser.add_argument('output', metavar='OUTPUT',
        type=argparse.FileType('w'),
        help='name of the output file (default=<stdout>)')
    parser.add_argument('email', metavar='EMAIL', type=str,
        help='your e-mail address (%(type)s default=%(default)s)')
    parser.add_argument('-v', action='version', version=version(parser.prog))

    arguments = parser.parse_args()

    try:
        fetch_ref(arguments.input, arguments.output, arguments.email)
    except ValueError, error:
        parser.error(error)


if __name__ == '__main__':
    main()
