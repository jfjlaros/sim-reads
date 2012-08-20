#!/usr/bin/env python

"""
Fetch a reference sequence from the NCBI.


The e-mail parameter is used to connect to the NCBI.
"""

from Bio import Entrez
import argparse
import sys

def fetchRef(accno, output, email):
    """
    Fetch a reference sequence.

    @arg accno: Accession number of a reference sequence.
    @type accno: str
    @arg outputName: Open handle to the output file.
    @type outputName: stream
    @arg email: E-mail address used to connect to the NCBI.
    @type email: str
    """
    Entrez.email = email

    handle = Entrez.efetch(db="nuccore", rettype="fasta", retmode="text",
        id=accno)

    output.write(handle.read())
    handle.close()
#fetchRef

def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])

    parser.add_argument(dest="input", metavar="INPUT", type=str,
        help="accession number of a reference sequence")
    parser.add_argument("-o", dest="output", type=argparse.FileType('w'),
        default=sys.stdout, help="name of the output file (default=<stdout>)")
    parser.add_argument("-e", dest="email", required=True, type=str,
        help="your e-mail address (%(type)s default=%(default)s)")

    arguments = parser.parse_args()

    try:
        fetchRef(arguments.input, arguments.output, arguments.email)
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
