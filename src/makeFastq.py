#!/usr/local/bin/python

"""
Generate simulated reads.


See the help of the positional arguments for more info.
"""
# TODO Make the quality score a command line parameter.

import SOAPpy
import random
import argparse

from Bio import Seq
from Bio import SeqIO

# Location of the Mutalyzer webservice.
mutalyzerServiceDescription = "https://mutalyzer.nl/services/?wsdl"

def getVariant(handle):
    """
    Read one line from the input file, strip the newline and return it.

    @arg handle: Open handle to the input file.
    @type handle: file handle
    """
    return handle.readline().strip("\n")
#getVariant

def writeFastq(results, sequence, numberOfReads, insertSize, variation,
    readLength):
    """
    Make the simulated reads.
    """
    outputHandle1 = open("%s_1.fq" % results, "w")
    outputHandle2 = open("%s_2.fq" % results, "w")
    readNumber = 1
    for i in range(numberOfReads):
        while True:
            position = random.randint(0, len(sequence) - insertSize)
            myInsertSize = int(random.normalvariate(insertSize, variation))
            if position + myInsertSize <= len(sequence):
                break
        #while

        outputHandle1.write("@%i/1\n%s\n+\n%s\n" % (readNumber, 
            sequence[position:position + readLength], "b" * readLength))
        outputHandle2.write("@%i/2\n%s\n+\n%s\n" % (
            readNumber, Seq.reverse_complement(
                sequence[position + myInsertSize - readLength:
                    position + myInsertSize]), "b" * readLength))
        readNumber += 1
    #for
    outputHandle2.close()
    outputHandle1.close()
#writeFastq

def mutate(args):
    """
    Use a file containing variants to obtain a mutated reference sequence from
    Mutalyzer. Then make paired end reads out of the mutated sequence.


    If a chromosomal accession number (option -r) is used, the options -b and
      -e are used to make a slice of this chromosome.
    To retrieve the reference sequence of a gene, use the -u option. The
      options -r, -b and -e are ignored.
    With the -o OUTPUT option the prefix for three output files is given:
      OUTPUT.txt, OUTPUT_1.fq and OUTPUT_2.fq


    @arg args: Argparse argument list.
    @type args: object
    """
    # Read the input file.
    allelicVariant = getVariant(args.input)
    line = getVariant(args.input)
    while line:
        allelicVariant = "%s;%s" % (allelicVariant, line)
        line = getVariant(args.input)
    #while
    
    # Set up the SOAP interface to Mutalyzer.
    mutalyzerService = SOAPpy.WSDL.Proxy(mutalyzerServiceDescription)

    # Retrieve the reference sequence.
    if not args.accno:
        args.accno = mutalyzerService.sliceChromosome(
            chromAccNo=args.reference, start=args.start, end=args.end,
            orientation=args.orientation)

    # Mutate the reference sequence.
    mutalyzerOutput = mutalyzerService.runMutalyzer(
        variant = "%s:g.[%s]" % (args.accno, allelicVariant))
    sequence = mutalyzerOutput.mutated

    # Write the chromosomal description to the results file.
    resultsHandle = open("%s.txt" % args.output, "w")
    resultsHandle.write("%s: %s:%i_%i\n\n" % (args.accno, args.reference,
        args.start, args.end))
    
    if int(mutalyzerOutput.errors) > 0:
        for i in mutalyzerOutput.messages[0]:
            raise ValueError("(%s): %s" % (i.errorcode, i.message))
    #if
    if "chromDescription" in mutalyzerOutput:
        description = mutalyzerOutput.chromDescription
    else:
        description = mutalyzerOutput.genomicDescription

    for i in description.split("[")[1][:-1].split(';'):
        resultsHandle.write("%s\n" % i)
    resultsHandle.close()

    writeFastq(args.output, sequence, args.number, args.insert, args.var,
        args.length)
#mutate

def local(args):
    """
    Use a local fasta file to make paired end reads.

    @arg args: Argparse argument list.
    @type args: object
    """
    for record in SeqIO.parse(args.refFile, "fasta"):
        writeFastq(args.output, str(record.seq), args.number, args.insert,
            args.var,args.length)
#local

def main():
    """
    Main entry point.
    """

    parent_parser = argparse.ArgumentParser('parent', add_help=False)
    parent_parser.add_argument("-o", dest="output", type=str, required=True,
        help="prefix of the names of the output files")
    parent_parser.add_argument("-s", dest="insert", type=int, default=300,
        help="mean insert size (default=%(default)s)")
    parent_parser.add_argument("-v", dest="var", type=int, default=25,
        help="standard deviation of the insert size (default=%(default)s)")
    parent_parser.add_argument("-l", dest="length", type=int, default=50,
        help="read length (default=%(default)s)")
    parent_parser.add_argument("-n", dest="number", type=int, default=1000000,
        help="number of reads (default=%(default)s)")

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    subparsers = parser.add_subparsers()

    usage = mutate.__doc__.split("\n\n\n")
    parser_mutate = subparsers.add_parser("mutate", parents=[parent_parser],
        description=usage[0], epilog=usage[1])
    parser_mutate.add_argument("-r", dest="reference", default="NC_000008.10",
        type=str, help="chromosomal accession number (default=%(default)s)")
    parser_mutate.add_argument("-b", dest="start", type=int, default=136800000,
         help="begin of the slice (default=%(default)s)")
    parser_mutate.add_argument("-e", dest="end", type=int, default=139000000,
        help="end of the slice (default=%(default)s)")
    parser_mutate.add_argument("-u", dest="accno", type=str, default="",
         help="accession number of a referencence sequence")
    parser_mutate.add_argument("-t", dest="orientation", default=1, const=2,
        action="store_const", help="reverse the orientation of the slice")
    parser_mutate.add_argument("-i", dest="input", type=argparse.FileType('r'),
        required=True, help="name of the input file")
    parser_mutate.set_defaults(func=mutate)

    parser_local = subparsers.add_parser("local", parents=[parent_parser],
        description=local.__doc__.split("\n\n")[0])
    parser_local.add_argument("-r", dest="refFile", required=True,
        type=argparse.FileType('r'), help="name of a local reference sequence")
    parser_local.set_defaults(func=local)

    arguments = parser.parse_args()

    try:
        arguments.func(arguments)
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
