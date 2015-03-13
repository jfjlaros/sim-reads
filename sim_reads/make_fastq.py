#!/usr/local/bin/python

"""
Generate simulated reads.


See the help of the positional arguments for more info.
"""
# TODO Make the quality score a command line parameter.

import argparse
import random
import SOAPpy

from Bio import Seq, SeqIO

# Location of the Mutalyzer webservice.
mutalyzer_service_description = "https://mutalyzer.nl/services/?wsdl"
read_number = 1

def get_variant(handle):
    """
    Read one line from the input file, strip the newline and return it.

    @arg handle: Open handle to the input file.
    @type handle: file handle
    """
    return handle.readline().strip("\n")
#get_variant

def write_fastq(results, sequence, number_of_reads, insert_size, variation,
    read_length):
    """
    Make the simulated reads.
    """
    global read_number

    for i in range(number_of_reads):
        while True:
            position = random.randint(0, len(sequence) - insert_size)
            this_insert_size = int(random.normalvariate(insert_size,
                variation))
            if position + this_insert_size <= len(sequence):
                break
        #while

        results[0].write("@%i/1\n%s\n+\n%s\n" % (read_number, 
            sequence[position:position + read_length], "b" * read_length))
        results[1].write("@%i/2\n%s\n+\n%s\n" % (
            read_number, Seq.reverse_complement(
                sequence[position + this_insert_size - read_length:
                    position + this_insert_size]), "b" * read_length))
        read_number += 1
    #for
#write_fastq

def mutate(input_handle, output, insert, var, length, number, reference,
        start, end, accno, orientation, heterozygous):
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
    allelic_variant = get_variant(input_handle)
    line = get_variant(input_handle)
    while line:
        allelic_variant = "%s;%s" % (allelic_variant, line)
        line = get_variant(input_handle)
    #while
    
    # Set up the SOAP interface to Mutalyzer.
    mutalyzer_service = SOAPpy.WSDL.Proxy(mutalyzer_service_description)

    # Retrieve the reference sequence.
    if not accno:
        accno = mutalyzer_service.sliceChromosome(chromAccNo=reference,
            start=start, end=end, orientation=orientation)

    # Mutate the reference sequence.
    mutalyzer_output = mutalyzer_service.runMutalyzer(
        variant = "%s:g.[%s]" % (accno, allelic_variant))
    sequence = mutalyzer_output.mutated

    # Write the chromosomal description to the results file.
    results_handle = open("%s.txt" % output, "w")
    results_handle.write("%s: %s:%i_%i\n\n" % (accno, reference, start, end))
    
    if int(mutalyzer_output.errors) > 0:
        for i in mutalyzer_output.messages[0]:
            raise ValueError("(%s): %s" % (i.errorcode, i.message))
    #if
    if "chromDescription" in mutalyzer_output:
        description = mutalyzer_output.chromDescription
    else:
        description = mutalyzer_output.genomicDescription

    for i in description.split("[")[1][:-1].split(';'):
        results_handle.write("%s\n" % i)
    results_handle.close()

    #outputHandle1 = open("%s_1.fq" % results, "w")
    #outputHandle2 = open("%s_2.fq" % results, "w")
    results = [open("%s_1.fq" % output, "w"),
        open("%s_2.fq" % output, "w")] 
    if heterozygous:
        write_fastq(results, sequence, number / 2, insert, var, length)
        write_fastq(results, mutalyzer_output.original, number / 2, insert,
            var, length)
    #if
    else:
        write_fastq(results, sequence, number, insert, var, length)
#mutate

def local(output, reference_handle, insert, var, length, number):
    """
    Use a local fasta file to make paired end reads.

    @arg args: Argparse argument list.
    @type args: object
    """
    results = [open("%s_1.fq" % output, "w"),
        open("%s_2.fq" % output, "w")] 

    for record in SeqIO.parse(reference_handle, "fasta"):
        write_fastq(results, str(record.seq), number, insert, var, length)
#local

def main():
    """
    Main entry point.
    """
    parent_parser = argparse.ArgumentParser('parent', add_help=False)
    parent_parser.add_argument("output", metavar="OUTPUT", type=str,
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
    parser_mutate.add_argument("input_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="name of the input file")
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
    parser_mutate.add_argument("-d", dest="heterozygous", default=False,
        action="store_true", help="make heterozygous variants")
    parser_mutate.set_defaults(func=mutate)

    parser_local = subparsers.add_parser("local", parents=[parent_parser],
        description=local.__doc__.split("\n\n")[0])
    parser_local.add_argument("reference_handle", metavar="REFERENCE",
        type=argparse.FileType('r'), help="name of a local reference sequence")
    parser_local.set_defaults(func=local)

    arguments = parser.parse_args()

    try:
        arguments.func(arguments)
    except ValueError, error:
        parser.error(error)

    try:
        arguments.func(**dict((k, v) for k, v in vars(arguments).items()
            if k not in ('func', 'subcommand')))
    except ValueError as error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
