#!/usr/local/bin/python

"""
Generate simulated reads from a reference sequence and a file with variants in 
HGVS notation.

Options
=======
    -r NAME  --reference=NAME
    Name of the reference sequence.

    -i NAME  --input=NAME
    Name of the input file.

    -o NAME  --output=NAME
    Name of the output files, will be appended with "_1.fq" and "_2.fq".

    -s NUM  --insert=NUM
    Mean insert size.

    -v NUM  --var=NUM
    Standard deviation of the insert size.

    -l NUM  --length=NUM
    Read length.

    -n NUM  --number=NUM
    Number of mates per left read.

@requires: SOAPpy
@requires: random
@requires: Bio.Seq
@requires: getopt
@requires: sys
"""

import SOAPpy
import random
import Bio.Seq
import getopt
import sys

# Location of the Mutalyzer webservice.
mutalyzerServiceDescription = "http://www.mutalyzer.nl/2.0/service.wsdl"

def getVariant(handle) :
    """
    Read one line from the input file, strip the newline and return it.

    @arg handle: Open handle to the input file.
    @type handle: file handle
    """

    return handle.readline().strip("\n")
#getVariant

def main() :
    """
    - Read the input file and make one big variant of the description in the
      input file.
    - Feed the variant to the Mutalyzer webservice to obtain a mutated
      reference sequence.
    - Make paired end reads out of the mutated sequence.
    """

    # Some default parameters.
    referenceSequence = "UD_129077402355"
    inputFile = "input.txt"
    outputFile1 = None
    outputFile2 = None
    insertSize = 300
    variation = 25
    readLength = 50
    numberOfMates = 5

    # Argument parsing.
    try :
        opts, args = getopt.getopt(sys.argv[1:], "r:i:o:s:v:l:n", [
            "reference=", "input=", "output=", "insert=", "var=", "length=",
            "number="])
    except getopt.GetoptError, err :
        print str(err)
        sys.exit(2)
    #except
    for option, argument in opts :
        if option in ("-r", "--reference") :
            referenceSequence = argument
        elif option in ("-i", "--input") :
            inputFile = argument
        elif option in ("-o", "--output") :
            outputFile1 = "%s_1.fq" % argument
            outputFile2 = "%s_2.fq" % argument
        #elif
        elif option in ("-s", "--insert") :
            insertSize = argument
        elif option in ("-v", "--var") :
            variation = argument
        elif option in ("-l", "--length") :
            readLength = argument
        elif option in ("-n", "--number") :
            numberOfMates = argument
        else :
            assert False, "Unhandled option."
    #for
    if not outputFile1 :
        print "No output given, use -o NAME"
        sys.exit(2)
    #if

    # Read the input file.
    inputHandle = open(inputFile, "r")
    allelicVariant = getVariant(inputHandle)
    line = getVariant(inputHandle)
    while line :
        allelicVariant = "%s;%s" % (allelicVariant, line)
        line = getVariant(inputHandle)
    #while
    inputHandle.close()
    
    # Make a mutated sequence.
    mutalyzerService = SOAPpy.WSDL.Proxy(mutalyzerServiceDescription)
    mutalyzerOutput = mutalyzerService.runMutalyzer(
        variant = "%s:g.[%s]" % (referenceSequence, allelicVariant))
    sequence = mutalyzerOutput.mutated

    # Make the simulated reads.
    outputHandle1 = open(outputFile1, "w")
    outputHandle2 = open(outputFile2, "w")
    readNumber = 1
    for i in range(len(sequence) - readLength - insertSize) :
        for j in range(numberOfMates) :
            myInsertSize = int(random.normalvariate(insertSize, variation))
            outputHandle1.write("@%i_forward/1\n%s\n+\n%s\n" % (
                readNumber, sequence[i:i + readLength], "I" * readLength))
            outputHandle2.write("@%i_reverse/2\n%s\n+\n%s\n" % (
                readNumber, Bio.Seq.reverse_complement(
                    sequence[i + myInsertSize - readLength:i + myInsertSize]),
                "I" * readLength))
            readNumber += 1
        #for
    #for
    outputHandle2.close()
    outputHandle1.close()
#main

if __name__ == "__main__" :
    main()
