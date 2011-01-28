#!/usr/local/bin/python

"""
Generate simulated reads from a reference sequence and a file with variants in 
HGVS notation.

Options
=======
    -r NAME  --reference=NAME
    Name of the reference sequence.

    -b NUM  --start=NUM
    Begin of the slice of the reference sequence.

    -e NUM  --end=NUM
    End of the slice of the reference sequence.

    -t NUM  --orientation=NUM
    Orientation of the slice: 1=forward, 2=reverse.

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
    Number of reads.

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
    referenceSequence = "NC_000008.10"
    referenceStart = 136800000
    referenceEnd = 139000000
    referenceOrientation = 1 # 1=forward, 2=reverse

    inputFile = "input.txt"
    outputFile1 = None
    outputFile2 = None
    resultsFile = None
    insertSize = 300
    variation = 25
    readLength = 50
    numberOfReads = 1000000

    # Argument parsing.
    try :
        opts, args = getopt.getopt(sys.argv[1:], "r:b:e:t:i:o:s:v:l:n:", [
            "reference=", "start", "end", "orientation", "input=", "output=",
            "insert=", "var=", "length=", "number="])
    except getopt.GetoptError, err :
        print str(err)
        sys.exit(2)
    #except
    for option, argument in opts :
        if option in ("-r", "--reference") :
            referenceSequence = argument
        elif option in ("-b", "--input") :
            referenceStart = int(argument)
        elif option in ("-e", "--input") :
            referenceEnd = int(argument)
        elif option in ("-t", "--input") :
            referenceOrientation = argument
        elif option in ("-i", "--input") :
            inputFile = argument
        elif option in ("-o", "--output") :
            outputFile1 = "%s_1.fq" % argument
            outputFile2 = "%s_2.fq" % argument
            resultsFile = "%s.txt" % argument
        #elif
        elif option in ("-s", "--insert") :
            insertSize = int(argument)
        elif option in ("-v", "--var") :
            variation = int(argument)
        elif option in ("-l", "--length") :
            readLength = int(argument)
        elif option in ("-n", "--number") :
            numberOfReads = int(argument)
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
    
    # Set up the SOAP interface to Mutalyzer.
    mutalyzerService = SOAPpy.WSDL.Proxy(mutalyzerServiceDescription)

    # Retrieve the reference sequence.
    UD = mutalyzerService.sliceChromosome(
        chromAccNo = referenceSequence, start = referenceStart, 
        end = referenceEnd, orientation = referenceOrientation)

    # Mutate the reference sequence.
    mutalyzerOutput = mutalyzerService.runMutalyzer(
        variant = "%s:g.[%s]" % (UD, allelicVariant))
    sequence = mutalyzerOutput.mutated

    # Write the chromosomal description to the results file.
    resultsHandle = open(resultsFile, "w")
    resultsHandle.write("%s: %s:%i_%i\n\n" % (
        UD, referenceSequence, referenceStart, referenceEnd))
    for i in mutalyzerOutput.chromDescription.split("[")[1][:-1].split(';') :
        resultsHandle.write("%s\n" % i)
    #for
    resultsHandle.close()

    # Make the simulated reads.
    outputHandle1 = open(outputFile1, "w")
    outputHandle2 = open(outputFile2, "w")
    readNumber = 1
    for i in range(numberOfReads) :
        position = random.randint(0, len(sequence) - insertSize)
        myInsertSize = int(random.normalvariate(insertSize, variation))
        outputHandle1.write("@%i_forward/1\n%s\n+\n%s\n" % (readNumber, 
            sequence[position:position + readLength], "I" * readLength))
        outputHandle2.write("@%i_reverse/2\n%s\n+\n%s\n" % (
            readNumber, Bio.Seq.reverse_complement(
                sequence[position + myInsertSize - readLength:
                    position + myInsertSize]), "I" * readLength))
        readNumber += 1
    #for
    outputHandle2.close()
    outputHandle1.close()
#main

if __name__ == "__main__" :
    main()
