#!/usr/bin/env python

from collections import OrderedDict, Counter
import sys

infiles = sys.argv[1:]


def convert_library_type(library_type, prog, paired):
    """
    Convert TopHat library type to other program naming convention.
    """
    new = None

    if prog == "TopHat2":
        new == library_type


    elif prog == "RSeQC":
        if library_type == 'fr-firststrand' and paired==True:
            new = '"1+-,1-+,2++,2--"'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = '"1++,1--,2+-,2-+"'
        elif library_type == 'fr-firststrand' and paired==False:
            new = '"+-,-+"'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = '"++,--"'
        elif library_type == 'fr-unstranded':
            new = ''
        else:
            new = ''


    elif prog == "HISAT2":
        if library_type == 'fr-firststrand' and paired==True:
            new = 'RF'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = 'FR'
        elif library_type == 'fr-firststrand' and paired==False:
            new = 'R'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = 'F'
        elif library_type == 'fr-unstranded':
            new = 'unstranded'


    elif prog == "htseq-count":
        if library_type == 'fr-firststrand' and paired==True:
            new = 'reverse'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = 'yes'
        elif library_type == 'fr-firststrand' and paired==False:
            new = 'reverse'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = 'yes'
        elif library_type == 'fr-unstranded':
            new = 'no'


    elif prog == "featureCounts":
        if library_type == 'fr-firststrand' and paired==True:
            new = '2'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = '1'
        elif library_type == 'fr-firststrand' and paired==False:
            new = '2'
        elif library_type == 'fr-secondstrand' and paired==False:
            new = '1'
        elif library_type == 'fr-unstranded':
            new = '0'

    elif prog == "Bowtie2":
        if library_type == 'fr-firststrand' and paired==True:
            new = '--fr'
        elif library_type == 'fr-secondstrand' and paired==True:
            new = '--rf'
        else:
            new = ''

    return new


overall_library_type = []
overall_paired = []

for infile in infiles:
    paired = False
    with open(infile, "r") as f:
        strands = OrderedDict()
        for line in f.readlines():
            line = line.strip()
            if line.startswith('This is PairEnd Data'):
                paired = True

            if line.startswith('Fraction of reads explained by "'):
                values = line.replace('Fraction of reads explained by "','').split('": ')
                strands[values[0]] = float(values[1])
            if line.startswith('Fraction of reads explained by other combinations: '):
                value = float(line.replace('Fraction of reads explained by other combinations: ',''))
                if value >= 0.2:
                    print("ERROR: A larger fraction %s of reads is explained by uncommon strand combinations!".format(value))
                    print (infile)
                    exit(1)

        if len(strands.keys()) != 2:
            print("ERROR: Unclear strand-specificity in:")
            print(infile)
            exit(1)

        threshold = 0.6      #min quotient threshold

        k = strands.keys()
        v = strands.values()

        specificity = "fr-unstranded"
        if '++,--' in strands.keys() and '+-,-+' in strands.keys():
            if strands['++,--'] >= threshold and strands['+-,-+'] <= threshold:
                specificity = "fr-secondstrand"
            elif strands['++,--'] <= threshold and strands['+-,-+'] >= threshold:
                specificity = "fr-firststrand"
        if '1++,1--,2+-,2-+' in strands.keys() and '1+-,1-+,2++,2--' in strands.keys():
            if strands['1++,1--,2+-,2-+'] >= threshold and strands['1+-,1-+,2++,2--'] <= threshold:
                specificity = "fr-secondstrand"
            elif strands['1++,1--,2+-,2-+'] <= threshold and strands['1+-,1-+,2++,2--'] >= threshold:
                specificity = "fr-firststrand"

        overall_library_type.append(specificity)
        overall_paired.append(paired)


## majority library_type
paired = Counter(overall_paired).most_common()[0][0]
library_type = Counter(overall_library_type).most_common()[0][0]


##print library_type[0][0]

if paired:
    print("library_type\tpaired-end")
else:
    print("library_type\tsingle-end")

##print("TopHat2\t{}".format(library_type))
print("RSeQC\t{}".format(convert_library_type( library_type, 'RSeQC', paired ) ))
##print("HISAT2\t{}".format(convert_library_type( library_type, 'HISAT2', paired ) ))
##print("htseq-count\t{}".format(convert_library_type( library_type, 'htseq-count', paired ) ))
print("Bowtie2\t{}".format(convert_library_type( library_type, 'Bowtie2', paired ) ))
print("featureCounts\t{}".format(convert_library_type( library_type, 'featureCounts', paired ) ))
