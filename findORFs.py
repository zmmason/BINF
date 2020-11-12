# !/usr/bin/env python3
# University of California, Santa Cruz - BME 160
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)

import sys

class FastAreader:
    """
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    """

    def __init__(self, fname=''):
        """contructor: saves attribute fname """
        self.fname = fname

    def doOpen(self):
        """ Handle file opens, allowing STDIN."""
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """ Read an entire FastA record and return the sequence header/sequence"""
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class OrfFinder:
    """
    after 'findORFs' parses a fastA file using 'fastAreader,' this program searches open reading frames
    returning putative genes determined and defined by the location of the start and stop codons. The program will
    output: ORF positions (1,2,3,-1,-2,-3), gene start location, gene stop location, and gene length, ordered by
    decreasing gene length:
    - program uses an optional set of start codons - ATG, GTG, TTG
    - program uses an optional set of stop codons - TAG, TGA, TAA
    - returns an ordered list of every putative gene (>100nt) in an ORF rather than only the largest.

    """

    def __init__(self):
        """
        initializing start/stop codons as well as the list that will hold all ORF/data information
        """
        self.startCodon = ['ATG', 'GTG', 'TTG']  # initializing start codons
        self.stopCodon = ['TAG', 'TAA', 'TGA']  # initializing stop codons
        self.allGenes = []  # initializing list for ORF/gene data

    def ORFgen(self, seq, compSeq):
        """
        object that takes incomming sequences and searches their ORF's for specific data related
        to start/stop locations. those start and stop locations are used to determine gene length.
        """
        startIndex = []
        stopIndex = []
        revStartIndex = []
        revStopIndex = []
        for i in range(0, 3):  # setting up ORFs for the forward strand (+1,+2,+3)
            """
            FORWARD STRAND calculations and data:
            This block of code is the same for both the forward and  compliment strands
            **** The comments for the forward strand will pertain to both forward and compliment strands****
            """
            orfPosition = i + 1
            for position in range(i, len(seq[i:]), 3):  # takes position every 3 nts
                codon = seq[position:position + 3]  # creates a codon (str) every 3 nts
                if codon in self.startCodon:  # adds all start codon occurrences to one list
                    startIndex.append(position + 1)
                if codon in self.stopCodon:  # adds all stop codon occurrences to one list
                    stopIndex.append(position + 3)
            for i in range(len(stopIndex)):
                for j in range(len(startIndex)):
                    if startIndex[j] < stopIndex[i]:
                        geneLen = (stopIndex[i] - startIndex[j] + 1)  # gets length of the gene
                        self.allGenes.append([orfPosition, startIndex[j], stopIndex[i], geneLen])
                        startIndex = [elem for elem in startIndex if elem > stopIndex[
                            i]]  # updates 'startIndex' by removing all starts up to the used stop
                    break
            startIndex.clear()
            stopIndex.clear()

        for i in range(0, 3):
            revOrfPosition = -1 * (i + 1)
            for position in range(i, len(compSeq[i:]), 3):
                codon = compSeq[position:position + 3]
                if codon in self.startCodon:
                    revStartIndex.append(position + 1)
                if codon in self.stopCodon:
                    revStopIndex.append(position + 3)

            for i in range(len(revStopIndex)):
                for j in range(len(revStartIndex)):
                    if revStartIndex[j] < revStopIndex[i]:
                        revGeneLen = (revStopIndex[i] - revStartIndex[j] + 1)
                        self.allGenes.append([revOrfPosition, revStartIndex[j], revStopIndex[i], revGeneLen])
                        revStartIndex = [elem for elem in revStartIndex if elem > revStopIndex[i]]
                    break
            revStartIndex.clear()
            revStopIndex.clear()

        self.allGenes.sort(key=lambda x: (x[3], x[0]),
                           reverse=True)  # sorting by increasing length, then by left position of gene
        return self.allGenes


def main ():
    """
    This program takes the sequences present in a fastA file  and searches open reading frames. Putative genes can then
    be determined and defined in those ORF's by finding the location of a start codon at the 5' end of the ORF and a
    stop codon at the 3' end. The program will output:
    + ORF positions (1,2,3,-1,-2,-3), gene start location, gene stop location, and gene length.
    - uses STDIN for input and STDOUT for output
    - print every putative gene (>100nt) in an ORF instead of only the largest.
    - program can vary in ORF size minimums based on a single inequality that can be altered
    """
    compDict = {'A':'t', 'T': 'a', 'G':'c', 'C':'g'} # nuc compliments for reverse strand
    myReader = FastAreader()  # make sure to change this to use stdin if using a filename for testing
    myOrf = OrfFinder() 
    for name, seq in myReader.readFasta():
        compSeq = seq[::-1]  # reverses the strand to be prepped for nt compliments
        for nuc in compDict:
            compSeq = compSeq.replace(nuc, compDict[nuc])  # replaces nt's with their compliments
        compSeq = compSeq.upper()
        data = myOrf.ORFgen(seq, compSeq) 

    print(name)
    for i in data:
        if i[3] > 100: # sets the length minimum for putative genes to be printed
            print('{:+d} {:>5d}..{:>5d} {:>5d}' .format(i[0], i[1], i[2], i[3]))


if __name__ == "__main__":
    main()
