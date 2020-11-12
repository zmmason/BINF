# !/usr/bin/env python3
# Rank motifs based on statistical under-representation
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys
import math
import argparse
import itertools


class CommandLine:
    """Handle the command line, usage and help requests"""

    def __init__(self, inOpts=None):
        """CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse."""
        self.parser = argparse.ArgumentParser(description='proper style BME205 program that reads a fasta file from '
                                                          'STDIN and ranks motifs based on how statistically '
                                                          'underrepresented the specific motif is.'
                                                          'Considered motif sequences are from 3 to 8 in length, '
                                                          'specified by minMotif and maxMotif. Specifing '
                                                          'the statistical cutoff will be done using negative '
                                                          'z-scores because the program is looking for '
                                                          'underrepresented sequences. Outputs are printed to a file '
                                                          'using STDOUT.',
                                              add_help=True,
                                              prefix_chars='-',
                                              usage='%(prog)s --minMotif int --maxMotif int --cutoff int < input.fa > '
                                                    'output.out')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        self.parser.add_argument('-mMin', '--minMotif', type=int, default=3, action='store',
                                 help='minimum motif size')
        self.parser.add_argument('-mMax', '--maxMotif', type=int, default=8, action='store',
                                 help='maximum motif size')
        self.parser.add_argument('-cut', '--cutoff', type=float, default=-5, action='store',
                                 help='z-score bases statistical cutoff')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class FastAreader:
    """
    (Class obtained from University of California, Santa Cruz - Biomolecular Engineering and Bioinformatics, BME 205)
    Reads a file of fasta sequences from STDIN and returns the header and sequence
    """

    def __init__(self, fname=''):
        """contructor: saves attribute fname"""
        self.fname = fname

    def doOpen(self):
        """Handle file opens, allowing STDIN"""
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """Read an entire FastA record and return the sequence header/sequence"""
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            line = fileH.readline()
            while not line.startswith('>'):
                if not line:  # we are at EOF
                    return header, sequence
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


class MissingMotif:
    """
    This program is designed to read a fasta file from STDIN and ranks motifs based on how statistically
    underrepresented that specific motif is. Using Z-scores for statistical cutoffs, Markov approximations are
    calculated to produce the expected count (E) of each K-mer.
    """

    def __init__(self, fullSeq, minMotif, maxMotif, cutoff):
        """
        Create and hold constructor information:
        Saves full sequence, min motif size, max motif size, and the z-score cutoff
        """
        self.fullSeq = fullSeq
        self.minMotif = minMotif
        self.maxMotif = maxMotif
        self.cutoff = cutoff
        self.motifList = []  # initializing the list of kmer dictionaries
        self.combined = {}  # initializing a dictionary of combined motif:rSeq

    def createKmer(self):
        """Create a list of all possible kmers from 1mers-maxkmers."""
        for kmer in range(1, self.maxMotif + 1, 1):
            kmerSet = {}  # creates dictionary to hold each kmer set created
            for itr in list(
                    itertools.product('ATCG', repeat=kmer)):  # starts each kmer at each nucleotide position in sequence
                motif = ''.join(itr)  # grabs kmer substring
                kmerSet[motif] = [0]  # initializing the motif in dict
                kmerSet[self.reverseComp(motif)] = kmerSet[motif]  # sets the reverse = fwd strand to be counted as 1
            self.motifList.append(kmerSet)

    def reverseComp(self, kmer):
        """Return the reverse compliment of the given Kmer sequence."""
        compDict = {'A': 't', 'T': 'a', 'G': 'c', 'C': 'g'}  # nuc compliments for reverse strand
        revMotif = kmer[::-1]  # reverses the strand to be prepped for nt compliments
        for nuc in compDict:
            revMotif = revMotif.replace(nuc, compDict[nuc])  # replaces nt's with their compliments
        revMotif = revMotif.upper()
        return revMotif

    def countKmer(self):
        """Search the .fa file for k-mers of 'k' size and count their frequencies."""
        canonical = ['A', 'T', 'C', 'G', 'a', 'c', 'g', 't']  # list of canonical bases
        for x in range(0, len(self.fullSeq), 1):
            for location in range(0, self.maxMotif, 1):
                kmerLoc = x + location + 1  # sets size of the kmer being located
                seqKmer = self.fullSeq[x:kmerLoc]  # gets the motif from the seq for comparison
                if any(base not in canonical for base in seqKmer):  # if the motif contains a non-canonical base, skips
                    continue
                if kmerLoc <= len(self.fullSeq):
                    self.motifList[location][self.fullSeq[x:kmerLoc]][0] += 1

    def kmerParser(self):
        """Parse the complete dict of kmers ('kmers':'count') and combine compliments ('fwd-kmer:rev-kmer':'count'). """
        usedMotif = []  # holds the used motifs to prevent reverse combinations
        for kmerSet in self.motifList:
            for motif in kmerSet:
                if motif not in usedMotif:  # makes sure only new (not including compliment strands) are added to list
                    rSeq = self.reverseComp(motif)  # get reverse compliment
                    usedMotif.append(motif)  # adds motif to used motif list
                    usedMotif.append(rSeq)  # adds motif to used motif list
                    motifCounts = kmerSet[motif][0]  # grabs the shared count for the motif/rSeq
                    self.combined['{0:8}:{1:8}'.format(motif, rSeq)] = motifCounts
                if motif in usedMotif:
                    pass
        return self.combined

    def kmerCalcs(self, combined):
        """calculate expected value and Z-score for each kmer."""
        combined = combined.replace(' ', '')
        splitMotifs = combined.split(':')  # getting only the motif of the motif:rSeq pair
        motif = splitMotifs[0]
        rSeq = splitMotifs[1]
        if len(motif) >= self.minMotif:  # can only use kmers >= 3
            outStart = self.motifList[len(motif)-2][motif[:-1]][0]  # gets the count of the motif[0:]
            outStop = self.motifList[len(motif)-2][motif[1:]][0]  # gets the count of the motif[1:-1]
            mid = self.motifList[len(motif)-3][motif[1:-1]][0]  # gets the count of the motif[1:-1]
            if mid == 0:  # handles first division by 0 possibility
                return None  # returns none in place of kmer data because None is handled in the following method
            else:
                e = (outStart * outStop) / mid  # expected value given the null distribution
                if e == 0:  # handles second division by 0 possibility
                    return None
                else:
                    s = self.motifList[len(motif)-1][motif][0]  # the specific kmer count being converted
                    n = len(self.fullSeq)  # approx genome length
                    p = e / n  # probability of success given the null distribution
                    m = n * p  # mean
                    sd = math.sqrt(m * (1 - p))  # standard deviation
                    zScore = (s - m) / sd  # zscore
                    if zScore <= self.cutoff:  # sets the zscore cutoff for printing
                        return motif, rSeq, s, e, zScore


def main(command=None):
    """
    This program is designed to search a fasta file for all possible motifs and outputs
    the total approx genome size, motif and its compliment strand, its occurrences in the fasta file,
    its expected occurrences, and its zscore. This program takes STDIN and grabs info from the user
    to create an output file of desired information relating to kmer range, and zscore cutoff
    """

    # input is a fasta file
    if command is None:
        myCommand = CommandLine()  # read options from the command line
    else:
        myCommand = CommandLine(command)  # interpret the list passed from the caller of main
    myReader = FastAreader()  # make sure to change this to use stdin if using a filename for testing
    minMotif = myCommand.args.minMotif
    maxMotif = myCommand.args.maxMotif
    cutoff = myCommand.args.cutoff
    allData = []  # need to create a list to hold all data to be sorted
    toPrint = []  # creates a list to hold all sorted kmers or similar size
    fullseq = ''  # initializes a string that will hold the full genome
    for header, sequence in myReader.readFasta():
        fullseq = fullseq + sequence  # adds each sequence in FASTA to create a single genome string
    print('N = {}'.format(len(fullseq)))  # prints the full length of the full genome
    myMotif = MissingMotif(fullseq, minMotif, maxMotif, cutoff)
    myMotif.createKmer()  # creates all kmers (1mer - maxMotif kmer)
    myMotif.countKmer()  # counts all kmers in the sequence
    myCount = myMotif.kmerParser()  # combines the motif and compliment motif with one count
    for combined in myCount:
        myData = myMotif.kmerCalcs(combined)
        if myData is not None:  # takes care of instances where None is present for kmer data
            allData.append(myData)
    for i in range(minMotif, maxMotif + 1):
        sizeList = []  # list to sorts motifs by size
        for j in allData:
            if len(j[0]) == i:
                sizeList.append(j)
        toPrint.append(sizeList)
    toPrint.reverse()  # reverses the order of the list to show largest kmers first
    for kmerSet in toPrint:
        kmerSet.sort(key=lambda x: x[4])  # sorts each kmer set by increasing zscores
        for i in kmerSet:
            print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(i[0], i[1], i[2], i[3], i[4]))


if __name__ == "__main__":
    main()
