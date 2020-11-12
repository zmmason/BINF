# !/usr/bin/env python3
# University of California, Santa Cruz - BME 160
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)

import sys


class FastAreader:
    """ 
    (Previously used class) Reads a file of fasta sequences from STDIN and returns the 
    header and sequence
    """

    def __init__(self, fname=''):
        """contructor: saves attribute fname"""
        self.fname = fname

    def doOpen(self):
        """Handle file opens, allowing STDIN"""
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """Read an entire FastA record and return the sequence header/sequence"""
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
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


class SetFinder:
    """
    finds all unique and essential subsequences that occur in a tRNA sequence.
    Each set of unique subsequences are minimized to be specific only to that tRNA set
    Prins unique elements in ordered by their starting position in the tRNA sequence.
    """

    def __init__(self):
        """ constructor: saves important frequently used lists"""
        self.tRNAsets = []
        self.headers = []
        self.sequences = []

    def computetRNAset(self, header, seq):
        """
        computes the set of subsequences for each tRNA sequence
        saves the sets along with their original header and origin sequence 
        """
        subStrings = []
        self.headers.append(header)  # adds each header to its the constructor list
        self.sequences.append(seq)  # adds each sequence to its the constructor list
        for i in range(len(seq)):
            for j in range(i + 1, len(seq) + 1):  # creating each individual subsequence
                subStrings.append((seq[i:j]))
        self.tRNAsets.append(set(subStrings))  # adds each tRNA subsequence set to its the constructor list

    def findUniques(self):
        """ finds all unique subsequences in each tRNA set"""
        uniqueSets = []
        for i in range(len(self.tRNAsets)):  # compares each tRNA set against every other tRNA set
            otherSet = []
            for j in range(len(self.tRNAsets)):  # saves all 'other' tRNA sets to a specific list
                if i != j:  # takes all sets that dont include the specific target set
                    otherSet.append(self.tRNAsets[j])
            otherSet = set.union(*otherSet)  # takes the union of all the 'other' tRNA sets to make one set
            uniqueSets.append(set(self.tRNAsets[i].difference(otherSet)))
        return uniqueSets

    def essentialFinder(self, uniqueSets):
        """ 
        finds essential uniques in each tRNA set where no member of that 
        unique subsequence set is a substring of any other member of that set 
        """
        uniqueEssentials = []
        for uniquetRNAset in uniqueSets:
            nonEssential = set()  # creating a set of nonEssential subsequences that will come from tRNA unique sets
            for subString in uniquetRNAset:  # finding all nonEssentials by searching the tRNA set for larger variations of a subsequence
                if subString[1:] in uniquetRNAset:
                    nonEssential.add(subString)
                if subString[:-1] in uniquetRNAset:
                    nonEssential.add(subString)
                    uniqueEssential = uniquetRNAset - nonEssential  # removes all non-essentials from original unique tRNA set
            uniqueEssentials.append(uniqueEssential)
        return uniqueEssentials

    def printtRNAset(self, uniqueEssentials):
        "prints unique tRNA elements ordered and alligned by their starting position"

        for index in range(len(self.tRNAsets)):  # prints each data component of all tRNA's
            essentialUniquetoPrint = []
            print(self.headers[index])
            print(self.sequences[index])
            for substring in uniqueEssentials[index]:  # creating format printing each unique
                locations = [i for i in range(len(self.sequences[index])) if self.sequences[index].startswith(substring, i)]
                for locationIndex in locations:
                    essentialUniquetoPrint.append(f"{locationIndex * '.'}{substring}")  # formatting
            essentialUniquetoPrint.sort(key=len)  # sort output to be alligned in order by location
            for item in essentialUniquetoPrint:
                print(item)


def main():
    """
    Program reads a file of fasta sequences from STDIN, finds all unique subsequences 
    that occur in each tRNA. Each set is minimized such that no member of that tRNA's
    unique subsequence set is a substring of any other member of this set. 
    The printed unique elements are ordered and alligned by their starting position in 
    the host tRNA sequence.
    """
    myReader = FastAreader()  # takes STDIN
    mySetFinder = SetFinder()
    headSeq = {}  # creates a dictionary to hold name, sequence combinations
    for name, sequence in myReader.readFasta():  # gets name, sequence from file
        name = name.replace(' ', '')
        sequence = sequence.replace('-', '')
        headSeq.update({name: sequence})  # adds name:sequence to dictionary
    for key in sorted(headSeq):  # alphabetically sorts dictionary by key
        header = key
        seq = headSeq[key]
        mySetFinder.computetRNAset(header, seq)  # sending sorted header, sequence to SetFinder class to be analyzed
    myUniques = mySetFinder.findUniques()
    myUniqueEssentials = mySetFinder.essentialFinder(
        myUniques)  # taking the calculated unique substrings to be further analyzed (essential)
    mySetFinder.printtRNAset(myUniqueEssentials)  # prints data


if __name__ == "__main__":
    main()
