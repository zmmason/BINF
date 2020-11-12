# !/usr/bin/env python3
# Reconstruct a String from its Genome Path
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys


class geneReconstruction:
    """A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1 symbols of Patterni are equal to the
    first k - 1 symbols of Patterni+1 for i from 1 to n-1. and returns A string Text of length k+n-1 where the i-th
    k-mer in Text is equal to Patterni for all i."""
    
    def __init__(self, contents):
        """constructor: saves kmer size, string"""
        self.contents = contents
        self.contentsLen = len(contents)

    def newString(self):
        """grabs kmers from the string"""
        string = self.contents[0]  # starts the return string with the first full kmer
        for i in range(1, self.contentsLen, 1):
            string += self.contents[i][-1]  # adds the last nt from each of the following kmers to build the string
        return string


def main():
    """A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1 symbols of Patterni are equal to the
    first k - 1 symbols of Patterni+1 for i from 1 to n-1. and returns A string Text of length k+n-1 where the i-th
    k-mer in Text is equal to Patterni for all i."""

    # contents = ['ACCGA', 'CCGAA', 'CGAAG', 'GAAGC', 'AAGCT']  # sample input
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.replace('\n', ''))  # adds each line in file to a list
    geneRecon = geneReconstruction(contents)
    print(geneRecon.newString())


if __name__ == '__main__':
    main()
