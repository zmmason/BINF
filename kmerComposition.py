# !/usr/bin/env python3
# Generate the k-mer Composition of a String
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys


class KmerGen:
    """reads a txt file and returns all kmers in seq of specified size"""
    def __init__(self, k, string):
        """constructor: saves kmer size, string"""
        self.k = int(k)
        self.string = string

    def getkmers(self):
        """grabs kmers from the string"""

        for i in range(0, len(self.string) - self.k, 1):
            kmer = self.string[i:i + self.k]
            if len(kmer) == self.k:  # makes sure kmer is k length
                print(kmer)


def main():
    """reads a txt file and returns all kmers in seq of specified size"""

    # contents = ['5', 'CAATCCAAC']  # sample input
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line)  # adds each line in file to a list
    k = contents[0].replace('\n', '')  # first item should be k
    string = contents[1]  # second item should be the seq
    myKmer = KmerGen(k, string)
    myKmer.getkmers()


if __name__ == '__main__':
    main()
