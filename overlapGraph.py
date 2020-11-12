# !/usr/bin/env python3
# Construct the Overlap Graph of a Collection of k-mers
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE
import sys


class OverlapCollection:
    """Given an arbitrary collection of k-mers Patterns, we form a graph having a node for each k-mer in Patterns and
     connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern'). The
     resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns)."""

    def __init__(self, seqs):
        """constructor: saves the sequences list"""
        self.seqs = seqs

    def collection(self):
        """ finds overlapping sequences and prints the pairs"""
        overlap = []  # list to hold pairs of overlapping sequences
        for i in range(0, len(self.seqs), 1):
            first = self.seqs[i]  # saves the primary sequence
            for j in range(0, len(self.seqs), 1):
                if first[1:] == self.seqs[j][0:-1]:  # gets the overlapping seq only if overlap is off by 1nt
                    second = self.seqs[j]
                    overlap.append(first + " -> " + second)  # formatting for printing
        overlap.sort()  # sorts by alphabetical by the primary seq
        for i in overlap:
            print(i)


def main():
    """Given an arbitrary collection of k-mers Patterns, we form a graph having a node for each k-mer in Patterns and
     connect k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern) is equal to Prefix(Pattern'). The
     resulting graph is called the overlap graph on these k-mers, denoted Overlap(Patterns)."""

    # contents = ['ATGCG', 'GCATG', 'CATGC', 'AGGCA', 'GGCAT']
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.replace('\n', ''))  # adds each line in file to a list
    construct = OverlapCollection(contents)
    construct.collection()


if __name__ == '__main__':
    main()
