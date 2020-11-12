# !/usr/bin/env python3
# Construct the De Bruijn Graph of a Collection of k-mers
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name:(zmmason)
# Group Members: NONE

import sys


class DeBruilinFromCol:
    """
    Given an arbitrary collection of k-mers Patterns, we define
    CompositionGraph(Patterns) as a graph with |Patterns| isolated edges. Every edge is labeled by a k-mer from
    Patterns, and the starting and ending nodes of an edge are labeled by the prefix and suffix of the k-mer labeling
    that edge. We then define the de Bruijn graph of Patterns, denoted DeBruijn(Patterns), by gluing identically
    labeled nodes in CompositionGraph(Patterns), which yields the following algorithm.
    minor influence from https://github.com/TatyanaV/Genome_Sequencing_Bioinformatics_II/blob/master/5.
    DeBruijnGraphfromk-mers.py
    """

    def __init__(self, kmers):
        """constructor: saves the kmer list"""
        self.kmers = kmers

    def deBruijnGraph(self):
        """gets identically labeled nodes from kmers with isolated edge between prefix and suffix"""
        patterns = {}  # holds the set of the prefix and suffixes
        adjacencyList = []  # holds the adjacency list to be printed
        for i in range(len(self.kmers)):
            try:
                patterns[self.kmers[i][:-1]].append(self.kmers[i][1:])  # linking kmer prefix to suffix
            except:
                patterns[self.kmers[i][:-1]] = [self.kmers[i][1:]]  # handles the case where prefix = suffix
        for key, val in patterns.items():
            adjacencyList.append(key + ' -> ' + ','.join(sorted([v for v in val])))  # sorts the internal edges (suffix)
        adjacencyList.sort()  # sorting in lexicographical order
        return adjacencyList

def main():
    """
    Given an arbitrary collection of k-mers Patterns (where some k-mers may appear multiple times), we define
    CompositionGraph(Patterns) as a graph with |Patterns| isolated edges. Every edge is labeled by a k-mer from
    Patterns, and the starting and ending nodes of an edge are labeled by the prefix and suffix of the k-mer labeling
    that edge. We then define the de Bruijn graph of Patterns, denoted DeBruijn(Patterns), by gluing identically
    labeled nodes in CompositionGraph(Patterns), which yields the following algorithm.
    """
    # kmers = ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']  # sample input
    kmers = []
    for line in sys.stdin:  # takes STDIN on#ly
        kmers.append(line.replace('\n', ''))
    construct = DeBruilinFromCol(kmers)
    adjacencyList = construct.deBruijnGraph()
    for i in adjacencyList:
        print(i)


if __name__ == '__main__':
    main()
