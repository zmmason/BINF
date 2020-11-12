# !/usr/bin/env python3
# Construct the De Bruijn Graph of a String
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name:(zmmason)
# Group Members: NONE

import sys


class ConstructDeBruijin:
    """
    Given a genome Text, PathGraphk(Text) is the path consisting of |Text| - k + 1 edges, where the i-th edge of 
    this path is labeled by the i-th k-mer in Text and the i-th node of the path is labeled by the i-th (k - 1)-mer in 
    Text. The de Bruijn graph DeBruijnk(Text) is formed by gluing identically labeled nodes in PathGraphk(Text).
    """
    def __init__(self, k, seq):
        """constructor: saving k, sequence, nodes & edges sets"""
        self.k = int(k)
        self.seq = seq
        self.nodes = {}  # set to hold nodes
        self.edges = {}  # set to hold adjacent edges

    def createDeBuijin(self):
        """finds nodes and edges by finding overlapping kmers-1 in the original seq"""
        kmers = []
        for i in range(0, len(self.seq) - self.k + 2, 1):
            kmers.append(self.seq[i:i + self.k - 1])
        kmers = sorted(list(set(kmers)))
        for i, j in enumerate(kmers):  # creates the nodes for the nodes set
            self.nodes[i] = j
        nodeDict = {j: i for i, j in self.nodes.items()}  # creates a dict to hold each node for getting adjacent kmers
        for i in range(0, len(self.seq) - self.k + 1, 1):  # gets the edges based on overlapping kmers-1
            if nodeDict[self.seq[i:i + self.k - 1]] in self.edges:
                self.edges[nodeDict[self.seq[i:i + self.k - 1]]].append(nodeDict[self.seq[i + 1:i + self.k]])
            else:
                self.edges[nodeDict[self.seq[i:i + self.k - 1]]] = [nodeDict[self.seq[i + 1:i + self.k]]]

    def adjacencyList(self):
        """uses the node and edge lists to set up the adjacency list to be printed"""
        adjacencyList = []
        for node, edge in self.edges.items():
            self.nodes[node], '->', ','.join(sorted([self.nodes[val] for val in edge]))
            adjacencyList.append(self.nodes[node] + ' -> ' + ','.join(sorted([self.nodes[val] for val in edge])))
        adjacencyList.sort()  # sorts by lexicographical order
        return adjacencyList


def main():
    """
    Given a genome Text, PathGraphk(Text) is the path consisting of |Text| - k + 1 edges, where the i-th edge of 
    this path is labeled by the i-th k-mer in Text and the i-th node of the path is labeled by the i-th (k - 1)-mer in 
    Text. The de Bruijn graph DeBruijnk(Text) is formed by gluing identically labeled nodes in PathGraphk(Text).
    """
    # contents = ['4', 'AAGATTCTCTAC']  # sample input
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.replace('\n', ''))
    k = contents[0]
    seq = contents[1]
    construct = ConstructDeBruijin(k, seq)
    construct.createDeBuijin()
    for i in construct.adjacencyList():
        print(i)


if __name__ == '__main__':
    main()
