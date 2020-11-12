# !/usr/bin/env python3
# Reconstruct a String from its k-mer Composition
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys


class constructString:
    """Reconstruct a string from its k-mer composition."""

    def getDebrujin(self, oriKmers):
        """Get deBruijin graph information from the given k-mers."""
        kmers = []  # list to hold set of each kmer to be used to create set of suffixes
        for pattern in oriKmers:
            kmers = kmers + [pattern[:-1]]  # gets the prefix of each kmer
        kmers = set(kmers)  # creates set of all prefix's
        preSufdict = {}  # dict to hold prefix-suffix for each kmer
        for kmer in kmers:
            preSufdict[kmer] = []  # creates empty list to prep for suffix
        for kmer in oriKmers:
            preSufdict[kmer[:-1]].append(kmer[1:])  # adds prefix-suffix pairs to dict
        return preSufdict

    def exploreNodes(self, adjacentList):
        """Find the first and last kmers of the string by counting overlap in prefix-suffix"""
        nodePathLocation = dict.fromkeys(adjacentList.keys(), 0)
        for node in adjacentList.keys():
            for out in adjacentList[node]:
                nodePathLocation[node] -= 1
                try:
                    nodePathLocation[out] += 1
                except:
                    nodePathLocation[out] = 1
        return nodePathLocation

    def eulerianPath(self, preSufdict):
        """Begin referencing nodes for possible eulerian paths."""
        stack = []  # holds possible pathways to stack kmers
        nodePathLocation = self.exploreNodes(preSufdict)
        stack.append([k for k, v in nodePathLocation.items() if v == -1][0])  # gets the first prefix in path
        path = []  # list to hold building path
        while stack:  # tries all prefix-suffix for overlap
            new = stack[-1]
            try:
                w = preSufdict[new][0]  # if unused, adds prefix-suffix to stack
                stack.append(w)
                preSufdict[new].remove(w)
            except:
                path.append(stack.pop())  # if none, nothing to stack
        return path[::-1]

    def createPath(self, kmers):
        """Start building the string using kmers that follow eularian path."""
        buildSeq = ''  # holds the string to be built

        for kmer in kmers:
            buildSeq += kmer[0]  # ads the first index of each kmer to build entire string
        buildSeq += kmer[1:]  # adds the last suffix to the string
        return buildSeq


def main():
    """Reconstruct a string from its k-mer composition."""
    # contents = ['4', 'CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']  # sample input
    contents = []
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.replace('\n', ''))
    oriKmers = contents[1:]  # sends only the kmers to the functions
    construct = constructString()
    preSufdict = construct.getDebrujin(oriKmers)
    eulerian = construct.eulerianPath(preSufdict)
    finalPath = construct.createPath(eulerian)
    print(finalPath)


if __name__ == "__main__":
    main()
