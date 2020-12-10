# !/usr/bin/env python3
# HMM Set 2
# Solve the Soft Decoding Problem
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys
import numpy as np  # numpy version 1.19.3


class softDecode:
    """
    Compute the The conditional probability Pr(x|π) that will be emitted by the HMM.
    Input:  A string x, followed by the alphabet Σ from which x was constructed, followed by the states States,
            transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
    """

    def __init__(self, emission, alphabet, states, transMatrix, emisMatrix):
        """Constructor: saves attributes from the input file."""
        self.emission = [em for em in emission]  # list of individual emissions in emission seq
        self.alphabet = alphabet
        x = 0  # counter to hold value for alphabet dict
        self.alphabetDict = {i: 0 for i in alphabet}  # initialize dict holding alphabet
        for i in self.alphabetDict:
            self.alphabetDict[i] = x
            x += 1
        self.states = states
        y = 0  # counter to hold value for states dict
        self.statesDict = {i: 0 for i in states}  # initialize dict holding states
        for i in self.statesDict:
            self.statesDict[i] = x
            y += 1
        self.transMatrix = transMatrix  # takes the log of each item in the transition matrix
        self.emisMatrix = emisMatrix  # takes the log of each item in the emission matrix
        self.endProb = 1 / len(self.emission)  # probability of ending node
        self.startProb = 1 / len(self.states)  # probability of starting node

    def forwardAlgo(self):
        """Calculate forward probabilities."""
        fwdNode = np.zeros((len(self.states), len(self.emission)))
        for i, alphabetVal in enumerate(self.emission):
            for j in range(len(self.states)):
                if i == 0:  # if first alphabet value then do this
                    fwdNode[j, i] = self.startProb * self.emisMatrix[j, self.alphabetDict[alphabetVal]]
                else:
                    values = [fwdNode[k, i - 1] * self.emisMatrix[j, self.alphabetDict[alphabetVal]] *
                              self.transMatrix[k, j] for k in range(len(self.states))]
                    fwdNode[j, i] = sum(values)
        finalState = np.multiply(fwdNode[:, -1], self.endProb)
        sinkVal = sum(finalState)
        return fwdNode, sinkVal

    def backwardAlgo(self):
        """Calculate backward probabilities."""
        bwdNode = np.zeros((len(self.states), len(self.emission)))
        for i in range(1, len(self.emission)+1):
            for j in range(len(self.states)):
                if i == 1:  # if first alphabet value then do this
                    bwdNode[j, -i] = self.endProb
                else:
                    values = [bwdNode[k, -i + 1] * self.emisMatrix[k, self.alphabetDict[self.emission[-i + 1]]] *
                              self.transMatrix[j, k] for k in range(len(self.states))]
                    bwdNode[j, -i] = sum(values)
        return bwdNode

    def forwardBackwardAlgo(self, forwarProbs, forwardSink, backwardProbs):
        """Calculate conditional probability Pr(πi = k|x) for each emission state"""
        condProb = np.zeros((len(self.states), len(self.emission)))
        for i in range(len(self.emission)):
            for j in range(len(self.states)):
                condProb[j, i] = (forwarProbs[j, i] * backwardProbs[j, i]) / forwardSink
        return condProb

    def printFormat(self, condProbs):
        """Format data for desired printing format"""
        condProbs = np.round(condProbs, 4)
        condProbs = condProbs.tolist()
        # creates list of probabilities in correspondence to their origin state
        print(*self.states, sep='\t')
        for i in range(len(self.emission)):
            indexProb = []
            for j in range(len(condProbs)):
                indexProb.append(condProbs[j][i])
            print(*indexProb, sep='\t')


def main():
    """Solve the Soft Decoding Problem."""
    contents = []  # list to hold the contents of the dataset
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.strip())
    statesLen = len(contents[6].split())  # counts amount of stated for to help parse input data
    transMatrix = np.array([line.split()[1:] for line in contents[7: statesLen + 7]], float)  # create array t-matrix
    emisMatrix = np.array([line.split()[1:] for line in contents[len(contents) - statesLen:]], float)  # array e-matrix
    decode = softDecode(contents[0], contents[2].split(), contents[4].split(), transMatrix, emisMatrix)
    forwarProbs, forwardSink = decode.forwardAlgo()
    backwardProbs = decode.backwardAlgo()
    condProbs = decode.forwardBackwardAlgo(forwarProbs, forwardSink, backwardProbs)
    decode.printFormat(condProbs)


if __name__ == '__main__':
    main()
