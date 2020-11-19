# !/usr/bin/env python3
# HMM Set 1
# Problem 20 - Compute the Probability of an Outcome Given a Hidden Path
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys
import numpy as np  # numpy version 1.19.3


class ViterbiAlgorithm:
    """Compute the The conditional probability Pr(x|π) that will be emitted by the HMM given the hidden path π, states,
    and transition matrix of an HMM."""

    def __init__(self, emission, alphabet, pi, states, matrix):
        """ Constructor: saves attributes from the input file. """
        self.emission = emission
        self.alphabet = alphabet
        self.pi = pi
        self.states = states
        self.matrix = matrix

    def emissionProb(self):
        """ Calculate probability of this emission. """
        pathPi = [self.states.index(state) for state in self.pi]  # initializing states as index for reference/counting
        emissions = [self.alphabet.index(em) for em in self.emission]  # initializing em as index for reference/counting
        prob = 1  # initializing the probability 
        for i in range(len(emissions)):  # iteration to use product rule to get full emissioin prob
            prob = prob * self.matrix[pathPi[i]][emissions[i]]
        return prob


def main():
    """Compute the The conditional probability Pr(x|π) that will be emitted by the HMM given the hidden path π, states,
    and transition matrix of an HMM."""
    # sample input
    # contents = ['xxyzyxzzxzxyxyyzxxzzxxyyxxyxyzzxxyzyzxzxxyxyyzxxzx',
    # '--------',
    # 'x   y   z',
    # '--------',
    # 'BBBAAABABABBBBBBAAAAAABAAAABABABBBBBABAABABABABBBB',
    # '--------',
    # 'A   B',
    # '--------',
    # '    x   y   z',
    # 'A   0.612   0.314   0.074',
    # 'B   0.346   0.317   0.336',]

    contents = []  # list to hold the contents of the dataset
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.strip())
    matrix = np.array([[float(i[1:]) for i in contents[9][1:].split()],  # create matrix of the emission and alphabet
                       [float(i[1:]) for i in contents[10][1:].split()]])
    viterbi = ViterbiAlgorithm(contents[0], contents[2].split(), contents[4], contents[6].split(), matrix)
    print(viterbi.emissionProb())


if __name__ == '__main__':
    main()
