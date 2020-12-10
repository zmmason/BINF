# !/usr/bin/env python3
# HMM Set 1
# Compute the Probability of a String Emitted by an HMM
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys
import numpy as np  # numpy version 1.19.3


class ViterbiAlgorithm:
    """
    Compute the Probability of a String Emitted by an HMM.
    Input:  A string x, followed by the alphabet Σ from which x was constructed, followed by the states States,
            transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: The probability Pr(x) that the HMM emits x.
    """

    def __init__(self, emission, alphabet, states, transMatrix, emisMatrix):
        """Constructor: saves attributes from the input file."""
        self.emission = [int(alphabet.index(em)) for em in emission]  # gets the index val for each character in em
        self.states = states
        self.transMatrix = np.exp(np.log(transMatrix))  # gets log of each item in the transition matrix
        self.emisMatrix = np.exp(np.log(emisMatrix))  # takes the log of each item in the emission matrix

    def outcome_likelihood(self):  # outcome-likelihood of HMM emitting emissions (sum of all hidden paths)
        """Create Viterbi, calculate and return built likelihood of emission string, given HMM information."""
        viterbi = [[0 for j in range(len(self.emission))] for i in range(len(self.states))]  # initialize viterbi
        viterbi = np.array(viterbi, float)  # creates array from initialized viterbi
        for state in range(len(self.states)):  # gets initial column of viterbi
            viterbi[state][0] = (1/len(self.states)) * self.emisMatrix[state][self.emission[0]]  # Pr emission
        for i in range(1, len(self.emission)):  # iterates through emission Pr and sums new edge for each node
            for j in range(len(self.states)):
                em = self.emisMatrix[j][self.emission[i]]  # Pr of individual emission
                print(em)
                for k in range(len(self.states)):
                    viterbi[j][i] += self.transMatrix[j][k] * em * viterbi[k][i-1]  # summing edge for node
        # probability Pr(x) that the HMM emits x
        total = sum(viterbi[s][len(self.emission)-1] for s in range(len(self.states)))
        return total


def main():
    """Compute the Probability of a String Emitted by an HMM given a (string, alphabet, states, t-matrix, e-matrix)."""
    # sample input
    # contents = ['xzyyzzyzyy', '--------', 'x   y   z', '--------', 'A   B', '--------', '    A   B', 'A   0.303   0.697', 'B   0.831   0.169', '--------', '    x   y   z', 'A   0.533   0.065   0.402', 'B   0.342   0.334   0.324']

    contents = []  # list to hold the contents of the dataset
    for line in open('sample.txt'):  # takes STDIN only
        contents.append(line.strip())
    statesLen = len(contents[4].split())  # counts amount of stated for to help parse input data
    transMatrix = np.array([line.split()[1:] for line in contents[7: statesLen + 7]], float)  # create array t-matrix
    emisMatrix = np.array([line.split()[1:] for line in contents[len(contents) - statesLen:]], float)  # array e-matrix
    viterbi = ViterbiAlgorithm(contents[0], contents[2].split(), contents[4].split(), transMatrix, emisMatrix)
    print(viterbi.outcome_likelihood())


if __name__ == '__main__':
    main()
