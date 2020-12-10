# !/usr/bin/env python3
# HMM Set 1
# Compute the Probability of a Hidden Path
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys


class ProbHiddenPath:
    """Compute the Probability of a Hidden Path."""

    def __init__(self, pi, matrixA, matrixB):
        """Constructor: saves attributes from the input file."""
        self.piPath = pi  # saves the sequence
        self.states = [matrixA[0], matrixB[0]]  # saves the states of the sequence
        matrixA = matrixA[1:].split()  # parsing states matrix
        matrixB = matrixB[1:].split()  # parsing states matrix
        #  creating the transition matrix from each state matrix
        self.transition = [[float(matrixA[0]), float(matrixA[1])], [float(matrixB[0]), float(matrixB[1])]]

    def hiddenPath(self):
        """
        Calculate the probability of a hidden path with the given HMM transition states assuming original
        have equal probabilities (Pr(A) = Pr(B)).
        """
        path = [self.states.index(state) for state in self.piPath]  # creates list of states in path
        probPi = 0.5  # initial probability - prob of obtaining a specific character (A,B -> 1/2, A,B,C -> 1/3, etc.)
        prev = path[0]  # holds the previous state in path for calculation of probabilities
        for state in path[1:]:
            trans_pr = self.transition[prev][state]
            probPi = probPi * trans_pr  # updates pi prob by getting the probability of the state * previous prob
            prev = state  # resets the 'prev' path index
        return probPi


def main():
    """
    Compute the Probability of a Hidden Path given a hidden path (pi), states of the path (states), and its
    transition matrix of an HMM (transition).
    """
    # sample input
    # contents = ['AABBBAABABAAAABBBBAABBABABBBAABBAAAABABAABBABABBAB',
    # '--------',
    # 'A   B',
    # '--------',
    # '    A   B',
    # 'A   0.194   0.806',
    # 'B   0.273   0.727']

    contents = []  # list to hold the contents of the dataset
    for line in open('sample.txt'):  # takes STDIN only
        contents.append(line.strip())
    path = ProbHiddenPath(contents[0], contents[5], contents[6])
    pathProbability = path.hiddenPath()
    print(pathProbability)


if __name__ == '__main__':
    main()
