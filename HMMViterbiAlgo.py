# !/usr/bin/env python3
# HMM Set 1
# Implement the Viterbi Algorithm
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys
import numpy as np  # numpy version 1.19.3


class ViterbiAlgorithm:
    """
    Compute the The conditional probability Pr(x|π) that will be emitted by the HMM.
    Input:  A string x, followed by the alphabet Σ from which x was constructed, followed by the states States,
            transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
    """

    def __init__(self, emission, alphabet, states, transMatrix, emisMatrix):
        """Constructor: saves attributes from the input file."""
        self.emission = [int(alphabet.index(em)) for em in emission]  # list of index codes for each em in emission seq
        self.alphabet = alphabet
        self.states = states
        self.transMatrix = np.log(transMatrix)  # takes the log of each item in the transition matrix
        self.emisMatrix = np.log(emisMatrix)  # takes the log of each item in the emission matrix
        # initializing viterbi path with theo max lowest bound & list to hold backpointers -> keeping track of the paths
        self.viterbi = np.array([[float(-9999) for j in range(len(self.emission))] for i in range(len(self.states))])
        self.backPointers = [[0 for j in range(len(self.emission))] for i in range(len(self.states))]

    def viterbiAlgo(self):
        """Get the path that maximizes Pr(x, π) over all possible paths π."""
        for state in range(len(self.states)):  # initializing start of viterbi with Pr(emisison)*(1/#states)
            self.viterbi[state][0] = np.log(1/len(self.states))+self.emisMatrix[state][self.emission[0]]
            self.backPointers[state][0] = -1  # sets starting nodes to -1 key
        for i in range(1, len(self.emission)):  # creating full viterbi graph after the initial node
            for state in range(len(self.states)):
                for prev in range(len(self.states)):  # getting total probability for
                    totalProb = self.emisMatrix[state][self.emission[i]] + \
                              self.transMatrix[prev][state] + self.viterbi[prev][i - 1]
                    if totalProb > self.viterbi[state][i]:  # find max-weight path to current node
                        self.viterbi[state][i] = totalProb
                        self.backPointers[state][i] = prev
        score = float(-9999)  # start backtrack from max-likelihood using max lowest value
        for state in range(len(self.states)):
            if self.viterbi[state][len(self.emission) - 1] > score:
                last = state
                score = self.viterbi[state][len(self.emission) - 1]
        path = [last]
        i = len(self.emission) - 1  # create max likelihood path backwards
        while i > 0:  # verification to go to next node
            next = self.backPointers[last][i]
            path.append(next)
            last = next
            i -= 1
        result = ''.join(str(self.states[state]) for state in path[::-1])
        return result


def main():
    """Compute the The conditional probability Pr(x|π) that will be emitted by the HMM."""
    # sample input
    # contents = ['xyxzzxyxyy', '--------', 'x   y   z', '--------', 'A   B', '--------', '    A   B', 'A   0.641   0.359', 'B   0.729   0.271', '--------', '    x   y   z', 'A   0.117   0.691   0.192', 'B   0.097   0.42    0.483']

    contents = []  # list to hold the contents of the dataset
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.strip())
    statesLen = len(contents[4].split())
    transMatrix = np.array([line.split()[1:] for line in contents[7: statesLen + 7]], float)
    emisMatrix = np.array([line.split()[1:] for line in contents[len(contents) - statesLen:]], float)
    viterbi = ViterbiAlgorithm(contents[0], contents[2].split(), contents[4].split(), transMatrix, emisMatrix)
    print(viterbi.viterbiAlgo())


if __name__ == '__main__':
    main()
