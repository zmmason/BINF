# !/usr/bin/env python3
# HMM Set 2
# Implement Baum-Welch Learning
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)
# Group Members: NONE

import sys
import numpy as np  # numpy version 1.19.3


class BaumWelch:
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
        self.transMatrix = transMatrix
        self.emisMatrix = emisMatrix
        self.endProb = 1 / len(self.emission)  # probability of ending node
        self.startProb = 1 / len(self.states)  # probability of starting node

    def forwardAlgo(self):
        """Calculate forward probabilities."""
        fwdNode = np.zeros((len(self.states), len(self.emission)))  # initializing matrix for fwd probs
        for i, alphabetVal in enumerate(self.emission):
            for j in range(len(self.states)):
                if i == 0:  # handles first index (last in bwd)
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
        bwdNode = np.zeros((len(self.states), len(self.emission)))  # initializing matrix for backward probs
        for i in range(1, len(self.emission) + 1):
            for j in range(len(self.states)):
                if i == 1:  # handles first index (last in fwd)
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
            for j in range(len(self.states)):  # calculate forward-backward conditional probability
                condProb[j, i] = (forwarProbs[j, i] * backwardProbs[j, i]) / forwardSink
        return condProb

    def siProb(self, forwarProbs, forwardSink, backwardProbs):
        """Calculate si probabilities for transition matrix."""
        siMatrix = np.zeros((len(self.states), len(self.emission) - 1, len(self.states)))  # initializing si matrix
        for i in range(len(self.emission) - 1):
            for j in range(len(self.states)):
                for k in range(len(self.states)): # get si prob for matrix
                    siMatrix[j, i, k] = (forwarProbs[j, i] * backwardProbs[k, i + 1] * self.transMatrix[j, k] *
                                         self.emisMatrix[k, self.alphabetDict[self.emission[i + 1]]]) / forwardSink
        return self.convertTransMatrix(siMatrix)

    def gammaProb(self, forwarProbs, forwardSink, backwardProbs):
        """Calculate gamma probabilities for emission matrix."""
        gammaMatrix = np.zeros((len(self.states), len(self.emission)))  # initializing gamma matrix
        for i in range(len(self.emission)):
            for j in range(len(self.states)):  # get gamma prob for matrix
                gammaMatrix[j, i] = (forwarProbs[j, i] * backwardProbs[j, i]) / forwardSink
        return self.convertEmisMatrix(gammaMatrix)

    def convertTransMatrix(self, siMatrix):
        """Convert si probabilities into transition matrix for printing."""
        with np.errstate(divide='ignore'):  # handles division by 0 errors
            finalTransMatrix = np.zeros((len(self.states), len(self.states)))  # initializing transition prob matrix
            for i in range(len(self.states)):
                for j in range(len(self.states)):
                    for k in range(len(self.emission) - 1):
                        finalTransMatrix[i, j] = finalTransMatrix[i, j] + siMatrix[i, k, j]
                    denomenator = [siMatrix[i, t_x, i_x] for t_x in range(len(self.emission) - 1) for i_x in
                                   range(len(self.states))]  # gets denominator for calc
                    denomenator = sum(denomenator)
                    finalTransMatrix[i, j] = finalTransMatrix[i, j] / denomenator  # calc final transition probs
            return finalTransMatrix

    def convertEmisMatrix(self, gammaMatrix):
        """Convert gamma probabilities into emission matrix for printing."""
        with np.errstate(divide='ignore'):  # handles division by 0 errors
            finalEmMatrix = np.zeros((len(self.states), len(self.alphabetDict)))  # initializing emission prob matrix
            for i in range(len(self.states)):
                for j in range(len(self.alphabet)):
                    indices = [idex for idex, value in enumerate(self.emission) if value == self.alphabet[j]]
                    numerator = sum(gammaMatrix[i, indices])  # gets numerator for calc
                    denomenator = sum(gammaMatrix[i, :])  # gets denominator for calc
                    finalEmMatrix[i, j] = numerator / denomenator  # calc final emission probs
            return finalEmMatrix

    def dataPrint(self, transMat, emissMat):
        """Using the data from transition and emission matrices , print to specified formatting."""
        roundTrans = np.round(transMat, 3)  # rounding array to 3rd dec place for formatting
        roundEmis = np.round(emissMat, 3)
        roundTransList = roundTrans.tolist()  # turns matrix into list
        roundEmisList = roundEmis.tolist()
        maxMatrix = ['\t'.join([str(x) for x in self.states])]
        for i in range(len(roundEmisList)):  # printing emission matrix info
            roundTransList[i].insert(0, self.states[i])  # adding the state to the specific probabilities
            maxMatrix.append('\t'.join([str(x) for x in roundTransList[i]]))
        maxMatrix.append('--------')
        maxMatrix.append('\t' + '\t'.join([str(x) for x in self.alphabet]))
        for i in range(len(roundTransList)):  # printing transition matrix info
            roundEmisList[i].insert(0, self.states[i])  # adding the state to the specific probabilities
            maxMatrix.append('\t'.join([str(x) for x in roundEmisList[i]]))
        return maxMatrix


def main():
    """Solve the Soft Decoding Problem."""
    contents = []  # list to hold the contents of the dataset
    for line in sys.stdin:  # takes STDIN only
        contents.append(line.strip())
    iteration = int(contents[0])  # number of iterations
    statesLen = len(contents[6].split())  # counts amount of stated for to help parse input data
    emission = contents[2]
    alphabet = contents[4].split()
    states = contents[6].split()
    transMatrix = np.array([line.split()[1:] for line in contents[9: statesLen + 9]], float)  # create array t-matrix
    emisMatrix = np.array([line.split()[1:] for line in contents[len(contents) - statesLen:]], float)  # array e-matrix
    lastMatrixSet = [[transMatrix, emisMatrix]]
    for i in range(iteration):  # iterate over desired iteration count specified in file
        decode = BaumWelch(emission, alphabet, states, lastMatrixSet[0][0], lastMatrixSet[0][1])
        forwarProbs, forwardSink = decode.forwardAlgo()
        backwardProbs = decode.backwardAlgo()
        transMat = decode.siProb(forwarProbs, forwardSink, backwardProbs)
        emissMat = decode.gammaProb(forwarProbs, forwardSink, backwardProbs)
        lastMatrixSet[0] = [transMat, emissMat]  # replaces last matrix data with iteration generated matrices
    matrix = decode.dataPrint(transMat, emissMat)  # send last iteration to be formatted
    for data in matrix:
        print(data)


if __name__ == '__main__':
    main()
