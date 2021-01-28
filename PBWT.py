# !/usr/bin/env python3
# PBWT
# University of California, Santa Cruz - BME 230A
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import numpy as np


class PBWTfunctions:
    """
    Algorithm for finding matches between a set of aligned strings.
    INPUT: list of aligned strings or MxN matrix
    OUTPUT OPTIONS: Mx(N+1) matrix, MxN matrix, list of aligned strings, suffix matrix, long matches
    """
    def constructReversePrefixSortMatrix(X):
        """
        Compute Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the
        index in X of the ith string ordered by jth reverse prefix
        """
        # Creates the Mx(N+1) matrix
        A = np.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0]) + 1], dtype=int)

        # Code to write - you're free to define extra functions
        # (inline or outside of this function) if you like.

        A[:, 0] = list(range(len(X)))  # appending first column of array with index in len(X)
        for i in range(len(X[0])):
            a = []  # list to hold intermediates
            b = []
            for j in range(len(X)):  # updating intermediates
                if int(X[A[j, i]][i:i + 1]) == 0:
                    a.append(A[j][i])
                else:
                    b.append(A[j][i])
            A[:, i + 1] = a + b  # construct the new positional prefix element and appending to array
        return A

    def constructYFromX(X):
        """
        Y is a transformation of X. Using Y to construct X from Y,
        """
        # Creates the MxN matrix
        Y = np.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0])], dtype=int)

        # Code to write - you're free to define extra functions
        # (inline or outside of this function) if you like.

        A = PBWTfunctions.constructReversePrefixSortMatrix(X)  # getting Mx(N+1) matrix for X
        for i in range(len(Y)):
            for j in range(len(Y[i])):
                Y[i, j] = X[A[i, j]][j]  # Constructing the MxN matrix
        return Y

    def constructXFromY(Y):
        """
        Construct Y for X.
        """
        Z = np.flip(np.rot90(np.flip(Y)))  # re-orientating Y to prep for forward-iterations
        ref = [[i] for i in Z[0]]  # initializing list to hold sequence and reference elements
        for i in range(1, len(Z)):
            totalZero = list(Z[i]).count(0)  # total zeros in Y columns (Z rows)
            currentZeros = 0
            currentOnes = 0
            hold = []  # initializing list to hold intermediates
            for j in range(len(Z[i])):
                if Z[i, j] == 0:
                    currentZeros += 1
                    hold.append([0] + ref[currentZeros - 1])  # using previous seq to get new intermediates
                else:
                    currentOnes += 1
                    hold.append([1] + ref[totalZero + currentOnes - 1])  # using previous seq to build new intermediates
            ref = hold  # Re-initializing reference data to build strings
        X = np.array(ref)  # Creates the MxN matrix from the reference list
        print(X)
        return list(map(lambda i: "".join(map(str, i)), X))  # Convert back to a list of strings


    def constructCommonSuffixMatrix(A, X):
        """
        Define the common suffix of two strings to be the maximum length suffix shared
        by both strings.
        """
        D = np.zeros(shape=A.shape, dtype=int)  # Creates the Mx(N+1) D matrix
        for i in range(1, len(D)):
            for j in range(1, len(D[i])):
                maxSuffix = 0  # initializes count to hold number of matching suffix for each comparison
                for k in range(len(X[A[i, j]][:j])):  # iterate through all possible suffix matches of 1 comparison
                    if X[A[i, j]][:j][k] == X[A[i-1, j]][:j][k]:
                        maxSuffix += len(X[A[i, j]][:j][k])  # if suffix matches, add length to 'maxSuffix'
                    elif X[A[i, j]][k:] != X[A[i-1, j]][k:]:  # handle a reset if the last index in a comparison != match
                        maxSuffix = 0  # reset the suffix count when last index in a comparison != match
                D[i, j] = maxSuffix
        return D


    def getLongMatches(X, minLength):
        """
        For a pair of strings X[x], X[y], find long matches ending at j is a common substring
        of X[x] and X[y] that ends at j.
        """
        assert minLength > 0
        A = PBWTfunctions.constructReversePrefixSortMatrix(X)
        D = PBWTfunctions.constructCommonSuffixMatrix(A, X)
        # For each column, in ascending order of column index
        for j in range(1, 0 if len(X) == 0 else len(X[0])):
            b, c = [], []
            # Iterate over the aligned symbols in column j in reverse prefix order
            for i in range(len(X)):
                # For each string in the order check if there is a long match
                if D[i, j] < minLength:
                    for x in b:
                        for y in c:
                            # return the match as tuple of two sequence
                            yield (x, y, j) if x < y else (y, x, j)
                    b, c = [], []

                # Partition the sequences by if they have '0' or '1' at position j.
                if X[A[i, j]][j] == '0':
                    b.append(A[i, j])
                else:
                    c.append(A[i, j])

            # Report any leftover long matches for the column
            for x in b:
                for y in c:
                    yield (x, y, j) if x < y else (y, x, j)