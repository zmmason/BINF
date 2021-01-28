# !/usr/bin/env python3
# Burrows–Wheeler transform
# Estimate the Parameters of an HMM
# University of California, Santa Cruz - BME 230A
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

class BWT:
    def makeBwt(t):
        """Create the BWT for the string t$"""

        t = t + '$'
        words = list(t)
        suff = []
        for i in range(len(words)):
            word = t[-1] + t[:-1]
            new = ''.join(word)
            t = new
            suff.append(new)
        sortSuff = sorted(suff)
        final = ''
        for i in range(len(words)):
            element = sortSuff[i]
            last = element[-1]
            final += last
        return final


    def invertBwt(bwt):
        """Inverts the Burrows-Wheeler Transform, returning the original string using
        inefficent algorithm"""

        matrix = [''] * len(bwt)
        for i in bwt:
            matrix = sorted(j + k for j, k in zip(bwt, matrix))
        for i in matrix:
            if i[-1] == '$':
                # print(i)
                return i
                break


class FmIndex(object):
    def __init__(self, t, alphabet):
        """ Create FM-index for t in naive manner """

        self.bwt = BWT.makeBwt(t)
        s = sorted(self.bwt)
        self.C = {}
        for i in range(len(s) - 1, -1, -1):
            self.C[s[i]] = i
        self.Occ = [{} for i in range(len(self.bwt))]
        for i in range(len(self.bwt)):
            for j in alphabet + "$":
                p = self.Occ[i - 1][j] if i > 0 else 0
                self.Occ[i][j] = p + (1 if self.bwt[i] == j else 0)

    def lf(self, i):
        """ Return the last-to-first mapping for index i of the bwt """

        x = list(self.C)
        lf = self.C.get(x[i - 1]) + self.Occ[i - 1][x[i - 1]]
        return lf

    def invertBwtUsingFmIndex(fmIndex):
        """ Returns t by using repeated lf search to reconstruct t$ backwards"""
        bwt = fmIndex.bwt
        first = fmIndex.bwt.index('$')
        t = ['$']
        matrix = [''] * len(bwt)
        for i in bwt:
            matrix = sorted(j + k for j, k in zip(bwt, matrix))
        for i in matrix:
            if i[-1] == '$':
                return i
                break

