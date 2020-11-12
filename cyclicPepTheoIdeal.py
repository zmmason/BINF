# !/usr/bin/env python3
# Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys


class TheoPepCyclic:
    """
    Program to find a Cyclic Peptide with a Theoretical Spectrum Matching an Ideal Spectrum
    given A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
    """
    def __init__(self, spectrum):
        """Create a constructor to hold important data and tables."""
        self.spectrum = spectrum  # original input spectrum
        self.possiblePeptides = []
        self.aaMass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99,  # dict to hold the AA masses
                       'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                       'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                       'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
        self.aaList = list(self.aaMass)  # list of the AAs

    def matchTheoretical(self):
        """Create possible peptide strings given the input spectrum"""
        sortedSpectrum = sorted(self.spectrum)  # sorting input spectrum to be for verification of compatible peptides
        parentMass = sortedSpectrum[len(self.spectrum) - 1]  # gets the total mass of the input spectrum to verify
        peptides = []  # list to hold valid peptides
        for aa in self.aaList:  # Add new amino acid and corresponding masses
            newPair = []  # list to hold each sub peptide
            newAA = []  # list to hold new AA
            newMass = []  # list to hold each sub peptide mass
            newAA.append(aa)
            newMass.append(self.aaMass[aa])
            newPair.append(newAA)
            newPair.append(newMass)
            peptides.append(newPair)
        while len(peptides) > 0:
            newPeptides = []  # holds all peptides that match the verification qualities
            for peptideInfo in peptides:
                pepList = peptideInfo[0]
                stringPep = pepList[0]
                pepMassList = peptideInfo[1]
                pepMass = pepMassList[0]
                if pepMass in self.spectrum:
                    newPeptides.append(peptideInfo)
                    if pepMass == parentMass:  # mass verification
                        if self.getCyclicSpec(stringPep) == sortedSpectrum:  # contents verification
                            return stringPep
            peptides = self.addPeps(newPeptides)  # loop to get all possible variations of the peptide spectrum

    def getCyclicSpec(self, stringPep):
        """Build the mass spectrum for each possible test peptide to be used for comparison against input spectrum."""
        allSubsMass = [0]  # list to hold each mass for sub peptide
        for i in range(1, len(stringPep) + 1):
            for j in range(len(self.aaList)):
                aa = self.aaList[j]
                aaMass = self.aaMass[aa]  # gets the mass of aa
                if aa == stringPep[i-1]:
                    subPeptideMass = allSubsMass[i-1] + aaMass  # gets the mass of each sub peptide
                    allSubsMass.append(subPeptideMass)
        peptideMass = allSubsMass[len(stringPep)]
        cyclicSpectrum = [0]  # list to hold generated spectrum based on matching sub peptides
        for i in range(len(stringPep)):
            for j in range(i + 1, len(stringPep) + 1):
                massChange = allSubsMass[j] - allSubsMass[i]  # gets change in mass for each sub peptide pair (i,j)
                cyclicSpectrum.append(massChange)
                if i > 0 and j < len(stringPep):
                    subPeptideMass = peptideMass - massChange  # gets the mass of combined sub peptides for calculations
                    cyclicSpectrum.append(subPeptideMass)
        return sorted(cyclicSpectrum)  # sorted for comparative verification

    def addPeps(self, builtPeptides):
        """Build possible peptide sequences to be cycled for verification."""
        newPeptideList = []  # list to hold new peptides
        for peptideInfo in builtPeptides:
            pepList = peptideInfo[0]  # list to hold temporary possible peptides
            stringPep = pepList[0]  # list to hold temporary peptide string
            pepMassList = peptideInfo[1]  # list to hold temporary possible peptide masses
            pepMass = pepMassList[0]  # list to hold temporary peptide mass

            for aa in self.aaList:
                newPair = []  # list to hold new peptide
                newPep = []
                newStringPep = ''  # string to build peptide
                newStringPep += stringPep
                newStringPep += aa
                newPep.append(newStringPep)
                newMass = []  # list to hold current peptide mass
                newCalcMass = pepMass + self.aaMass[aa]
                newMass.append(newCalcMass)
                newPair.append(newPep)
                newPair.append(newMass)
                newPeptideList.append(newPair)
        return newPeptideList
    
    def getPeptideSeq(self, peptide):
        """Build list of possible peptide sequences based on the peptides built from spectrum match for final spec."""
        for i in range(len(peptide)):
            if i == 0:
                self.possiblePeptides.append(peptide)
                self.possiblePeptides.append(peptide[::-1])  # adds the reverse peptide
            else:
                startString = peptide[i:]
                endString = peptide[:i]
                newFullPep = startString + endString
                self.possiblePeptides.append(newFullPep)
                self.possiblePeptides.append(newFullPep[::-1])  # adds the reverse peptide
        return

    def formatMatchSpectrum(self):
        """Get the spectrum of the matched peptide string and prepare it for output formatting."""
        allPepMasses = []  # list to hold the mass of each substring peptide for each peptide
        for massSeq in self.possiblePeptides:
            cirSeqMass = []  # list to hold the mass  of each substring peptide for a peptide
            for aminoAcid in massSeq:
                mass = self.aaMass[aminoAcid]  # getting mass of the substring peptides
                cirSeqMass.append(mass)
            allPepMasses.append(cirSeqMass)
        massSets = []  # list to hold strings of each peptide spectrum calculated
        for i in allPepMasses:  # creating desired output formatting as a string
            massString = ''  # holds the output string for peptide spectrum
            for j in range(len(i)):
                if j < len(i) - 1:
                    massString += str(i[j]) + '-'
                else:
                    massString += str(i[j])
            massSets.append(massString)
        return massSets

def main():
    """
    Program to find a Cyclic Peptide with a Theoretical Spectrum Matching an Ideal Spectrum
    given A collection of (possibly repeated) integers Spectrum corresponding to an ideal experimental spectrum.
    """

    # contents = ['0 113 128 186 241 299 314 427']  # sample input
    contents = sys.stdin.read().splitlines()
    inputSpectrum = contents[0]
    inputSpectrum = inputSpectrum.split()
    spectrum = list(int(item) for item in inputSpectrum)  # gets the integer value for each spectrum string
    myPeptides = TheoPepCyclic(spectrum)
    TheoPeptide = myPeptides.matchTheoretical()
    myPeptides.getPeptideSeq(TheoPeptide)
    allTheoSpec = myPeptides.formatMatchSpectrum()
    print(*allTheoSpec, sep=' ')  # prints the formatted output on one line with spaces between each spectrum


if __name__ == '__main__':
    main()
