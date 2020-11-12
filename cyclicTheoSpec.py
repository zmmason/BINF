# !/usr/bin/env python3
# Generate the Theoretical Spectrum of a Cyclic Peptide
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE

import sys
from itertools import islice, cycle


class TheoreticalSpectrum:
    """
    Given an amino acid string, return the collection of all of the masses of its subpeptides
    in addition to the mass 0 and the mass of the entire peptide.
    """

    def __init__(self, peptide):
        """Build consructor to hold peptide sequence and the integer mass table for reference."""
        self.peptide = peptide
        self.aaMass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101,  # holds the integer mass table for AA
                       'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'W': 186,
                       'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163}

    def cyclospectrum(self):
        """Build a theoretical cyclospectrum list to hold is the collection of masses (zero->subpeptide->complete)."""
        cycSpec = [0]  # creates and starts the cyclospectrum list at zero
        aaMasses = list(self.aaMass[aa] for aa in self.peptide)  # creates a list of masses fr each aa in peptide
        for i in range(1, len(aaMasses), 1):  # getting all possible sub-peptides and their masses
            for j in range(0, len(aaMasses), 1):
                aaGroup = islice(cycle(aaMasses), j, j + i)  # gets each sub-peptide slice (ordered combo)
                cycSpec.append(sum(aaGroup))
        cycSpec.append(sum(aaMasses))
        return cycSpec


def main():
    """
    Given an amino acid string, return the collection of all of the masses of its subpeptides
    in addition to the mass 0 and the mass of the entire peptide.
    """
    # contents = ['LEQN']  # sample input
    contents = []
    for line in sys.stdin:
        contents.append(line.strip())
    theo = TheoreticalSpectrum(contents[0])
    spectrum = theo.cyclospectrum()
    spectrum.sort()  # sorts the spectrum list in numerical order
    print(*spectrum, end=' ')  # prints the spectrum list on one line without commas


if __name__ == '__main__':
    main()
