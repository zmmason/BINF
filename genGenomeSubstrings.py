# !/usr/bin/env python3
# Find Substrings of a Genome Encoding a Given Amino Acid String
# University of California, Santa Cruz - BME 205
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)
# Group Members: NONE
import itertools
import sys


class GenomeEncoding:
    """Return all substrings of Text encoding Peptide (if any such substrings exist) given a DNA and AA string."""

    def __init__(self, seq, peptide):
        """Build a constructor to hold original seq and AA along with crucial data."""
        self.seq = seq  # original DNA sequence
        self.peptide = peptide  # original peptide sequence
        self.allPepSeqs = []  # list to hold all possible nuc sequences based on the peptide sequence
        self.codonTable = {  # holds all amino acids and their associated codons
           'F': ['TTT', 'TTC'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
           'Y': ['TAT', 'TAC'], 'C': ['TGT', 'TGC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
           '-': ['TAA', 'TGA', 'TAG'], 'W': ['TGG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
           'H': ['CAT', 'CAC'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Q': ['CAA', 'CAG'],
           'I': ['ATT', 'ATC', 'ATA'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'N': ['AAT', 'AAC'],
           'K': ['AAA', 'AAG'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
           'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
           'E': ['GAA', 'GAG']
        }

    def getCodonSeqs(self):
        """Add all possible forward stranded peptides to list of all possible peptide sequences."""
        combinations = list(self.codonTable[aa] for aa in self.peptide)  # creates a list of possible codons based on AA
        self.allPepSeqs = list(''.join(codon) for codon in itertools.product(*combinations))  # creates list of peptides
        return

    def getRevCodonSeqs(self):
        """Add all possible compliment stranded peptides to list of all possible peptide sequences."""
        compDict = {'A': 't', 'T': 'a', 'G': 'c', 'C': 'g'}  # nuc compliments for reverse strand
        revPep = []  # list to hold the temporary reverse peptides before incorporation into the complete list
        for seq in self.allPepSeqs:
            revSeq = seq[::-1]  # reverses the strand to be prepped for nt compliments
            for nuc in compDict:
                revSeq = revSeq.replace(nuc, compDict[nuc])  # replaces nt's with their compliments
            revSeq = revSeq.upper()
            revPep.append(revSeq)
        for i in revPep:
            self.allPepSeqs.append(i)  # adds the reverse strand peptide to the list of possible peptide seqs
        return

    def printEncodePep(self):
        """Search the original DNA sequence for desired peptide and print the peptide sequence."""
        for i in range(len(self.seq)):
            if self.seq[i:i+(len(self.peptide*3))] in self.allPepSeqs:
                print(self.seq[i:i+(len(self.peptide*3))])  # print peptide if DNA region matches a possible peptide
        return


def main():
    """Return all substrings of Text encoding Peptide (if any such substrings exist) given a DNA and AA string."""

    # contents = ['ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA']  # sample input
    contents = []
    for line in sys.stdin:
        contents.append(line.strip())
    myPeptide = GenomeEncoding(contents[0], contents[1])
    myPeptide.getCodonSeqs()
    myPeptide.getRevCodonSeqs()
    myPeptide.printEncodePep()


if __name__ == '__main__':
    main()
