# !/usr/bin/env python3
# University of California, Santa Cruz - BME 160
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)

import sys

class FastAreader:
    """
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    """

    def __init__(self, fname=''):
        """contructor: saves attribute fname """
        self.fname = fname

    def doOpen(self):
        """ Handle file opens, allowing STDIN."""
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        """ Read an entire FastA record and return the sequence header/sequence"""
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class NucParams:
    """Program to calculate the physical-chemical properties of a protein sequence."""
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inString=''):
        """
        initializing atributes of the dna/rna sequence provided from the fastAreader.
        - created dictionaries: nucleotide composition, amino acid composition, codon composition
        - all dictionaries hold a zero value to be updated when 'addSequence' adds to them
        """

        self.nucComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0}

        self.codonComp = {}  # creates dict to hold codons from the rna & dna codon tables with '0' set values
        for i in range(len(NucParams.rnaCodonTable)):
            newCodon = (NucParams.rnaCodonTable.keys() or NucParams.dnaCodonTable.keys())
            newCodon = dict.fromkeys(newCodon, 0)
            self.codonComp.update(newCodon)

        self.aaComp = {  # empty dict for aa
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
        }

        self.addSequence(inString)  # gives the sequences to 'addSequence'

    def addSequence(self, inSeq):
        """
         Method accepts new sequences in the form of nucleotides (ACGTUN), and is presumed to start in frame 1.
         This data is added to the dictionary data in the '__init__' method
        """
        for nuc in inSeq:  # for loop to cyle and count individual nucs in the input string to add to 'self.nucComp'"""
            if nuc in self.nucComp:  # count valid nuc
                self.nucComp[nuc] += 1  # adds 1 to dict if found in sequence

        invalBase = 'N'  # invalid base to search and remove codon that contains it
        self.codonList = [(inSeq[i:i + 3]) for i in
                          range(0, len(inSeq), 3)]  # splits the string sequence every 3 nuc's and puts into list
        for codon in self.codonList:
            if invalBase not in codon:  # removing codon with invalid base
                if codon in self.codonComp:  # count valid codon
                    self.codonComp[codon] += 1  # adds 1 to dict if found in sequence

        aaList = []  # list to hold aa's from rna & dna dicts
        for aa in NucParams.rnaCodonTable or NucParams.dnaCodonTable:
            aaCount = NucParams.rnaCodonTable[aa] or NucParams.dnaCodonTable[aa], self.codonComp[
                aa]  # sets format '(aa, 0)'
            aaList.append(aaCount)  # adds each aa to list
        for aa, count in aaList:
            total = self.aaComp.get(aa, 0) + count  # uses dict of aa's and adds for each occurrance
            self.aaComp[aa] = total  # adds total of each aa to dict

    def aaComposition(self):
        """
        Method will return a dictionary of counts over the 20 amino acids and stop codons
        in the RNA/DNA codon tables.
        """
        return (self.aaComp)

    def nucComposition(self):
        """
        Method returns a dictionary of counts of valid nucleotides (ACGTNU) found in
        the analysis.
        """
        return (self.nucComp)

    def codonComposition(self):
        """
        This returns a dict that counts # of specific codons (ACGTU). Codons that contain 'N' bases are
        discarded in this method. All codon counts are stored as RNA codons.
        """
        return (self.codonComp)

    def nucCount(self):
        """
        This returns an integer value which is the sum of every valid nucleotide (ACGTUN) found.
        """
        return sum(self.nucComp.values())


def main ():
    """
    program takes a fatsa file translated by the 'FastAreader' and runs through 'nucParams'
    in 'sequenceAnalysis.py'. It is able to then analyze genetic sequences (DNA/RNA) to produce &
    print computational data for that fasta file relating to: amino acids, codons, nucleotides.
    """
    myReader = FastAreader()  # make sure to change this to use stdin if using a filename for testing
    myNuc = NucParams()
    for head, inSeq in myReader.readFasta():
        myNuc.addSequence(inSeq)
    gcCount = (myNuc.nucComposition().get('G')+myNuc.nucComposition().get('C'))
    print("Sequence length = {:.2f} Mb\n" .format(myNuc.nucCount()/1e6))  # prints seq length
    print("GC content = {:.1f}%\n" .format(gcCount/myNuc.nucCount()*100))  # prints gc content in seq length
        
    nucs = sorted(myNuc.rnaCodonTable.items(), key=lambda aa: aa[1])  # sort codons in alpha order, by Amino Acid
    
    for nuc, aa in nucs:  # calculate relative codon usage for each codon and print
        val = myNuc.codonComp[nuc]/(sum(myNuc.codonComp.values()))  # ratio of specific codon to total codons
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val*100, myNuc.codonComp[nuc]))
        

if __name__ == "__main__":
    main()







