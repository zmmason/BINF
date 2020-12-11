# !/usr/bin/env python3
# Bioinformatics toolbox classes 
# University of California, Santa Cruz
# Biomolecular Engineering and Bioinformatics
# Name: Zachary Mason (zmmason)


# NOTE: THIS PROGRAM CONTAINS MULTIPLE CLASSES THAT ARE TO BE USED SEPARATELY OR IN CONJUNCTION WITH ONE ANOTHER.
#       THESE CLASSES ARE DESIGNED TO BE IMPORTED FROM THIS BINF TOOLBOX FOR UTILIZATION

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


class OrfFinder:
    """
    after 'findORFs' parses a fastA file using 'fastAreader,' this program searches open reading frames
    returning putative genes determined and defined by the location of the start and stop codons. The program will
    output: ORF positions (1,2,3,-1,-2,-3), gene start location, gene stop location, and gene length, ordered by
    decreasing gene length:
    - program uses an optional set of start codons - ATG, GTG, TTG
    - program uses an optional set of stop codons - TAG, TGA, TAA
    - returns an ordered list of every putative gene (>100nt) in an ORF rather than only the largest.

    """

    def __init__(self):
        """
        initializing start/stop codons as well as the list that will hold all ORF/data information
        """
        self.startCodon = ['ATG', 'GTG', 'TTG']  # initializing start codons
        self.stopCodon = ['TAG', 'TAA', 'TGA']  # initializing stop codons
        self.allGenes = []  # initializing list for ORF/gene data

    def ORFgen(self, seq, compSeq):
        """
        object that takes incomming sequences and searches their ORF's for specific data related
        to start/stop locations. those start and stop locations are used to determine gene length.
        """
        startIndex = []
        stopIndex = []
        revStartIndex = []
        revStopIndex = []
        for i in range(0, 3):  # setting up ORFs for the forward strand (+1,+2,+3)
            """
            FORWARD STRAND calculations and data:
            This block of code is the same for both the forward and  compliment strands
            **** The comments for the forward strand will pertain to both forward and compliment strands****
            """
            orfPosition = i + 1
            for position in range(i, len(seq[i:]), 3):  # takes position every 3 nts
                codon = seq[position:position + 3]  # creates a codon (str) every 3 nts
                if codon in self.startCodon:  # adds all start codon occurrences to one list
                    startIndex.append(position + 1)
                if codon in self.stopCodon:  # adds all stop codon occurrences to one list
                    stopIndex.append(position + 3)
            for i in range(len(stopIndex)):
                for j in range(len(startIndex)):
                    if startIndex[j] < stopIndex[i]:
                        geneLen = (stopIndex[i] - startIndex[j] + 1)  # gets length of the gene
                        self.allGenes.append([orfPosition, startIndex[j], stopIndex[i], geneLen])
                        startIndex = [elem for elem in startIndex if elem > stopIndex[
                            i]]  # updates 'startIndex' by removing all starts up to the used stop
                    break
            startIndex.clear()
            stopIndex.clear()

        for i in range(0, 3):
            revOrfPosition = -1 * (i + 1)
            for position in range(i, len(compSeq[i:]), 3):
                codon = compSeq[position:position + 3]
                if codon in self.startCodon:
                    revStartIndex.append(position + 1)
                if codon in self.stopCodon:
                    revStopIndex.append(position + 3)

            for i in range(len(revStopIndex)):
                for j in range(len(revStartIndex)):
                    if revStartIndex[j] < revStopIndex[i]:
                        revGeneLen = (revStopIndex[i] - revStartIndex[j] + 1)
                        self.allGenes.append([revOrfPosition, revStartIndex[j], revStopIndex[i], revGeneLen])
                        revStartIndex = [elem for elem in revStartIndex if elem > revStopIndex[i]]
                    break
            revStartIndex.clear()
            revStopIndex.clear()

        self.allGenes.sort(key=lambda x: (x[3], x[0]),
                           reverse=True)  # sorting by increasing length, then by left position of gene
        return self.allGenes


class ProteinParam:
    """
    Program to calculate the physical-chemical properties of a protein sequence.

    INPUT:  - a protein sequence (example: VLSPADKTNVKAAW)

    OUTPUT: - number of amino acids,
            - total molecular weight,
            - Molar extinction coefficient,
            - Mass extinction coefficient,
            - theoretical isoelectric point (pI),
            - amino acid composition
    """

    def __init__(self, protein):
        """ 
        - initializes and sets up the input string to be manipulated by methods
        - sets up the amino acid composition dictionary (aaComp) to be used by methods

        """
        self.aa2mw = {  # dictionary containing molecular weights for each amino acid key
            'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
            'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
            'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
            'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
        self.mwH2O = 18.015  # water released with peptide bond formation
        self.aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}
        self.aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
        self.aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        self.aaNterm = 9.69
        self.aaCterm = 2.34
        protein = protein.upper()  # converst user input to uppercase
        protein = protein.replace(" ", "")  # removes spaces from user input
        self.protein = protein

        self.aaComp = {}  # creates the dictionary
        for aa in ProteinParam.aa2mw:  # for loop to cycle and count the amino acids in the input string against the 'aa2mw' dictionary
            myAAcomp = dict.fromkeys(aa, self.protein.count(aa))
            self.aaComp.update(myAAcomp)  # adds the aa and its comp too the dictionary

    def aaCount(self):
        """ method that counts and returns the number of characters (length) in the input string"""
        self.protein
        return len(self.protein)  # returns the length of the protein sequence when method is called

    def pI(self):
        """ method to estimate the theoretical isolelectric point using the particular pH that yields a neutral net Charge"""
        self.aa2chargePos.update({'aaNterm': 9.69})
        self.aaComp.update({'aaNterm': 1})
        self.aa2chargeNeg.update({'aaCterm': 2.34})
        self.aaComp.update({'aaCterm': 1})

        for pH in range(0, 1400 + 1):
            """ for loop that itterates throught pH's to find the pH that produces a net charge that is approx. 0.00"""
            pH1 = pH / 100
            netPos = sum(
                self.aaComp[aa] * (10 ** self.aa2chargePos[aa]) / (10 ** self.aa2chargePos[aa] + (10 ** pH1)) for aa in
                self.aa2chargePos)
            netNeg = sum(self.aaComp[aa] * ((10 ** pH1) / (10 ** self.aa2chargeNeg[aa] + (10 ** pH1))) for aa in
                         self.aa2chargeNeg)
            netCharge = 100 * (
                        netPos - netNeg)  # this makes the charges larger than decimals to be grabbed without extra conversions
            if 0 < netCharge < 1:  # grabbing the pH that produces a charge closest to zero
                self.aaComp.pop('aaNterm')  # deletes term from the dictionary
                self.aaComp.pop('aaCterm')
                return pH1

    def aaComposition(self):
        """method that returns the dictionary created in the '__init__' method; Shows the composition of each aa in the string"""
        return self.aaComp

    def molarExtinction(self):
        """ method calculates the extinction coefficient which indicates the amount of light a protein absorbs at a certain wavelength"""
        extinction = sum(self.aaComp[aa] * ProteinParam.aa2abs280[aa] for aa in ProteinParam.aa2abs280)
        return extinction

    def massExtinction(self):
        """ method to calculate the Mass Extinction coefficient from the Molar Extinction coefficient"""
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        """ This method calculates the molecular weight (MW) of the protein sequence. 
            Uses the protein calculated composition and sums the the weights of the individual Amino acids 
            (excluding the waters that are released with peptide bond formation). """

        molecWeight = ProteinParam.mwH2O + sum(
            self.aaComp[aa] * (ProteinParam.aa2mw[aa] - ProteinParam.mwH2O) for aa in ProteinParam.aa2mw)
        return molecWeight
    
    
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
        codonList = [(inSeq[i:i + 3]) for i in
                          range(0, len(inSeq), 3)]  # splits the string sequence every 3 nuc's and puts into list
        for codon in codonList:
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
        return self.aaComp

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

