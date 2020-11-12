# !/usr/bin/env python3
# University of California, Santa Cruz - BME 160
# Biomolecular Engineering and Bioinformatics
# Name: (zmmason)

import sys


class ProteinParam:
    """
    Program to calculate the physical-chemical properties of a protein sequence.
    INPUT:  - a protein sequence
    OUTPUT: - number of amino acids,
            - total molecular weight,
            - Molar extinction coefficient,
            - Mass extinction coefficient,
            - theoretical isoelectric point (pI),
            - amino acid composition
    """

    def __init__(self, protein, ):
        """ 
        - initializes and sets up the input string to be manipulated by methods
        - sets up the amino acid composition dictionary (aaComp) to be used by methods
        """
        protein = protein.upper()  # converts user input to uppercase
        protein = protein.replace(" ", "")  # removes spaces from user input
        self.protein = protein
        self.aaComp = {}  # creates dict to hold aa composition
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
        for aa in self.aa2mw:  # for loop to cycle and count the amino acids in the input string
            myAAcomp = dict.fromkeys(aa, self.protein.count(aa))
            self.aaComp.update(myAAcomp)

    def aaCount(self):
        """ method that counts and returns the number of characters (length) in the input string"""
        return len(self.protein)  # returns the length of the protein sequence when method is called

    def pI(self):
        """ Estimate the theoretical isolelectric point using the particular pH that yields a neutral net Charge """

        self.aa2chargePos.update({'aaNterm': self.aaNterm})
        self.aaComp.update({'aaNterm': 1})
        self.aa2chargeNeg.update({'aaCterm': self.aaCterm})
        self.aaComp.update({'aaCterm': 1})

        for pH in range(0, 1400 + 1):
            """ Iterate through pH's to find the pH that produces a net charge that is approx. 0.00."""
            pH1 = pH / 100
            netPos = sum(
                self.aaComp[aa] * (10 ** self.aa2chargePos[aa] / 10 ** self.aa2chargePos[aa] + (10 ** pH1)) for aa in
                self.aa2chargePos)
            netNeg = sum(self.aaComp[aa] * ((10 ** pH1) / (10 ** self.aa2chargeNeg[aa] + (10 ** pH1))) for aa in
                         self.aa2chargeNeg)
            netCharge = 100 * (netPos - netNeg)  # this makes the charges larger than decimals to be grabbed
            if 0 < netCharge < 1:  # grabbing the pH that produces a charge closest to zero
                self.aaComp.pop('aaNterm')  # deletes term from the dictionary
                self.aaComp.pop('aaCterm')
                return pH1

    def aaComposition(self):
        """Return the composition of each aa in the string."""
        return self.aaComp

    def _charge_(self):
        """Calculate the net charge on the protein at a specific pH (specified as a parameter of this method).
           The method is used by the pI method. """
        netPos = sum(
            self.aaComp[aa] * (10 ** self.a2chargePos[aa] / 10 ** self.aa2chargePos[aa] + (10 ** pH1)) for aa in
            self.aa2chargePos)
        netNeg = sum(
            self.aaComp[aa] * ((10 ** pH1) / (10 ** self.aa2chargeNeg[aa] + (10 ** pH1))) for aa in self.aa2chargeNeg)
        netCharge = 100 * (netPos - netNeg)  # used for calculations in pI()
        return netCharge

    def molarExtinction(self):
        """ Calculate the extinction coefficient specific wavelength"""
        extinction = sum(self.aaComp[aa] * self.aa2abs280[aa] for aa in self.aa2abs280)
        return extinction

    def massExtinction(self):
        """ method to calculate the Mass Extinction coefficient from the Molar Extinction coefficient"""
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        """ This method calculates the molecular weight (MW) of the protein sequence. 
            Uses the protein calculated composition and sums the the weights of the individual Amino acids 
            (excluding the waters that are released with peptide bond formation). """

        molecWeight = self.mwH2O + sum(self.aaComp[aa] * (self.aa2mw[aa] - self.mwH2O) for aa in self.aa2mw)
        return molecWeight


def main():
    """ Main method asks for user input (protein sequence) and call the appropriate methods to return the desired output.
        This method prints all outputs and repeats the process (asks for new input) untill the program termination key 
        (Ctrl-D) is hit to end the process.
    """

    inString = input('protein sequence?')  # user input
    while inString:  # while loop to make sure that all desired methods are called
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())  # creates a list of the keys found in the aaComp dictionary
        keys.sort()  # sorts the keys alphabetically
        if myAAnumber == 0:
            myAAnumber = 1  # handles the case where no AA are present
        for key in keys:  # divides the amount of aa in the string by the total aa's and prints the new LIST
            print("\t{} = {:.2%}".format(key, myAAcomposition[key] / myAAnumber))
        inString = input('protein sequence?')


if __name__ == "__main__":
    main()

# In[ ]:
