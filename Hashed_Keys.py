import random
import numpy as np

import os
import sys



def readDatabase():
    """read SMILES from SMilES coconut database"""
    with open("COCONUT_DB.smi", "r") as db:
        lines = db.readlines()
        
    coconut_smiles = []
    for line in lines:
        SMILES, ID = line.split(" ")
        coconut_smiles.append(SMILES)

    print("database read!")
    print("database length: ", len(coconut_smiles))
    return coconut_smiles


def databaseBitStrings(coconut_smiles):
    """calculate bit strings for every coconut smiles entry"""
    coconut_bitStrings = []
    for SMILES in coconut_smiles:
        h = HashedKey(SMILES)

        uniqueSubstring = h.createUniqueSubstrings()
        hashKeys = h.generateHashKeys()
        bit_positions = h.pickRandomBit()
        bitString = h.manipulateBitString()

        coconut_bitStrings.append(bitString)
        
    return coconut_bitStrings
        

class HashedKey(object):
    def __init__(self,SMILES):
        """reads a SMILES and initializes a bit string."""
        self.SMILES = SMILES

        self.bitString128 = np.zeros(1024)

    def createUniqueSubstrings(self):
        """creates a unique subset of symbols based on a given SMILES."""
        # only up to 8 characters should be implemented as well
        self.uniqueSubstrings = []
        queue = []
        
        init_index = 0
        running_index = 1
        
        queue.append(self.SMILES[init_index:running_index])
    
        while queue:
            string = queue.pop()

            if string not in self.uniqueSubstrings:
                self.uniqueSubstrings.append(string)


            if (running_index - init_index) > 10: # number as input
                init_index += 1
                running_index = init_index + 1
                queue.append(self.SMILES[init_index:running_index])
            elif running_index < len(self.SMILES):
                running_index += 1
                queue.append(self.SMILES[init_index:running_index])
            elif init_index < len(self.SMILES)-1:
                init_index += 1
                running_index = init_index + 1
                queue.append(self.SMILES[init_index:running_index])

        return self.uniqueSubstrings
                        

    def generateHashKeys(self):
        """generates a hash of a given list of substrings"""
        self.hashedSubstrings = []
        for string in self.uniqueSubstrings:
            hashed_string = hash(string)
            self.hashedSubstrings.append(hashed_string)

        return self.hashedSubstrings


    def pickRandomBit(self):
        """chooses a random bit position on the bit string"""
        self.bit_positions = []
        for hashed_string in self.hashedSubstrings:
            random.seed(hashed_string)
            value = random.randrange(0, 1024)
            self.bit_positions.append(value)

        return self.bit_positions
        

    def manipulateBitString(self):
        """manipulates the 1024 bit string in XOR style"""
        for value in self.bit_positions:
            self.bitString128[value] = 1

        return self.bitString128

    
                
    def bitStringComparison(self, coconut_bitStrings): 
        """calculates the on bits appearances on bit strings from the coconut bits strings and a given bit string"""
        
        self.on_bits = []
        
        for other_bitString128 in coconut_bitStrings:
            on_both = 0
            on_target = 0
            on_query = 0
            for target, query in zip(other_bitString128, self.bitString128):
                if target == 1 and query == 1:
                    on_both += 1
                elif target == 1:
                    on_target += 1
                elif query == 1:
                    on_query += 1

            self.on_bits.append( (on_both, on_target, on_query) )

        return self.on_bits


    def calcTanimotoSimilarity(self):
        """calculates the tanimoto similarity between bit strings from coconut and the given input string"""
        self.tanimoto_on_coconut = []
        
        for on_both, on_target, on_query in self.on_bits:
            tanimoto = on_both / (on_target + on_query + on_both)
            self.tanimoto_on_coconut.append(tanimoto)

        return self.tanimoto_on_coconut

    
    def calcTverskySimilarity(self, alpha=1, beta=1):
         """calculates the tversky similarity between bit strings from coconut and the given input string
         Alpha and Beta are non negative values/weights"""
         self.tversky_on_coconut = []
        
         for on_both, on_target, on_query in self.on_bits:
             tversky = on_both / (alpha * (on_target - on_both) + beta * (on_query - on_both) + on_both)
             self.tversky_on_coconut.append(tversky)
             
         return self.tversky_on_coconut

    
    def calcDiceSimilarity(self):
         """calculates the dice similarity between bit strings from coconut and the given input string"""
         self.dice_on_coconut = []
         
         for on_both, on_target, on_query in self.on_bits:
             dice = 2 * on_both / (on_target + on_query)
             self.dice_on_coconut.append(dice)

         return self.dice_on_coconut


    
    def similarityThreshold(self, threshold, coconut_smiles, metric_on_coconut):
        """filters molecules with a certain threshold and return similarity and index"""
        above_threshold_entries = []
        for index, similarity in enumerate(metric_on_coconut):
            if similarity >= threshold:
                above_threshold_entries.append( (index, coconut_smiles[index], similarity) )

        return above_threshold_entries
                
    

if __name__ == "__main__":
    hashseed = os.getenv('PYTHONHASHSEED')
    if not hashseed:
        os.environ['PYTHONHASHSEED'] = '73'
        os.execv(sys.executable, ['python3'] + sys.argv)

    coconut_smiles = readDatabase()
    coconut_bitStrings = databaseBitStrings(coconut_smiles[:500])
    
    #SMILES from Penicillin G
    h = HashedKey("CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C")
    
    uniqueSubstring = h.createUniqueSubstrings()
    hashKeys = h.generateHashKeys()
    bit_positions = h.pickRandomBit()
    bitString = h.manipulateBitString()
    on_bits = h.bitStringComparison(coconut_bitStrings)
    
    tanimotoSimilarity = h.calcTanimotoSimilarity()
    #print(sorted(tanimotoSimilarity)[-2])
    tverskySimilarity = h.calcTverskySimilarity(alpha=1, beta=1)
    #print(sorted(tverskySimilarity)[-2])
    diceSimilarity = h.calcDiceSimilarity()
    #print(sorted(diceSimilarity)[-2])
    above_threshold_entries = h.similarityThreshold(threshold=0.89, coconut_smiles=coconut_smiles, metric_on_coconut=diceSimilarity)
    print(above_threshold_entries)
    
    # things to do:
    # - save all coconut bit strings to a file (takes several hours, test it with a small batch np.save("name.npy", name, allow_pickle=True)
    # - dice similarity, tversky similarity ... 3.4.1.1 Molecular Descriptors and Fingerprints (Machine Learning in Chemistry)
    # - siamese neural network for similarity check ??!!??!!
    # - similarity check for bitstring matches, how to calculate the similarity? RDKit, better by yourself
    #DONE - String size should be 8 !!!!!!!!!!!!!!!!!!!!!!!!
    #-> now you have less matchings so you should do folding
    #--> automated system when folding is neccesarry
    # - is 1024 the right size for the bitstring???
    # - graphical evaluation:
    #     - compare different similarity metrics (histogram)
    #     - compare bit similarity to SMILES similarity
    # - WHAT ELSE???
