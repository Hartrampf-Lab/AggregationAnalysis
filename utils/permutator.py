###########################################################################
#
#
#   AFPS Deletion Permutator Program  - permutator.py
#   A robust data analytical method to investigate sequence dependence in flow-based peptide synthesis.
#   Bálint Tamás, Pietro Luigi Willi, Héloïse Bürgisser, Nina Hartrampf
#   DOI: 10.1039/D3RE00494E
#   UZH, Institute of Chemistry.
#
#
###########################################################################

import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import rdkit.Chem as Chem
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

class DeletionPermutator():
    """
    The DeletionPermutator class is used to generate all possible deletions 
    and calculate the respective probability of each deletion. This tool is 
    particularly useful for determining the side products that are likely to
    by found in the LC-MS spectra

    """
    def __init__(self, df: pd.DataFrame, sn: int):
        self.df = df.query(f"serial_n == {sn}").reset_index(drop=True)
        
    @staticmethod    
    def convert_to_binary_with_leading_zeros(n: int, num_bits: int):
        """
        The convert_to_binary_with_leading_zeros converts an integer 
        to binary with leading zeros and a one at the end. the binary 
        string is then converted to a numpy array of integers 1 and 0.
        """
        binary_str = bin(n)[2:]
        binary_str = binary_str.zfill(num_bits) # add leading zeros
        return np.array(list('1' + binary_str)).astype(int) # add 1 at the end and convert to array

    def permute_sequence(self):
        """
        The permute_sequence method uses the property of binary numbers
        to generate all possible combinations of deletions. The method
        then calculates the probability of each combination by using the
        same binary number that was used to generate the combination.
        The probability is calculated by multiplying the areas with the 
        binary array and add this to the inverse areas multiplied with 
        the complement of the binary array. The probability is then
        calculated by multiplying all the probabilities of the amino acids
        in the sequence. The method returns a dataframe with the permutated 
        sequence and the probability of the sequence.
        """
        permutated_dict = {}
        num_bits = len(self.df)-1
        n_combinations = 2**(len(self.df)-1)-1
        area = self.df["first_area"]
        seq = np.array(list(self.df["AA Name"]))
        inverse_area = 1-area
        print("Calculating permutation probabilities...")
        for i in tqdm(range(n_combinations)):
            n = n_combinations-i
            binary = self.convert_to_binary_with_leading_zeros(n, num_bits) # convert index to binary
            complement = abs(binary-1) # calculate complement of binary array
            complement[0] = 1 # set first element to 1
            prob = area*binary + inverse_area*complement # calculate probability for each amino acid
            prob = prob.prod() # calculate probability of sequence
            seq_i = np.argwhere(binary==1).flatten()
            comb = "".join(seq[seq_i])
            permutated_dict[comb] = float(prob)
        permute_array = np.array(list(permutated_dict.items()))
        permute_array_norm = permute_array[:,1].astype(float) / max(permute_array[:,1].astype(float)) #normalise probabilities to max probability
        permute_array = np.concatenate((permute_array, permute_array_norm.reshape(-1,1)), axis=1)
        permute_array = permute_array[np.argsort(permute_array[:,1].astype(float))]  
        permute_array = np.flip(permute_array,axis=0)
        self.permute_df = pd.DataFrame(permute_array, columns=["Sequence", "Probability", "Probability Normalised"])
        self.permute_df = self.permute_df.drop(self.permute_df.query(f"Sequence=='{''.join(self.df['AA Name'])}'").index) 
        self.permute_df["Rank"] = list(range(1,len(self.permute_df)+1))
        self.permute_df = self.permute_df.reset_index(drop=True)
        return self.permute_df
    
    def compute_mass(self):
        """
        The compute_mass method calculates the mass of each deletion permutation.
        The mass is calculated by using the rdkit function CalcExactMolWt.
        """
        print("Calculating masses...")
        for seq in tqdm(self.permute_df["Sequence"]):
            mol = Chem.MolFromFASTA(seq)
            self.permute_df.loc[self.permute_df["Sequence"]==seq, "Mass"] = rdMolDescriptors.CalcExactMolWt(mol)-1
        return self.permute_df
