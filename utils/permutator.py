import pandas as pd
import numpy as np
from tqdm import tqdm

class DeletionPermutator():
    """


    """
    def __init__(self, df: pd.DataFrame, sn: int):
        self.df = df.query(f"serial_n == {sn}").reset_index(drop=True)
        
    @staticmethod    
    def convert_to_binary_with_leading_zeros(n: int, num_bits: int) -> str:
        binary_str = bin(n)[2:]
        binary_str = binary_str.zfill(num_bits)
        return np.array(list('1' + binary_str)).astype(int)

    def permute_sequence(self):
        permutated_dict = {}
        num_bits = len(self.df)-1
        n_comb = 2**(len(self.df)-1)-1
        area = self.df["first_area"]
        seq = np.array(list(self.df["AA Name"]))
        inverse_area = 1-area
        for i in tqdm(range(n_comb)):
            n = n_comb-i
            binary = self.convert_to_binary_with_leading_zeros(n, num_bits)
            complement = abs(binary-1)
            complement[0] = 1
            prob = area*binary + inverse_area*complement
            prob = prob.prod()
            seq_i = np.argwhere(binary==1).flatten()
            comb = "".join(seq[seq_i])
            permutated_dict[comb] = float(prob)
        permute_array = np.array(list(permutated_dict.items()))
        permute_array_norm = permute_array[:,1].astype(float) / max(permute_array[:,1].astype(float))
        permute_array = np.concatenate((permute_array, permute_array_norm.reshape(-1,1)), axis=1)
        permute_array = permute_array[np.argsort(permute_array[:,1].astype(float))]  
        permute_array = np.flip(permute_array,axis=0)
        permute_df = pd.DataFrame(permute_array, columns=["Sequence", "Probability", "Probability Normalised"])
        return permute_df