"""
Created on Sun April 9 12:35:31 2023
@author: Murat Cem Köse
"""
import numpy as np
import random
import pandas as pd

######## Common functions

def transform_Huffman_base3(Huffman_codec,string,to_b3=True):
    """Function to convert strings to Huffman base 3 codes and vica versa.  
    Parameters
    ----------
    Huffman_codec : DataFrame
        Huffman_codec.
    string : String
        The sequence to be converted.
    to_b3 : Bool
        The direction of the conversion. True from string to b3 and False from b3 to string.
    """
    if to_b3 == True:
        base3 = []
        [base3.append(Huffman_codec.loc[i,3]) for i in string];
        return base3
    else:
        word = ""
        i_start = 0
        while i_start+4 < len(string):
            try:
                b3 = Huffman_codec.iloc[np.where(Huffman_codec[3] == string[i_start:i_start+6])[0][0]].name
                word=word+b3
                i_start = i_start+6
            except:
                try:
                    b3 = Huffman_codec.iloc[np.where(Huffman_codec[3] == string[i_start:i_start+5])[0][0]].name
                    word=word+b3
                    i_start = i_start+5
                except:
                    raise Exception("Base 3 can not be found in the Huffman codec.")
        return word

def base3_DNA_transformation(DNA_convertion_table,string,direction=True,prior=None):
    """Function to convert base 3 sequence to DNA sequence vica versa.  
    Parameters
    ----------
    DNA_convertion_table : DataFrame
        DNA_convertion_table.
    string : String
        The sequence to be converted.
    direction : Bool
        The direction of the conversion. True from b3 to seq and False from seq to b3.
    prior: String
        The prior sequence to start with.
    """
    if direction == True:
        DNA = ""
        if prior == None:
            DNA = DNA+DNA_convertion_table.loc["A",string[0]]
        else:
            DNA = DNA+DNA_convertion_table.loc[prior,string[0]]   ## The prior can be checked if it is a single character.
        for i in range(1,len(string)):
                DNA = DNA+DNA_convertion_table.loc[DNA[i-1],string[i]]
        return DNA
    else:
        base3 = ""
        for i in range(len(string)-2,-1,-1):
            base3 = str(np.where(DNA_convertion_table.loc[string[i]] == string[i+1])[0][0])+base3
        if prior == None:
            base3 = str(np.where(DNA_convertion_table.loc["A"] == string[0])[0][0])+base3
        else:
            base3 = str(np.where(DNA_convertion_table.loc[prior] == string[0])[0][0])+base3
        return base3
    
def complement_DNA(DNA):
    conversion_dict = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join([conversion_dict.get(i) for i in DNA])

def reverse_completement_DNA(DNA):
    return complement_DNA(DNA)[::-1]

def encrypt_decript(DNA, keystream, add):
    """Function to perform keystream encryption or decryption.  
    Parameters
    ----------
    DNA : String
        DNA sequence to be work with.
    keystream : Integer
        The choice of the keystream.
    add : Bool
        The direction of the conversion. True for encryption and False for decryption.
    """
    bases = []
    kstrits = []

    DNA = DNA.lower().translate(str.maketrans('acgt', '0123'))

    bases = [int(base) for base in DNA]
    kstrits = [int(trit) for trit in keystream]

    L = len(keystream)
    if len(DNA) != L:
        raise ValueError("encrypt error")

    # compute differences
    for i in range(L-1, 0, -1):
        bases[i] = (bases[i] - bases[i - 1]) % 4

    # subtract 1 to recover base-3 trits
    for i in range(1, L):
        bases[i] = (bases[i] - 1) % 3

    # add/subtract keystream trits as appropriate
    for i in range(1, L):
        bases[i] = (bases[i] + kstrits[i]) % 3 if add else (bases[i] - kstrits[i]) % 3

    # add 1 to generate base-4 differences
    for i in range(1, L):
        bases[i] = (bases[i] + 1) % 4

    # add successively to generate 0,1,2,3 coding
    for i in range(1, L):
        bases[i] = (bases[i] + bases[i - 1]) % 4

    DNA = ''.join(str(base) for base in bases)
    DNA = DNA.translate(str.maketrans('0123', 'acgt'))

    # return result
    return DNA.upper()

######## Encryption functions
    
def get_S1(base3):
    return "".join(base3)

def get_S1_base(S1):
    return np.base_repr(len(S1),base=3)

def get_S2(S1_base3):
    return "0"*(25-len(S1_base3))+ S1_base3

def get_S3(S1,S2):
    limit = 0
    while (limit+len(S1+S2))%25!=0:
        limit = limit + 1
    return "0"*(limit+len(S1+S2)-len(S1+S2))

def get_IX(ID,n_F):
    base3_n = np.base_repr(n_F,base=3)
    i3 = "0"*(12-len(str(base3_n)))+str(base3_n)
    P = sum([int(i) for i in i3[0::2]])+sum([int(i) for i in ID[0::2]])
    P = P%3
    return ID+i3+str(P)

def get_index_border_begining(string):
    if string[0] == "A":
        return "T"
    elif string[0] == "T":
        return "A"
    else:
        return str(random.randint(0, 1)).replace("0","A").replace("1","T")

def get_index_border_end(string):
    if string[-1] == "C":
        return "G"
    elif string[-1] == "G":
        return "C"
    else:
        return str(random.randint(0, 1)).replace("0","C").replace("1","G")
    
def add_border_NTs(string):
    return get_index_border_begining(string)+string+get_index_border_end(string)

######## Decryption functions

def remove_border_nt(DNA):
    return DNA[1:-1]

def partition_F(F):
    F_IX = remove_border_nt(F)
    IX_DNA = F_IX[-15:]
    F_randomized = F_IX[:100]
    return IX_DNA,F_randomized

def get_ID_and_n_F(IX):
    return IX[:2],int(IX[2:-1])

def revert_S4(S4):
    S2 = S4[0:25]
    len_S1 = int(S2,3)
    S1 = S4[25:25+len_S1]
    S3 = S4[25+len_S1:]
    return S1, S2, S3