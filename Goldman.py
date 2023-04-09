"""
Created on Sun April 9 12:35:31 2023
@author: Murat Cem Köse
"""
import numpy as np
import pandas as pd
import random
import utils


class GoldmanObject:
    def __init__(self,sequence):
        """Contructor function for Object class.
        Parameters
        ----------
        sequence : String
            String variable to be be encrypted or decrypted.
        """
        self.sequence = sequence
        self.Huffman_codec = self._read_Huffman_codec()
        self.DNA_convertion_table = self._generate_DNA_conv_table()
        self.keystreams = ["0002000010110102111112122210011122221010102121222022000221201020221002121121000212222021211121122221",
"0202020122121210120001200210222112020222022222220220001221012111022121120202022211221112202002121022",
"0221221101200221120220011002222100000020200021121021020122100021201010210202002000101020022121100100",
"0100122100011112100120210020011102201122122100100120122212000021220022012202201100010212222110222020"]
            
    def encrypt(self,ID):
        """Function to encrypt the given sequence.
        Parameters
        ----------
        ID : String
            ID variable. It needs to be a 2 digit string in base 3.
        """
        if len(ID) != 2 or not all(int(i) < 3 for i in ID):
            raise Exception("Please provide a 2 digit ID in base 3")
        else:
            self.ID = ID
            
        base3 = utils.transform_Huffman_base3(self.Huffman_codec,self.sequence,to_b3=True) ## check if all characters are valid.
        S1 = utils.get_S1(base3)
        S1_base3 = utils.get_S1_base(S1)
        S2 = utils.get_S2(S1_base3)
        S3 = utils.get_S3(S1,S2)
        S4 = S2+S1+S3
        S5 = utils.base3_DNA_transformation(self.DNA_convertion_table,S4,direction=True,prior=None)

        n25_chunks = [S5[i:i+25] for i in range(0, len(S5), 25)]
        n_Fragments = (len(n25_chunks)-4)+1
        
        F_list = []
        for i in range(0,n_Fragments):
            F_list.append("".join(n25_chunks[i:i+4]))
        F_list = [F_list[i] if i%2 == 0 else utils.reverse_completement_DNA(F_list[i]) for i in range(len(F_list))]
        
        F_list_randomized = []
        for i in range(len(F_list)):
            F_list_randomized.append(utils.encrypt_decript(F_list[i], self.keystreams[i%4], True))
        
        IX_list = []
        for i in range(len(F_list_randomized)):
            IX_list.append(utils.get_IX(ID,i))
        
        IX_DNA_list = []
        for i in range(len(IX_list)):
            IX_DNA_list.append(utils.base3_DNA_transformation(self.DNA_convertion_table,IX_list[i],prior = F_list_randomized[i][-1]))
        
        F_IX_list = []
        for i in range(len(IX_DNA_list)):
            F_IX_list.append(F_list_randomized[i]+IX_DNA_list[i])
        
        F_complete_list = []
        for i in range(len(F_IX_list)):
            F_complete_list.append(utils.add_border_NTs(F_IX_list[i]))

        self.fragments = F_complete_list
        
        return "".join(F_complete_list)
            
            
            
    def decrypt(self):
        """Function to decrypt the given DNA sequence.
        """
        DNA_seq = self.sequence.upper()
        if self._valid_dna_sequence(string = DNA_seq) == False:
            raise Exception("Please provide a valid DNA sequence")
            
        F_complete_list = [DNA_seq[i:i+117] for i in range(0, len(DNA_seq), 117)]
        
        IX_DNA_list = []
        F_list_randomized = []
        for i in range(len(F_complete_list)):
            IX_DNA,F_randomized = utils.partition_F(F_complete_list[i])
            IX_DNA_list.append(IX_DNA)
            F_list_randomized.append(F_randomized)
            
        IX_list = []
        for i in range(len(IX_DNA_list)):
            IX_list.append(utils.base3_DNA_transformation(self.DNA_convertion_table,IX_DNA_list[i],direction=False,prior = F_list_randomized[i][-1]))
            
        ID_list = []  ## Need to check if these IDs all match
        order_list_base3 = []
        for i in range(len(IX_list)):
            ID, n_F = utils.get_ID_and_n_F(IX_list[i])
            ID_list.append(ID)
            order_list_base3.append(n_F)
        order_list = [int(str(i),3) for i in order_list_base3]
            
        F_list = []
        for i in range(len(F_list_randomized)):
            F_list.append(utils.encrypt_decript(F_list_randomized[i], self.keystreams[order_list[i]%4], False))
            
        F_list = [F_list[i] if order_list[i]%2 == 0 else utils.reverse_completement_DNA(F_list[i]) for i in range(len(F_list))]
        
        F_chunks_dict = {}
        for i in range(len(F_list)):
            F_chunks_dict.update({order_list[i]:[F_list[i][j:j+25] for j in range(0, len(F_list[i]), 25)]})
        
        n_Fragments = max(F_chunks_dict.keys())+1
        
        n_n25_chunks = n_Fragments+3
        
        n25_chunks_dict = {}
        for i in range(0,n_n25_chunks):
            n25_chunks_dict.update({i:[]})
        
        for key in F_chunks_dict.keys():
            F_chunks_list = F_chunks_dict.get(key)
            for i in range(0,len(F_chunks_list)):
                indexer = key+i
                data = n25_chunks_dict.get(indexer)
                data.append(F_chunks_list[i])
                n25_chunks_dict.update({indexer:data})


        S5 = "".join([np.unique(n25_chunks_dict.get(i))[0]for i in range(0,n_n25_chunks)]) ## if there is a mutation or a gap, this will give error, need to add a clause for mutations

        S4 = utils.base3_DNA_transformation(self.DNA_convertion_table,S5,direction=False,prior=None)

        S1,S2,S3 = utils.revert_S4(S4)

        return utils.transform_Huffman_base3(self.Huffman_codec,S1,to_b3=False)
        
        
    def _valid_dna_sequence(self,string):
        """Function for checking if a string is a DNA sequence.
        Parameters
        ----------
        string : String
            String to be checked.
        """
        dna_bases = ['A', 'C', 'G', 'T']
        for base in string:
            if base not in dna_bases:
                return False
        return True
    
    def _read_Huffman_codec(self):
        """Reading Huffman codec 
        """
        try:
            Huffman_codec = pd.read_csv("./View_huff3.cd.new",sep=" ",encoding= 'unicode_escape')["0\t"].str.split("\t",expand=True)
            Huffman_codec[1] = Huffman_codec[1].str.replace("space"," ")
            Huffman_codec.index = Huffman_codec[1]
            del Huffman_codec[0]
            del Huffman_codec[1]
            return Huffman_codec
        except:
            "Please make sure View_huff3.cd.new file is in the same folder."
            
    def _generate_DNA_conv_table(self):
        """Generating the DNA to base 3 conversion table 
        """
        DNA_convertion_table = pd.DataFrame({"A":["C","G","T"],"C":["G","T","A"],"G":["T","A","C"],"T":["A","C","G"]}).T
        DNA_convertion_table.columns = DNA_convertion_table.columns.astype(str)
        return DNA_convertion_table