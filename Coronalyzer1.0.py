import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
#from tabulate import tabulate
import warnings
warnings.filterwarnings("ignore")
pd.plotting.register_matplotlib_converters()
#import seaborn as sns
from Bio import pairwise2
from Bio import SeqIO
from Bio import Entrez
from tabulate import *
pd.plotting.register_matplotlib_converters()


def Mutaciones():
    
    # This dictionary is used for reshaping the numpy array for each of the genes. 
    # For example, the gene ORF1ab numpy array will be reshaped with rows=115 and cols=115.
    # The second element in the list represent 'n' number of dummy 'N' to be added at the end of
    # each gene nucletode seq to make to compactible with the rows and cols of the array.
    numpy_image_dict = {'gene=ORF1ab':[(115,115),7],
                        'gene=S':[(62,62),22],
                        'gene=ORF3a':[(28,30),12],
                        'gene=E':[(15,16),12], 
                        'gene=M':[(26,27),33],
                        'gene=ORF6':[(14,14),10],    
                        'gene=ORF7a':[(19,20),14],
                        'gene=ORF7b':[(12,12),12],
                        'gene=ORF8':[(19,20),14],
                        'gene=N':[(36,36),36],
                        'gene=ORF10':[(11,11),4]}
    
    
    # This dictionary has the codons for the amino acids from mRNA.
    amacid_dict = {'UUU':('F','PHE','Phenylalanine'),'UUC':('F','PHE','Phenylalanine'),
                   'UUA':('L','LEU','Leucine'),'UUG':('L','LEU','Leucine'),
                   'UCU':('S','SER','Serine'),'UCC':('S','SER','Serine'),
                   'UCA':('S','SER','Serine'),'UCG':('S','SER','Serine'), 
                   'UAU':('Y','TYR','Tyrosine'),'UAC':('Y','TYR','Tyrosine'),
                   'UAA':('STOP'),'UAG':('STOP'),
                   'UGU':('C','CYS','Cysteine'),'UGC':('C','CYS','Cysteine'),
                   'UGA':('STOP'),'UGG':('W','TRP','Tryptophan'),
                   'CUU':('L','LEU','Leucine'),'CUC':('L','LEU','Leucine'),
                   'CUA':('L','LEU','Leucine'),'CUG':('L','LEU','Leucine'),
                   'CCU':('P','PRO','Proline'),'CCC':('P','PRO','Proline'),
                   'CCA':('P','PRO','Proline'),'CCG':('P','PRO','Proline'),
                   'CAU':('H','HIS','Histidine'),'CAC':('H','HIS','Histidine'),
                   'CAA':('Q','GLU','Glutamine'),'CAG':('Q','GLU','Glutamine'),
                   'CGU':('R','ARG','Arginine'),'CGA':('R','ARG','Arginine'),
                   'CGG':('R','ARG','Arginine'),'CGC':('R','ARG','Arginine'),
                   'AUU':('I','ILE','Isoleucine'),'AUC':('I','ILE','Isoleucine'),
                   'AUA':('I','ILE','Isoleucine'),'AUG':('M','MET','Methionine'),
                   'ACU':('T','THR','Threonine'),'ACC':('T','THR','Threonine'),
                   'ACA':('T','THR','Threonine'),'ACG':('T','THR','Threonine'),
                   'AAU':('N','ASN','Asparagine'),'AAC':('N','ASN','Asparagine'),
                   'AAA':('K','LYS','Lysine'),'AAG':('K','LYS','Lysine'),
                   'AGU':('S','SER','Serine'),'AGC':('S','SER','Serine'),
                   'AGG':('R','ARG','Arginine'),'AGA':('R','ARG','Arginine'),
                   'GUU':('V','VAL','Valine'),'GUC':('V','VAL','Valine'),
                   'GUA':('V','VAL','Valine'),'GUG':('V','VAL','Valine'),
                   'GCU':('A','ALA','Alanine'),'GCC':('A','ALA','Alanine'),
                   'GCA':('A','ALA','Alanine'),'GCG':('A','ALA','Alanine'),
                   'GAU':('D','ASP','Aspartate'),'GAC':('D','ASP','Aspartate'),
                   'GAA':('E','GLU','Glutamate'),'GAG':('E','GLU','Glutamate'),
                   'GGU':('G','GLY','Glycine'),'GGC':('G','GLY','Glycine'),
                   'GGA':('G','GLY','Glycine'),'GGG':('G','GLY','Glycine')
    }
    #This is the core module developed as part of the effort.
    class dna:
    
        # Constructor method
        def __init__(self,dna_seq):
            dna_seq = dna_seq.upper() # Convert the nucleotide bases to Upper Case 
            for seq in dna_seq:
                # Valid nucleotide bases. If not a valid sequence raise an Error
                if seq not in ['A','T','G','C',' ','N']: 
                    error = 'Wrong DNA Sequence {}!!'.format(seq)
                    raise ValueError(error)
            # Remove all of the empty characters in the nucleotide sequence.
            dna_seq = dna_seq.replace(' ','') 
            self.dir_3_5=dna_seq
            self.dir_5_3=self.dir_5_3_strand()
            self.mRna = None
            self.amino_acid = None
            self.num_array = None
            self.nucl_len = len(dna_seq)
    
        def __repr__(self):
            return "DNA has {} nucleotide and they are {} :".format(self.nucl_len,self.dir_3_5)
    
        def __eq__(self, other):
            if other is None:
                return False
            return self.seq == other.seq
    
        #def replicate(self):
        #    return
    
        def transcription(self):
            # This is a method that imitates the transcription of a gene to mRNA for Protein transalation.
            # This is mostly of the future use.
            trans=''
            for nuc in self.dir_5_3:
                if nuc == 'A':
                    trans += 'U'
                if nuc == 'T':
                    trans += 'A'
                if nuc == 'C':
                    trans += 'G'
                if nuc == 'G':
                    trans += 'C'
                if nuc == 'N':
                    trans += 'N'
            self.mRna = trans
            return self.mRna
    
        def translation(self):
            # This is the method where the transcripted mRNA gets translated into Amino Acid. Each 3
            # base in the mRNA codes for an amino acid.
            begin = 'No'
            ac = ''
            for i in range(0,len(self.mRna)-3,3):
                if self.mRna[i:3] == 'AUG':
                    begin = 'Yes'
                if self.mRna[i:3] in ('UAA','UAG','UGA'):
                    being = 'No'
                if begin == 'Yes':
                    ac+= amacid_dict[self.mRna[i:3+i]][0]
            self.amino_acid = ac
            return self.amino_acid
    
        def dir_5_3_strand(self):
            dir_5_3 = ''
            # This is a method which reads the 3 - 5 prime sequence and creates the 5 - 3 prime sequence.
            for nuc in self.dir_3_5:
                if nuc == 'A':
                    dir_5_3 += 'T'
                if nuc == 'T':
                    dir_5_3 += 'A'
                if nuc == 'C':
                    dir_5_3 += 'G'
                if nuc == 'G':
                    dir_5_3 += 'C'
                if nuc == 'N':
                    dir_5_3 += 'N'
            return dir_5_3
    
        def numpfy(self):
            # This method takes in a dna sequence and convert them into numpy array.
            # Each of the nucleotide sequence is converted into one of the below numbers 
            # which then can be used in for analysis and comparison.
            arr = ''
            for i in self.dir_3_5:
                if i == 'A':
                    arr += '0 '
                if i == 'T':
                    arr += '255 '
                if i == 'C':
                    arr += '100 '
                if i == 'G':
                    arr += '200 '
                if i == 'N':
                    arr += '75 '   
            arr_np = np.fromstring(arr,dtype=np.uint8,sep=' ')        
            self.num_array = arr_np
            return self.num_array
        
        
        
    def read_dna_seq(file_name):
        # This method reads the dna sequence from the file downloaded from NCBI and crates a python dictionary.
        fil = open(file_name,'r')
        fil_list = fil.readlines()
        fil.close
        
        genome = {}
        gene_name = ''
        protein_name = ''
        gene_seq = ''
        for i in fil_list:
            if i[0] == '>':
                # Reads each line from the file and creates a dictionary with the following information for each
                # gene. {<'gene_name-1'>:[<protein_name>,nucleotide sequence],
                #        <'gene_name-2'>:[<protein_name>,nucleotide sequence],
                #        <'gene_name-2'>:[<protein_name>,nucleotide sequence]}
                if list(genome.keys()) != []:
                    gene_seq = gene_seq.replace('\n','')
                    genome[gene_name].append(gene_seq)
                gene_seq = ''
                g_st = i.find('[gene=')
                g_end = i[g_st:].find(']')
                p_st = i.find('[protein=')
                p_end = i[p_st:].find(']') 
    
                if g_st > 0 and g_end > 0:
                    gene_name = i[g_st+1:g_st+g_end]
                    genome[gene_name] = []
                
                if p_st > 0 and p_end > 0:
                    protein_name = i[p_st+1:p_st+p_end]
                    genome[gene_name].append(protein_name)
            else:
                gene_seq += i
        gene_seq = gene_seq.replace('\n','')
        genome[gene_name].append(gene_seq)    
        return genome
    
    def gene_mod(genome):
        # This method modifies each of the sequence with dummy nucleotide 'N' so that for the shape of the numpy array.
        genome_keys = list(genome.keys())
        for k in genome_keys:
            if len(numpy_image_dict[k]) > 1:
                N = numpy_image_dict[k][1]
                seq = add_N(N,genome[k][1])
                genome[k][1] = seq
        return genome
    
    def add_N(n,seq):
        # This method is called from gene_mod() method, for creating dummy nucleotide 'N'.
        for i in range(0,n):
            seq += 'N'
        return seq
    
    nombre1=input("ingrese nombre del archivo (con el .txt por favor): ")
    nombre2=input("ingrese nombre del archivo (con el .txt por favor): ")
    # Read the dna sequence file-1 previously downloaded from NCBI.
    dict_seq_1 = read_dna_seq(nombre1)
    # Modify the sequence with dummy 'N' nucleotide.
    dict_seq_1 = gene_mod(dict_seq_1)
    
    # Read the dna sequence file-2 previously downloaded from NCBI.
    dict_seq_2 = read_dna_seq(nombre2)
    # Modify the sequence with dummy 'N' nucleotide.
    dict_seq_2 = gene_mod(dict_seq_2)
    
    # Create matplotlib subplots for each gene. 
    f,ax = plt.subplots(nrows=11,ncols=3,figsize=(25,30))
    gene_name = list(numpy_image_dict.keys())
    row = 0
    col = 0
    mut_dict={}
    for i in gene_name:
        G = i[5:]
        # Loop thru each gene in the Cornona Virus nucleotide sequence.
        gene_us = dna(dict_seq_1['gene='+G][1])
        # Invoke the transcription method of the class dna 
        gene_us.transcription()
        # Invoke the mothod that converts the gene sequence into a numpy array.
        numpfy_usa = gene_us.numpfy()
        # Reshape the numpy array with a predeifned shape from the numpy_image_dict dictionary.
        numpfy_usa = numpfy_usa.reshape(numpy_image_dict['gene='+G][0])
        # sub-plot the numpy array with matplotlib pcolor method.
        ax[row][col].pcolor(numpfy_usa)
        ax[row][col].set_title(G+' Gene - USA')
        col+=1
        gene_china = dna(dict_seq_2['gene='+G][1])
        # Invoke the transcription method of the class dna 
        gene_china.transcription()
        # Invoke the mothod that converts the gene sequence into a numpy array.
        numpfy_china = gene_china.numpfy()
        # Reshape the numpy array with a predeifned shape from the numpy_image_dict dictionary.
        numpfy_china = numpfy_china.reshape(numpy_image_dict['gene='+G][0])
        # sub-plot the numpy array with matplotlib pcolor method.
        ax[row][col].pcolor(numpfy_china)
        ax[row][col].set_title(G+' Gene - CHINA')
        col+=1
    
        # To find the gene mutation subtract the numpy array from base sequence with the newer sequence. Here the 
        # the Chinese sequence is the base sequence and the USA sequence is a newer sequence.
        mut = numpfy_china - numpfy_usa
        if mut.any():
            # Here we are looking for a non zero value in the mutated numpy array (result of the subtracting the 2 numpy arrays).
            # Presence of non-zero value means that there is difference between the 2 numpy arrays and the gene has 
            # mutataions. If there are mutations in the gene create a python dictionary "mut_dict" with details as below.
            # {'<Gene_Name-1>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]], '<Gene_Name-2>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]]}
            mut_nec = np.nonzero(mut)
            x=mut_nec[0]
            y=mut_nec[1]
            l=0
            mut_dict[G]=[]
            for i in x:
                us_base = numpfy_usa[i][y[l]]
                ch_base = numpfy_china[i][y[l]]
                mut_base = mut[i][y[l]]
                info_list = [ch_base,us_base,mut_base,(i,y[l])]
                mut_dict[G].append(info_list)
                print("Mutated DNA Base {} in China and Base {} in USA at position {} For the Gene {}".format(ch_base,us_base,(i,y[l]),G))
                l+= 1
        # Giving a title to the matplotlib subplot
        ax[row][col].pcolor(mut)
        ax[row][col].set_title(G+' Gene - Mutataion')
        row+= 1
        col=0
    
    f.tight_layout()
    # Saving the matplotlib subplot as a jpg.
    f.savefig('Sars_Cov-2_Gene_Mutation.jpg')

def AnalisiSeq():
    
    nombre3=input("ingrese nombre del archivo (con el .fna por favor): ")
    coronavirus = open(nombre3, "r")
    print('HEADER:',coronavirus.readline())
    coronavirus = coronavirus.readlines()
    COVID_seq = ''
    for line in coronavirus:
        line = line.strip()
        COVID_seq += line
    print(COVID_seq[0:1000])
    
    def basic_properties(DNAseq):
        total_base = len(DNAseq)
        num_Adenine = DNAseq.count('A')
        num_Guanine = DNAseq.count('G')
        num_Thymine = DNAseq.count('T')
        num_Cytosine = DNAseq.count('C')
        
        if total_base != num_Adenine + num_Guanine + num_Thymine + num_Cytosine:
            print('Something is not right')
        else : pass
        
        A_percent = num_Adenine / total_base
        G_percent = num_Guanine / total_base
        T_percent = num_Thymine / total_base
        C_percent = num_Cytosine / total_base
        
        #visualization
        x = np.arange(4)
        bases = ['Adenine', 'Guanine', 'Thymine' ,'Cytosine']
        values = [num_Adenine, num_Guanine, num_Thymine, num_Cytosine]
        plt.bar(x,values)
        plt.xticks(x, bases)
        plt.show()
        table = [['total base',total_base,'Percentage',str('100%')],
                 ['Adenine:',num_Adenine, 'Percentage:',str(round(A_percent*100,2))+'%'],
                ['Guanine:',num_Guanine, 'Percentage:',str(round(G_percent*100,2))+'%'],
                 ['Thynime:',num_Thymine, 'Percentage:',str(round(T_percent*100,2))+'%'],
                 ['Cytosine:',num_Cytosine, 'Percentage:',str(round(C_percent*100,2))+'%']]
        print(tabulate(table))
        print('GC content:', round((((num_Guanine + num_Cytosine) / total_base)*100),2),'%')
        
    print(basic_properties(COVID_seq))
    
    def transcription(DNAseq):
        mRNAseq = DNAseq.replace('T','U')
        print(len(mRNAseq))
        return mRNAseq
    COVID_mRNA = transcription(COVID_seq)
    print("--------------------------mRNA sequence------------------------------")
    print(COVID_mRNA[0:1000])
    print("---------------------------Protein Sequence-----------------------------")
    def translate(seq):
        table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
        protein =""
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table:
                protein+= table[codon]
        return protein
    
    protein_seq = translate(COVID_seq)
    print(protein_seq[0:1000])
    
    def visualization(protein_seq):
        # composition of Amino Acids
        plt.figure(figsize=(35,50))
        x = np.arange(22)
        AA = ['Arginine', 'Histidie','Lysine',
              'Aspartic Acid','Glutamic Acid',
              'Serine','Threonine','Asparagine','Glutamine',
              'Cysteine','Selenocysteine','Glycine','Prolien',
              'Alanine','Valine','Isoleucine','Leucine',
              'Methionine','Phenylalanine','Tyrosine','Tryptophan',
              'Stop Codon']
        values = [protein_seq.count('R'),protein_seq.count('H'),protein_seq.count('K'),
                  protein_seq.count('D'),protein_seq.count('E'),
                  protein_seq.count('S'),protein_seq.count('T'),protein_seq.count('N'),protein_seq.count('Q'),
                  protein_seq.count('C'),protein_seq.count('U'),protein_seq.count('G'),protein_seq.count('P'),
                  protein_seq.count('A'),protein_seq.count('V'),protein_seq.count('I'),protein_seq.count('L'),
                  protein_seq.count('M'),protein_seq.count('F'),protein_seq.count('Y'),protein_seq.count('W'),
                  protein_seq.count('_')]
        plt.subplot(2,2,1)
        plt.rc('font',size = 20)
        plt.barh(AA,values,height=0.6)
        plt.title('AA in general')
        
        # visualization by groups as shown in figure
        x = np.arange(4)
        Electric =  protein_seq.count('R')+protein_seq.count('H')+protein_seq.count('K')+protein_seq.count('D')+protein_seq.count('E')
        Uncharged = protein_seq.count('S')+protein_seq.count('T')+protein_seq.count('N')+protein_seq.count('Q')
        Special =   protein_seq.count('C')+protein_seq.count('U')+protein_seq.count('G')+protein_seq.count('P')
        Hydrophobic = protein_seq.count('A')+protein_seq.count('V')+protein_seq.count('I')+protein_seq.count('L')+protein_seq.count('M')+protein_seq.count('F')+protein_seq.count('Y')+protein_seq.count('W')
        
        plt.subplot(2,2,2)
        types = ['Elecrically charged','Polar uncharged',
                 'Special case', 'Hydrophobic Side Chain']
        values = [Electric,Uncharged,Special,Hydrophobic]
        plt.barh(types, values, height = 0.6)
        plt.title('AA in groups')
        
        # Visualization of positive and negative amino acid in Hydrophobic side chain group
        plt.subplot(2,2,3)
        x = np.arange(2)
        positive = protein_seq.count('R')+protein_seq.count('H')+protein_seq.count('K')
        negative = protein_seq.count('D')+protein_seq.count('E')
        types = ['Positive','Negative']
        values = [positive, negative]
        plt.barh(types,values, height=0.6)
        plt.title('Charge diff in Electrically charged side chain')
        plt.show()
        
        # Visualization abundance of Amino Acid by relative frequenceies of amino aicd residue in secondary structure
        alpha_helix = protein_seq.count('A')+protein_seq.count('C')+protein_seq.count('L')+protein_seq.count('M')+protein_seq.count('E')+protein_seq.count('Q')+protein_seq.count('H')+protein_seq.count('K')
        beta_sheet =protein_seq.count('V')+protein_seq.count('I')+protein_seq.count('F')+protein_seq.count('Y')+protein_seq.count('W')+protein_seq.count('T')
        turn = protein_seq.count('G')+protein_seq.count('S')+protein_seq.count('D')+protein_seq.count('N')+protein_seq.count('P')
        x = np.arange(3)
        plt.subplot(2,2,4)
        types = ['Alpha helix','beta sheet','Turn']
        values = [alpha_helix, beta_sheet, turn]
        plt.barh(types,values, height=0.6)
        plt.title('AA residues in secondary structure')
        plt.show()
   
    print(visualization(protein_seq))
    print("---------------------------Protein Splicing-----------------------------")
    def proteinseq_splic(protein_sequence):
        protein_group = protein_sequence.split('_')
        return protein_group
    protein_list = proteinseq_splic(protein_seq)
    print(proteinseq_splic(protein_seq)[0:100])
    print('Number of proteins:',len(protein_list))
    def Functional_protein(protein_list):
        funct_group = []
        for protein in protein_list:
            if len(protein) > 20:
                funct_group.append(protein)
        return funct_group
    function_protein = Functional_protein(protein_list)
    print("---------------------------Functional Protein-----------------------------")
    print(function_protein[0:52])
    print('Number of Functional proteins:',len(function_protein))
    
    dna_ORF = []
    
    def rf(seq):
        #getting all reading frames
        i=0
        while i+2<len(seq):
            if seq[i:i+3]=='ATG':
                j=i+3
                while j+2<len(seq):
                    if seq[j:j+3] in ['TGA', 'TAA', 'TAG']:
                        dna_ORF.append(seq[i:j+3])
                    j=j+3
            i=i+3
        return dna_ORF
    
    def RC(seq):
        #reverse compliment
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        bases = list(seq) 
        bases = reversed([complement.get(base,base) for base in bases])
        bases = ''.join(bases)
        return bases
       
    def ORF(seq):
        #open reading frame
        #num=0
        seqRC=RC(seq)
        for i in [0,1,2]:
            rf(seq[i:])
            rf(seqRC[i:])
            
    RC(COVID_seq)[0:1000]
    
    print("-----------------------------DNA ORF---------------------------------")
    ORF(COVID_seq)
    print(dna_ORF[0:10])
    print(len(dna_ORF))
    
    # converting dnaORF to proteinORF
    candidate_pro_ORF = []
    def dna_pro_ORF(dnaORF):
        for ORF in dnaORF:
            candidate_pro_ORF.append(translate(ORF))
        return candidate_pro_ORF

    print(dna_pro_ORF(dna_ORF)[0:10])
    
    valid_group = []
    def valid_ORF(protein_seq):
        for protein in protein_seq:
            if protein.count('_') == 1:
                valid_group.append(protein)
        return valid_group
    
    print(valid_group[0:10])
    len(valid_group)   
    
    functional_ORF = Functional_protein(valid_group)
    print(functional_ORF[0:20])
    print('number of functional protein:',len(functional_ORF))
    
    sort_functional_ORF = list(sorted(functional_ORF, key = len, reverse=True))
    largest_FP = sort_functional_ORF[0]
    print(largest_FP[0:1000])
    print(len(largest_FP))
    visualization(largest_FP)
    
    second_FP = sort_functional_ORF[1]
    print(second_FP[0:1000])
    print('protein size:',len(second_FP))
    visualization(second_FP)
    
    third_FP = sort_functional_ORF[2]
    print(third_FP[0:1000])
    print('protein size:',len(third_FP))
    visualization(third_FP)
    
def Comparar():
    nombre4=input("ingrese nombre del archivo (con el .fasta por favor): ")
    nombre5=input("ingrese nombre del archivo (con el .fasta por favor): ")
    nombre6=input("ingrese nombre del archivo (con el .fasta por favor): ")
    SARS = SeqIO.read(nombre4, 'fasta')
    MERS = SeqIO.read(nombre5, 'fasta')
    COV2 = SeqIO.read(nombre6, 'fasta')
    
    #Squiggle wont work in Kaggle kernels, you need to run this command on terminal.
    
    
    # Alignments using pairwise2 alghoritm
    
    SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
    print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)
    
    MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
    print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)
    
    MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
    print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)
    
    # Plot the data
    X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
    Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
    plt.title('Sequence identity (%)')
    plt.bar(X,Y,color=(0.2, 0.4, 0.6, 0.6))



def descargar_fasta():
    archivo=input("Ingrese el nombre del archivo.fasta:")

    Entrez.email = "mingyan24@126.com"
    hd1 = Entrez.efetch(db="nucleotide",id=[archivo],rettype='fasta')
    seq = SeqIO.read(hd1,'fasta')
    fw = open(archivo,'w')
    SeqIO.write(seq,fw,'fasta')
    fw.close()
    os.getcwd()
    
def Menu():
 
    correcto=False
    num=0
    while(not correcto):
        try:
            num = int(input("Elige metodo de análisis: "))
            correcto=True
        except ValueError:
            print('Error, vuelve a elegir')
     
    return num
 
salir = False
opcion = 0
 
while not salir:
 
    print ("1. Mutaciones entre variantes")
    print ("2. Análisis de Secuencia Covid")
    print ("3. Comparacion secuencia variantes")
    print ("4. Salir")
     
    print ("Elige una opcion")
 
    opcion = Menu()
 
    if opcion == 1:
        Mutaciones()
    elif opcion == 2:
        AnalisiSeq()
    elif opcion == 3:
        Comparar()
    elif opcion == 4:
        salir = True
    else:
        print ("Elige una opción")
 
print ("Fin")



