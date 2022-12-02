'''
Versión 1.1.0 de Coronalyzer. Programa desarrollado para analizar las secuencias de las diferentes variantes del Sars-Cov2. 
Tipo de análisis: Análisis de Secuencia, Análisis de mutaciones, Análisis de Similitud entre secuencias.
Función añadida: Descarga de archvivos PDB, Visualización 
Funciones extra: Descarga de Fasta y Fasta cds, preparación de fastas para clustl (Align Tree)
Autor: Pedro Manuel López Zarzuela 
Fecha: 03/12/2021
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
from tabulate import tabulate
import warnings
warnings.filterwarnings("ignore")
pd.plotting.register_matplotlib_converters()
import seaborn as sns
import Bio
from Bio import pairwise2
from Bio import ExPASy
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBList
from Bio.Align import substitution_matrices
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo
from io import StringIO
#from Bio import Phylo
from tkinter import *
import tkinter as tk
from prody import *
from numpy import *
from matplotlib.pyplot import *
import matplotlib.image as mpimg
import pylab
import pandas as pd
import nglview as nv

text1=""
text2=""
text3=""
text4=""
text5=""
text6=""
text7=""
text8=""

'''Función do_nothing() para ejectutar el menú de tkinter sin llamar a las funciones del programa y que pueda ejecutarse'''

def donothing():
   x = 0

'''Función limpiarCampos(), función que limpia los campos de las entradas (Entry) de nuestro menu contextual de Tkinter'''

def limpiarCampos():
    
    miNombre1.set("")
    miNombre2.set("")
    miNombre3.set("")
    miNombre4.set("")
    miNombre5.set("")
    
'''Función descargar_fasta(), como su propio nombre indica nos descarga los fastas introduciendo el codigo de del NCBI
descarga los fasta y los fasta_cds'''

def descargar_fasta():
    global text4
    text4=nombreID4.get()
    print(text4)
    
    Entrez.email = "8tato8@gmail.com"
    hd1 = Entrez.efetch(db="nucleotide",id=[text4],rettype='fasta')
    seq = SeqIO.read(hd1,'fasta')
    fw = open(text4,'w')
    SeqIO.write(seq,fw,'fasta')
    
    
    handle = Entrez.efetch(db="nucleotide", id=[text4], rettype="fasta_cds_na", retmode="text")
    seq=handle.read()
    text_file = open(text4+'.txt', "w")
    n = text_file.write(seq)
    
    text_file.close()
    
    fw.close()
    os.getcwd()
    
    return(text4)

'''Función descargar_PDB, como su propio nombre indica descarga archivos del PDB introduciendo el codigo de 4 letras'''

def descargar_PDB():
    global text5
    text5=nombreID5.get()
    pdbl=PDBList()
    pdbl1=[text5]
    for i in pdbl1:
        pdbl.download_pdb_files(pdbl1, pdir='PDB', file_format='mmCif')
    
    for i in pdbl1:
        pdbl.download_pdb_files(pdbl1, pdir='PDB', file_format='pdb')
    
    for i in pdbl1:
        pdbl.download_pdb_files(pdbl1, pdir='PDB', file_format='fasta')
        
    return(text5)

'''Función Analizar_Seq(), hace un analisis de la secuencia obtenida a partir de un archivo .fasta: 
Imprime secuencia, la traduce a mRNA, transcribe a proteina, reliza el splicin, la RC y ORF'''

def Analizar_Seq():
    global text1

    text1=nombreID1.get()

    print(text1)
    
    coronavirus = open(text1, "r")
 
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
        
    basic_properties(COVID_seq)
    
    def transcription(DNAseq):
        mRNAseq = DNAseq.replace('T','U')
        print(len(mRNAseq))
        return mRNAseq
    COVID_mRNA = transcription(COVID_seq)
    print("------------------------------------------------------mRNA Sequence-------------------------------------------------------")
    print(COVID_mRNA[0:1000])
    print("------------------------------------------------------Protein Sequence-------------------------------------------------------")
    
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
    
    protein_seq = translate(COVID_seq)
    protein_seq[0:1000]
    
    def visualization(protein_seq):
        # composición de los Aminoácidos
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
        
        # Visualización de grupos segun su carga 
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
        
        # Visualización de los aminoacidos positivos y negativos en grupos de cadena laterales hidrofóbicas
        plt.subplot(2,2,3)
        x = np.arange(2)
        positive = protein_seq.count('R')+protein_seq.count('H')+protein_seq.count('K')
        negative = protein_seq.count('D')+protein_seq.count('E')
        types = ['Positive','Negative']
        values = [positive, negative]
        plt.barh(types,values, height=0.6)
        plt.title('Charge diff in Electrically charged side chain')
        plt.show()
        
        # Visualizacion de la abundancia de Aminoacidos segun su estructura secundaria en la cadena proteica 
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
        
    visualization(protein_seq)
    
    print("------------------------------------------------------Protein Splicing-------------------------------------------------------")
    def proteinseq_splic(protein_sequence):
        protein_group = protein_sequence.split('_')
        return protein_group
    protein_list = proteinseq_splic(protein_seq)
    print(proteinseq_splic(protein_seq)[0:100])
    print('Number of proteins:',len(protein_list))
    
    print("-----------------------------------------------------Functional Protein-------------------------------------------------------")
    def Functional_protein(protein_list):
        funct_group = []
        for protein in protein_list:
            if len(protein) > 100:
                funct_group.append(protein)
        return funct_group
    function_protein = Functional_protein(protein_list)
    print(function_protein[0:52])
    print('Number of Functional proteins:',len(function_protein))
    
    print("--------------------------------------------------DNA OPEN READING FRAME------------------------------------------------------")

    record = SeqIO.read(text1, "fasta")
    table = 11
    min_pro_len = 100
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
         for frame in range(3):
             length = 3 * ((len(record)-frame) // 3) #Multiple of three
             for pro in nuc[frame:frame+length].translate(table).split("*"):
                 if len(pro) >= min_pro_len:
                     print("%s...%s - length %i, strand %i, frame %i" \
                           % (pro[:30], pro[-3:], len(pro), strand, frame))
 
    return (text1,text2,text3)
     
def Analizar_Mut():
    global text1
    global text2
    text1=nombreID1.get()
    text2=nombreID2.get()
    print(text1,text2)
  
    # Este diccionario se usa para restablecer el tamaño de la numpy array por cada uno de los genes. 
    # La numpy array del gen  ORF1ab se establecera con rows=115 y cols=115.
    # El segundo elemento de la lista hace referencia a el numero 'n' dummy Nucleotides que seran añadidos al final
    # de cada gen para hacer las secuencias compactas.
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
    #Este es el modulo principal .
    class dna:
    
        # metodo cosntructor
        def __init__(self,dna_seq):
            dna_seq = dna_seq.upper() # Convirte bases en mayusculas 
            for seq in dna_seq:
                # Valida las bases. Si la secuencia no es valida salta un error
                if seq not in ['A','T','G','C',' ','N']: 
                    error = 'Wrong DNA Sequence {}!!'.format(seq)
                    raise ValueError(error)
            # Quita todos los caracteres vacios de la secuencia.
            dna_seq = dna_seq.replace(' ','') 
            self.dir_3_5=dna_seq
    
        def __repr__(self):
            return "DNA has {} nucleotide and they are {} :".format(self.nucl_len,self.dir_3_5)
    
        def __eq__(self, other):
            if other is None:
                return False
            return self.seq == other.seq
    
        def numpfy(self):
            # Este metodo comnvierte el Dna en un numpy array.
            # Cada nucleotido de la secuencia se traduce a un numero que vemos a continuación 
            # Esto puede ser usado para analisis o comparación, como el One-Hot-Encoding.
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
        global fil
        # Este metodo lee la secuencia fasta descargada del NCBI y crea un diccionario con ella.
        try:
            fil = open(file_name,'r')
            fil_list = fil.readlines()
            fil.close
        except: FileNotFoundError()
        
        genome = {}
        gene_name = ''
        protein_name = ''
        gene_seq = ''
        for i in fil_list:
            if i[0] == '>':
                # Leemos cada linea del archivo y creamos un diccionario con la siguiente información
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
        # Substitución por dummy nucleotide 'N' para compactar la numpy array.
        genome_keys = list(genome.keys())
        for k in genome_keys:
            if len(numpy_image_dict[k]) > 1:
                N = numpy_image_dict[k][1]
                seq = add_N(N,genome[k][1])
                genome[k][1] = seq
        return genome
    
    def add_N(n,seq):
        # Esta funcion la usa gene_mod para añadir dummy 'N'.
        for i in range(0,n):
            seq += 'N'
        return seq
   
    # Lee la secuencia 1 de DNA descargada anteriormente del NCBI.
    dict_seq_1 = read_dna_seq(text1+'.txt')
    # Modifica la secuencia con dummy 'N' nucleotide.
    dict_seq_1 = gene_mod(dict_seq_1)
    
    # Lee la secuencia 2 de DNA descargada anteriormente del NCBI.
    dict_seq_2 = read_dna_seq(text2+'.txt')
    # Modifica la secuencia con dummy 'N' nucleotide.
    dict_seq_2 = gene_mod(dict_seq_2)
    
    # Creamos plots con Matplotlib para cada gen. 
    f,ax = plt.subplots(nrows=11,ncols=3,figsize=(25,30))
    gene_name = list(numpy_image_dict.keys())
    row = 0
    col = 0
    mut_dict={}
    for i in gene_name:
        G = i[5:]
        # Loopeamos a traves de cada gen de la secuencia de Cornona Virus.
        gene_us = dna(dict_seq_1['gene='+G][1])
        # Invocamos el metodo que convierte gene la secuencia de genes en un numpy array.
        numpfy_usa = gene_us.numpfy()
        # Reestablecemos las dimensiones del numpy array.
        numpfy_usa = numpfy_usa.reshape(numpy_image_dict['gene='+G][0])
        # Hacemos plot del array.
        ax[row][col].pcolor(numpfy_usa)
        ax[row][col].set_title(G+' Gene '+text1)
        col+=1
        gene_china = dna(dict_seq_2['gene='+G][1])
        # Invocamos el metodo que convierte gene la secuencia de genes en un numpy array.
        numpfy_china = gene_china.numpfy()
        # Reestablecemos las dimensiones del numpy array.
        numpfy_china = numpfy_china.reshape(numpy_image_dict['gene='+G][0])
        # Hacemos el plot de la array
        ax[row][col].pcolor(numpfy_china)
        ax[row][col].set_title(G+' Gene '+text2)
        col+=1
    
        # Para encontrar las mutaciones lo que hacemos es restar la nueva secuencia a la  secuencia base.
        # Chinese sequence is the base sequence and the USA sequence is a newer sequence.
        mut = numpfy_china - numpfy_usa
        if mut.any():
            # Here we are looking for a non zero value in the mutated numpy array (result of the subtracting the 2 numpy arrays).
            # La presencia de un valor diferente de zero nos indica que en esa localizacion existe una mutación para ese nucleotido . 
            #{'<Gene_Name-1>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]], '<Gene_Name-2>': [[<value_of_base_seq>, <value_of_newer_seq>, <value_in_mutated_numpy>, (x_value,y_value)]]}
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
                cbase = 0
                usbase = 0
                dims = 0
                pos = 0
                posp = 0
                posh = 0
                
                if ch_base == 0:
                    cbase = 'A'
                elif ch_base == 255:
                    cbase = 'T'
                elif ch_base == 100:
                    cbase = 'C'
                elif ch_base == 200:
                    cbase = 'G'
                    
                if us_base == 0:
                    usbase = 'A'
                elif us_base == 255:
                    usbase = 'T'
                elif us_base == 100:
                    usbase = 'C'
                elif us_base == 200:
                    usbase = 'G'
                    
                
                if G == 'ORF1ab':
                    dims == 115
                elif G == 'S':
                    dims = 62
                elif G == 'ORF3a':
                    dims = 28
                elif G == 'E':
                    dims = 15
                elif G == 'M':
                    dims = 26
                elif G == 'ORF6':
                    dims = 14
                elif G == 'ORF7a':
                    dims = 19
                elif G == 'ORF7b':
                    dims = 12
                elif G == 'ORF8':
                    dims = 19
                elif G == 'N':
                    dims = 36
                elif G == 'ORF10':
                    dims = 11
                      
                pos = dims*i+y[l]
                posh = pos/3+1
                posp = round(posh)
                
                print("Mutated DNA Base {} in {} and Base {} in {} at position {} / {} / {} For the Gene {} \n".format(cbase,text2,usbase,text1,(i,y[l]),pos,posp,G))
                l+= 1
                
        # Poniendole título al plot
        ax[row][col].pcolor(mut)
        ax[row][col].set_title(G+' Gene - Mutation')
        row+= 1
        col=0

    f.tight_layout()
    # Guardamos el matplotlib subplot como un jpg.
    f.savefig('Sars_Cov-2_Gene_Mutation.jpg')
    img = mpimg.imread('Sars_Cov-2_Gene_Mutation.jpg')
    imgplot = plt.imshow(img)
    #plt.show()
    
    return (text1,text2,text3)

def Analizar_Sim():
    global text1
    global text2
    
    text1=nombreID1.get()
    text2=nombreID2.get()
    print(text1,text2)
  
    Seq1 = SeqIO.read(text1, 'fasta')
    Seq2 = SeqIO.read(text2, 'fasta') 
    
    Seq2_Seq1 = pairwise2.align.globalxx(Seq2.seq, Seq1.seq, one_alignment_only=True, score_only=True)
    print('Seq2/Seq1 Similarity (%):', Seq2_Seq1 / len(Seq1.seq) * 100)
    
    
    alignments = pairwise2.align.localms(Seq1.seq, Seq2.seq, 1, -2, -2, -0.5)
    print(pairwise2.format_alignment(*alignments[0])) 
     
      # Printeamos los plots con la data
    X = ['Seq2/Seq1']
    Y = [Seq2_Seq1/len(Seq1.seq)*100]
    plt.title('Sequence identity (%)')
    plt.bar(X,Y,color=(0.2, 0.4, 0.6, 0.6))
    plt.show()
    return (text1,text2)
     
def Analizar_Mult():
    global text1
    global text2
    global text3
    
    text1=nombreID1.get()
    text2=nombreID2.get()
    text3=nombreID3.get()
    print(text1,text2,text3)
  
    Seq1 = SeqIO.read(text1, 'fasta')
    Seq2 = SeqIO.read(text2, 'fasta')
    Seq3 = SeqIO.read(text3, 'fasta')
    
    Seq1_Seq3 = pairwise2.align.globalxx(Seq1.seq, Seq3.seq, one_alignment_only=True, score_only=True)
    print('Seq1/Seq3 Similarity (%):', Seq1_Seq3 / len(Seq1.seq) * 100)
    
    Seq2_Seq3 = pairwise2.align.globalxx(Seq2.seq, Seq3.seq, one_alignment_only=True, score_only=True)
    print('Seq2/Seq3 Similarity (%):', Seq2_Seq3 / len(Seq2.seq) * 100)
    
    Seq2_Seq1 = pairwise2.align.globalxx(Seq2.seq, Seq1.seq, one_alignment_only=True, score_only=True)
    print('Seq2/Seq1 Similarity (%):', Seq2_Seq1 / len(Seq1.seq) * 100)
    
    # Printeamos la data
    X = ['Seq1/Seq3', 'Seq2/Seq3', 'Seq2/Seq1']
    Y = [Seq1_Seq3/ len(Seq1.seq) * 100, Seq2_Seq3/ len(Seq2.seq)*100, Seq2_Seq1/len(Seq1.seq)*100]
    plt.title('Sequence identity (%)')
    plt.bar(X,Y,color=(0.2, 0.4, 0.6, 0.6))
    plt.show()
    
    return (text1,text2,text3)

def adn_aln_local():
    
    global text1
    global text2
    
    text1=nombreID1.get()
    text2=nombreID2.get()
    print(text1,text2)
    
    seq1 = SeqIO.read(text1, "fasta")
    seq2 = SeqIO.read(text2, "fasta")
    

    print("--------------------------------------------------------ADN ALIGNMENT-------------------------------------------------------------\n")

    alignments = pairwise2.align.localms(seq1.seq, seq2.seq, 1, -2, -2, -0.5)
    print(pairwise2.format_alignment(*alignments[0])) 
    
    return(text1,text2)

def protein_aln_local():
    
    global text1
    global text2
    
    text1=nombreID1.get()
    text2=nombreID2.get()
    print(text1,text2)
    
    seq1 = SeqIO.read(text1, 'fasta')
    seq2 = SeqIO.read(text2, 'fasta')


    print("--------------------------------------------------------PROTEIN ALIGNMENT-------------------------------------------------------------\n")

    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.localds(seq1.seq, seq2.seq, blosum62, -10, -1)
    print(pairwise2.format_alignment(*alignments[0]))  
    
    return(text1,text2)

def protein_aln_global():
    
    global text1
    global text2
    
    text1=nombreID1.get()
    text2=nombreID2.get()
    print(text1,text2)

    seq1 = SeqIO.read(text1, 'fasta')
    seq2 = SeqIO.read(text2, 'fasta')

    print("--------------------------------------------------------PROTEIN ALIGNMENT-------------------------------------------------------------\n")
    
    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    print(pairwise2.format_alignment(*alignments[0]))  
    
        
    return(text1,text2)

    
def ramachadran():
    
    global text1
       
    text1=nombreID1.get()
      
    confProDy(auto_show=False)
    confProDy(auto_secondary=True)
    p38 = parsePDB(text1,compressed=False)
    
    showContactMap(p38.chain_A)
    
    betas = p38.ca.getBetas()
    plot(betas);
    xlabel('Residue index')
    ylabel('B-factor')
    
    chain = p38['A']
    Phi = []; Psi = []; c = []
    for res in chain.iterResidues():
        try:
            phi = calcPhi(res)
            psi = calcPsi(res)
        except:
            continue
        else:
            Phi.append(phi)
            Psi.append(psi)
            if res.getResname() == 'GLY':
                c.append('black')
            else:
                secstr = res.getSecstrs()[0]
                if secstr == 'H':
                    c.append('red')
                elif secstr == 'G':
                    c.append('darkred')
                elif secstr == 'E':
                    c.append('blue')
                else:
                    c.append('grey')
    
    return(text1)

def phylotree(): 
    
    # Upload the new alignment file to your folder or working directory 
    # Open the alignment file as a MultipleSeqAlignment object 
    with open("variantes.aln","r") as aln: 
        alignment = AlignIO.read(aln,"clustal")
    print(type(alignment))

    # Open and initiate the Distance Calculator using the Identity model 
    from Bio.Phylo.TreeConstruction import DistanceCalculator 
    calculator = DistanceCalculator('identity')
    # Write the Distance Matrix 
    distance_matrix = calculator.get_distance(alignment)
    print(distance_matrix)

    # Open and initiate the Tree Constructor 
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    constructor = DistanceTreeConstructor(calculator)
    # Build the tree 
    variant_tree = constructor.build_tree(alignment)
    variant_tree.rooted = True
    print(variant_tree)

    # Save the tree to a new file 
    Phylo.write(variant_tree, "Cov_tree.xml", "phyloxml")

    # Import matplotlib and create a basic tree 
    import matplotlib
    import matplotlib.pyplot as plt
    variant_tree.clade[0, 0].color = "blue"
    variant_tree.root.color = "red"

    fig = Phylo.draw(variant_tree)

    # Make a better looking tree using the features of matplotlib 

    fig = plt.figure(figsize=(13, 5), dpi=100) # create figure & set the size 
    matplotlib.rc('font', size=12)              # fontsize of the leaf and node labels 
    matplotlib.rc('xtick', labelsize=10)       # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=10)       # fontsize of the tick labels
    #turtle_tree.ladderize()
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(variant_tree, axes=axes)
    fig.savefig("Covid_cladogram")
        
def MUSCLE():
    
    in_file = "C:/Users/hades/Coronalyzer/variantes.fasta"
    out_file = "C:/Users/hades/Coronalyzer/aligned.fasta"
    muscle_exe = "C:/Users/hades/Coronalyzer/muscle3.8.31_i86win32.exe"
    cline = MuscleCommandline(muscle_exe, input=in_file, phyiout=out_file)
    stdout, stderr = cline(in_file)
    align = AlignIO.read(out_file, "phylip")
    print(align)
    
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\AQUÍ EMPIEZ EL LAUNCHER TKINTER/////////////////////////////////////////////////
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\CORONALYZER 1.1.0////////////////////////////////////////////////////// 
root =tk.Tk()
root.iconbitmap("icon.ico")
root.title("Coronalyzer1.0")
root.geometry("700x650")
root.resizable(False,False)

text = tk.Text(root,bg="yellow", height=1)
text.pack()

text.insert('4.0', '------------------------BIENVENIDO A LA DEMO DE CORONALYZER---------------------')

menubar = tk.Menu(root)
root.config(bg="#0059b3", menu=menubar, width=300, height=300)
imagen = PhotoImage(file="corona.gif")
my_label=tk.Label(root, image=imagen)
my_label.place(x=0, y=0, relwidth=1, relheight=1,)
filemenu=tk.Menu(menubar, tearoff=0)
filemenu1=tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Sequence",command=Analizar_Seq)
filemenu.add_command(label="Compare 3 Seq.",command=Analizar_Mult)
filemenu.add_command(label="Compare 2 Seq.",command=Analizar_Sim)
filemenu.add_command(label="ADN Alignment",command=adn_aln_local)
filemenu.add_command(label="Mutations",command=Analizar_Mut)
filemenu.add_command(label="Protein Alignment local",command=protein_aln_local)
filemenu.add_command(label="Protein Alignment global",command=protein_aln_global)
filemenu.add_command(label="Protein Plot",command=ramachadran)
filemenu.add_command(label="MUSCLE",command=MUSCLE)
filemenu.add_command(label="Phylotree",command=phylotree)
filemenu.add_command(label="3D View",command=ramachadran)
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="Analysis", menu=filemenu)


miFrame=tk.Frame(root)
miFrame.pack()
miFrame.config(bg="grey")

miNombre1= tk.StringVar()
miNombre2= tk.StringVar()
miNombre3= tk.StringVar()
miNombre4=tk.StringVar()
miNombre5= tk.StringVar()

nombreID1=tk.Entry(miFrame, textvariable=miNombre1)
nombreID1.grid(row=1, column=1, padx=10, pady=10)

nombreID2=tk.Entry(miFrame, textvariable=miNombre2)
nombreID2.grid(row=2, column=1, padx=10, pady=10)

nombreID3=tk.Entry(miFrame, textvariable=miNombre3)
nombreID3.grid(row=3, column=1, padx=10, pady=10)

nombreID4=tk.Entry(miFrame, textvariable=miNombre4)
nombreID4.grid(row=4, column=1, padx=10, pady=10)

nombreID5=tk.Entry(miFrame, textvariable=miNombre5)
nombreID5.grid(row=5, column=1, padx=10, pady=10)

idLabel=tk.Label(miFrame, text="File 1:", width=10, height=1, borderwidth=3, relief="raised")
idLabel.grid(row=1, column=0, sticky="e", padx=10, pady=10)

idLabel=tk.Label(miFrame, text="File 2:", width=10, height=1, borderwidth=3, relief="raised")
idLabel.grid(row=2, column=0, sticky="e", padx=10, pady=10)

idLabel=tk.Label(miFrame, text="File 3:", width=10, height=1, borderwidth=3, relief="raised")
idLabel.grid(row=3, column=0, sticky="e", padx=10, pady=10)

idLabel=tk.Label(miFrame, text="Fasta File:", width=14, height=1, borderwidth=3, relief="raised")
idLabel.grid(row=4, column=0, sticky="e", padx=10, pady=10)

idLabel=tk.Label(miFrame, text="PDB File:", width=14, height=1, borderwidth=3, relief="raised")
idLabel.grid(row=5, column=0, sticky="e", padx=10, pady=10)

boton4=tk.Button(root, text="Clean all fields", command=limpiarCampos).pack()
boton3=tk.Button(root, text="Download Fasta", command=descargar_fasta).pack()
boton3=tk.Button(root, text="Download PDB", command=descargar_PDB).pack()


root.config(menu=menubar)

root.mainloop()