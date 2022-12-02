#import Biopython
import Bio as Bio
#import our three important packages
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from pylab import *
#import the sequences we will use. These are 16s sequences from GenBank
#example: https://www.ncbi.nlm.nih.gov/nuccore/FJ039971.1?report=genbank 
t1 = Bio.SeqIO.read("sars.fasta", "fasta")
t2 = Bio.SeqIO.read("cov2.fasta", "fasta")
t3 = Bio.SeqIO.read("mers.fasta", "fasta")
t4 = Bio.SeqIO.read("alpha.fasta", "fasta")
t5 = Bio.SeqIO.read("Delta.fasta", "fasta")
t6 = Bio.SeqIO.read("HKU1.fasta", "fasta")
t7 = Bio.SeqIO.read("Gamma.fasta", "fasta")
t8 = Bio.SeqIO.read("Beta.fasta", "fasta")
t9 = Bio.SeqIO.read("229E.fasta", "fasta")


#rename each of the sequences 
#this step is not required, it will just make the tree easier to understand 
print(t1.description)
print(t2.description)
print(t3.description)


t1.id = 'Sars-CoV1'
t2.id = 'Sars-CoV2'
t3.id = 'Mers'
t4.id = 'Alpha'
t5.id = 'Delta'
t6.id = 'HKU1'
t7.id = 'Gamma'
t8.id = 'Beta'
t9.id = '229E'



# Combine all of the individual sequences into a new file 
variantes = SeqIO.write([t1,t2,t3,t4,t5,t6,t7,t8,t9], "variantes.fasta", "fasta")
# Load the turtles sequences into MUSCLE 
#https://www.ebi.ac.uk/Tools/msa/muscle/
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
matplotlib.rc('font', size=15)              # fontsize of the leaf and node labels 
matplotlib.rc('xtick', labelsize=10)       # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=10)       # fontsize of the tick labels
turtle_tree.ladderize()
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(variant_tree, axes=axes)
fig.savefig("Covid_cladogram")

