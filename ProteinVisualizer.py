'''
Solo funciona en Jupyter Notebooks 
import Bio
import pylab
import urllib
import pandas as pd
import nglview as nv
from Bio.PDB import PDBParser,MMCIFParser

p = PDBParser(PERMISSIVE=1) 
#find PDB file
urllib.request.urlretrieve('https://files.rcsb.org/download/6YYT.pdb', '6vxx.pdb')
parser = PDBParser()
structure = parser.get_structure('6VXX', '6vxx.pdb')
structure
view = nv.show_biopython(structure)
view
http://prody.csb.pitt.edu/_static/ipynb/workshop2020/prody_basic.html
https://medium.com/mlearning-ai/3d-sars-cov-2-protein-visualization-with-biopython-7c4f1955b1db
'''
#En el programa tendremos que implementar esto o abrirlo desde JN directamente
from prody import *
from numpy import *
from matplotlib.pyplot import *
%matplotlib inline
confProDy(auto_show=False)
confProDy(auto_secondary=True)
p38 = parsePDB('6vxx',compressed=False)
showProtein(p38)
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
                
scatter(Phi, Psi, c=c, s=10);
xlabel('Phi (degree)');
ylabel('Psi (degree)');