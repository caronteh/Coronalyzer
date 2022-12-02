# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 00:55:55 2021

@author: hades
"""

#En el programa tendremos que implementar esto o abrirlo desde JN directamente
from prody import *
from numpy import *
from matplotlib.pyplot import *
%matplotlib inline
confProDy(auto_show=False)
confProDy(auto_secondary=True)
p38 = parsePDB('6vxx',compressed=False)
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