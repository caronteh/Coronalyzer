# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 14:23:54 2022

@author: hades
"""


from Bio.Align.Applications import MuscleCommandline
in_file = "C:/Users/hades/Coronalyzer/variantes.fasta"
out_file = "C:/Users/hades/Coronalyzer/aligned.fasta"
muscle_exe = "C:/Users/hades/Coronalyzer/muscle3.8.31_i86win32.exe"
cline = MuscleCommandline(muscle_exe, input=in_file, phyiout=out_file)
stdout, stderr = cline(in_file)
from io import StringIO
from Bio import AlignIO
align = AlignIO.read(out_file, "phylip")
print(align)





