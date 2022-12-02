# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 00:00:54 2022

@author: hades
"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

sequence_data = open("blast_example.fasta").read() 
result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data) 
with open('results.xml', 'w') as save_file: 
  blast_results = result_handle.read() 
  save_file.write(blast_results) 
  
blast_records = NCBIXML.parse(result_handle)

E_VALUE_THRESH = 1e-20 

for record in NCBIXML.parse(open("results.xml")): 
    if record.alignments: 
       print("\n") 
       print("query: %s" % record.query[:100]) 
       for align in record.alignments: 
          for hsp in align.hsps: 
             if hsp.expect < E_VALUE_THRESH: 
                print("****Alignment****")
                print("sequence:", align.title)
                print("length:", align .length)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                print()