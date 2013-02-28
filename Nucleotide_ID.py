#!/usr/bin/env python
#Nucleotide ID, Amino acid ID, substitutions per codon site, Dn/Ds
from Bio import SeqIO, Seq, Entrez,AlignIO
import csv,sys
from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
import difflib
import os

NT_Percent_diff=[]
AA_Percent_same=[]
Ts_Tv=[]


path = raw_input("Enter path to files: ")
listing = os.listdir(path)

for infile in listing:
    if "fasta" in infile:
        handle=AlignIO.read(infile,"fasta")
        for j in range(1,3):
##Nucleotide Diff
            count=0
            for i in range(0, len(handle[0].seq)):
                if handle[0].seq[i] != handle[j].seq[i]:
                    if handle[j].seq[i]!='N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-':
                        count=count +1
                    
            NT_Percent_diff.append(float(count)/len(handle[0].seq))

#Ts/Tv
            Ts_Tv.append(infile)
            Ts=0
            Tv=0
            for i in range(0, len(handle[0].seq)):
                if handle[0].seq[i] != handle[j].seq[i]:
                    if handle[j].seq[i]!='N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-':
                        if handle[0].seq[i]=='A' and handle[j].seq[i]=='G':
                            Ts=Ts+1
                            print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                        elif handle[0].seq[i]=='G' and handle[j].seq[i]=='A':
                            Ts=Ts+1
                            print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                        elif handle[0].seq[i]=='C' and handle[j].seq[i]=='T':
                            Ts=Ts+1
                            print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                        elif handle[0].seq[i]=='T' and handle[j].seq[i]=='C':
                            Ts=Ts+1
                            print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                        else:
                            Tv=Tv+1
            if Tv==0:
                rat_sv="%i:0" %Ts
            else:
                rat_sv=float(Ts)/float(Tv)        
            Ts_Tv.append(rat_sv)       


