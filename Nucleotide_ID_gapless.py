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

sub_site=[]
Dn_Ds=[]
path = '/Volumes/home/Users/stsmall/Komodo_scripts/seqs/gapless_sub_site'
listing = os.listdir(path)

for infile in listing:
    if "fasta" in infile:
        #print infile
        handle=AlignIO.read(infile,"fasta")
        for j in range(1,3):
#codon site substitutions
            sub_site.append(infile)
            a=handle[0].seq
            b=handle[j].seq
            for i in range(0,3):
                count=0
                x=a[i::3]
                y=b[i::3]
                for k in range(0, len(x)):
                    if x[k] != y[k]:
                        if y[k]!='N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-':
                            count=count +1                        
                sub_site.append(count)

#Dn/Ds
            Dn_Ds.append(infile)
            c=handle[0].seq
            d=handle[j].seq
            Dn=0
            Ds=0
            m=0
            p=3
            while p<len(c):
                w=c[m:p]
                z=d[m:p]
                #print "%s:%s" %(w,z)
                if str(w)!=str(z):  #triplet comparison
                    #print "%s:%s" %(w.translate(table=5),z.translate(table=5))
                    if str(w.translate(table=5))==str(z.translate(table=5)):   #translated triplet
                        Ds=Ds+1
                     #   print "Ds=%i" %Ds
                    else:       
                        Dn=Dn+1
                     #   print "Dn=%i" %Dn
                m=m+3
                p=p+3
            ratio=float(Dn)/float(Ds)    
            Dn_Ds.append("%i:%i=%f" %(Dn,Ds,ratio))
            