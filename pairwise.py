#!/usr/bin/env python
from Bio import SeqIO, Seq, AlignIO
import csv,sys
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
import difflib

def pairwise(seqfile):
    pair_diff=[]
    sample=open("/Volumes/home/Users/stsmall/Desktop/snn_sample/sample.dat","w")
    line1="OTUs\t"
    sample.writelines(line1)
    for seq in seqfile:
        sample.write("%s " %seq.id)
    for i in range(0,(len(seqfile)-1)):
        for j in range(i+1,len(seqfile)):
            seq1=seqfile[i].seq
            seq2=seqfile[j].seq
            dist = sum((cSeq1 != cSeq2 for cSeq1, cSeq2 in zip(seq1,seq2)))
            pair_diff.append(dist)
        line="\n%s\t" %seqfile[i].id
        sample.writelines(line)
        for item in pair_diff:
            sample.write("%i " %item)
        pair_diff=[]         
    sample.close()

pairwise(AlignIO.read("/Volumes/home/Users/stsmall/Desktop/snn_sample/seqfile1","phylip"))
