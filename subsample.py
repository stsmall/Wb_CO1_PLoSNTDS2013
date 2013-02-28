#!/usr/bin/env python
from Bio import SeqIO, Seq, AlignIO
import csv,sys
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
import random

def subsample(sequences):
    #PN,A2,A1,Y1,Y2,MO,ML=69,134,149,23,42,59,13
    PN_data={}
    A2_data={}
    A1_data={}
    Y1_data={}
    Y2_data={}
    MO_data={}
    ML_data={}
    my_data={}
    ids=[]
    for seq in sequences:
        ids.append(seq.id)
        my_data[seq.id]=seq.seq
        if "PN" in seq.id:
            PN_data[seq.id]=seq.seq
        if "A1" in seq.id:
            A1_data[seq.id]=seq.seq
        if "A2" in seq.id:
            A2_data[seq.id]=seq.seq            
        if "Y1" in seq.id:
            Y1_data[seq.id]=seq.seq
        if "Y2" in seq.id:
            Y2_data[seq.id]=seq.seq            
        if "ML" in seq.id:
            ML_data[seq.id]=seq.seq            
        if "MO" in seq.id:
            MO_data[seq.id]=seq.seq            
    POP1=PN_data
    POP2=ML_data
    pop1=random.sample(POP1,10) #without replacement
    pop2=random.sample(POP2,10) #without replacement
    POPtotal=pop1+pop2
    my_alignment=[]
    for item in POPtotal:
        for seq in sequences:
            if item==seq.id:
                my_alignment.append(seq)
                break
    alignment=MultipleSeqAlignment(my_alignment)            
    
    #def sample_wr(population, k):
    #    "Chooses k random elements (with replacement) from a population"
    #    n = len(population)
    #    _random, _int = random.random, int  # speed hack 
    #    result = [None] * k
    #    for i in xrange(k):
    #        j = _int(_random() * n)
    #        result[i] = population[j]
    #    return result

    return AlignIO.write(alignment, "/Volumes/home/Users/stsmall/Desktop/snn_sample/seqfile", "phylip")


subsample(AlignIO.read("/Volumes/home/Users/stsmall/Desktop/snn_sample/Wb_assembly_All_patients.phy","phylip"))
