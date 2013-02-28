#!/usr/bin/env python
"""All modules are available in the python package from www.enthought.com"""
from __future__ import division
from Bio import SeqIO, Seq, AlignIO
import csv,sys
#from Bio.Blast import NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
import difflib
import os
import numpy as np
from itertools import izip
from pandas import *

#1calculates pi per site -->Pi=[]
def pi(alignment,seqs,bases):
    Pi=[]
    for k in range(0,bases):
        m=0
        for i in range(0,seqs):
            for j in range(i+1,seqs):
                if alignment[i,k]!=alignment[j,k]:
                    m=m+1
        Pi.append(m)        
    X=Series(Pi)
    pi=X/(seqs*(seqs-1)/2.)
    return (pi)
    
#2calculates %bp freq per site --> freq_A=[] freq_T=[] freq_C=[] freq_G=[]
def freqbp(alignment,seqs,bases):
    freq_A=[]
    freq_T=[]
    freq_C=[]
    freq_G=[]
    for k in range(0,bases):
        freq_A.append(alignment[:,k].count('A'))
        freq_T.append(alignment[:,k].count('T'))
        freq_G.append(alignment[:,k].count('G'))
        freq_C.append(alignment[:,k].count('C'))
    freq={'freqA':freq_A,'freqT':freq_T,'freqG':freq_G,'freqC':freq_C}
    df_freq=DataFrame(freq)
    print "**%s**\nA=%s\nT=%s\nG=%s\nC=%s\n" %(infile,sum(freq_A)/(bases*seqs),sum(freq_T)/(bases*seqs),sum(freq_G)/(bases*seqs),sum(freq_C)/(bases*seqs))
    return (df_freq)
    
#3indexes sites where there are >2 mutations, violates infinite sites -->infinite_sites=[]
def infinite_sites(bases,df_freq):
    infinite_sites=[]
    for k in range(0,bases):
        if list(df_freq.ix[k]).count(0)<2:
            infinite_sites.append(k)
    print "Infinite site violations=%i\n" %(len(infinite_sites))
    return (infinite_sites)

#4indexes where there are sinlgton, doubletons...freq distribution -->single_SNP=[]
def SNP_freq(bases,df_freq):
    snp_freq=[]
    for S in range(1,int(round(seqs/2))):
        single_SNP=0
        for k in range(0,bases):
            for j in df_freq.ix[k]:
                if j == S:
                    single_SNP=single_SNP+1
        snp_freq.append(single_SNP)
    #plot
    return snp_freq

#5coding position subs-->sub_site=[]
def CDS_subs(alignment,seqs,bases):
    ref=alignment[0].seq #make first seq reference
    sub_site=[]
    for k in range(0,3):
        site_count=0
        while k < bases:
            if 'N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-' not in alignment[1:,k]: 
                site_count=site_count + ((seqs-1)-(alignment[1:,k].count(ref[k])))
            k=k+3
        sub_site.append(site_count)
    print "CDS substitutions for %s=%s\n\n" %(infile,sub_site)
    return (sub_site)

#6Ratio of non-synonomous substitution to synonomous substitutions per gene-->Dn_Ds=[]
def Dn__Ds(alignment,seqs,bases):
    Dn_Ds=[]
    ref=alignment[0].seq
    dn=0
    ds=0
    m=0
    p=3
    while p<bases:
        for seq in alignment:
            if str(ref[m:p])!=str(seq.seq[m:p]):  #triplet comparison
                if str(ref[m:p].translate(table=5))==str(seq.seq[m:p].translate(table=5)):   #translated triplet
                    ds=ds+1
                else:       
                    dn=dn+1
        m=m+3
        p=p+3
    ratio=dn/ds
    Dn_Ds.append("%i:%i=%f" %(dn,ds,ratio))
    print "Dn:Ds ratio for %s = %s\n\n" %(infile,Dn_Ds)
    return (Dn_Ds)


#7Transition : transversion rate per site for each gene-->Ts_Tv=[]
def TS_TV(alignment,seqs,bases):
    Ts_Tv=[]
    ref=alignment[0].seq
    Ts=0
    Tv=0
    for k in range(0,bases): 
        for seq in alignment:
            if ref[k] != seq.seq[k]:
                if seq.seq[k]!='N' or 'W' or 'Y' or 'R' or 'K' or 'D' or '-':
                    if ref[k]=='A' and seq.seq[k]=='G':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='G' and seq.seq[k]=='A':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='C' and seq.seq[k]=='T':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    elif ref[k]=='T' and seq.seq[k]=='C':
                        Ts=Ts+1
                        #print "%s,%s" %(handle[0].seq[i],handle[j].seq[i])
                    else:
                        Tv=Tv+1
    if Tv==0:
        rat_sv="-->%i:0" %(Ts)
    else:
        rat_sv=(Ts)/(Tv)        
    Ts_Tv.append("%s: Ts:%s,Tv:%s,Ts/Tv:%s" %(infile,Ts,Tv,rat_sv))       
    print "%s\n\n" %Ts_Tv
    return (Ts_Tv)

#8filtering out non-synonymous changes
def DNfilter(alignment,seqs,bases):
    DN_site=[]
    posit=[]
    ref=alignment[0].seq
    numb_changes=0
    m=0
    p=3
    while p<bases+1:
        for seq in alignment:
            if str(ref[m:p])!=str(seq.seq[m:p]): #if the seq is the same it doesnt matter
                if str(ref[m:p].translate(table=5))!=str(seq.seq[m:p].translate(table=5)):  #if the seq is diff does the AminoAcid change?
                    numb_changes=numb_changes+1 #freq
                    posit=[i for i, (a1,a2) in enumerate(izip(str(ref[m:p]),str(seq.seq[m:p]))) if a1!=a2] #where is the change?
        if numb_changes>0:
            DN_site.append("%i,%i,%i" %(m+1,posit[0]+1,numb_changes))  #x is bp position, number changes is freq out of 489, posit is 1st, 2nd, or 3rd
        if len(posit)>1:
            print m+1
            print posit
        numb_changes=0
        posit=[]
        m=m+3
        p=p+3
    return (DN_site)
    
#9number of diff per site
def site_freq(alignment,seqs,bases):
    ref=alignment[0].seq
    site_freq=[]
    site=0
    for k in range(0,bases): 
        for seq in alignment:
            if ref[k] != seq.seq[k]:
                site=site+1
        site_freq.append(site)
        site=0
    return (site_freq)

############IN PROGRESS###################################

#Polymerase error correction
    #corrects snps in XX Freq back to concensus
    #corrects pairwise distances at other sites by randomly increasing or decreasing the hamming distance
    #########

#corrects for infinite sites violations (non recombining only)
    ######### 4 gamete test, DNAdist
"""assumes non-recombining loci"""

#Nucleotide difference
#difflib.SequenceMatcher(None,string1,string2).ration() #s=str(seq.seq[0:200])
def levenshtein(a,b):
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]

######################XXXXXXXXXXXXXXX###########################

check=raw_input("Please check the following:\n 1) all files are in fasta format\n 2) contain an alignment\n 3) are all the same length with no gaps\n When ready to start enter 'Y': ")

if check=="Y" or "y":
    in_path = raw_input("Enter path to input files:")
    out_path = raw_input("Enter where to save out files:")
    if out_path=='same':
        out_path=in_path    
    listing = os.listdir(in_path)
    #try/except if not fasta or aligned
    for infile in listing:
        if ".phy" in infile and ".txt" not in infile:
            alignment=AlignIO.read(in_path+'/'+infile,"phylip")
            bases=alignment.get_alignment_length()
            seqs=len(alignment)
            analysis=raw_input("\n\nPlease specify which analyses to run:\n 1:pi per site\n 2:percent nucleotide freq\n 3:infinite site violations\n 4:SNP frequency spectrum\n 5:Coding site substitutions\n 6:Dn/Ds\n 7:Ts:Tv\n 8:DNfiltering\n 9:SNPfreq\n 10:quit\n **Enter as a continous number, e.g. 12345678 will run all**\n run:")
            print "\n\n gene: %s\n bases: %i\n Sequences: %i\n" %(infile,bases,seqs)
            if '1' in analysis:
                x=pi(alignment,seqs,bases)
                x.to_csv('%s%s_Pi.txt' %(out_path,infile))
            if '2' in analysis:
                freqbp(alignment,seqs,bases).to_csv('%s%s_freq_bases.txt' %(out_path,infile))
            if '3' in analysis:
                infinite_site=infinite_sites(bases,freqbp(alignment,seqs,bases))
                if len(infinite_site)!=0:
                    f_i=open('%s%s_Infinite_sites.txt' %(out_path,infile),'w')
                    for item in infinite_site:
                        f_i.write("%s\n" % item)
                    f_i.close()
            if '4' in analysis:
                snp_freq=SNP_freq(bases,freqbp(alignment,seqs,bases))
                f_s=open('%s%s_SNPs_freq.txt' %(out_path,infile),'w')
                for item in snp_freq:
                    f_s.write("%s\n" % item)
                f_s.close()
            if '5' in analysis:
                CDS=CDS_subs(alignment,seqs,bases)
                f_s=open('%s%s_CDS_subs.txt' %(out_path,infile),'w')
                for item in CDS:
                    f_s.write("%s\n" % item)
                f_s.close()
            if '6' in analysis:
                dnds=Dn__Ds(alignment,seqs,bases)
                f_s=open('%s%s_Dn__Ds.txt' %(out_path,infile),'w')
                for item in dnds:
                    f_s.write("%s\n" % item)
                f_s.close()
            if '7' in analysis:
                tstv=TS_TV(alignment,seqs,bases)
                f_s=open('%s%s_TS_TV.txt' %(out_path,infile),'w')
                for item in tstv:
                    f_s.write("%s\n" % item)
                f_s.close()
            if '8' in analysis:
                Dn_site=DNfilter(alignment,seqs,bases)
                if len(Dn_site)!=0:
                    f_i=open('%s%s_DNfiltering.txt' %(out_path,infile),'w')
                    for item in Dn_site:
                        f_i.write("%s\n" % item)
                    f_i.close()
            if '9' in analysis:
                snp_freq=site_freq(alignment,seqs,bases)
                f_s=open('%s%s_SNPs_freq.txt' %(out_path,infile),'w')
                for item in snp_freq:
                    f_s.write("%s\n" % item)
                f_s.close()
            if '10' in analysis:
                break