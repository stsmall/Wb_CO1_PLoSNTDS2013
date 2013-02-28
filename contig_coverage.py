#!/usr/bin/env python
contig_name=[]
contig_len=[]
A=[]
G=[]
T=[]
C=[]

newline=open("/Volumes/home/Users/Desktop/wuchereria_bancrofti_1_supercontigs.fasta","r")
contig_coverage=open("/Volumes/home/Users/Desktop/contig_len.txt","w")
for line in newline:
    if '>' in line:
        contig_name.append(line.split()[0])
        contig_coverage.write("%i\n" %(sum(A)+sum(G)+sum(C)+sum(A)))
        #contig_len.append(sum(A)+sum(G)+sum(C)+sum(A))
        A=[]
        G=[]
        T=[]
        C=[]
    else:
        A.append(line.count('A'))
        T.append(line.count('T'))
        G.append(line.count('G'))
        C.append(line.count('C'))
