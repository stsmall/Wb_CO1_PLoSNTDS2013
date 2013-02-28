#!/usr/bin/env python
from pandas import *
import numpy

newline = [line.strip() for line in open("/Volumes/home/Users/stsmall/Desktop/Wuchereria bancrofti/Wb_CO1/Analysis/Pairwise_GD/outfile_All",'r')]
dna_dist={}
index=[]

for line in newline:
    if 'T' in line:
        index.append(line.split()[0])

for line in newline:
    for name in line.split():
        if 'T' in name:
            key=name
            list=[]
        else:
            list.append(name)
    dna_dist[key]=list

within=[]
among=[]
same=[]
for keys in dna_dist:
    for k in range(0, len(dna_dist.keys())):
        if keys==index[k]:
            print 'same %s=%s' %(keys,index[k])
            same.append(dna_dist[keys][k])
        elif keys[1:5] in index[k]:
            within.append(dna_dist[keys][k])
            #print 'within %s=%s' %(keys,index[k])
        elif keys[1:5] not in index[k]:
            among.append(dna_dist[keys][k])
            #print 'among %s=%s' %(keys,index[k])
                
within=Series(within)
among=Series(among)
among.to_csv('Anderson_among.txt')
within.to_csv('Anderson_within.txt')

#np.asarray(within).T  #transpose
