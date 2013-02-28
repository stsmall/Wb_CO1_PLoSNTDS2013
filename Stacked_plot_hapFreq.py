#!/usr/bin/env python
from pandas import *
import numpy

x=open("/Volumes/home/Users/stsmall/Desktop/hap_freq_stacked.csv",'r')
newline=x.readline().split('\r')
Haps=[]
for line in newline:
    enteries=line.split(",")
    enteries[0]   #hap name
    PN=0
    A2=0
    A1=0
    Y1=0
    Y2=0
    ML=0
    MO=0
    for items in enteries:
        if "PN" in items:
            PN=PN+1
        if "A2" in items:
            A2=A2+1
        if "A1" in items:
            A1=A1+1
        if "Y1" in items:
            Y1=Y1+1
        if "Y2" in items:
            Y2=Y2+1
        if "MO" in items:
            MO=MO+1
        if "ML" in items:
            ML=ML+1
    Haps.append('%s,%i,%i,%i,%i,%i,%i,%i' %(enteries[0],PN,A2,A1,Y1,Y2,MO,ML))
f_s=open('Haps_freq.txt','w')
for item in Haps:
    f_s.write("%s\n" % item)
f_s.close()
x.close()