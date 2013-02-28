#!/bin/bash
#resamples from individual file for 1,2,5,10 times

times=1 #change j, sed, snn
n=0
while [ $n -lt 10 ]; do
    for i in T0059PN.phy T0083PN.phy T0097PN.phy T0150A1.phy T0186A1.phy T0388A1.phy T0145A2.phy T0142A2.phy T0346A2.phy T0363Y1.phy T0557Y2.phy T0582Y2.phy T0609ML.phy T1602ML.phy T1384MO.phy T1358MO.phy T1019MS.phy; do
        for j in 1; do
            s=$( wc -l <$i )  
            rnumber=$((RANDOM%$s+1))  
            cat $i | head -n $rnumber | tail -n 1  >>seqfile  
        done
    done
    sed -e "1i\\
        17 450" seqfile >seqfile1.phy
    
    python pairwise.py

PN=$((3*times))
A1=$((3*times))
A2=$((3*times))
Y1=$((1*times)) 
Y2=$((2*times))
ML=$((2*times))
MO=$((2*times))
MS=$((1*times))
./snn 8 1000 3 3 3 1 2 2 2 1 <sample.dat >sample.out
cat sample.out | grep Snn | sed 's:(::' | sed 's:)::' >>sample_patient_original_$times.txt
    rm seqfile
    rm seqfile1.phy
    rm sample.out
    rm sample.dat
    let n=n+1
done



