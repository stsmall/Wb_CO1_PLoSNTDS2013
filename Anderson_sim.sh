#!/bin/sh
#sims
for m in 1 100;do
    n=0
    while [ $n -lt 1 ]; do
        ms 140 1 -T -I 7 20 20 20 20 20 20 20 $m | tail +4 | grep -v // >treefile | seq-gen -mHKY -l 1000 -s .01 <treefile >seqfile
        head -1 seqfile >seqfile1;tail -n +2 seqfile | sort -n >>seqfile1_$m
		#python pairwise.py
		#Dnadist>>andersoncode
        let n=n+1 
    done
done

##data

#n=0
#while [ $n -lt 1000 ]; do
#    python subsample.py | python pairwise.py | ./snn n 1000 10 10 <sample.dat >>PN_ML_sample.out
#    let n=n+1
#    echo $n
#done
