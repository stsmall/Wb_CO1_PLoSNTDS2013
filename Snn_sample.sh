#!/bin/sh
#sims
for m in .01 .1 1 10 20;do
    n=0
    while [ $n -lt 1000 ]; do
        ms 40 1 -T -I 2 20 20 $m | tail +4 | grep -v // >treefile | seq-gen -mHKY -l 1000 -s .015 <treefile >seqfile
        head -1 seqfile >seqfile1;tail -n +2 seqfile | sort -n >>seqfile1 | python pairwise.py
        ./snn n 1000 20 20 <sample.dat >>sample$m.out
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
