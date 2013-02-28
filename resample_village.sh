#!/bin/sh
#resamples from pooled village file for 5,10 times

times=5 #change j,sed,snn
n=0
while [ $n -lt 1000 ]; do
    for i in Peneng.phy Albulum1.phy Albulum2.phy Yautong1.phy Yautong2.phy Molienge.phy Moihauk.phy Musunguwah.phy; do
        for j in 1 2 3 4 5; do
            s=$( wc -l <$i )  
            rnumber=$((RANDOM%$s+2))  
            cat $i | head -n $rnumber | tail -n 1  >>seqfile
        done
    done
    #phylip_len=$( wc -l <seqfile ) # add header for phylip 
    sed -e "1i\\
        40 450" seqfile >seqfile1.phy

    python pairwise.py

    ./snn 8 1000 5 5 5 5 5 5 5 5 <sample.dat >>sample.out

    cat sample.out | grep Snn | sed 's:(::' | sed 's:)::' >sample_village_unique_$times.txt
    rm seqfile
    rm seqfile1.phy
    rm sample.out
    rm sample.dat
    let n=n+1
done
