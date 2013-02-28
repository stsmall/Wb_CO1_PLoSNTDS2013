#!/bin/bash

for nsam in 100 150 200 250 300;do



ms $nsam 1000 -t 12.40398 | msstats | cut -f 15 > sum2.txt

count=0
function Haplotype_Prob()
{
	number=$1
	reps=$2
	my_file=$3
	
	for line in $(cat $my_file); do
		if [ $line -ge $number ]; then
			count=$(($count+1))
		fi
	done
		
#prob=$(($count/1000))
echo $count >> file.txt

}

Haplotype_Prob 34 1000 sum2.txt

done