#!/bin/bash

cd ./../bin
SQM=./SQM

M=50
N=50
p=10
k=3
l=6
m=3
s=500.0


for SS in {01..07} do
for M in 50 100 150; do
    for N in 30 50 75; do
	for p in 7 10 15 20; do
	    for i in {01..05}; do
		prefix=Test_${M}_${N}_${p}_${i}
		options="-f$prefix -M$M -N$N -p$p -k$k -l$l -m$m -s$s"
		echo "$SQM $options --brief Scatter_Search_${SS}"
		$SQM $options --brief Scatter_Search_${SS}
	    done
	done
    done
done

exit