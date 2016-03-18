#!/bin/bash

cd ./../bin
SQM=./SQM

k=3
l=6
m=3
s=500.0


# Model
for M in 50 100 150; do
    for N in 30 50 75; do
	for p in 7 10 15 20; do
	    echo "### Start test cases with ###"
	    echo "${M} demand points"
	    echo "${N} potential sites"
	    echo "${p} servers"
	    for i in {01..10}; do
		echo "## Start Test case ${i} ##"
		prefix=Test_${M}_${N}_${p}_${i}
		options="-f$prefix -M$M -N$N -p$p -k$k -l$l -m$m -s$s"
		$SQM $options --brief Local_Search
		if [ $? -eq 139 ]; then
		    echo "It crashed!"
		    exit 1
		fi
	    done
	done
    done
done
