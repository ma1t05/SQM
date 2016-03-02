#!/bin/bash

cd ./../bin
SQM=./SQM

prefix=Test_01
M=50
N=50
p=10
k=3
l=6
m=3
s=500.0


# Model
for M in {50,100,150}; do
    for N in {30,50,75}; do
	for p in {7,10,15,20}; do
	    for i in {01..05}; do
		prefix="Test_$M_$N_$p_$i"
		options="-f$prefix -M$M -N$N -p$p -k$k -l$l -m$m -s$s"
		$SQM $options --superbrief Path_Relinking
	    done
	done
    done
done

exit

echo "Test Model"
#$SQM $options model
# Berman Heuristic
echo "Test Berman Heuristic"
$SQM $options heuristic
# Local Search
echo "Test Local Search"
$SQM $options Local_Search
# Path Relinking
echo "Test Path Relinking"
$SQM $options --verbose Path_Relinking
# Random
echo "Test Random"
$SQM $options --verbose random
# GRASP
echo "Test GRASP"
$SQM $options --verbose GRASP
