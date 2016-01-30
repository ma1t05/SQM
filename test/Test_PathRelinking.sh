#!/bin/bash

SQM=./../git/SQM/bin/SQM

prefix=Test_01
M=50
N=50
p=10
k=3
l=6
m=3
s=500.0

options="-f$prefix -M$M -N$N -p$p -k$k -l$l -m$m -s$s"

# Model
for command in Local_Search; do
    echo "Test $command"
    $SQM $options --brief $command
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
