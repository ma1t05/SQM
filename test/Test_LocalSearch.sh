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


for prefix in Test_01 Test_02; do
    options="-f$prefix -M$M -N$N -p$p -k$k -l$l -m$m -s$s"
    $SQM $options --brief Local_Search
done

