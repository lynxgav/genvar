#!/bin/bash
n=1000;
R0=$1;

for (( i=1;i<=$n;i+=1))
do
        ../genvar $i $R0 > out$i
done

