#!/bin/bash
R0=`./logscale 1 100 10`

i_r=0 # iterator for R0

for r in $R0
do
i_r=`expr $i_r + 1`

        folder=R0_$i_r
        mkdir -p $folder
        cd $folder

        qsub -cwd -V -b n ../single.sh $r
        cd ..
done

