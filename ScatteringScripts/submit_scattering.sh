#!/bin/bash

CC1NTDIR=/lustre/or-scratch/cades-bsd/9cq/CC1/CC1_LDOpenMMDask
while read line ; 
do 
    echo $line;
    IFS='/' read -r -a linearray <<< $line; 
    swarmdir="${linearray[0]}/${linearray[1]}/${linearray[2]}"; 
    cd $swarmdir/Scattering
    #cc1ntfname=$CC1NTDIR/$line
    #echo $cc1ntfname
    sbatch $CC1NTDIR/ScatteringScripts/run_ScatteringBundl.sh
    cd ../../../../
 done < $1
