#!/bin/bash

CC1NTDIR=/lustre/or-scratch/cades-bsd/9cq/CC1/CC1_LDOpenMMDask
while read line ; 
do 
    echo $line;
    IFS='/' read -r -a linearray <<< $line; 
    swarmdir="${linearray[0]}/${linearray[1]}/${linearray[2]}"; 
    cd $swarmdir
    cc1ntfname=$CC1NTDIR/$line
    echo $cc1ntfname
    mkdir -p Scattering
    cd Scattering
    vmd -dispdev text -eofexit -args $cc1ntfname  < ${CC1NTDIR}/ScatteringScripts/psfgen_13x4_orderfrag.tcl >> vmd1.txt
    cd ../../../../
 done < $1
