#!/bin/bash 

for n in {1..16}
do
mkdir -p "RandomConformersScattering/RandomSep_DotYZ/R${n}"
echo $n 
vmd -dispdev text -eofexit -args "CrossLinkREST_E2EVectorScoring_TYRMinDistLt75_perpY_parZ_Scoregt15_pathlist.csv" "RandomConformersScattering/RandomSep_DotYZ/R${n}"< ScatteringScripts/psfgen_13x1_orderfrag_randomConf.tcl > RandomConformersScattering/RandomSep_DotYZ/R${n}/vmd.out
done
   
