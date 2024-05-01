#!/bin/bash


for n in {1..16}
do
#mkdir -p RandomConformersScattering/RandSep_ShiftCC1/R${n}
cd RandomConformersScattering/RandomSep_DotYZ/R${n}
sbatch ../../../ScatteringScripts/run_ScatteringBundl_RandConf_RandSep_13x1.sh
cd ../../../
done 
