#!/bin/bash


for n in {1..16}
do
cd RandomConformersScattering/RandSep_ShiftCC1/R${n}
sbatch ../../../ScatteringScripts/run_ScatteringBundl_RandConf_RandSep_13x1.sh
cd ../../../
done 
