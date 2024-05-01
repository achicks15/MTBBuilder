#!/bin/bash

#SBATCH -J MicroBundl_Scatt
#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --gres=gpu:0
#SBATCH --mem=0G
#SBATCH -t 48:00:00
#SBATCH -o ./%j-output.txt
#SBATCH -e ./%j-output.txt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=9cq@ornl.gov

module load vmd
module load anaconda3

source activate ~/anaconda_SassConMatch/SassConMatch 

SSDIR=/lustre/or-scratch/cades-bsd/9cq/CC1/CC1_LDOpenMMDask/ScatteringScripts
cd $SLURM_SUBMIT_DIR

/home/9cq/Pepsi-SAXS 3jalfull_13x4_gmoltype_4Xcc1nt.pdb -n 75 -ms 0.5 -au 1 -ns 250 -o 3jalfull_13x4_gmoltype_4Xcc1nt_SAXS.dat
/home/9cq/Pepsi-SANS 3jalfull_13x4_gmoltype_4Xcc1nt.pdb -o 3jall_13x4_gmoltype_4Xcc1nt_85DSans.dat -n 75 -ms 0.5 -au 1 -ns 250 --d2o 0.85 --deut 0.51 --deuterated 'Z'
/home/9cq/Pepsi-SANS 3jalfull_13x4_gmoltype_4Xcc1nt.pdb -o 3jalfull_13x4_gmoltype_4Xcc1nt_42DSans.dat -n 75 -ms 0.5 -au 1 -ns 250 --d2o 0.42 --deut 0.51 --deuterated 'Z'

for BUNDL in Hexagon
#Linear Triangle Square Diamond Trapezoid Hexagon
do
mkdir -p $BUNDL 
cd $BUNDL
vmd -dispdev text -args ../3jalfull_13x4_gmoltype_4Xcc1nt.psf ../3jalfull_13x4_gmoltype_4Xcc1nt.pdb < ${SSDIR}/make_${BUNDL}_randomSep.tcl > hex.out
python3 ${SSDIR}/split_dcd.py
/home/9cq/Pepsi-SAXS 3jalfull_${BUNDL}_mdtraj.pdb -o 3jalfull_${BUNDL}_SAXS_mdtraj.dat -n 75 -ms 0.5 -au 1 -ns 250
/home/9cq/Pepsi-SANS 3jalfull_${BUNDL}_mdtraj.pdb -o 3jalfull_${BUNDL}_SANS_mdtraj_85D_CC1NTd.dat -n 75 -ms 0.5 -au 1 -ns 250 --d2o 0.85 --deut 0.51 --deuterated 'C' 
/home/9cq/Pepsi-SANS 3jalfull_${BUNDL}_mdtraj.pdb -o 3jalfull_${BUNDL}_SANS_mdtraj_42D_CC1NTd.dat -n 75 -ms 0.5 -au 1 -ns 250 --d2o 0.42 --deut 0.51 --deuterated 'C'
cd ../
done
