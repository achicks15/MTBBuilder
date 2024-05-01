#!/bin/python

import numpy as np
import pandas as pd
import mdtraj as md
from mdtraj.formats import PDBTrajectoryFile
import glob 

psffile = glob.glob("../*.psf")[0] ## Only works if there is only one psf file per directory

dcd = md.load('3jalfull_Hexagon.dcd', top=psffile)

frame1 = dcd.xyz[0]
frame2 = dcd.xyz[1]
frame3 = dcd.xyz[2]
frame4 = dcd.xyz[3]
frame5 = dcd.xyz[4]
frame6 = dcd.xyz[5]
frame7 = dcd.xyz[6]

combined_frames = np.concatenate([frame1, frame2, frame3, frame4, frame5, frame6, frame7])
#print(combined_frames.shape)

topology, bonds_micro = dcd.topology.to_dataframe()

ATubIndx = topology[~topology['segmentID'].str.find('AT').astype(bool)].index
BTubIndx = topology[~topology['segmentID'].str.find('BT').astype(bool)].index
CC1Indx = topology[~topology['segmentID'].str.find('C1').astype(bool)].index
GTPIndx = topology[~topology['segmentID'].str.find('GT').astype(bool)].index
MGIndx = topology[~topology['segmentID'].str.find('MG').astype(bool)].index

topology.loc[ATubIndx, 'chainID']='A'
topology.loc[BTubIndx, 'chainID']='B'
topology.loc[CC1Indx, 'chainID']='C'
topology.loc[GTPIndx, 'chainID']='G'
topology.loc[MGIndx, 'chainID']='M'

top2 = topology.copy()
top3 = topology.copy()
top4 = topology.copy()
top5 = topology.copy()
top6 = topology.copy()
top7 = topology.copy()

topology['segmentID'] = 'MT1'
top2['segmentID'] = 'MT2'
top3['segmentID'] = 'MT3'
top4['segmentID'] = 'MT4'
top5['segmentID'] = 'MT5'
top6['segmentID'] = 'MT6'
top7['segmentID'] = 'MT7'

PDBTrajectoryFile.set_chain_names(['A','B','C','G','M']*7)

top_microbundl = pd.concat([topology, top2, top3, top4, top5, top6, top7]).reset_index(drop=True)
bundltraj = md.Trajectory(combined_frames, topology=md.Topology.from_dataframe(top_microbundl))
bundltraj.save_pdb('3jalfull_Hexagon_mdtraj.pdb')
quit
