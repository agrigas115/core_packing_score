#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:35:43 2019

@author: agrigas115
"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Core Packing Decoy Detector') 
parser.add_argument('-i', '--input', help='input pdb_file', required=True)  
args = parser.parse_args()

hydrophobic_residues = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE']

def KL(a, b):
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)

    a += 0.000001
    b += 0.000001
    return np.sum(a * np.log(a / b))

dun18_phobic_n = np.asarray([3.62562554e+00, 1.35715579e-01, 2.25915067e-01, 3.18384519e-01,
       3.83488309e-01, 4.92414990e-01, 6.23505334e-01, 7.84231301e-01,
       1.03739519e+00, 1.40941685e+00, 2.10367863e+00, 3.66541294e+00,
       6.33410086e+00, 2.74910468e+00, 1.15632714e-01, 8.73490667e-04,
       3.06000000e-04, 3.06000000e-04, 3.06000000e-04, 3.06000000e-04,
       3.06000000e-04, 3.06000000e-04, 3.06000000e-04, 3.06000000e-04])

dun18_phobic_bins = np.asarray([0.        , 0.04166667, 0.08333333, 0.125     , 0.16666667,
       0.20833333, 0.25      , 0.29166667, 0.33333333, 0.375     ,
       0.41666667, 0.45833333, 0.5       , 0.54166667, 0.58333333,
       0.625     , 0.66666667, 0.70833333, 0.75      , 0.79166667,
       0.83333333, 0.875     , 0.91666667, 0.95833333, 1.        ])


decoy_name = args.input
file = './output/'+decoy_name+'_packinglist_all.txt'
target_hydrophobic_phi = []
with open(file, 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
sourcelines.pop(0)
for line in sourcelines:
    restype = line.split()[1]
    if restype in hydrophobic_residues:
        packing_fraction = float(line.split()[2])
        target_hydrophobic_phi.append(packing_fraction)

(n, bins, patches) = plt.hist(target_hydrophobic_phi, bins=dun18_phobic_bins, density=True)
plt.close()
 	    
relative_S_decoy = (KL(n, dun18_phobic_n))
	           
with open(decoy_name+'_dkl.txt', 'w') as f:
    f.write(str(relative_S_decoy))
        
