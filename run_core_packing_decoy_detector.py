#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:32:20 2020

@author: agrigas115
"""
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Core Packing Decoy Detector') 
parser.add_argument('-i', '--input', help='input pdb_file', required=True)  
args = parser.parse_args()

matlab_root = '/Applications/MATLAB_R2019a.app/bin/matlab'

pdb_file = args.input

print('Beginning scoring on...'+pdb_file)

### Preprocessing ###
print('\t Preprocessing...')
with open('./preprocessing_code/pdb_list.txt', 'w') as f:
    f.write(pdb_file[:-4])


stream = os.popen('pwd')
path = stream.read()
path = path[:-1]

preprocess_step1 = 'python ./preprocessing_code/download_preprocess_pdb.py '+path+'/ ./preprocessing_code/pdb_list.txt 1 '+path+'/preprocessing_code/'
os.system(preprocess_step1)
os.chdir('./preprocessing_code')
preprocess_step2 = 'bash  ./bash_edge_code_script_local.sh pdb_list.txt '+path+'/ 1 ./vol_code ./vol_code'
os.system(preprocess_step2)
print('\t Calculating rSASA...')
os.chdir('../rSASA_code/')
rsasa_command = matlab_root+' -nosplash -nodisplay -nojvm -r "database_sasa"'
os.system(rsasa_command)

print('\t Calculating Voronoi Volumes...')

print('\t Calculating Residue Volumes...')

print('\t Calculating D_kl...')

print('\t Running SNN...')

print('\t Cleaning up...')