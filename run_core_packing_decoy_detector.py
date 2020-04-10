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

pdb_file = args.input

print('Beginning scoring on...'+pdb_file)

### Preprocessing ###
print('\t Preprocessing...')
with open('./preprocessing_code/pdb_list.txt', 'w') as f:
    f.write(pdb_file[:-4])


stream = os.popen('pwd')
pwd = stream.read()

os.system('python ./preprocessing_code/download_preprocess_pdb.py ./ ./preprocessing_code/pdb_list.txt 1 '+pwd+'/preprocessing_code/')

