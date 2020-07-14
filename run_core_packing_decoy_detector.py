#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 15:32:20 2020

@author: agrigas115
"""
import os
import subprocess
import argparse
import sys

import tensorflow as tf
tf.get_logger().setLevel('INFO')
tf.autograph.set_verbosity(1)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
from tensorflow.keras.models import model_from_json

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

try:
    import numpy as np
except ModuleNotFoundError:
    print('We need numpy. If you use anaconda, try: conda install numpy')
    raise
try:
    import matlab.engine
except ModuleNotFoundError or OSError:
    print('To run matlab scripts from python we need to install the matlab python enginge \n cd matlabroot/extern/engines/python \n python setup.py install')
    raise

parser = argparse.ArgumentParser(description='Core Packing Decoy Detector') 
parser.add_argument('-i', '--input', help='input pdb_file', required=True)  
args = parser.parse_args()

matlab_root = '/Applications/MATLAB_R2019a.app/bin/matlab'

pdb_file = args.input
struc_name = pdb_file[:-4].lower()


print('Beginning scoring on... '+pdb_file)

### Preprocessing ###
print('\t Preprocessing...')
with open('./preprocessing_code/pdb_list.txt', 'w') as f:
    f.write(pdb_file[:-4])

stream = os.popen('pwd')
path = stream.read()
path = path[:-1]

command = ['python', './preprocessing_code/download_preprocess_pdb.py', path+'/', './preprocessing_code/pdb_list.txt', '1', path+'/preprocessing_code/']
with open(os.devnull, "w") as fnull:
    result = subprocess.call(command, stdout = fnull, stderr = fnull)

os.chdir('./preprocessing_code')
command = ['bash', './bash_edge_code_script_local.sh', 'pdb_list.txt', path+'/', '1', '../', '../']
with open(os.devnull, "w") as fnull:
    result = subprocess.call(command, stdout = fnull, stderr = fnull)

#########################################
print('\t Calculating rSASA...')
os.chdir('../rSASA_code/')
eng = matlab.engine.start_matlab()
eng.database_sasa(nargout=0)
#########################################
print('\t Calculating Voronoi Volumes...')
os.chdir('../voro_code/')
# Generate Surface
eng = matlab.engine.start_matlab()
eng.generate_surface('../'+struc_name,nargout=0)

# Run pomelo
command = ['./bin/pomelo', '-mode', 'GENERIC', '-i', '../'+struc_name+'_param.lua', '-o', './result']
with open(os.devnull, "w") as fnull:
    result = subprocess.call(command, stdout = fnull, stderr = fnull)

# Save output
copy_command = 'cp ./result/setVoronoiVolumes.dat ./'+struc_name+'_vor.txt' 
os.system(copy_command)

print('\t Finding Open Voronoi Cells...')
# Generate surface with larger box size
eng = matlab.engine.start_matlab()
eng.generate_surface_10('../'+struc_name,nargout=0)

# Run pomelo
command = ['./bin/pomelo', '-mode', 'GENERIC', '-i', '../'+struc_name+'_param_10.lua', '-o', './result_10']
with open(os.devnull, "w") as fnull:
    result = subprocess.call(command, stdout = fnull, stderr = fnull)

# Save output
copy_command = 'cp ./result_10/setVoronoiVolumes.dat ./'+struc_name+'_vor_10.txt' 
os.system(copy_command)
#########################################
print('\t Calculating Residue Volumes... this is the longest step :(')

os.chdir('../vol_code/')
# =============================================================================
# with open('../preprocessing_code/tasklist.sh', 'r') as f:
#     info = f.read()
# sourcelines = info.splitlines()
# vol_command = sourcelines[0].split()
# with open(os.devnull, "w") as fnull:
#     result = subprocess.call(vol_command, stdout = fnull, stderr = fnull)
# =============================================================================

#########################################
print('\t Overlap Energy...')
os.chdir('../')
command  = 'python calc_single_overlap_energy.py -i '+struc_name
os.system(command)    

#########################################
print('\t Calculating Packing Fraction...')
command  = 'python calc_packing_fraction.py -i '+struc_name
os.system(command)

#########################################
print('\t Calculating D_kl...')
command  = 'python calc_kl_div_packing_hydro.py -i '+struc_name
os.system(command)

#########################################
print('\t Loading Data...')
features = ['number of residues', 'fc 10**-3', 'fc 10**-2', 'fc 10**-1' ,'<phi> 10**-3','sigma(phi) 10**-3', 
 '<phi> 10**-2', 'sigma(phi) 10**-2', '<phi> 10**-1', 'sigma(phi) 10**-1', '<U> 10**-3',
 '<U> 10**-2', '<U> 10**-1', 'D_kl', '<hydro> 10**-3', 'sigma(hydro) 10**-3', '<hydro> 10**-2', 
 'sigma(hydro) 10**-2', '<hydro> 10**-1', 'sigma(hydro) 10**-1']

hydrophobicity_dict = {'ARG':0, 'ASP':0.09, 'GLU':0.16, 'LYS':0.16, 'ASN':0.25,
                       'GLN':0.29, 'PRO':0.39, 'HIS':0.4, 'SER':0.42, 'THR':0.48,
                       'GLY':0.52, 'TYR':0.64, 'ALA':0.67, 'CYS':0.74, 'MET':0.84,
                       'TRP':0.85, 'VAL':0.89, 'PHE':0.96, 'LEU':0.97, 'ILE':1.0}

decoy_features = np.zeros((1, 20))
with open(struc_name+'_dkl.txt', 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
dkl = float(sourcelines[0])

decoy_features[0,13] = dkl

with open('./output/'+struc_name+'_packinglist_all.txt', 'r') as f:
    info = f.read()
sourcelines = info.splitlines()

rsasa_list = [0, 10**-3, 10**-2, 10**-1]
avg_packing_rsasa = []
std_packing_rsasa = []
avg_overlap_rsasa = []
fraction_core_rsasa = []
avg_hydro_rsasa = []
std_hydro_rsasa = []
for i in range(0, len(rsasa_list)-1):
    rsasa_lower = rsasa_list[i]
    rsasa_upper = rsasa_list[i+1]
    packing_list = []
    overlap_list = []
    hydro_list = []
    core_count = 0
    total_res_count = 0
    for line in sourcelines:
        rsasa = float(line.split()[-1])
        if rsasa_lower <= rsasa < rsasa_upper:
            packing_fraction = float(line.split()[2])
            overlap = float(line.split()[3])
            packing_list.append(packing_fraction)
            overlap_list.append(overlap)
            core_count += 1
            restype = line.split()[1]
            hydro = hydrophobicity_dict[restype]
            hydro_list.append(hydro)
        total_res_count += 1
    
    avg_packing = np.average(packing_list)
    std_packing = np.std(packing_list)
    avg_packing_rsasa.append(avg_packing)
    std_packing_rsasa.append(std_packing)
    
    avg_overlap = np.average(overlap_list)
    avg_overlap_rsasa.append(avg_overlap)
    
    fraction_core = core_count / total_res_count
    fraction_core_rsasa.append(fraction_core)
    
    avg_hydro = np.average(hydro_list)
    std_hydro = np.std(hydro_list)
    avg_hydro_rsasa.append(avg_hydro)
    std_hydro_rsasa.append(std_hydro)
    
decoy_features[0,0] = total_res_count

decoy_features[0,1] = fraction_core_rsasa[0]
decoy_features[0,2] = fraction_core_rsasa[1]
decoy_features[0,3] = fraction_core_rsasa[2]

decoy_features[0,4] = avg_packing_rsasa[0]
decoy_features[0,5] = std_packing_rsasa[0]
decoy_features[0,6] = avg_packing_rsasa[1]
decoy_features[0,7] = std_packing_rsasa[1]
decoy_features[0,8] = avg_packing_rsasa[2]
decoy_features[0,9] = std_packing_rsasa[2]

decoy_features[0,10] = avg_overlap_rsasa[0]
decoy_features[0,11] = avg_overlap_rsasa[1]
decoy_features[0,12] = avg_overlap_rsasa[2]

decoy_features[0,14] = avg_hydro_rsasa[0]
decoy_features[0,15] = std_hydro_rsasa[0]
decoy_features[0,16] = avg_hydro_rsasa[1]
decoy_features[0,17] = std_hydro_rsasa[1]
decoy_features[0,18] = avg_hydro_rsasa[2]
decoy_features[0,19] = std_hydro_rsasa[2]

# Saving Raw Features #
with open('./output/'+struc_name+'_features.txt', 'w') as f:
    for i in range(0, np.shape(decoy_features)[1]):
        f.write("{: <20} {: <20}".format(*[features[i]+':',str(decoy_features[0,i])]))
        f.write('\n')

# Scaling Features #
decoy_features[0,10] = np.log(decoy_features[0,10])
decoy_features[0,11] = np.log(decoy_features[0,11])
decoy_features[0,12] = np.log(decoy_features[0,12])

mu_list = []
with open('./SNN_wbs/mu.txt', 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
for line in sourcelines:
    mu_list.append(float(line.split()[0]))
mu_list = np.asarray(mu_list)
mu_list = np.reshape(mu_list, (1,20))


sigma_list = []
with open('./SNN_wbs/sigma.txt', 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
for line in sourcelines:
    sigma_list.append(float(line.split()[0]))
sigma_list = np.asarray(sigma_list)
sigma_list = np.reshape(sigma_list, (1,20))

decoy_features = (decoy_features - mu_list) / sigma_list

#########################################
print('\t Running SNN...')
# load json and create model
json_file = open('./SNN_wbs/SNN_GDT_trained.json', 'r')
loaded_model_json = json_file.read()
json_file.close()
model = model_from_json(loaded_model_json)
# load weights into new model
model.load_weights('./SNN_wbs/SNN_GDT_trained.h5')

optimizer = tf.keras.optimizers.Adam(10**-3)
model.compile(loss='logcosh',
                optimizer=optimizer,
                metrics=['mae', 'mse', 'logcosh'])

score = model.predict(decoy_features).flatten()[0]

print('Packing score = '+str(score))
with open('./output/'+struc_name+'_score.txt', 'w') as f:
    f.write(struc_name+': '+str(score))


#########################################
print('\t Cleaning up...')
command = 'rm '+struc_name+'_dkl.txt'
os.system(command)

command = 'rm '+struc_name+'_H.pdb'
os.system(command)

command = 'rm '+struc_name+'_noH.pdb'
os.system(command)

command = 'rm -r ./rSASA_code/'+struc_name+'_H_tempdir'
os.system(command)

command = 'rm '+struc_name+'_ordered*.pdb'
os.system(command)

command = 'rm '+struc_name+'_overlaplist.txt'
os.system(command)

command = 'rm '+struc_name+'_param*'
os.system(command)

command = 'rm '+struc_name+'_total_position*'
os.system(command)

command = 'rm '+struc_name+'1;_vol.txt'
os.system(command)

command = 'rm -r ./voro_code/result'
os.system(command)

command = 'rm -r ./voro_code/result_10'
os.system(command)

command = 'rm ./voro_code/'+struc_name+'_vor.txt'
os.system(command)

command = 'rm ./voro_code/'+struc_name+'_vor_10.txt'
os.system(command)
#########################################
print('Done \o\ \o/ /o/')

