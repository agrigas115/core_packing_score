#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 16:07:48 2019
Calculating the packing fraction
-i : structure name
@author: agrigas115
"""
import argparse
import scipy.io as sio

parser = argparse.ArgumentParser(description='Core Packing Decoy Detector') 
parser.add_argument('-i', '--input', help='input pdb_file', required=True)  
args = parser.parse_args()

decoy_name = args.input
vol_txt = './'+decoy_name+'1;_vol.txt'
#volume data file names
vor_txt = './voro_code/'+decoy_name+'_vor.txt'
#voronoi volume data file names
vor_txt_10 = './voro_code/'+decoy_name+'_vor_10.txt'
#second set of voronoi volume data to find the closed voronoi cells
rSASA_txt = './rSASA_code/rSASA_data/'+decoy_name+'_H_sasa_data'
#rSASA data file names
resid_txt = './'+decoy_name+'.txt'
#pdb data file name

overlap_file = './'+decoy_name+'_overlaplist.txt'
overlap_data = []
with open(overlap_file, 'r') as f:
    info = f.read()
sourcelines = info.splitlines()
for line in sourcelines:
    overlap = line.split()[1]
    if overlap == '0':
        overlap = '0                 '
    overlap_data.append(line.split()[1])


with open(vor_txt, 'r') as f:
    vors_list = f.read() 
vor_lines = vors_list.splitlines()
vor_lines.pop(0)
starting_index = 0
if float(vor_lines[0].split()[1]) == 0:
    vor_lines.pop(0)
    starting_index = 1
if float(vor_lines[0].split()[1]) == 0:
    vor_lines.pop(0)
    starting_index = 2
#The first pop removes the labeling of the columns in the vor file
#The second check is looking to see if there is volume for the 0th particle
#Sometimes there is volume their because my .txt files which i ran the 
#Pomelo calulation on start at 0, sometimes at 1
#Starting index is then used to calculate the volume per residue
vor_data = []
for item in vor_lines:
    line_split = line_split = item.split()
    vor_data.append([float(line_split[0]), float(line_split[1])])
#reading in vor data and creating a least where each entry is a list
#with [resid, voronoi_volume]

with open(vor_txt_10, 'r') as f:
    vors_list = f.read() 
vor_lines = vors_list.splitlines()
vor_lines.pop(0)
starting_index = 0
if float(vor_lines[0].split()[1]) == 0:
    vor_lines.pop(0)
    starting_index = 1
if float(vor_lines[0].split()[1]) == 0:
    vor_lines.pop(0)
    starting_index = 2
vor_data_10 = []
for item in vor_lines:
    line_split = line_split = item.split()
    vor_data_10.append([float(line_split[0]), float(line_split[1])])
#reading in vor data and creating a least where each entry is a list
#with [resid, voronoi_volume]
closed_vor_list = []
for k in range(0, len(vor_data)):
    if vor_data[k][1] == vor_data_10[k][1]:
        closed_vor_list.append(1)
    elif abs(vor_data[k][1]-vor_data_10[k][1])/(vor_data[k][1]) <= 0.01:
        closed_vor_list.append(1)
    else:
        closed_vor_list.append(0)
#making a list of the closed voronoi cells
#if it is closed: 1
#if it is not closed: 0

with open(vol_txt, "r") as f:
    volumes_list = f.read()
vol_lines = volumes_list.splitlines()
vol_data = []
for item in vol_lines:
    line_split = item.split()
    vol_data.append([float(line_split[1]), float(line_split[-1])])
#reading in vol data and creating a least where each entry is a list
#with [solvent_exposed, volume]


with open(resid_txt, 'r') as f:
    txt_list = f.read()
txt_lines = txt_list.splitlines()
resid_data = []
prev_res = 0
dict_resid = {}
for item in txt_lines:
    line_split = item.split()
    resid_data.append(float(line_split[2]))
    current_res = line_split[2]
    if prev_res != current_res:
        dict_resid[float(line_split[2])] = line_split[1]
        prev_res = line_split[2]
#readinging in the original pdb in order to associate each atom number
#with the right resid
mat_load = (sio.loadmat(rSASA_txt))
rSASA_data = mat_load['each_res_data']
#reading in the rSASA data, which is stored as a .mat file so it is opened
#using scipy.io

vol_per_residue = []
for i in range(starting_index, int(resid_data[-1])+1):
    vol_list = []
    skip = 0
    for j in range(0, len(resid_data)):
        if resid_data[j] == i:
            vol_list.append(vol_data[j][1])
            if vol_data[j][1] == 0:
                skip = 1
    if skip == 0:
        vol_per_residue.append([float(i), sum(vol_list), 'yes_count'])
    else:
        vol_per_residue.append([float(i), sum(vol_list), 'dont_count'])
#calculating the total volume of each residue
#looping over the total number of residues, checking to see if that atom
#is part of that residue and then adding those together saving to a list
#added a check to remove any residues that have a 0 volume entry
this_structure_data = []
for i in range(0, len(vol_per_residue)):
    if closed_vor_list[i] == 1:
        packing = str(vol_per_residue[i][1]/vor_data[i][1])
        this_structure_data.append([i+starting_index, dict_resid[i+starting_index], packing, overlap_data[i], rSASA_data[i][1]])
    if closed_vor_list[i] == 0:
        packing = '0'
        this_structure_data.append([i+starting_index, dict_resid[i+starting_index], packing, overlap_data[i], rSASA_data[i][1]])
#calculating packing fraction and making a list to save all the info to
#the packinglist file


with open('./output/'+decoy_name+'_packinglist_all.txt', 'w') as f:
        for row in this_structure_data:
            f.write("{: <8} {: <8} {: <23} {: <23} {: <23}".format(*row))
            f.write('\n')