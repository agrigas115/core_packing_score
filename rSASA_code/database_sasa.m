%% Get sasa for dataset of pdb files
function [] = database_sasa()
clear;
close all;
clc;


% get pdb dir
parentdir = ['./'];
pdbdir = [parentdir '../'];

% get savedir
savedir = [parentdir 'rSASA_data/'];
if ~exist(savedir,'dir')
    mkdir(savedir);
end
% get mdir
mdir = [parentdir 'pdb_mat/'];

% get sasa for all proteins in pdb
flist = dir(fullfile([pdbdir '*_H.pdb']));
NF = length(flist);

for ff = 1:NF
    get_single_pdb_sasa(mdir,pdbdir,flist(ff).name,savedir);
end
