# Core Packing Score
Implementation of the core packing score published here: https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3914

Contact me with any questions: alex.grigas@yale.edu

This method will first calculate the relative Solvent Accessibile Surface Area of each residue, the packing fraction of each residue, the overlap energy of each residue and the Kullback-Leibler divergence of the distribution of hydrophobic packing fraction compared to an x-ray crystal structure reference distribution. Then, these inputs are used in a feed-forward neural network to predict GDT_TS.

## Requirements
Requirments to run code

- Python 3.7
- Matlab w/ Bioinformatics Toolbox
- Reduce
- Naccess
- Pomelo

Python Modules
- Numpy
- Biopython
- Scipy
- Tensorflow

# Installation

## Matlab
This method calls Matlab from Python. In order to do this, you need to first the Matlab Engine API for Python. From the Matlab command window, type:
```
cd (fullfile(matlabroot,'extern','engines','python'))
system('python setup.py install')
```
This can also be accomplished by navigating to the Matlab root folder on the commandline and entering 'python setup.py install'
See: https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
Make sure to also install the Matlab Bioinformatics Toolbox.

## Reduce
Reduce is a commonly used software to add hydrogen atoms to PDB structures. Please install the Reduce executable to the preproccessing_code directory. Rename the executable to 'reduce'. See: http://kinemage.biochem.duke.edu/software/reduce.php

## Naccess
Naccess calculates the solvent accessible surface area of proteins. Please install all of the Naccess files to the rSASA_code directory.
See: http://wolf.bms.umist.ac.uk/naccess/

## Pomelo
Pomelo calculates the Voronoi volumes of the residues using the surface Voronoi Tessellation. Make sure to compile Pomelo in the 'Generic' mode, as this is what allows Pomelo to use generic surface triangulations as the input defining the surface of the particles.
See: https://github.com/spatialfruitsalad/pomelo

## Volume Estimation
In the vol_code directory is 'grid_method_packing_fraction.cxx'. To compile:
```
g++ grid_method_packing_fraction.cxx -o vol
```
If you don't have it already, you can install GCC by:
```
brew install gcc
```
If you, like me, use anaconda, Hombrew and Anaconda don't always get along. If you have trouble calling g++, try temproarily changing you PATH, see: https://hashrocket.com/blog/posts/keep-anaconda-from-constricting-your-homebrew-installs

# Usage
This code will calculate the following for an input PDB structure:

- relative Solvent Accessible Surface Area (rSASA)
- residue packing fraction
- residue overlap energy
- Kullback-Liebler Divergence of the hydrophobic packing fraction

It will save this information for you and then predict the Global Distance Test score of the structure.
```
python run_core_packing_decoy_detector.py -i abcd.pdb
```
All of the results will be saved to the output directory.
