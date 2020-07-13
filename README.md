# Core Packing Score
Implementation of the core packing score published here: https://arxiv.org/abs/2001.01161

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
