#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 21:50:48 2018

@author: agrigas115
"""
import glob

a = glob.glob('T*')

with open('pdb.txt', 'w') as f:
    for item in a:
        f.write("%s\n" % item)

print(len(a))