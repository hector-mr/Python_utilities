# Merge early spectra (104 years) and spectra for older ages (30 columns) 


# coding: utf-8


import numpy as np
from astropy.table import Table, hstack, vstack, Column
import os




path_files = os.getcwd() + '/'
path_write = os.getcwd() + '/Merged/'



#Sort files numerically
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):     
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


# Sort the files numerically, given that Python reads them in a random way
expl_names_early=[]
for file in os.listdir(path_files):
	if ((file.endswith('.fits')) & ('early' in file)):
   		expl_names_early.append(file)
    	expl_names_early=sorted(expl_names_early,key=numericalSort) 


expl_files_early = []
for i in range(len(expl_names_early)):
	expl_files_early.append(i)
	expl_files_early[i] = Table(Table.read(path_files + expl_names_early[i],format='fits')[6])
expl_files_early = expl_files_early




expl_names=[]
for file in os.listdir(path_files):
	if ((file.endswith('.fits')) and ('early' not in file)):
		expl_names.append(file)
		expl_names=sorted(expl_names,key=numericalSort) 


expl_files = []
for i in range(len(expl_names)):
	expl_files.append(i)
	expl_files[i] = Table.read(path_files + expl_names[i],format='fits')
expl_files = expl_files



# Merge tables for each model
n1 = len(expl_files)

for i in range(n1):

	vstack([expl_files_early[i], expl_files[i]]).write(path_write + expl_names[i], overwrite=True)