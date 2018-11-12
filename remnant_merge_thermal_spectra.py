# coding: utf-8


import numpy as np
from astropy.table import Table, hstack, vstack, Column
import os




path_files = os.getcwd() + '/'
path_write = os.getcwd() + '/../'



#Sort files numerically
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):     
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


# Sort the files numerically, given that Python reads them in a random way
expl_names=[]
for file in os.listdir(path_files):
	if file.endswith('.fits'):
		expl_names.append(file)
		expl_names=sorted(expl_names,key=numericalSort) 


expl_files = []
for i in range(len(expl_names)):
	expl_files.append(i)

	# Take first column (energy)
	if i == 0:
		Energies = Table.read(path_files + expl_names[i],format='fits').columns[0]

	expl_files[i] = Table.read(path_files + expl_names[i],format='fits').columns[1]

	# For the purposes of this single script. Use same column name as DEM spectra
	expl_names[i] = expl_names[i][:-6]
	expl_names[i] = expl_names[i].split("_")[-3] + '_' + expl_names[i].split("_")[-1]\
	 + '_' + expl_names[i].split("_")[-2] + '_' + expl_names[i].split("_")[-4]


expl_files = expl_files


lenfi = len(expl_files)

# Merge all tables (energy, spectra)
t = np.zeros((lenfi+1), dtype = 'object')
t[0] = Column(Energies, name = 'Energy', unit = 'keV')

for l in range(lenfi):

    t[l+1] = Column(expl_files[l], name = expl_names[l], unit = 'ph+1s-1keV-1cm-2')



fsrs = path_files.split("/")[-3] + '_'
expl = path_files.split("/")[-2].split("_")[1]
remnant = path_files.split("/")[-2].split("_")[0]
explremnant = expl + '_' + remnant

Table_end = Table([x for x in t])
Table_end.write(path_write + 'spectra_ph_thermal_%s.fits' % (fsrs + explremnant), format = 'fits', overwrite = 'True') 