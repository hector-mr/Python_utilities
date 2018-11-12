# Merge early spectra (104 years) and spectra for older ages (30 columns) 


# coding: utf-8


import numpy as np
from astropy.table import Table, hstack, vstack, Column
import os, sys, time
import astropy.io.fits as fits




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
	expl_files_early[i] = fits.open(path_files + expl_names_early[i], format='fits')
expl_files_early = expl_files_early




expl_names = []
for file in os.listdir(path_files):
	if ((file.endswith('.fits')) and ('early' not in file)):
		expl_names.append(file)
		expl_names=sorted(expl_names,key=numericalSort) 


expl_files = []
for i in range(len(expl_names)):
	expl_files.append(i)
	expl_files[i] = fits.open(path_files + expl_names[i], format='fits')
expl_files = expl_files



# Merge tables for each model
n1 = len(expl_files)

for i in range(n1):


	prihdr = fits.Header()
	prihdr['DATE'] = time.strftime("%m/%d/%Y")
	prihdr['COMMENT'] = '= Physical profile of a modeled SNR for several expansion ages. Includes shocked ejecta and ambient medium'
	prihdu = fits.PrimaryHDU(header=prihdr)



	colsmrgd = []
	agesmrgd = []


	for fdx in range(np.shape(expl_files_early[i])[0] -1):

		# Append only if extension ages are lower than the ones from the second file
		if int(expl_files_early[i][fdx+1].header['EXTNAME'][:-2]) < int(expl_files[i][1].header['EXTNAME'][:-2]):
				
			colsmrgd.append(fits.BinTableHDU.from_columns(expl_files_early[i][fdx+1]))
			agesmrgd.append(int(expl_files_early[i][fdx+1].header['EXTNAME'][:-2]))
			

	# Append all the second file
	for gdx in range(np.shape(expl_files[i])[0] -1):
		
		colsmrgd.append(fits.BinTableHDU.from_columns(expl_files[i][gdx+1]))
		agesmrgd.append(int(expl_files[i][gdx+1].header['EXTNAME'][:-2]))
			

	tbhdu = [0]* (1 + len(colsmrgd)) # Header + early merged cols + cols
	tbhdu[0] = prihdu


	for tdx in range(len(tbhdu)-1):
		
		tbhdu[tdx + 1] = fits.BinTableHDU.from_columns(colsmrgd[tdx].columns, name = str(int(agesmrgd[tdx])) + 'yr', header = None)
		
		
	tbhdulist = fits.HDUList(tbhdu)
	tbhdulist.writeto(path_write + expl_names[i], overwrite=True)