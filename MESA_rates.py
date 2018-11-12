# I created this script in order to retrieve the rates from the MESA binary files located in ~/MESA/data/rates_data/cache. There is an executable file in ~/MESA/rates/test/ called show_rates which can only be executed from there,
# e.g. ./show_rates ~/MESA/data/rates_data/cache/r_na22_ap_mg25_1.bin > r_na22_ap_mg25_1.txt, because it needs to access to ../../data/version_number

# PLACE THIS SCRIPT IN THE DIRECTORY ~/MESA/rates/test  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# coding: utf-8


import numpy as np
import scipy as scipy
import os
from os import walk
import subprocess
import tempfile
import glob
import atpy



pathexe = './show_rates' # ~/MESA/rates/test/show_rates
homedir = '/home/hector'
pathfiles = homedir + '/MESA/data/rates_data/cache/'
save_path = homedir + '/MESA_projects/WD_ignition/MESA_rates'


# Find all the .bin files in the "pathfiles" directory


filelist = []
for file in os.listdir(pathfiles):
	if file.endswith(".bin"):
		filelist.append(file)
		filelist = sorted(filelist) # Sorts the files by name, given that Python reads them in a random way



# Write the complete path of these .bin files so that show_rates works

files = []
for i in range(len(filelist)):
	files.append(i)
	files[i] = pathfiles + filelist[i]


 # Create the .txt tables

tables = []
names = []
output = []
for i,item in enumerate(filelist):
	names.append(i)
	names[i] = os.path.join(save_path, item.replace(".bin", "")+".txt") # Save the tables, which are named according to "item", in a different directory (save_path) 
	# Also, item.replace(".bin", "") removes .bin from the strings in order to have more beautiful file names, e.g. r_na22_ap_mg25_1.txt instead of r_na22_ap_mg25_1.bin.txt
for i,item in enumerate(files):
	tables.append(i)
	output.append(i)
	tables[i] = open(names[i],"w")  
	output[i] = subprocess.Popen(pathexe + " " + item,shell=True, stdout=subprocess.PIPE).stdout # Execute command on the terminal and save the output
	output[i] = output[i].read() # Execute command on the terminal and save the output
	tables[i].write(output[i].replace("D", "E"))  # ESSENTIAL! Replace in all the tables, e.g., 1D+3 by 1E+3 in order to be able to use the data in Python!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	tables[i].close()