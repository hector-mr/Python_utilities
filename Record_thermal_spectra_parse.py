# coding: utf-8

# Create thermal spectra for the FS and the RS and for each time step (RS_evo, FS_evo)


import numpy as np
from astropy.table import Table, vstack, Column
import os
from scipy import interpolate
import pyatomdb as pyat



#####################################################################################################################################

# IMPORTANT: READ THIS
# Definition of filename as an optional argument to run the script multiple times in
# the terminal with all the different values via a simple loop in a shell script:
# python crspectra_resp.py filename
# where filename is any subdirectory in the present location
# See https://docs.python.org/3/library/argparse.html


import argparse
pathw = os.getcwd() + '/'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'index', metavar='str', type=str, choices=[x[1] for x in os.walk(pathw)][0],
         nargs='+', help='introduce any subdirectory in the present location')

    args = parser.parse_args()


filename = args.index[0] # Subdirectory selected to create the synthetic spectra with PyAtomDB. 
# e.g. DDTa.0.10, DDTa.0.50, DDTa.1.00, DDTa.5.00


# Note: usage of os.walk(pathw)
#for (path, dirs, files) in os.walk(path):
#    print path
#    print dirs
#    print files



# So, if the script is to be run from the directory itself (not the upper level suggested at the beginning of the script),
# instead of [x[2] for x in os.walk(pathw)], just comment put this parsing section and redefine path_files:

#path_files = pathw


#####################################################################################################################################


path_files = pathw + filename + '/'


age_sk = np.genfromtxt(path_files + 'data_sk_vs_rad.dat', skip_footer=2)[:,2] # age recorded for each shock
n_sk = len(age_sk)


explmodel = path_files.split(os.sep)[-2].split("_")[-1]
rhoISM = path_files.split(os.sep)[-3].split("_")[-1]


t_RS = np.zeros((n_sk+1), dtype='object')
t_FS = np.zeros((n_sk+1), dtype='object')




# IMPORTANT: READ THIS
# PyAtomDB calculates emissivities in [ct cm^3 s^-1 bin^-1] 
# for given energy bins, e.g., E = [1, 2, 3, 4, 5...] eV,
# and then retrieves a photon spectrum whose length is
# len(ebins) -1, corresponding to the mean energy for
# each bin, in this case E = [1.5, 2.5, 3.5, 4.5...] eV
# (see dummyfirst flag in make_ion_spectrum)

# By default, crhydronei calculates thermal RS/FS photon
# spectra in 10000 energy bins with a step of 1.2 eV, from 
# 0.95 keV to 12.0938 keV. To summarize,
# Input: 10000 energy bins
# Output: 9999 mean energy values, 9999 spectral values


En_start = 9.50E-2
En_step = 1.20E-3
En_nelem = 10000
# Input energies
En_th_bins = np.asarray([En_start + En_step*x for x in range(En_nelem)]) # Energy bin limits from thermal spectra
# Ouput energies that will be saved in the final files
En_th = np.asarray([En_start + En_step/2. + En_step*x for x in range(En_nelem - 1)]) # Average values for each bin ("Ebin")

ebins = En_th_bins 
Ebin = En_th



# These have the information for each shock
thermal_photons_xspec_unabs_RS_evo = np.genfromtxt(path_files + 'thermal_photons_xspec_unabs_RS_evo.dat')
thermal_photons_xspec_unabs_RS_evo = np.reshape(thermal_photons_xspec_unabs_RS_evo,(n_sk, En_nelem-1, 3))

thermal_photons_xspec_unabs_FS_evo = np.genfromtxt(path_files + 'thermal_photons_xspec_unabs_FS_evo.dat')
thermal_photons_xspec_unabs_FS_evo = np.reshape(thermal_photons_xspec_unabs_FS_evo,(n_sk, En_nelem-1, 3))


# The units for these thermal spectra are [ph cm^-2 s^-1 bin^-1]
# So, given that the energy step is 1.2 eV,
# [ph cm^-2 s^-1 bin^-1] * 1 bin / (1.2E-3 keV) = [ph cm^-2 s^-1 keV^-1]


##############################################################################

# thermal.f90

#!!        file= 'thermal_photons_xspec_unabs_FS(RS)_evo.dat'          
#          write(i_write_unit_2,"(A7, I5, 1p1e12.3, 2I5)") "#Epoch ", &
#                                                 i_sk_ct, time_yr_run, &
#                                                 i_skip_Doppler_1, &
#                                                 i_skip_Doppler_2


#!!            file= 'thermal_photons_xspec_unabs_FS(RS)_evo.dat'
#			  write(i_write_unit_2,"(1p3e15.5)") &
#					 binlo, binhi, &           
#					 totfluxk/bin_kev  ![N/cm^2/s] per E bin        
#		  enddo
#		  
#		  write(*,*) "Finished summing for shell", i_sk_ct, i_FS_RS

##############################################################################




t_RS[0] = Column(Ebin, name='Energy', unit='keV')
t_FS[0] = Column(Ebin, name='Energy', unit='keV')

for n in range(n_sk):

    t_RS[n+1] = Column(thermal_photons_xspec_unabs_RS_evo[n][:,2] / En_step, name = str(int(age_sk[n])) + 'yr', unit = 'ph+1cm-2s-1keV-1') 
    t_FS[n+1] = Column(thermal_photons_xspec_unabs_FS_evo[n][:,2] / En_step, name = str(int(age_sk[n])) + 'yr', unit = 'ph+1cm-2s-1keV-1') 
    # See https://astropy.readthedocs.io/en/v0.1/wcs/units.html for all the units (including counts and photons)




#####################################################################################################################################
# Finally, save a FITS table with the energy bins in the first column and the different photon spectra (ages) in the remaining ones


Table_end_RS = Table([x for x in t_RS])
Table_end_RS.write('spectra_ph_thermal_RS_' + explmodel + '_' + ISMrho + '.fits', overwrite='True')  

Table_end_FS = Table([x for x in t_FS])
Table_end_FS.write('spectra_ph_thermal_FS_' + explmodel + '_' + ISMrho + '.fits', overwrite='True') 