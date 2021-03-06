# Create synthetic X-ray spectra with CR-NEI_hydro DEMs and PyAtomDB
# This script does NOT convolve the photon spectra with RMF and ARF files


# coding: utf-8


import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table, vstack, Column

import os
import astropy.io.fits as pyfits
import pyatomdb as pyat
#import time
#t0 = time.time()



#####################################################################################################################################

# IMPORTANT: READ THIS
# Definition of filename as an optional argument to run the script multiple times in
# the terminal with all the different values via a simple loop in a shell script:
# python create_spectra_resp.py filename
# where filename is any subdirectory in the present location
# See https://docs.python.org/3/library/argparse.html


import argparse
pathw = os.getcwd() + '/'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'index', metavar='str', type=str, choices=[x for x in os.listdir(pathw) if 'fits' in x],
         nargs='+', help='introduce any FITS file in the present location')

    args = parser.parse_args()


filename = args.index[0] # File selected to create an apec_nei_line.fits file which can be read by PyAtomDB. It is retrieved as a list

#####################################################################################################################################


path_files = pathw #+ filename + '/'
expl_file = pyfits.open(path_files + filename, format='fits')




# Load the AtomDB files and the responses

linefile = os.path.expandvars("$ATOMDB/apec_nei_line.fits")
cocofile = os.path.expandvars("$ATOMDB/apec_nei_comp.fits")

#myrmf  = os.path.expandvars('$ATOMDB/N103B_xi03.rmf')
#myarf  = os.path.expandvars('$ATOMDB/N103B_xi03_1p5.arf')




# Define energy range for the spectra


# IMPORTANT: READ THIS
# PyAtomDB calculates emissivities in [ct cm^3 s^-1 bin^-1]
# for given energy bins, e.g., E = [1, 2, 3, 4, 5...] eV,
# and then retrieves a photon spectrum whose length is
# len(ebins) -1, corresponding to the mean energy for
# each bin, in this case E = [1.5, 2.5, 3.5, 4.5...] eV
# (see dummyfirst flag in make_ion_spectrum)

# By default, crhydronei calculates thermal RS/FS photon
# spectra in 10000 energy bins with a step of 1.2 eV, from 
# 0.95 keV to 12.0938 keV. The same array will be used here
# for the DEM-RS/FS photon spectra. To summarize,
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




# Define a broadening, in keV, for the lines
#de = 0.01




# Define the temperature at which to plot (keV)
te = np.loadtxt(os.path.expandvars('$ATOMDB/Te_Adam.txt'))




# find the indices which are closest to the various temperatures
ite = []

for s in range(len(te)):
    ite.append(s)
    ite[s] = pyat.spectrum.get_index( te[s], teunits='keV', logscale=False)
    
ite = np.asarray(ite)




#####################################################################################################################################
# DEM(T)
# Each file will have a variable number of columns depending on the corresponding ionization states
# For example, each expl_file_DEM[i][j][2] will have 7 columns (CI, CII, CIII, CIV, CV, CVI, CVII), i.e., C+0+6 - 7 columns
# DEMs for each element and ion


# CR-NEI-hydro structure:
#       19 elements( H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe, Ni, Na, Al, P, Ti, Cr, Mn )
myelem = np.array([1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28, 11, 13, 15, 22, 24, 25])
mymetal = myelem[2:] # In the end, I will not care about H/He




##### CREATION OF THE SYNTHETIC SPECTRA #####




lente = len(te)
lenmetal = len(mymetal)

Spectra_em = np.zeros((1, lente, lenmetal), dtype=object)
Spectra_em_sumions = np.zeros((1, lente, lenmetal), dtype=object) 
Spectra_em_sumelem = np.zeros((1, lente), dtype=object)
Spectra_em_sumTe = np.zeros(1, dtype=object)
Spectra_ph = np.zeros(1, dtype=object)
# Unless dtype=object is included, setting array element within a sequence error will raise
# If...elif...else, NEVER USE If...if...else or Python will get stuck and just apply the final "else"

    
Spectra_em_sumTe[0] = [0.]*len(Ebin)#*len(ebins)


for j, it in enumerate(ite):
    
    Spectra_em_sumelem[0][j] = [0.]*len(Ebin)#*len(ebins)

    
    for k, el in enumerate(mymetal): # I do not care about H or He (see the [k+4] for expl_file)

        Spectra_em[0][j][k] = np.zeros(el+1, dtype=object) # Note the dtype=object again
            
        # So make_ion_spectrum(ebins,ite,6,7) is fully ionized carbon
        # And make_ion_spectrum(ebins,ite,6,1) is neutral carbon
        # Whereas make_ion_spectrum(ebins,ite,6,0/8) are null arrays

            
        for ion in range(1, el+2): #ion is the ion charge plus one (e.g. 3 for C III)

            Spectra_em[0][j][k][ion-1] = [0.]*len(ebins)

            # Generate the spectra multiplying the emissivities and the differential emission measures:
            # epsilon [ph cm+3 s-1] * DEM [cm-5 eV-1] = [ph cm-2 s-1 eV-1]
            # See the expression for the DEMs in the commented-out lines before this loop

            Spectra_em[0][j][k][ion-1] = expl_file[k+4].data[j][ion-1]*\
            pyat.spectrum.make_ion_spectrum(ebins, it, el, ion, linefile = linefile,\
             cocofile = cocofile)
         
                
                
                
        # Now I start to stack all the spectra. Over the various ions...  
        #Spectra_em_sumions[i][j][k] = [0.]*len(ebins)
        Spectra_em_sumions[0][j][k] = np.sum(Spectra_em[0][j][k])
        
    
        # ...plus the different elements...
        Spectra_em_sumelem[0][j] += Spectra_em_sumions[0][j][k]
    

    # ...and finally, the electron temperatures (many layers of each explosion model)
    Spectra_em_sumTe[0] += Spectra_em_sumelem[0][j]


# Last step, APEC normalization: necessary to multiply the convolved spectra by 1E14
# Also, to go from bin-1 to keV-1, a 1/En_step factor is necessary
Spectra_ph[0] = 1.E14 / En_step * Spectra_em_sumTe[0]
    



#####################################################################################################################################
# Finally, save a FITS table with the energy bins in the first column and the spectrum in the second one

u = np.zeros((2), dtype='object')
u[0] = Column(Ebin, name='Energy', unit='keV')

#explmodel = path_files.split(os.sep)[-2][6:]
#ISMrho = path_files.split(os.sep)[-3][-3:]

for el in filename.split("_"): # Split string based on underscores. Pretty convenient for my purposes!

    if "yr" in el:

        #age = int(el.replace("yr", ""))
        age = el



# CHANGE COLUMN NAME ACCORDINGLY. THOUGHT TO BE, E.G., 108yr
# See https://astropy.readthedocs.io/en/v0.1/wcs/units.html for all the units (including counts and photons)
u[1] = Column(Spectra_ph[0], name = age, unit='count+1cm-2s-1keV-1') 
    



# Save table    
Table_end_ph = Table([y for y in u])
Table_end_ph.write('spectra_ph_Enth_%s' % (filename), format='fits', overwrite='True') 