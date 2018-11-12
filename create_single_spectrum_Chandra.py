# Create synthetic X-ray spectra with CR-NEI_hydro DEMs and PyAtomDB
# This script DOES convolve the photon spectra with RMF and ARF files
# Telescope: Chandra. Change accordingly (Suzaku, XARM, Athena,...)


# coding: utf-8


import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table, vstack, Column

import os
import astropy.io.fits as pyfits
import pyatomdb as pyat

import warnings
warnings.filterwarnings("ignore") 
# Ignore (XARM) -- RuntimeWarning: divide by zero encountered in divide lammax = (const.HC_IN_KEV_A/ebins[0])

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

# Chandra
myrmf  = os.path.expandvars('$ATOMDB/acisi_aimpt_cy20.rmf')
myarf  = os.path.expandvars('$ATOMDB/acisi_aimpt_cy20.arf')
# Linear ebins with En_step = 0.01147 keV = 11.47 eV

# Suzaku
#myrmf  = os.path.expandvars('$ATOMDB/N103B_xi03.rmf')
#myarf  = os.path.expandvars('$ATOMDB/N103B_xi03_1p5.arf')
# Aside from ~10 bins, linear ebins with En_step = 0.002 keV = 2 eV

# XARM
#myrmf  = os.path.expandvars('$ATOMDB/xarm_res_h5ev_20170818.rmf')
#myarf  = os.path.expandvars('$ATOMDB/xarm_res_pnt_pa_20170818.arf')
# Linear ebins with En_step = 0.0005 keV = 0.5 eV




# Define energy range for the spectra

ebins =  pyat.spectrum.get_response_ebins(myrmf)

En_step = np.zeros( len(ebins)-1 ) # Energy step [keV] (assuming it is not constant through the whole range)
Ebin = np.zeros( len(ebins)-1 ) # Average energy value in each bin

for r in range(len(ebins)-1):
    
    En_step[r] = ebins[r+1] - ebins[r]

    Ebin[r] = ebins[r] +  En_step[r] / 2. 



# IMPORTANT: READ THIS
# PyAtomDB calculates emissivities in [ct cm^3 s^-1 bin^-1]
# for given energy bins, e.g., E = [1, 2, 3, 4, 5...] eV,
# and then retrieves a photon spectrum whose length is
# len(ebins) -1, corresponding to the mean energy for
# each bin, in this case E = [1.5, 2.5, 3.5, 4.5...] eV
# (see dummyfirst flag in make_ion_spectrum)

# To summarize,

# Chandra
# Input: 1025 energy bins from 0.255 to 12 keV
# Output: 1024 mean energy values, 1024 spectral values

# Suzaku
# Input: 7901 energy bins from 0.2 to 16 keV
# Output: 7900 mean energy values, 7900 spectral values

# XARM
# Input: 32769 energy bins from 0.0 to 16.384 keV
# Output: 32768 mean energy values, 32768 spectral values

# Then, I will convolve the output with the spectral responses
# (see the end of the script)



# Define a broadening, in keV, for the lines
#de = 0.01




# Define the temperature at which to plot (keV)
te = np.loadtxt(os.path.expandvars('$ATOMDB/Te_Adam.txt'))




# find the indices which are closest to the various temperatures
ite = []

for s in range(len(te)):
    ite.append(s)
    ite[s] = pyat.spectrum.get_index( te[s], teunits = 'keV', logscale = False)
    
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
Spectra_resp = np.zeros(1, dtype=object)
# Unless dtype=object is included, setting array element within a sequence error will raise
# If...elif...else, NEVER USE If...if...else or Python will get stuck and just apply the final "else"

    
Spectra_em_sumTe[0] = [0.]*len(Ebin)


for j, it in enumerate(ite):
    
    Spectra_em_sumelem[0][j] = [0.]*len(Ebin)

    
    for k, el in enumerate(mymetal): # I do not care about H or He (see the [k+4] for expl_file)

        Spectra_em[0][j][k] = np.zeros(el+1, dtype=object) # Note the dtype=object again
            
        # So make_ion_spectrum(ebins,ite,6,7) is fully ionized carbon
        # And make_ion_spectrum(ebins,ite,6,1) is neutral carbon
        # Whereas make_ion_spectrum(ebins,ite,6,0/8) are null arrays

            
        for ion in range(1, el+2): #ion is the ion charge plus one (e.g. 3 for C III)

            Spectra_em[0][j][k][ion-1] = [0.]*len(ebins)

            # Generate the spectra multiplying the emissivities and the differential emission measures:
            # epsilon [ph cm+3 s-1] * DEM [cm-5] = [ph cm-2 s-1 bin-1]

            Spectra_em[0][j][k][ion-1] = expl_file[k+4].data[j][ion-1]*\
            pyat.spectrum.make_ion_spectrum(ebins, it, el, ion, linefile = linefile,\
             cocofile = cocofile)
         
                
                
                
        # Now I start to stack all the spectra. Over the various ions...  
        Spectra_em_sumions[0][j][k] = np.sum(Spectra_em[0][j][k])
        
    
        # ...plus the different elements...
        Spectra_em_sumelem[0][j] += Spectra_em_sumions[0][j][k]
    

    # ...and finally, the electron temperatures (many layers of each explosion model)
    Spectra_em_sumTe[0] += Spectra_em_sumelem[0][j]




##### Fold spectrum with response #####




# APEC normalization: necessary to multiply the photon spectra by 1E14
Spectra_resp[0] = pyat.spectrum.apply_response(1.E14 * Spectra_em_sumTe[0], rmf = myrmf, arf = myarf)[1] # Fluxes


Energ_resp = pyat.spectrum.apply_response(1.E14 * Spectra_em_sumTe[0], rmf = myrmf, arf = myarf)[0] # Energy bins


Eresp_step = np.zeros( len(Energ_resp)-1 ) # Energy step [keV] (it is not constant through the whole range)
Eresp = np.zeros( len(Energ_resp)-1 ) # Average energy value in each bin


for rr in range(len(Energ_resp)-1):

    Eresp_step[rr] = Energ_resp[rr+1] - Energ_resp[rr]
    Eresp[rr] = Energ_resp[rr] +  (Eresp_step[rr]) / 2. 


# Now, AFTER the convolution, to go from bin-1 to keV-1, a 1/Eresp_step factor is necessary

# Chandra: Eresp_step ~0.0146 keV = 1.46 eV (NOT exactly the same for each bin)  
# Eresp[0] = 0.0109500002582 keV; Eresp[-1] = 14.9430999756 keV; len(Eresp) = 1024 (same as unconvolved)

# Suzaku: Eresp_step ~0.00365 keV = 3.65 eV (NOT exactly the same for each bin)  
# Eresp[0] = 0.00200750006479 keV; Eresp[-1] = 14.9485750198 keV; len(Eresp) = 4096 < len(unconcolved) = 7900

# XARM: Eresp_step ~0.00050 keV = 0.50 eV (NOT exactly the same for each bin)
# Eresp[0] = 0.000275000009424 keV; Eresp[-1] = 16.3837499619 keV; len(Eresp) = 32768 (same as unconvolved)


for stp, stpel in enumerate(Eresp_step): 

    Spectra_resp[0][stp] = Spectra_resp[0][stp] / stpel 




#####################################################################################################################################
# Finally, save a FITS table with the energy bins in the first column and the spectrum in the second one

t = np.zeros((2), dtype = 'object')
t[0] = Column(Eresp, name = 'Energy', unit = 'keV')

#explmodel = path_files.split(os.sep)[-2][6:]
#ISMrho = path_files.split(os.sep)[-3][-3:]

for el in filename.split("_"): # Split string based on underscores. Pretty convenient for my purposes!

    if "yr" in el:

        #age = int(el.replace("yr", ""))
        age = el



# CHANGE COLUMN NAME ACCORDINGLY. THOUGHT TO BE, E.G., 108yr
# See https://astropy.readthedocs.io/en/v0.1/wcs/units.html for all the units (including counts and photons)
t[1] = Column(Spectra_resp[0], name = age, unit = 'ph+1s-1keV-1') 
    



# Save table    
Table_end = Table([x for x in t])
Table_end.write('spectra_resp_Chandra_%s' % (filename), format = 'fits', overwrite = 'True') 