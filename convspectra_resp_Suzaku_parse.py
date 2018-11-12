# Convolve photon, synthetic X-ray spectra calculated with CR-NEI_hydro 
# DEMs and PyAtomDB using RMF and ARF files
# Telescope: Suzaku. Observation: Tycho. Change accordingly (Chandra, XARM, Athena,...)

# This script is meant to convolve multiple photon spectra grouped as a FITS file,
# being the first column the energy and the rest the spectra for different ages


# coding: utf-8


import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table, vstack, Column

import os
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


# Load pre-calculated spectrum and get energy bins

path_files = pathw #+ filename + '/'
sp_file = Table.read(path_files + filename)


Energy = sp_file.columns[0]
En_step = Energy[1] - Energy[0] # 1.20E-3
En_start = Energy[0] - En_step/2. # 9.50E-2
En_nelem = len(Energy) + 1 # 10000
ebins = np.asarray([En_start + En_step*x for x in range(En_nelem)]) # Energy bin limits from thermal spectra
bins = len(ebins)


# Note: this original unconvolved spectrum has units of keV-1,but PyAtomDB operates with bin-1, so need to unpack
# Will have to do the opposite later, for the convolved spectrum, to go from bin-1 to keV-1
sp_to_con = Table(sp_file.columns[1:])




# Load the responses

# Chandra
#myrmf  = os.path.expandvars('$ATOMDB/acisi_aimpt_cy20.rmf')
#myarf  = os.path.expandvars('$ATOMDB/acisi_aimpt_cy20.arf')
# Linear ebins with En_step = 0.01147 keV = 11.47 eV

# Suzaku
#myrmf  = os.path.expandvars('$ATOMDB/N103B_xi03.rmf')
#myarf  = os.path.expandvars('$ATOMDB/N103B_xi03_1p5.arf')
# Aside from ~10 bins, linear ebins with En_step = 0.002 keV = 2 eV

# XARM
#myrmf  = os.path.expandvars('$ATOMDB/xarm_res_h5ev_20170818.rmf')
#myarf  = os.path.expandvars('$ATOMDB/xarm_res_pnt_pa_20170818.arf')
# Linear ebins with En_step = 0.0005 keV = 0.5 eV


# Change accordingly, as well as final file name
myrmf  = os.path.expandvars('/home/hector/SNR_Data/Suzaku_Data/Tycho_Data/output/clevt/xi03.rmf')
myarf  = os.path.expandvars('/home/hector/SNR_Data/Suzaku_Data/Tycho_Data/output/clevt/xi03.arf')




# Get the response ebins
ebins_rmf =  pyat.spectrum.get_response_ebins(myrmf)

En_step_rmf = np.zeros( len(ebins_rmf)-1 ) # Energy step [keV] (assuming it is not constant through the whole range)
Ebin_rmf = np.zeros( len(ebins_rmf)-1 ) # Average energy value in each bin

for r in range(len(ebins_rmf)-1):
    
    En_step_rmf[r] = ebins_rmf[r+1] - ebins_rmf[r]

    Ebin_rmf[r] = ebins_rmf[r] +  En_step_rmf[r] / 2. 
    
nbins_rmf = len(ebins_rmf)




# Convolve rhe photon spectra
ages = sp_to_con.colnames
nages = len(ages)
spec_out = np.zeros(nages, dtype=object)


for c, cc in enumerate(sp_to_con.columns): 
      
    # Get cumulative sum of spectrum
    cumulspectrum = np.cumsum(sp_to_con.columns[c])

    # Add in a zero at the beginning
    cumulspectrum = np.append(0, cumulspectrum)

    # Interpolate onto ebins_rmf
    cumulspectrum_rmf = np.interp(ebins_rmf, ebins, cumulspectrum)

    # Get counts/bin on this rmf energy grid
    spectrum_rmf_in = cumulspectrum_rmf[1:]-cumulspectrum_rmf[:-1]

    # Apply the response

    if c == 0:

        ebins_out, spec_out[c] = pyat.spectrum.apply_response(spectrum_rmf_in, rmf = myrmf, arf = myarf)

    else:

        spec_out[c] = pyat.spectrum.apply_response(spectrum_rmf_in, rmf = myrmf, arf = myarf)[1]




Eresp_out_step = np.zeros( len(ebins_out)-1 ) # Energy step [keV] (it is not constant through the whole range)
Eresp_out = np.zeros( len(ebins_out)-1 ) # Average energy value in each bin


for rr in range(len(ebins_out)-1):

    Eresp_out_step[rr] = ebins_out[rr+1] - ebins_out[rr]
    Eresp_out[rr] = ebins_out[rr] +  (Eresp_out_step[rr]) / 2
    
    
    
    
# Now, AFTER the convolution, to go from bin-1 to keV-1, a 1/Eresp_step factor is necessary

# Chandra: Eresp_step ~0.0146 keV = 1.46 eV (NOT exactly the same for each bin)  
# Eresp[0] = 0.0109500002582 keV; Eresp[-1] = 14.9430999756 keV; len(Eresp) = 1024 (same as unconvolved)

# Suzaku: Eresp_step ~0.00365 keV = 3.65 eV (NOT exactly the same for each bin)  
# Eresp[0] = 0.00200750006479 keV; Eresp[-1] = 14.9485750198 keV; len(Eresp) = 4096 < len(unconcolved) = 7900

# XARM: Eresp_step ~0.00050 keV = 0.50 eV (NOT exactly the same for each bin)
# Eresp[0] = 0.000275000009424 keV; Eresp[-1] = 16.3837499619 keV; len(Eresp) = 32768 (same as unconvolved)




# create polaw continuum for high-energy bins (above 11 keV)

def polaw(E, gamma, norm):
    
    return E**gamma*norm


highval = 11.0 # keV

wh_10keV = np.where(Eresp_out >= 10.0)[0][0]
wh_11keV = np.where(Eresp_out >= 11.0)[0][0]

Eresp_out_10keV = Eresp_out[wh_10keV]
Eresp_out_11keV = Eresp_out[wh_11keV]




for s, sp in enumerate(spec_out):



    for stp, stpel in enumerate(Eresp_out_step):



        spec_out[s][stp] = spec_out[s][stp] / stpel * En_step # See previous comments to understand why En_step comes into play



        if stpel == Eresp_out_step[-1]: # For the power law

            spec_out_10keV = spec_out[s][wh_10keV]
            spec_out_11keV = spec_out[s][wh_11keV]


            Gamma = ( np.log(spec_out_11keV) - np.log(spec_out_10keV) ) / ( np.log(Eresp_out_11keV) - np.log(Eresp_out_10keV) )
            Norm  = spec_out_11keV / (Eresp_out_11keV**Gamma)



        if Eresp_out[stp] > highval:
        
            spec_out[s][stp] = polaw(Eresp_out[stp], Gamma, Norm)




#####################################################################################################################################
# Finally, save a FITS table with the energy bins in the first column and the spectra in the other ones

t = np.zeros(( nages + 1 ), dtype = 'object')
t[0] = Column(Eresp_out, name = 'Energy', unit = 'keV')


for l in range(nages):

    t[l+1] = Column(spec_out[l], name = ages[l], unit = 'ph+1s-1keV-1')

    # See https://astropy.readthedocs.io/en/v0.1/wcs/units.html for all the units (including counts and photo
    


# Save table
# filename assumed to be, e.g., spectra_ph_Enth_RS_sch115_2p0
# Change accordingly from Tycho to whichever SNR    
Table_end = Table([x for x in t])
Table_end.write('spectra_resp_Suzaku_Tycho%s' % (filename.split('Enth')[1]), format = 'fits', overwrite = 'True') 