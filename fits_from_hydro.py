# Create FITS file from .dat DEM(T), output from CR_NEI_hydro
# Assume separate DEM files for each species


# coding: utf-8


import numpy as np
from astropy.table import Table
import os
import astropy.io.fits as fits
import time


### DEMFILES

filename = '500yr'
path_demfiles=os.getcwd() + '/DEM_' + filename + '/'



myelem = np.array([1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28, 11, 13, 15, 22, 24, 25])
# Names for the final tables within the FITS table. Include Te
nam= np.array(['Te','H','He','C','N','O','Ne','Mg','Si','S','Ar','Ca','Fe','Ni','Na','Al','P','Ti','Cr','Mn'])


nion_H = 2 # H+0,+1
nion_befHe = 2
nion_He = 3 # He+0,...+2 
nion_befC = 5
nion_C = 7 # C+0,...+6 
nion_befN = 12
nion_N = 8 # N+0,...+7
nion_befO = 20 
nion_O = 9 # O+0,...+8 
nion_befNe = 29
nion_Ne = 11 # Ne+0,...+10 
nion_befMg = 40
nion_Mg = 13 # Mg+0,...+12 
nion_befSi = 53
nion_Si = 15 # Si+0,...+14 
nion_befS = 68
nion_S = 17 # S+0,...+16 
nion_befAr = 85
nion_Ar = 19 # Ar+0,...+18
nion_befCa = 104 
nion_Ca = 21 # Ca+0,...+20
nion_befFe = 125 
nion_Fe = 27 # Fe+0,...+26
nion_befNi = 152 
nion_Ni = 29 # Ni+0,...+28 
nion_befNa = 181
nion_Na = 12 # Na+0,...+11
nion_befAl = 193 
nion_Al = 14 # Al+0,...+13 
nion_befP = 207
nion_P = 16 # P+0,...+15 
nion_befTi = 223
nion_Ti = 23 # Ti+0,...+22 
nion_befCr = 246
nion_Cr = 25 # Cr+0,...+24
nion_befMn = 271 
nion_Mn = 26 # Mn+0,...+25 
# 297 ions in total



T_e = np.genfromtxt(path_demfiles + 'DEM_H.dat', skip_header=1)[:,0]
DEM__H = np.genfromtxt(path_demfiles + 'DEM_H.dat', skip_header=1)[:,1:]
DEM__He = np.genfromtxt(path_demfiles + 'DEM_He.dat', skip_header=1)[:,1:]
DEM__C = np.genfromtxt(path_demfiles + 'DEM_C.dat', skip_header=1)[:,1:]
DEM__N = np.genfromtxt(path_demfiles + 'DEM_N.dat', skip_header=1)[:,1:]
DEM__O = np.genfromtxt(path_demfiles + 'DEM_O.dat', skip_header=1)[:,1:]
DEM__Ne = np.genfromtxt(path_demfiles + 'DEM_Ne.dat', skip_header=1)[:,1:]
DEM__Mg = np.genfromtxt(path_demfiles + 'DEM_Mg.dat', skip_header=1)[:,1:]
DEM__Si = np.genfromtxt(path_demfiles + 'DEM_Si.dat', skip_header=1)[:,1:]
DEM__S = np.genfromtxt(path_demfiles + 'DEM_S.dat', skip_header=1)[:,1:]
DEM__Ar = np.genfromtxt(path_demfiles + 'DEM_Ar.dat', skip_header=1)[:,1:]
DEM__Ca = np.genfromtxt(path_demfiles + 'DEM_Ca.dat', skip_header=1)[:,1:]
DEM__Fe = np.genfromtxt(path_demfiles + 'DEM_Fe.dat', skip_header=1)[:,1:]
DEM__Ni = np.genfromtxt(path_demfiles + 'DEM_Ni.dat', skip_header=1)[:,1:]
DEM__Na = np.genfromtxt(path_demfiles + 'DEM_Na.dat', skip_header=1)[:,1:]
DEM__Al = np.genfromtxt(path_demfiles + 'DEM_Al.dat', skip_header=1)[:,1:]
DEM__P = np.genfromtxt(path_demfiles + 'DEM_P.dat', skip_header=1)[:,1:]
DEM__Ti = np.genfromtxt(path_demfiles + 'DEM_Ti.dat', skip_header=1)[:,1:]
DEM__Cr = np.genfromtxt(path_demfiles + 'DEM_Cr.dat', skip_header=1)[:,1:]
DEM__Mn = np.genfromtxt(path_demfiles + 'DEM_Mn.dat', skip_header=1)[:,1:]




Te = fits.Column(name='Te',format='1E', unit='K', array=T_e)


prihdr = fits.Header()
prihdr['DATE'] = time.strftime("%m/%d/%Y")
prihdr['COMMENT'] = "= DEMs for various ions and electron temperatures"
prihdu = fits.PrimaryHDU(header=prihdr)


cols = [0]*(len(myelem) + 1) # Te + Elements

tbhdu = [0]*(len(myelem) + 1 + 1) # Header + Te + Elements



DEM_H = np.zeros(nion_H, dtype=object)
DEM_He = np.zeros(nion_He, dtype=object)
DEM_C = np.zeros(nion_C, dtype=object)
DEM_N = np.zeros(nion_N, dtype=object)
DEM_O = np.zeros(nion_O, dtype=object)
DEM_Ne = np.zeros(nion_Ne, dtype=object)
DEM_Mg = np.zeros(nion_Mg, dtype=object)
DEM_Si = np.zeros(nion_Si, dtype=object)
DEM_S = np.zeros(nion_S, dtype=object)
DEM_Ar = np.zeros(nion_Ar, dtype=object)
DEM_Ca = np.zeros(nion_Ca, dtype=object)
DEM_Fe = np.zeros(nion_Fe, dtype=object)
DEM_Ni = np.zeros(nion_Ni, dtype=object)
DEM_Na = np.zeros(nion_Na, dtype=object)
DEM_Al = np.zeros(nion_Al, dtype=object)
DEM_P = np.zeros(nion_P, dtype=object)
DEM_Ti = np.zeros(nion_Ti, dtype=object)
DEM_Cr = np.zeros(nion_Cr, dtype=object)
DEM_Mn = np.zeros(nion_Mn, dtype=object)



for j in range(nion_H):
        DEM_H[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__H[:,j])
for j in range(nion_He):
        DEM_He[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__He[:,j])
for j in range(nion_C):
        DEM_C[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__C[:,j])
for j in range(nion_N):
        DEM_N[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__N[:,j])
for j in range(nion_O):
        DEM_O[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__O[:,j])
for j in range(nion_Ne):
        DEM_Ne[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ne[:,j])
for j in range(nion_Mg):
        DEM_Mg[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Mg[:,j])
for j in range(nion_Si):
        DEM_Si[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Si[:,j])
for j in range(nion_S):
        DEM_S[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__S[:,j])
for j in range(nion_Ar):
        DEM_Ar[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ar[:,j])
for j in range(nion_Ca):
        DEM_Ca[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ca[:,j])
for j in range(nion_Fe):
        DEM_Fe[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Fe[:,j])
for j in range(nion_Ni):
        DEM_Ni[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ni[:,j])
for j in range(nion_Na):
        DEM_Na[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Na[:,j])
for j in range(nion_Al):
        DEM_Al[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Al[:,j])
for j in range(nion_P):
        DEM_P[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__P[:,j])
for j in range(nion_Ti):
        DEM_Ti[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ti[:,j])
for j in range(nion_Cr):
        DEM_Cr[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Cr[:,j])
for j in range(nion_Mn):
        DEM_Mn[j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Mn[:,j])



# Final step: creation of the FITS table (in actuality, tables within a table)
    
cols[0] = fits.BinTableHDU.from_columns([Te], name='Te') 
cols[1] = fits.BinTableHDU.from_columns(DEM_H, name='H')
cols[2] = fits.BinTableHDU.from_columns(DEM_He, name='He')
cols[3] = fits.BinTableHDU.from_columns(DEM_C, name='C')
cols[4] = fits.BinTableHDU.from_columns(DEM_N, name='N')
cols[5] = fits.BinTableHDU.from_columns(DEM_O, name='O')
cols[6] = fits.BinTableHDU.from_columns(DEM_Ne, name='Ne')
cols[7] = fits.BinTableHDU.from_columns(DEM_Mg, name='Mg')
cols[8] = fits.BinTableHDU.from_columns(DEM_Si, name='Si')
cols[9] = fits.BinTableHDU.from_columns(DEM_S, name='S')
cols[10] = fits.BinTableHDU.from_columns(DEM_Ar, name='Ar')
cols[11] = fits.BinTableHDU.from_columns(DEM_Ca, name='Ca')
cols[12] = fits.BinTableHDU.from_columns(DEM_Fe, name='Fe')
cols[13] = fits.BinTableHDU.from_columns(DEM_Ni, name='Ni')
cols[14] = fits.BinTableHDU.from_columns(DEM_Na, name='Na')
cols[15] = fits.BinTableHDU.from_columns(DEM_Al, name='Al')
cols[16] = fits.BinTableHDU.from_columns(DEM_P, name='P')
cols[17] = fits.BinTableHDU.from_columns(DEM_Ti, name='Ti')
cols[18] = fits.BinTableHDU.from_columns(DEM_Cr, name='Cr')
cols[19] = fits.BinTableHDU.from_columns(DEM_Mn, name='Mn')


tbhdu[0] = prihdu
    

for c in range(len(myelem)+1): # Do not forget about Te!

    # Final table encompassing header, Te, and DEMs
    tbhdu[c+1] = fits.BinTableHDU.from_columns(cols[c].columns, name = nam[c], header = None) 
        
        
tbhdulist = fits.HDUList(tbhdu)
#tbhdulist.writeto('DDTa_' + filename + 'XSPEC.fits', clobber=True) # clobber forces rewriting, and remove ".dat" from the names 
tbhdulist.writeto('DDTa_' + filename + 'XSPEC.fits', overwrite=True) # clobber forces rewriting, and remove ".dat" from the names 