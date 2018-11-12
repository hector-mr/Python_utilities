# coding: utf-8

# Create DEM FITS files for all the recorded shocks in CR-NEI-hydro


import numpy as np
from astropy.table import Table
import os
import astropy.io.fits as fits
import time


path_files = os.getcwd() + '/'
filename = 'DEM_allages_FS.dat'


dem = []

with open(path_files + filename) as Dem:
    for n, line in enumerate(Dem):
            dem.append(line.split())




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


Nions = 297
Nelem = 19




# NOTATION: for instance, DDTa_500yr_2p0.dat refers to a DDTa model 
# with an expansion age of 500 years and an ISM density of 2g/cm3

# These DEM files have as many blocks as shocks recorded by crhydronei,
# comprise 15198 numbers (51 Te + 51*297 DEMs)
# and are preceded by a line with the shock number and the expansion age

# Then, len(dem) = 2*30 shocks = 60



sk_ct = []
age = []

for i in range(len(dem)):
    
    if len(dem[i]) == 2:

        sk_ct.append(int(dem[i][0]))
        age.append(int(round(float(dem[i][1])))) # Round to closest integer
        
sk_ct = np.asarray(sk_ct)
nshocks = sk_ct[-1]
age = np.asarray(age)



DEM__H = []
DEM__He = []
DEM__C = []
DEM__N =[]
DEM__O = []
DEM__Ne = []
DEM__Mg = []
DEM__Si = [] 
DEM__S = [] 
DEM__Ar = []
DEM__Ca = []
DEM__Fe = []
DEM__Ni = []
DEM__Na = []
DEM__Al = []
DEM__P = []
DEM__Ti = []
DEM__Cr = []
DEM__Mn = []


for i in np.arange(len(dem)):

    if len(dem[i]) > 2:
        
        DEM__H.append(dem[i][51 : 51*(nion_H + 1)])
        DEM__He.append(dem[i][51*(nion_befHe + 1) : 51*(nion_befC + 1)])
        DEM__C.append(dem[i][51*(nion_befC + 1) : 51*(nion_befN + 1)])
        DEM__N.append(dem[i][51*(nion_befN + 1) : 51*(nion_befO + 1)])
        DEM__O.append(dem[i][51*(nion_befO + 1) : 51*(nion_befNe + 1)])
        DEM__Ne.append(dem[i][51*(nion_befNe + 1) : 51*(nion_befMg + 1)])
        DEM__Mg.append(dem[i][51*(nion_befMg + 1) : 51*(nion_befSi + 1)])
        DEM__Si.append(dem[i][51*(nion_befSi + 1) : 51*(nion_befS + 1)])
        DEM__S.append(dem[i][51*(nion_befS + 1) : 51*(nion_befAr + 1)])
        DEM__Ar.append(dem[i][51*(nion_befAr + 1) : 51*(nion_befCa + 1)])
        DEM__Ca.append(dem[i][51*(nion_befCa + 1) : 51*(nion_befFe + 1)])
        DEM__Fe.append(dem[i][51*(nion_befFe + 1) : 51*(nion_befNi + 1)])
        DEM__Ni.append(dem[i][51*(nion_befNi + 1) : 51*(nion_befNa + 1)])
        DEM__Na.append(dem[i][51*(nion_befNa + 1) : 51*(nion_befAl + 1)])
        DEM__Al.append(dem[i][51*(nion_befAl + 1) : 51*(nion_befP + 1)])
        DEM__P.append(dem[i][51*(nion_befP + 1) : 51*(nion_befTi + 1)])
        DEM__Ti.append(dem[i][51*(nion_befTi + 1) : 51*(nion_befCr + 1)])
        DEM__Cr.append(dem[i][51*(nion_befCr + 1) : 51*(nion_befMn + 1)])
        DEM__Mn.append(dem[i][51*(nion_befMn + 1) : ])

        

# Some values not recorded properly (e.g. 1.045473-103)
try:
    DEM__H = np.asarray(DEM__H).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__H[i])):
            try:
                DEM__H[i][j] = float(DEM__H[i][j])
            except ValueError:
                DEM__H[i][j] = 0.0    
DEM__H = np.reshape(DEM__H,(len(age), nion_H, 51))

try:
    DEM__He = np.asarray(DEM__He).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__He[i])):
            try:
                DEM__He[i][j] = float(DEM__He[i][j])
            except ValueError:
                DEM__He[i][j] = 0.0    
DEM__He = np.reshape(DEM__He,(len(age), nion_He, 51))

try:
    DEM__C = np.asarray(DEM__C).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__C[i])):
            try:
                DEM__C[i][j] = float(DEM__C[i][j])
            except ValueError:
                DEM__C[i][j] = 0.0    
DEM__C = np.reshape(DEM__C,(len(age), nion_C, 51))

try:
    DEM__N = np.asarray(DEM__N).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__N[i])):
            try:
                DEM__N[i][j] = float(DEM__N[i][j])
            except ValueError:
                DEM__N[i][j] = 0.0    
DEM__N = np.reshape(DEM__N,(len(age), nion_N, 51))

try:
    DEM__O = np.asarray(DEM__O).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__O[i])):
            try:
                DEM__O[i][j] = float(DEM__O[i][j])
            except ValueError:
                DEM__O[i][j] = 0.0    
DEM__O = np.reshape(DEM__O,(len(age), nion_O, 51))

try:
    DEM__Ne = np.asarray(DEM__Ne).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Ne[i])):
            try:
                DEM__Ne[i][j] = float(DEM__Ne[i][j])
            except ValueError:
                DEM__Ne[i][j] = 0.0    
DEM__Ne = np.reshape(DEM__Ne,(len(age), nion_Ne, 51))

try:
    DEM__Mg = np.asarray(DEM__Mg).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Mg[i])):
            try:
                DEM__Mg[i][j] = float(DEM__Mg[i][j])
            except ValueError:
                DEM__Mg[i][j] = 0.0    
DEM__Mg = np.reshape(DEM__Mg,(len(age), nion_Mg, 51))

try:
    DEM__Si = np.asarray(DEM__Si).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Si[i])):
            try:
                DEM__Si[i][j] = float(DEM__Si[i][j])
            except ValueError:
                DEM__Si[i][j] = 0.0    
DEM__Si = np.reshape(DEM__Si,(len(age), nion_Si, 51))

try:
    DEM__S = np.asarray(DEM__S).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__S[i])):
            try:
                DEM__S[i][j] = float(DEM__S[i][j])
            except ValueError:
                DEM__S[i][j] = 0.0    
DEM__S = np.reshape(DEM__S,(len(age), nion_S, 51))

try:
    DEM__Ar = np.asarray(DEM__Ar).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Ar[i])):
            try:
                DEM__Ar[i][j] = float(DEM__Ar[i][j])
            except ValueError:
                DEM__Ar[i][j] = 0.0    
DEM__Ar = np.reshape(DEM__Ar,(len(age), nion_Ar, 51))

try:
    DEM__Ca = np.asarray(DEM__Ca).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Ca[i])):
            try:
                DEM__Ca[i][j] = float(DEM__Ca[i][j])
            except ValueError:
                DEM__Ca[i][j] = 0.0    
DEM__Ca = np.reshape(DEM__Ca,(len(age), nion_Ca, 51))

try:
    DEM__Fe = np.asarray(DEM__Fe).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Fe[i])):
            try:
                DEM__Fe[i][j] = float(DEM__Fe[i][j])
            except ValueError:
                DEM__Fe[i][j] = 0.0    
DEM__Fe = np.reshape(DEM__Fe,(len(age), nion_Fe, 51))

try:
    DEM__Ni = np.asarray(DEM__Ni).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Ni[i])):
            try:
                DEM__Ni[i][j] = float(DEM__Ni[i][j])
            except ValueError:
                DEM__Ni[i][j] = 0.0    
DEM__Ni = np.reshape(DEM__Ni,(len(age), nion_Ni, 51))

try:
    DEM__Na = np.asarray(DEM__Na).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Na[i])):
            try:
                DEM__Na[i][j] = float(DEM__Na[i][j])
            except ValueError:
                DEM__Na[i][j] = 0.0    
DEM__Na = np.reshape(DEM__Na,(len(age), nion_Na, 51))

try:
    DEM__Al = np.asarray(DEM__Al).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Al[i])):
            try:
                DEM__Al[i][j] = float(DEM__Al[i][j])
            except ValueError:
                DEM__Al[i][j] = 0.0    
DEM__Al = np.reshape(DEM__Al,(len(age), nion_Al, 51))

try:
    DEM__P = np.asarray(DEM__P).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__P[i])):
            try:
                DEM__P[i][j] = float(DEM__P[i][j])
            except ValueError:
                DEM__P[i][j] = 0.0    
DEM__P = np.reshape(DEM__P,(len(age), nion_P, 51))

try:
    DEM__Ti = np.asarray(DEM__Ti).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Ti[i])):
            try:
                DEM__Ti[i][j] = float(DEM__Ti[i][j])
            except ValueError:
                DEM__Ti[i][j] = 0.0    
DEM__Ti = np.reshape(DEM__Ti,(len(age), nion_Ti, 51))

try:
    DEM__Cr = np.asarray(DEM__Cr).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Cr[i])):
            try:
                DEM__Cr[i][j] = float(DEM__Cr[i][j])
            except ValueError:
                DEM__Cr[i][j] = 0.0    
DEM__Cr = np.reshape(DEM__Cr,(len(age), nion_Cr, 51))

try:
    DEM__Mn = np.asarray(DEM__Mn).astype(np.float)
except ValueError:
    for i in np.arange(nshocks):
        for j in range(len(DEM__Mn[i])):
            try:
                DEM__Mn[i][j] = float(DEM__Mn[i][j])
            except ValueError:
                DEM__Mn[i][j] = 0.0    
DEM__Mn = np.reshape(DEM__Mn,(len(age), nion_Mn, 51))





T_e = np.asarray(dem[1][0:51]).astype(np.float)

Te = fits.Column(name='Te',format='1E', unit='K', array=T_e)


prihdr = fits.Header()
prihdr['DATE'] = time.strftime("%m/%d/%Y")
prihdr['COMMENT'] = "= DEMs for various ions and electron temperatures"
prihdu = fits.PrimaryHDU(header=prihdr)


cols = []  
tbhdu = []
tbhdulist = []


arange = np.arange(len(age))


DEM_H = np.zeros((len(age),nion_H), dtype=object)
DEM_He = np.zeros((len(age),nion_He), dtype=object)
DEM_C = np.zeros((len(age),nion_C), dtype=object)
DEM_N = np.zeros((len(age),nion_N), dtype=object)
DEM_O = np.zeros((len(age),nion_O), dtype=object)
DEM_Ne = np.zeros((len(age),nion_Ne), dtype=object)
DEM_Mg = np.zeros((len(age),nion_Mg), dtype=object)
DEM_Si = np.zeros((len(age),nion_Si), dtype=object)
DEM_S = np.zeros((len(age),nion_S), dtype=object)
DEM_Ar = np.zeros((len(age),nion_Ar), dtype=object)
DEM_Ca = np.zeros((len(age),nion_Ca), dtype=object)
DEM_Fe = np.zeros((len(age),nion_Fe), dtype=object)
DEM_Ni = np.zeros((len(age),nion_Ni), dtype=object)
DEM_Na = np.zeros((len(age),nion_Na), dtype=object)
DEM_Al = np.zeros((len(age),nion_Al), dtype=object)
DEM_P = np.zeros((len(age),nion_P), dtype=object)
DEM_Ti = np.zeros((len(age),nion_Ti), dtype=object)
DEM_Cr = np.zeros((len(age),nion_Cr), dtype=object)
DEM_Mn = np.zeros((len(age),nion_Mn), dtype=object)




for i in arange:
    
    cols.append([i])
    cols[i] = cols[i]*(len(myelem) + 1) # Te + Elements for each age
    
    tbhdu.append([i])
    tbhdu[i] = tbhdu[i]*(len(myelem) + 1 + 1) # Header + Te + Elements for each age
    
    
    for j in range(nion_H):
        DEM_H[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__H[i][j])
    for j in range(nion_He):
        DEM_He[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__He[i][j])
    for j in range(nion_C):
        DEM_C[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__C[i][j])
    for j in range(nion_N):
        DEM_N[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__N[i][j])
    for j in range(nion_O):
        DEM_O[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__O[i][j])
    for j in range(nion_Ne):
        DEM_Ne[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ne[i][j])
    for j in range(nion_Mg):
        DEM_Mg[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Mg[i][j])
    for j in range(nion_Si):
        DEM_Si[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Si[i][j])
    for j in range(nion_S):
        DEM_S[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__S[i][j])
    for j in range(nion_Ar):
        DEM_Ar[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ar[i][j])
    for j in range(nion_Ca):
        DEM_Ca[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ca[i][j])
    for j in range(nion_Fe):
        DEM_Fe[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Fe[i][j])
    for j in range(nion_Ni):
        DEM_Ni[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ni[i][j])
    for j in range(nion_Na):
        DEM_Na[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Na[i][j])
    for j in range(nion_Al):
        DEM_Al[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Al[i][j])
    for j in range(nion_P):
        DEM_P[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__P[i][j])
    for j in range(nion_Ti):
        DEM_Ti[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Ti[i][j])
    for j in range(nion_Cr):
        DEM_Cr[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Cr[i][j])
    for j in range(nion_Mn):
        DEM_Mn[i][j] = fits.Column(name=str(j),format='1E', unit='cm-5 K-1', array=DEM__Mn[i][j])
        
        
        
    # Final step: creation of the FITS table (in actuality, tables within a table)
    
    cols[i][0] = fits.BinTableHDU.from_columns([Te], name='Te') 
    cols[i][1] = fits.BinTableHDU.from_columns(DEM_H[i], name='H')
    cols[i][2] = fits.BinTableHDU.from_columns(DEM_He[i], name='He')
    cols[i][3] = fits.BinTableHDU.from_columns(DEM_C[i], name='C')
    cols[i][4] = fits.BinTableHDU.from_columns(DEM_N[i], name='N')
    cols[i][5] = fits.BinTableHDU.from_columns(DEM_O[i], name='O')
    cols[i][6] = fits.BinTableHDU.from_columns(DEM_Ne[i], name='Ne')
    cols[i][7] = fits.BinTableHDU.from_columns(DEM_Mg[i], name='Mg')
    cols[i][8] = fits.BinTableHDU.from_columns(DEM_Si[i], name='Si')
    cols[i][9] = fits.BinTableHDU.from_columns(DEM_S[i], name='S')
    cols[i][10] = fits.BinTableHDU.from_columns(DEM_Ar[i], name='Ar')
    cols[i][11] = fits.BinTableHDU.from_columns(DEM_Ca[i], name='Ca')
    cols[i][12] = fits.BinTableHDU.from_columns(DEM_Fe[i], name='Fe')
    cols[i][13] = fits.BinTableHDU.from_columns(DEM_Ni[i], name='Ni')
    cols[i][14] = fits.BinTableHDU.from_columns(DEM_Na[i], name='Na')
    cols[i][15] = fits.BinTableHDU.from_columns(DEM_Al[i], name='Al')
    cols[i][16] = fits.BinTableHDU.from_columns(DEM_P[i], name='P')
    cols[i][17] = fits.BinTableHDU.from_columns(DEM_Ti[i], name='Ti')
    cols[i][18] = fits.BinTableHDU.from_columns(DEM_Cr[i], name='Cr')
    cols[i][19] = fits.BinTableHDU.from_columns(DEM_Mn[i], name='Mn')


    tbhdu[i][0] = prihdu
    
    
    for c in range(len(myelem)+1): # Do not forget about Te!

        # Final table encompassing header, Te, and DEMs
        tbhdu[i][c+1] = fits.BinTableHDU.from_columns(cols[i][c].columns, name = nam[c], header = None)
       
        
    tbhdulist.append(i)
    tbhdulist[i] = fits.HDUList(tbhdu[i])
    


    # Final DEM tables, e.g. SCH088_500yr_2p0.fits


    ## path_files.split(os.sep) will yield, for '/home/hector/CR_NEI_hydro/ChN_Ia_ddta_ifort_500yr_2p0/model_ddta/',
	## ['', 'home', 'hector', 'CR_NEI_hydro', 'ChN_Ia_ddta_ifort_500yr_2p0', 'model_ddta', '']

    explmodel = path_files.split(os.sep)[-2].split("_")[-1]
    rhoISM = path_files.split(os.sep)[-3].split("_")[-1]

    tbhdulist[i].writeto(explmodel + '_' + str(age[i]) + 'yr' + '_' + rhoISM + '_FS.fits', overwrite = True) 
    # clobber forces rewriting, and remove ".dat" from the names
