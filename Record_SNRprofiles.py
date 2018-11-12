# coding: utf-8

# Record SNR profiles for each crhydronei directory
# Information coming from data_vs_rad.dat 


import numpy as np
from astropy.table import Table, vstack, Column
import os, sys, time
import astropy.io.fits as fits



pathw=os.getcwd()  # Current work directory


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)



subdirect = []

for i, el in enumerate([x[1] for x in os.walk(pathw)][0]):

	if len(el) > 25: # Record elements like ChN_Ia_sch115_ifort_2000yr_2p0

		subdirect.append(el)

# List of all subdirectories. Change the parameter [0] depending on the needs (i.e., on the depth of the nested grid of subdirectories)
# It should always be x[1], which corresponds to the directories present in the current work path. x[2] comprises the files present in the current work path


subdirect = sorted(subdirect, key=str.lower) 	# Otherwise, capitalized strings go first and the vector becomes a mess


path_files = []

for i in range(len(subdirect)):
    path_files.append(i)
    if 'ddta' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_ddta/'
    if 'ddt12' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_ddt12/'
    if 'ddt16' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_ddt16/'
    if 'ddt24' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_ddt24/'
    if 'ddt40' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_ddt40/'
    elif 'sch088' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_sch088/'
    elif 'sch097' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_sch097/'
    elif 'sch106' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_sch106/'
    elif 'sch115' in subdirect[i]:
        path_files[i] = pathw + '/' + subdirect[i] + '/model_sch115/'




nion_befH = 0
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


elements = np.asarray(['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Fe', 'Ni', 'Na', 'Al', 'P', 'Ti', 'Cr', 'Mn'])
idx_ion_limit = [] # For instance, the last component will be [nion_befMn, nion_befMn + nion_Mn]. Simplify notation




tnion_befH = 0
tnion_H = 1 # H+1
tnion_befHe = 1
tnion_He = 2 # He+1,...+2 
tnion_befC = 3
tnion_C = 6 # C+1,...+6 
tnion_befN = 9
tnion_N = 7 # N+1,...+7
tnion_befO = 16 
tnion_O = 8 # O+1,...+8 
tnion_befNe = 24
tnion_Ne = 10 # Ne+1,...+10 
tnion_befMg = 34
tnion_Mg = 12 # Mg+1,...+12 
tnion_befSi = 46
tnion_Si = 14 # Si+1,...+14 
tnion_befS = 60
tnion_S = 16 # S+1,...+16 
tnion_befAr = 76
tnion_Ar = 18 # Ar+1,...+18
tnion_befCa = 94 
tnion_Ca = 20 # Ca+1,...+20
tnion_befFe = 114 
tnion_Fe = 26 # Fe+1,...+26
tnion_befNi = 140 
tnion_Ni = 28 # Ni+1,...+28 
tnion_befNa = 168
tnion_Na = 11 # Na+1,...+11
tnion_befAl = 179 
tnion_Al = 13 # Al+1,...+13 
tnion_befP = 192
tnion_P = 15 # P+1,...+15 
tnion_befTi = 207
tnion_Ti = 22 # Ti+1,...+22 
tnion_befCr = 229
tnion_Cr = 24 # Cr+1,...+24
tnion_befMn = 253 
tnion_Mn = 25 # Mn+1,...+25 
# (297 ions - 19 elements) ion temperatures in total


idx_tion_limit = [] # Different for the ion temperatures: neutrals do NOT belong to the array. See the explanation below




for index, el in enumerate(elements):
    
    idx_ion_limit.append([0, 0])
    idx_ion_limit[index] = [vars()['nion_bef' + el], vars()['nion_bef' + el] + vars()['nion_' + el]] # Index limits
    
    idx_tion_limit.append([0, 0])
    idx_tion_limit[index] = [vars()['tnion_bef' + el], vars()['tnion_bef' + el] + vars()['tnion_' + el]] # Index limits




def heaviside(x):
    """ Heaviside step function """
    if (x > 0):
        result = 1
    elif (x == 0):
        result = 0.5
    else:
        result = 0

    return result




# Assume nshocks is not equal for each model


pc_to_cm = 3.086e18
Msun = 1.98E33
dist = 10 #kpc



for idx, path in enumerate(path_files):

    #20 profiles, regardless of nshocks: 
    # num_div_print = 20    !Number of profile printouts
    data_vs_rad = np.loadtxt(path + 'data_vs_rad.dat')
    data_vs_rad_age = np.unique(data_vs_rad[:,19]) # time_yr_run
    
    
    if data_vs_rad_age[0] == 0:
        data_vs_rad_age = data_vs_rad_age[1:]
       
        
    data_vs_rad_iRS = np.zeros(len(data_vs_rad_age), dtype=int)
    data_vs_rad_iCD = np.zeros(len(data_vs_rad_age), dtype=int)
    data_vs_rad_iFS = np.zeros(len(data_vs_rad_age), dtype=int)

    
    
    r = [ []*i for i in range(len(data_vs_rad_age)) ]
    dm = [ []*i for i in range(len(data_vs_rad_age)) ]
    lagm = [ []*i for i in range(len(data_vs_rad_age)) ]
    
    r_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    dm_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    lagm_amb = [ []*i for i in range(len(data_vs_rad_age)) ]


    
    data_vs_rad_age_block = [ []*i for i in range(len(data_vs_rad_age)) ]



    uxx = [ []*i for i in range(len(data_vs_rad_age)) ]
    te = [ []*i for i in range(len(data_vs_rad_age)) ]
    rho = [ []*i for i in range(len(data_vs_rad_age)) ]
    xnee = [ []*i for i in range(len(data_vs_rad_age)) ]
    ti = [ []*i for i in range(len(data_vs_rad_age)) ]
    xni = [ []*i for i in range(len(data_vs_rad_age)) ]
    tau = [ []*i for i in range(len(data_vs_rad_age)) ]


    uxx_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    te_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    rho_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    xnee_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    ti_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    xni_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    tau_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
        
        
        
        
    for ag, datag in enumerate(data_vs_rad_age):
    
        data_vs_rad_age_block[ag] = np.where(data_vs_rad[:,19] == datag)[0]
    
        data_vs_rad_iRS[ag] = np.where(data_vs_rad[:,11][data_vs_rad_age_block[ag]] == 1)[0][0] # rev_shock
        data_vs_rad_iCD[ag] = np.where(data_vs_rad[:,18][data_vs_rad_age_block[ag]] == 1)[0][0] # rad_o_CD 
        data_vs_rad_iFS[ag] = np.where(data_vs_rad[:,17][data_vs_rad_age_block[ag]] == 1)[0][0] # rad_o_FS
    
    
        r[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,2] * pc_to_cm # radius in cm
        dm[ag] = 4*np.pi * r[ag]**2 * data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,3]
        lagm[ag] = np.zeros(len(r[ag]))
        
        r_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,2] * pc_to_cm # radius in cm
        dm_amb[ag] = 4*np.pi * r_amb[ag]**2 * data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,3]
        lagm_amb[ag] = np.zeros(len(r_amb[ag]))


        uxx[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,5] # v [km/s]
        te[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,34] # te
        rho[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,3] # rho
        xnee[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,28] # xnee
        # CAREFUL: Neutrals not heated by shock. Hence, the 19 neutral ions have NO ion temperature, which explains the "- Nelem" factor
        ti[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,35 : 35 + Nions - Nelem] # ti for all the ion species
        xni[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,35 + Nions - Nelem : 35 + Nions - Nelem + Nions] # ni for all the ion species
        tau[ag] = data_vs_rad[data_vs_rad_age_block[ag]][0:data_vs_rad_iCD[ag]][:,-1] # tau = ne*t
    
    
        uxx_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,5] # v [km/s]
        te_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,34] # te
        rho_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,3] # rho
        xnee_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,28] # xnee
        ti_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,35 : 35 + Nions - Nelem] # ti
        xni_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,35 + Nions - Nelem : 35 + Nions - Nelem + Nions] # ni
        tau_amb[ag] = data_vs_rad[data_vs_rad_age_block[ag]][data_vs_rad_iCD[ag]:data_vs_rad_iFS[ag]][:,-1] # tau = ne*t
    
    
    
    
        for j, el in enumerate(r[ag]): #RS
        
            if j == 0:
            
                lagm[ag][j] = dm[ag][j] * el / Msun
            
            else:
            
                # TRAPEZOIDAL RULE FOR INTEGRATING
                lagm[ag][j] = lagm[ag][j-1] + 0.5 * (dm[ag][j]+dm[ag][j-1]) * (r[ag][j]-r[ag][j-1]) / Msun
            
            
    
            
        for o, al in enumerate(r_amb[ag]): #FS
        
            if o == 0:
            
                lagm_amb[ag][o] = 0.
            
            else:
            
                # TRAPEZOIDAL RULE FOR INTEGRATING
                lagm_amb[ag][o] = lagm_amb[ag][o-1] + 0.5*(dm_amb[ag][o]+dm_amb[ag][o-1])*(r_amb[ag][o]-r_amb[ag][o-1]) / Msun




    lagm = np.asarray(lagm)
    r = np.asarray(r)
    uxx = np.asarray(uxx)        
    te = np.asarray(te)
    rho = np.asarray(rho)
    xnee = np.asarray(xnee)
    ti = np.asarray(ti)
    tau = np.asarray(tau)

    lagm_amb = np.asarray(lagm_amb)
    r_amb = np.asarray(r_amb)
    uxx_amb = np.asarray(uxx_amb)
    te_amb = np.asarray(te_amb)
    rho_amb = np.asarray(rho_amb)
    xnee_amb = np.asarray(xnee_amb)
    ti_amb = np.asarray(ti_amb)
    tau_amb = np.asarray(tau_amb)




    # Correct some nans

    for i in range(np.shape(ti)[0]): # np.shape(ti)[0] = np.shape(ti_amb)[0]  (ages x radial coordinate)[0] = ages
        
        for j in range(np.shape(ti)[1]):
            
            for a in range(Nions - Nelem):
            
            # Correct some nans
                if np.isfinite(ti[i][j][a]) == False:
                
                    if a > 0:
                
                        ti[i][j][a] = ti[i][j][a-1]
                    
                    else:
                    
                        ti[i][j][a] = ti[i][j][0]
                        
                        
            for aa in range(Nions):
                        
                if np.isfinite(xni[i][j][aa]) == False:
                    
                    xni[i][j][aa] = -99
                    
                    
                    
                    
        for k in range(np.shape(ti_amb[i])[0]):
                 
            for b in range(Nions - Nelem):
                    
            # Correct some nans
                if np.isfinite(ti_amb[i][k][b]) == False:
                
                    if b > 0:
                
                        ti_amb[i][k][b] = ti_amb[i][k][b-1]
                    
                    else:
                    
                        ti_amb[i][k][b] = ti_amb[i][k][0]
                        
                  
            for bb in range(Nions):
                           
                if np.isfinite(xni_amb[i][k][bb]) == False:
                    
                    xni_amb[i][k][bb] = -99




    # For ti_average calculations

    wh_neutral = np.zeros(Nelem, dtype=int)

    for ineu, neu in enumerate(idx_ion_limit):
        wh_neutral[ineu] = neu[0]
        
    wh_noneutral = np.arange(Nions, dtype=int)

    wh_neutral = list(wh_neutral)
    wh_noneutral = list(wh_noneutral)
        
    for inoneu in sorted(wh_neutral, reverse=True):

        del wh_noneutral[inoneu] # Remove opposite set




    # Mean and average ion temperature, whose indices are (age, radial/mass coordinate, ion)
    # But beware: ti are log-values!!!!!!!!

    ti_mean = np.zeros((np.shape(ti)[0],np.shape(ti)[1]))
    ti_mean_amb = []

    ti_aver = np.zeros((np.shape(ti)[0],np.shape(ti)[1]))
    ti_aver_amb = []

    


    # Also calculate ti, xni, xel = sum(xni), zel_aver = average charge state. But beware: xni are log-values!!!!!!!!


    for index, el in enumerate(elements):
        
        vars()['ti_' + el] = [] # Ion temperatures (ti_H, ti_He,...)
        vars()['ti_' + el + '_amb'] = [] # Ion temperatures (ti_H_amb, ti_He_amb,...)
        
        vars()['xni_' + el] = [] # Ion fractions (xni_H, xni_He,...)
        vars()['xni_' + el + '_amb'] = [] # Ion fractions (xni_H_amb, xni_He_amb,...)
            
        

        
        for i in range(np.shape(ti)[0]):
            
            vars()['ti_' + el].append([0.])
            vars()['ti_' + el][i] = (vars()['ti_' + el][i])*np.shape(ti[i])[0]
            
            vars()['xni_' + el].append([0.])
            vars()['xni_' + el][i] = (vars()['xni_' + el][i])*np.shape(xni[i])[0]
            
            
            
            
            for j in range(np.shape(ti)[1]):
            
                vars()['ti_' + el][i][j] = ti[i][j][idx_tion_limit[index][0] : idx_tion_limit[index][1]]
                
                vars()['xni_' + el][i][j] = xni[i][j][idx_ion_limit[index][0] : idx_ion_limit[index][1]]
                
                
                
                
                if index == 0: # Avoid looping unnecessarily
                    
                    ti_mean[i][j] = np.log10(np.mean(10**ti[i][j]))
            
                    # Beware of the wh_noneutral indices
                    ti_aver[i][j] = np.log10(np.sum(10**ti[i][j]*10**xni[i][j][wh_noneutral])/np.sum(10**xni[i][j][wh_noneutral]))
                
                
                
            
            vars()['ti_' + el + '_amb'].append([0.])
            vars()['ti_' + el + '_amb'][i] = (vars()['ti_' + el + '_amb'][i])*np.shape(ti_amb[i])[0]
            
            vars()['xni_' + el + '_amb'].append([0.])
            vars()['xni_' + el + '_amb'][i] = (vars()['xni_' + el + '_amb'][i])*np.shape(xni_amb[i])[0]
            
            
            
            
            if index == 0: # Avoid looping unnecessarily
                    
                ti_mean_amb.append([0.]) # NEITHER APPEND K NOR ZERO!!!!!!!!
                ti_mean_amb[i] = (ti_mean_amb[i])*(np.shape(ti_amb[i])[0])
                ti_mean_amb[i] = np.asarray(ti_mean_amb[i])
        
                ti_aver_amb.append([0.]) # NEITHER APPEND K NOR ZERO!!!!!!!!
                ti_aver_amb[i] = (ti_aver_amb[i])*(np.shape(ti_amb[i])[0])
                ti_aver_amb[i] = np.asarray(ti_aver_amb[i])
            
            
            
            
            for k in range(np.shape(ti_amb[i])[0]):
                
                vars()['ti_' + el + '_amb'][i][k] = ti_amb[i][k][idx_tion_limit[index][0] : idx_tion_limit[index][1]]
                
                vars()['xni_' + el + '_amb'][i][k] = xni_amb[i][k][idx_ion_limit[index][0] : idx_ion_limit[index][1]]
                
                
                
                
                if index == 0: # Avoid looping unnecessarily
                    
                    ti_mean_amb[i][k] = np.log10(np.mean(10**ti_amb[i][k]))
            
                    ti_aver_amb[i][k] = np.log10(np.sum(10**ti_amb[i][k]*10**xni_amb[i][k][wh_noneutral])\
                                                 /np.sum(10**xni_amb[i][k][wh_noneutral]))
                
                
                
                
            if index == 0: # Avoid looping unnecessarily
                
                ti_mean_amb[i] = np.asarray(ti_mean_amb[i])
        
                ti_aver_amb[i] = np.asarray(ti_aver_amb[i])
                
                
                
        
        # This is needed for the next calculations (below dotted line)
        vars()['ti_' + el] = np.asarray(vars()['ti_' + el])
        vars()['ti_' + el + '_amb'] = np.asarray(vars()['ti_' + el + '_amb'])
        
        vars()['xni_' + el] = np.asarray(vars()['xni_' + el])
        vars()['xni_' + el + '_amb'] = np.asarray(vars()['xni_' + el + '_amb'])
        
        
        
        
        if index == 0: # Avoid looping unnecessarily
            
            ti_mean_amb = np.asarray(ti_mean_amb)
            ti_aver_amb = np.asarray(ti_aver_amb)
        

        
        #------------------------------------------------------------------------------------------------------------------------------# 
        
        
        
        vars()['x' + el] = np.sum(10**vars()['xni_' + el], axis=2) # Number fractions (xH, xHe,...) [xH = np.sum(10**xni_H, axis=2)]
        vars()['x' + el + '_amb'] = []  # Number fractions (xH_amb, xHe_amb,...)
        
        
        # Herman says:

        #If you are only interested in the ion fractions within an element inside a particular grid, 
        #then you can use the cnc array in ion_cnc_evo.dat as well. That file only contains info in 
        #the shocked plasma (i.e. from RS to FS). However, to integrate and obtain the overall ejecta 
        #or CSM distribution, you still need to know the number densities, right? i.e., from the ion 
        #fractions, you need to multiply with the abundance as a function of radius (also available 
        #in the same file), and then with the volumes to weight the contribution from each shell. 
        #If you use densities from data_vs_rad.dat, you can skip the abundance part.   

        #For the plot I sent you, as mentioned above, you need to first multiply the volume of each 
        #radial grid (or shell) to the population and add them up, before you renormalize with the maximum ion fraction.

        vars()['z' + el + '_aver'] = np.zeros((np.shape(ti)[0],np.shape(ti)[1])) # Average charge states (zH_aver, zHe_aver,...)
        vars()['z' + el + '_aver_amb'] = [] # Average charge states (zH_aver_amb, zHe_aver_amb,...)
        
        
        
        
        vars()['EMi_' + el] = np.zeros(np.shape(vars()['xni_' + el])) # Integrated emission measure for each ion and mass coordinate
        vars()['EM_' + el] = np.zeros(np.shape(vars()['x' + el])) # Integrated emission measure for each element and mass coordinate
        vars()['EM_' + el + '_tot'] = np.zeros(len(data_vs_rad_age)) # Total emission measure for each element
        
        #vars()['EMi_' + el + '_amb'] = []
        #vars()['EM_' + el + '_amb'] = []
        #vars()['EM_' + el + '_tot_amb']= np.zeros(len(data_vs_rad_age))
        
        
        
        
        for i in range(np.shape(ti)[0]):
        
            for j in range(np.shape(ti)[1]):
                
                # Only the shocked ejecta contributes to the emission measure. Remove layers with r < RS
                hv_EM = heaviside(r[i][j] - r[i][data_vs_rad_iRS[i]])
            
                for a in range(vars()['nion_' + el]): # Number of ions for each element
                    
                    vars()['z' + el + '_aver'][i][j] += a*10**vars()['xni_' + el][i][j][a] # Charge * abundance
                    
                    
                    # Volume integral of (n_e*n_ion) normmalized to 10 kpc
                    vars()['EMi_' + el][i][j][a] = np.trapz(y = hv_EM * 10**vars()['xni_' + el][i][0:j][:,a] * 10**xnee[i][0:j] \
                                                            * 4*np.pi * r[i][0:j]**2 , x = r[i][0:j]) / (4*np.pi*(pc_to_cm*1000*dist)**2)
                
                
                vars()['z' + el + '_aver'][i][j] /= vars()['x' + el][i][j] # /= np.sum(10**xni_H[i][j]) # Normalizes with total abund.
                
                
                vars()['EM_' + el][i][j] = np.sum(vars()['EMi_' + el][i][j])
            
            
            vars()['EM_' + el + '_tot'][i] = vars()['EM_' + el][i][-1]
                
                
            

            vars()['x' + el + '_amb'].append([0.]) # NEITHER APPEND "i" NOR INTEGER ZERO!!!!!!!!
            vars()['x' + el + '_amb'][i] = (vars()['x' + el + '_amb'][i])*(np.shape(ti_amb[i])[0]) 
            vars()['x' + el + '_amb'][i] = np.asarray(vars()['x' + el + '_amb'][i])
            vars()['z' + el + '_aver_amb'].append([0.])
            vars()['z' + el + '_aver_amb'][i] = (vars()['z' + el + '_aver_amb'][i])*(np.shape(ti_amb[i])[0])
            
            #vars()['EMi_' + el + '_amb'].append([0.]) # NEITHER APPEND "i" NOR INTEGER ZERO!!!!!!!!
            #vars()['EMi_' + el + '_amb'][i] = np.zeros(np.shape(vars()['xni_' + el + '_amb'][i]))
            #vars()['EM_' + el + '_amb'].append([0.])
            #vars()['EM_' + el + '_amb'][i] = (vars()['EM_' + el + '_amb'][i])*(np.shape(ti_amb[i])[0])
            
            
            
            for k in range(np.shape(ti_amb[i])[0]):
            
                vars()['x' + el + '_amb'][i][k] = np.sum(10**vars()['xni_' + el + '_amb'][i][k])

                for b in range(vars()['nion_' + el]):

                    vars()['z' + el + '_aver_amb'][i][k] += b*10**vars()['xni_' + el + '_amb'][i][k][b]
                    
                    
                    #vars()['EMi_' + el + '_amb'][i][k][b] = \
                    #       np.trapz(y = np.asarray(10**np.asarray(vars()['xni_' + el + '_amb'][i]))[0:k][:,b] *\
                    #                np.asarray(10**np.asarray(xnee_amb[i]))[0:k] * 4*np.pi * r_amb[i][0:k]**2 , x = r_amb[i][0:k])
                
            
                vars()['z' + el + '_aver_amb'][i][k] /= vars()['x' + el + '_amb'][i][k] # /= np.sum(10**xni_amb_H[i][k]) 
                
                
                #vars()['EM_' + el + '_amb'][i][k] = np.sum(vars()['EMi_' + el + '_amb'][i][k])
                
            
            #vars()['EM_' + el + '_tot_amb'][i] = vars()['EM_' + el + '_amb'][i][-1]
     

            

            vars()['x' + el + '_amb'][i] = np.asarray(vars()['x' + el + '_amb'][i])
            vars()['z' + el + '_aver_amb'][i] = np.asarray(vars()['z' + el + '_aver_amb'][i])
        
        
        
        
        vars()['x' + el + '_amb'] = np.asarray(vars()['x' + el + '_amb'])
        vars()['z' + el + '_aver_amb'] = np.asarray(vars()['z' + el + '_aver_amb'])
    






    # Record final FITS table with model parameters for every mass coordinate
    # Each sub-table will contain information for a given age


    prihdr = fits.Header()
    prihdr['DATE'] = time.strftime("%m/%d/%Y")
    prihdr['COMMENT'] = '= Physical profile of a modeled SNR for several expansion ages. Includes shocked ejecta and ambient medium'
    prihdu = fits.PrimaryHDU(header=prihdr)


    cols = [0]*len(data_vs_rad_age)
    tbhdu = [0]*(len(data_vs_rad_age) + 1) # Header + cols
    tbhdu[0] = prihdu


    for ag, datag in enumerate(data_vs_rad_age):
        
        cols.append([0.])
        
        r_conc = np.concatenate((r[ag]/pc_to_cm, r[ag][-1]/pc_to_cm + r_amb[ag]/pc_to_cm), axis=0)
        r_Table = fits.Column(name = 'r', format='1E', unit='pc', array = r_conc)
        
        lagm_conc = np.concatenate((lagm[ag], lagm[ag][-1] + lagm_amb[ag]), axis=0)
        lagm_Table = fits.Column(name = 'lagm', format='1E', unit='Msun', array = lagm_conc)

        uxx_conc = np.concatenate((uxx[ag], uxx_amb[ag]), axis=0)
        uxx_Table = fits.Column(name = 'v', format='1E', unit='km+1s-1', array = uxx_conc)
        
        te_conc = np.concatenate((te[ag], te_amb[ag]), axis=0)
        te_Table = fits.Column(name = 'te', format='1E', unit='K', array = te_conc)
                                 
        rho_conc = np.concatenate((rho[ag], rho_amb[ag]), axis=0)
        rho_Table = fits.Column(name = 'rho', format='1E', unit='g+1cm-3', array = rho_conc)
                                  
        xnee_conc = np.concatenate((xnee[ag], xnee_amb[ag]), axis=0)
        xnee_Table = fits.Column(name = 'xnee', format='1E', unit='cm-3', array = xnee_conc)
                                   
        tau_conc = np.concatenate((tau[ag], tau_amb[ag]), axis=0)
        tau_Table = fits.Column(name = 'tau', format='1E', unit='cm-3s+1', array = tau_conc)
                                  
        ti_aver_conc = np.concatenate((ti_aver[ag], ti_aver_amb[ag]), axis=0)
        ti_aver_Table = fits.Column(name = 'ti_aver', format='1E', unit='K', array = ti_aver_conc)
        
        
        
        cols_Table = [r_Table, lagm_Table, uxx_Table, te_Table, rho_Table, xnee_Table, tau_Table, ti_aver_Table]
        
        
        
        for index, el in enumerate(elements):
        
            vars()['x' + el + '_conc'] = np.concatenate((vars()['x' + el][ag], vars()['x' + el + '_amb'][ag]), axis=0)
            vars()['x' + el + '_Table'] = fits.Column(name = 'x' + el, format='1E', unit='cm-3', \
                                                      array = vars()['x' + el + '_conc'])
                                                        
            vars()['z' + el + '_aver' + '_conc'] = np.concatenate((vars()['z' + el + '_aver'][ag],\
                                                                   vars()['z' + el + '_aver' + '_amb'][ag]), axis=0)
            vars()['z' + el + '_aver' + '_Table'] = fits.Column(name = 'z' + el + '_aver', format='1E', \
                                                                array = vars()['z' + el + '_aver' + '_conc'])
                                                        
            vars()['EM_' + el + '_conc'] = np.concatenate((vars()['EM_' + el][ag], np.zeros(len(r_amb[ag]))), axis=0)
            vars()['EM_' + el + '_Table'] = fits.Column(name = 'EM_' + el, format='1E', unit='cm-5', \
                                                        array = vars()['EM_' + el + '_conc'])
        
        
        
            cols_Table.append(vars()['x' + el + '_Table'])
            cols_Table.append(vars()['z' + el + '_aver' + '_Table'])
            cols_Table.append(vars()['EM_' + el + '_Table'])
            
            
        
        
        cols[ag] = fits.BinTableHDU.from_columns(cols_Table) 
        tbhdu[ag + 1] = fits.BinTableHDU.from_columns(cols[ag].columns, name = str(int(datag)) + 'yr', header = None)



    ## path_files.split(os.sep) will yield, for '/home/hector/CR_NEI_hydro/ChN_Ia_ddta_ifort_500yr_2p0/model_ddta/',
    ## ['', 'home', 'hector', 'CR_NEI_hydro', 'ChN_Ia_ddta_ifort_500yr_2p0', 'model_ddta', '']

    explmodel = path.split(os.sep)[-2].split("_")[-1]
    rhoISM = path.split(os.sep)[-3].split("_")[-1]


    # Write table
    tbhdulist = fits.HDUList(tbhdu)
    tbhdulist.writeto('SNRprofile_%s_%s.fits' % (explmodel, rhoISM), overwrite = True) 