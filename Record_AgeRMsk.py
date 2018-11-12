# coding: utf-8

# Record age, radius, and shocked masses (RS and FS) for each crhydronei directory


import numpy as np
from astropy.table import Table, vstack, Column
import os, sys, time
from scipy import interpolate



pathw = os.getcwd()  # Current work directory


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




pc_to_cm = 3.086e18
Msun = 1.98E33


nam = ['Age', 'R_FS', 'R_RS', 'M_sk', 'M_sk_amb']
un = ['yr', 'pc', 'pc', 'Msun', 'Msun']




for idx, path in enumerate(path_files):
    
    # Assume nshocks is not equal for each model
    data_sk_vs_rad = np.genfromtxt(path + 'data_sk_vs_rad.dat', skip_footer=2) # Skip last two rows
    sk_ct = data_sk_vs_rad[:,0]
    age = data_sk_vs_rad[:,2] # time_yr_run
    R_FS = data_sk_vs_rad[:,17] # FS_rad_pc
    R_RS = data_sk_vs_rad[:,21] # RS_rad_pc
    nshocks = sk_ct[-1]
    
    
    
    # Shocked mass for each age
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
    mshocked = np.zeros(len(data_vs_rad_age))
    
    r_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    dm_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    lagm_amb = [ []*i for i in range(len(data_vs_rad_age)) ]
    mshocked_amb = np.zeros(len(data_vs_rad_age))


    
    # Different from nshocks...for instance, for nshocks = 30, len(data_vs_rad_age) = 21, whereas len(data_vs_rad[:,19]) = 41874
    # But the first value for data_vs_rad_age was 0, so I dropped it
    data_vs_rad_age_block = [ []*i for i in range(len(data_vs_rad_age)) ]
        
        
        
        
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
    
    
    
    
        for j, el in enumerate(r[ag]): #RS
        
            if j == 0:
            
                lagm[ag][j] = dm[ag][j] * el / Msun
            
            else:
            
                # TRAPEZOIDAL RULE FOR INTEGRATING
                lagm[ag][j] = lagm[ag][j-1] + 0.5 * (dm[ag][j]+dm[ag][j-1]) * (r[ag][j]-r[ag][j-1]) / Msun
            
            
        mshocked[ag] = lagm[ag][-1] - lagm[ag][data_vs_rad_iRS[ag]]
 


        for i in range(1, len(mshocked)): # In case the RS rebounces against the center
    
            if mshocked[i] < mshocked[i-1]:
        
                mshocked[i] = mshocked[i-1]
            
            
            
            
            
        for o, al in enumerate(r_amb[ag]): #FS
        
            if o == 0:
            
                lagm_amb[ag][o] = 0.
            
            else:
            
                # TRAPEZOIDAL RULE FOR INTEGRATING
                lagm_amb[ag][o] = lagm_amb[ag][o-1] + 0.5*(dm_amb[ag][o]+dm_amb[ag][o-1])*(r_amb[ag][o]-r_amb[ag][o-1]) / Msun
            
            
        mshocked_amb[ag] = lagm_amb[ag][-1]
    





        # INTERPOLATE TO FIND MASS FOR EACH SHOCK
        
        interp_m_shocked = interpolate.interp1d(data_vs_rad_age, mshocked)
        mskinterp = np.zeros(len(age))
        
        interp_m_shocked_amb = interpolate.interp1d(data_vs_rad_age, mshocked_amb)
        mskinterp_amb = np.zeros(len(age))
        
        
    for n in range(len(age)):
        
        if n > 0:
            # data_vs_rad_age[0] might be bigger than age[0]
            mskinterp[n] = interp_m_shocked(age[n])
            mskinterp_amb[n] = interp_m_shocked_amb(age[n])
                    
        else: #Shocked mass for first recorded shock (age[n] = age[0])
            mskinterp[n] = interpolate.interp1d([0., age[n], age[n+1]],\
             [0., mshocked[n], mshocked[n+1]])(age[n])
            mskinterp_amb[n] = interpolate.interp1d([0., age[n], age[n+1]],\
             [0., mshocked_amb[n], mshocked_amb[n+1]])(age[n])
        
        
    
    

    ## path_files.split(os.sep) will yield, for '/home/hector/CR_NEI_hydro/ChN_Ia_ddta_ifort_500yr_2p0/model_ddta/',
    ## ['', 'home', 'hector', 'CR_NEI_hydro', 'ChN_Ia_ddta_ifort_500yr_2p0', 'model_ddta', '']

    explmodel = path.split(os.sep)[-2].split("_")[-1]
    rhoISM = path.split(os.sep)[-3].split("_")[-1]
    
    
    cols = [age, R_FS, R_RS, mskinterp, mskinterp_amb]
    
    
    # Record final tables
    tabl = Table(cols, names = nam)
    
    for p in range(len(nam)):
        
        tabl[nam[p]].unit = un[p]
        
    tabl.write('AgeRMsk_%s_%s.fits' % (explmodel, rhoISM), overwrite = True)