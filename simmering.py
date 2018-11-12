
# Information about the thermonuclear runaway in Type Ia supernovae models. Creates a FITS table with all the information


# coding: utf-8

import numpy as np
import scipy as scipy
import scipy.stats as stats
import astropy.stats as astats
import numpy.random as random
from astropy.table import Table
import mesa as ms
import os



path=os.getcwd()  # Current work directory
#pathw='/../' # Write the table in a different directory



# Creates a title which includes some parts of the current work directory. I am removing all the "/" and other characters
title='_'#+os.path.dirname(path)[-22:-18]+'_'+os.path.dirname(path)[-12:-8]+'_'+os.path.dirname(path)[-7:] 


# Start an instance of the history.data




dir1=path+'/LOGS1'
dir2=path+'/LOGS2'
#dir1='/../'
#dir2='/../'
#dir1 = raw_input ('Please, insert the name of the directory where the first history.data and profiles are:\n')
#dir2 = raw_input ('Please, insert the name of the directory where the second history.data and profiles are:\n')


for path, dirs, files in os.walk(dir1): #os.walk avoids [Errno21] Is a directory
    m1=ms.history_data(dir1,clean_starlog=False)
for path, dirs, files in os.walk(dir2): 
    m2=ms.history_data(dir2,clean_starlog=False)





log_center_Rho_1=m1.get('log_center_Rho')
log_center_T_1=m1.get('log_center_T')
log_center_P_1=m1.get('log_center_P')
star_age_1=m1.get('star_age')
star_mass_1=m1.get('star_mass')
model_number_1=m1.get('model_number')
#log_L_1=m1.get('log_L')
#log_R_1=m1.get('log_R')
#log_Teff_1=m1.get('log_Teff')
#log_g_1=m1.get('log_g')

#center_h1_1=m1.get('center_h1')
#center_he4_1=m1.get('center_he4')
#center_c12_1=m1.get('center_c12')
#center_c13_1=m1.get('center_c13')
#center_n14_1=m1.get('center_n14')
#center_o16_1=m1.get('center_o16')
#center_o18_1=m1.get('center_o18')
#center_ne20_1=m1.get('center_ne20')
#center_ne22_1=m1.get('center_ne22')
#center_ne23_1=m1.get('center_ne23')
#center_na23_1=m1.get('center_na23')
#center_mg24_1=m1.get('center_mg24')
#center_si28_1=m1.get('center_si28')
#center_ni56_1=m1.get('center_ni56')
#center_gamma_1=m1.get('center_gamma')
center_ye_1=m1.get('center_ye')

#surface_h1_1=m1.get('surface_h1')
#surface_he4_1=m1.get('surface_he4')
#surface_c12_1=m1.get('surface_c12')
#surface_o16_1=m1.get('surface_o16')
#log_average_h1_1=m1.get('log_average_h1')
#log_average_he4_1=m1.get('log_average_he4')
#log_average_c12_1=m1.get('log_average_c12')
#log_average_o16_1=m1.get('log_average_o16')
#log_average_ne20_1=m1.get('log_average_ne20')

#burn_c_1=m1.get('burn_c')
#burn_o_1=m1.get('burn_o')
#c12_c12_1=m1.get('c12_c12')
#o16_o16_1=m1.get('o16_o16')
#burn_ne_1=m1.get('burn_ne')
#burn_mg_1=m1.get('burn_mg')

#conv_mx1_top_1=m1.get('conv_mx1_top')
#conv_mx1_bot_1=m1.get('conv_mx1_bot')
#conv_mx2_top_1=m1.get('conv_mx2_top')
#conv_mx2_bot_1=m1.get('conv_mx2_bot')
#mx1_top_1=m1.get('mx1_top')
#mx1_bot_1=m1.get('mx1_bot')
#mx2_top_1=m1.get('mx2_top')
#mx2_bot_1=m1.get('mx2_bot')

#epsnuc_M_1_1=m1.get('epsnuc_M_1')
#epsnuc_M_2_1=m1.get('epsnuc_M_2')
#epsnuc_M_3_1=m1.get('epsnuc_M_3')
#epsnuc_M_4_1=m1.get('epsnuc_M_4')
#epsnuc_M_5_1=m1.get('epsnuc_M_5')
#epsnuc_M_6_1=m1.get('epsnuc_M_6')
#epsnuc_M_7_1=m1.get('epsnuc_M_7')
#epsnuc_M_8_1=m1.get('epsnuc_M_8')

#c_core_mass_1=m1.get('c_core_mass')
#o_core_mass_1=m1.get('o_core_mass')

tdyn_1=m1.get('dynamic_timescale')
#tkh_1=m1.get('kh_timescale')
#tnuc_1=m1.get('nuc_timescale')

#log_LC_1=m1.get('log_LC')
#log_LZ_1=m1.get('log_LZ')
#log_Lnuc_1=m1.get('log_Lnuc')
#log_Lneu_1=m1.get('log_Lneu')
#log_Lneu_nuc_1=m1.get('log_Lneu_nuc')
#log_Lneu_nonnuc_1=m1.get('log_Lneu_nonnuc')
#mass_loc_of_max_eps_nuc_1=m1.get('mass_loc_of_max_eps_nuc')
#log_abs_Lgrav_1=m1.get('log_abs_Lgrav')

mass_conv_core_1=m1.get('mass_conv_core')
#cz_bot_mass_1=m1.get('cz_bot_mass')
#cz_top_mass_1=m1.get('cz_top_mass')
#cz_log_eps_nuc_1=m1.get('cz_log_eps_nuc')
#neutron_rich_core_mass_1=m1.get('neutron_rich_core_mass')
#envelope_mass_1=m1.get('envelope_mass')

#center_eps_grav_1=m1.get('center_eps_grav')
#center_non_nuc_neu_1=m1.get('center_non_nuc_neu')
#center_eps_nuc_1=m1.get('center_eps_nuc')

#total_mass_h1_1=m1.get('total_mass_h1')
#total_mass_he4_1=m1.get('total_mass_he4')
#total_mass_c12_1=m1.get('total_mass_c12')
#total_mass_o16_1=m1.get('total_mass_o16')

#total_eps_grav_1=m1.get('total_eps_grav')
total_nuclear_heating_1=m1.get('total_nuclear_heating')
total_non_nuc_neu_cooling_1=m1.get('total_non_nuc_neu_cooling')

#max_eps_nuc_1=m1.get('max_eps_nuc')
#max_eps_nuc_lgT_1=m1.get('max_eps_nuc_lgT')
#max_eps_nuc_lgRho_1=m1.get('max_eps_nuc_lgRho')
#max_eps_nuc_m_1=m1.get('max_eps_nuc_m')
#max_eps_nuc_xm_1=m1.get('max_eps_nuc_xm')
#max_eps_nuc_lgP_1=m1.get('max_eps_nuc_lgP')
#max_eps_nuc_lgR_1=m1.get('max_eps_nuc_lgR')
#max_eps_nuc_opacity_1=m1.get('max_eps_nuc_opacity')
#max_eps_nuc_cp_1=m1.get('max_eps_nuc_cp')
max_eps_nuc_t_heat_1=m1.get('max_eps_nuc_t_heat')
#max_eps_nuc_csound_1=m1.get('max_eps_nuc_csound')



log_center_Rho_2=m2.get('log_center_Rho')
log_center_T_2=m2.get('log_center_T')
log_center_P_2=m2.get('log_center_P')
star_age_2=m2.get('star_age')
star_mass_2=m2.get('star_mass')
model_number_2=m2.get('model_number')
#log_L_2=m2.get('log_L')
#log_R_2=m2.get('log_R')
#log_Teff_2=m2.get('log_Teff')
#log_g_2=m2.get('log_g')

#center_h1_2=m2.get('center_h1')
#center_he4_2=m2.get('center_he4')
#center_c12_2=m2.get('center_c12')
#center_c13_2=m2.get('center_c13')
#center_n14_2=m2.get('center_n14')
#center_o16_2=m2.get('center_o16')
#center_o18_2=m2.get('center_o18')
#center_ne20_2=m2.get('center_ne20')
#center_ne22_2=m2.get('center_ne22')
#center_ne23_2=m2.get('center_ne23')
#center_na23_2=m2.get('center_na23')
#center_mg24_2=m2.get('center_mg24')
#center_si28_2=m2.get('center_si28')
#center_ni56_2=m2.get('center_ni56')
#center_gamma_2=m2.get('center_gamma')
center_ye_2=m2.get('center_ye')

#surface_h1_2=m2.get('surface_h1')
#surface_he4_2=m2.get('surface_he4')
#surface_c12_2=m2.get('surface_c12')
#surface_o16_2=m2.get('surface_o16')
#log_average_h1_2=m2.get('log_average_h1')
#log_average_he4_2=m2.get('log_average_he4')
#log_average_c12_2=m2.get('log_average_c12')
#log_average_o16_2=m2.get('log_average_o16')
#log_average_ne20_2=m2.get('log_average_ne20')

#burn_c_2=m2.get('burn_c')
#burn_o_2=m2.get('burn_o')
#c12_c12_2=m2.get('c12_c12')
#o16_o16_2=m2.get('o16_o16')
#burn_ne_2=m2.get('burn_ne')
#burn_mg_2=m2.get('burn_mg')

#conv_mx1_top_2=m2.get('conv_mx1_top')
#conv_mx1_bot_2=m2.get('conv_mx1_bot')
#conv_mx2_top_2=m2.get('conv_mx2_top')
#conv_mx2_bot_2=m2.get('conv_mx2_bot')
#mx1_top_2=m2.get('mx1_top')
#mx1_bot_2=m2.get('mx1_bot')
#mx2_top_2=m2.get('mx2_top')
#mx2_bot_2=m2.get('mx2_bot')

#epsnuc_M_1_2=m2.get('epsnuc_M_1')
#epsnuc_M_2_2=m2.get('epsnuc_M_2')
#epsnuc_M_3_2=m2.get('epsnuc_M_3')
#epsnuc_M_4_2=m2.get('epsnuc_M_4')
#epsnuc_M_5_2=m2.get('epsnuc_M_5')
#epsnuc_M_6_2=m2.get('epsnuc_M_6')
#epsnuc_M_7_2=m2.get('epsnuc_M_7')
#epsnuc_M_8_2=m2.get('epsnuc_M_8')

#c_core_mass_2=m2.get('c_core_mass')
#o_core_mass_2=m2.get('o_core_mass')

tdyn_2=m2.get('dynamic_timescale')
#tkh_2=m2.get('kh_timescale')
#tnuc_2=m2.get('nuc_timescale')

#log_LC_2=m2.get('log_LC')
#log_LZ_2=m2.get('log_LZ')
#log_Lnuc_2=m2.get('log_Lnuc')
#log_Lneu_2=m2.get('log_Lneu')
#log_Lneu_nuc_2=m2.get('log_Lneu_nuc')
#log_Lneu_nonnuc_2=m2.get('log_Lneu_nonnuc')
#mass_loc_of_max_eps_nuc_2=m2.get('mass_loc_of_max_eps_nuc')
#log_abs_Lgrav_2=m2.get('log_abs_Lgrav')

mass_conv_core_2=m2.get('mass_conv_core')
#cz_bot_mass_2=m2.get('cz_bot_mass')
#cz_top_mass_2=m2.get('cz_top_mass')
#cz_log_eps_nuc_2=m2.get('cz_log_eps_nuc')
#neutron_rich_core_mass_2=m2.get('neutron_rich_core_mass')
#envelope_mass_2=m2.get('envelope_mass')

#center_eps_grav_2=m2.get('center_eps_grav')
#center_non_nuc_neu_2=m2.get('center_non_nuc_neu')
#center_eps_nuc_2=m2.get('center_eps_nuc')

#total_mass_h1_2=m2.get('total_mass_h1')
#total_mass_he4_2=m2.get('total_mass_he4')
#total_mass_c12_2=m2.get('total_mass_c12')
#total_mass_o16_2=m2.get('total_mass_o16')

#total_eps_grav_2=m2.get('total_eps_grav')
total_nuclear_heating_2=m2.get('total_nuclear_heating')
total_non_nuc_neu_cooling_2=m2.get('total_non_nuc_neu_cooling')

#max_eps_nuc_2=m2.get('max_eps_nuc')
#max_eps_nuc_lgT_2=m2.get('max_eps_nuc_lgT')
#max_eps_nuc_lgRho_2=m2.get('max_eps_nuc_lgRho')
#max_eps_nuc_m_2=m2.get('max_eps_nuc_m')
#max_eps_nuc_xm_2=m2.get('max_eps_nuc_xm')
#max_eps_nuc_lgP_2=m2.get('max_eps_nuc_lgP')
#max_eps_nuc_lgR_2=m2.get('max_eps_nuc_lgR')
#max_eps_nuc_opacity_2=m2.get('max_eps_nuc_opacity')
#max_eps_nuc_cp_2=m2.get('max_eps_nuc_cp')
max_eps_nuc_t_heat_2=m2.get('max_eps_nuc_t_heat')
#max_eps_nuc_csound_2=m2.get('max_eps_nuc_csound')



log_center_Rho=np.concatenate((log_center_Rho_1,log_center_Rho_2),axis=0)
log_center_T=np.concatenate((log_center_T_1,log_center_T_2),axis=0)
log_center_P=np.concatenate((log_center_P_1,log_center_P_2),axis=0)
star_age=np.concatenate((star_age_1,star_age_2),axis=0)
star_mass=np.concatenate((star_mass_1,star_mass_2),axis=0)
#model_number=np.concatenate((model_number_1,model_number_2),axis=0)
#log_L=np.concatenate((log_L_1,log_L_2),axis=0)
#log_R=np.concatenate((log_R_1,log_R_2),axis=0)
#log_Teff=np.concatenate((log_Teff_1,log_Teff_2),axis=0)
#log_g=np.concatenate((log_g_1,log_g_2),axis=0)

#center_h1=np.concatenate((center_h1_1,center_h1_2),axis=0)
#center_he4=np.concatenate((center_he4_1,center_he4_2),axis=0)
#center_c12=np.concatenate((center_c12_1,center_c12_2),axis=0)
#center_c13=np.concatenate((center_c13_1,center_c13_2),axis=0)
#center_n14=np.concatenate((center_n14_1,center_n14_2),axis=0)
#center_o16=np.concatenate((center_o16_1,center_o16_2),axis=0)
#center_o18=np.concatenate((center_o18_1,center_o18_2),axis=0)
#center_ne20=np.concatenate((center_ne20_1,center_ne20_2),axis=0)
#center_ne22=np.concatenate((center_ne22_1,center_ne22_2),axis=0)
#center_ne23=np.concatenate((center_ne23_1,center_ne23_2),axis=0)
#center_na23=np.concatenate((center_na23_1,center_na23_2),axis=0)
#center_mg24=np.concatenate((center_mg24_1,center_mg24_2),axis=0)
#center_si28=np.concatenate((center_si28_1,center_si28_2),axis=0)
#center_gamma=np.concatenate((center_gamma_1,center_gamma_2),axis=0)
center_ye=np.concatenate((center_ye_1,center_ye_2),axis=0)

#surface_h1=np.concatenate((surface_h1_1,surface_h1_2),axis=0)
#surface_he4=np.concatenate((surface_he4_1,surface_he4_2),axis=0)
#surface_c12=np.concatenate((surface_c12_1,surface_c12_2),axis=0)
#surface_o16=np.concatenate((surface_o16_1,surface_o16_2),axis=0)
#log_average_h1=np.concatenate((log_average_h1_1,log_average_h1_2),axis=0)
#log_average_he4=np.concatenate((log_average_he4_1,log_average_he4_2),axis=0)
#log_average_c12=np.concatenate((log_average_c12_1,log_average_c12_2),axis=0)
#log_average_o16=np.concatenate((log_average_o16_1,log_average_o16_2),axis=0)
#log_average_ne20=np.concatenate((log_average_ne20_1,log_average_ne20_2),axis=0)

#burn_c=np.concatenate((burn_c_1,burn_c_2),axis=0)
#burn_o=np.concatenate((burn_o_1,burn_o_2),axis=0)
#c12_c12=np.concatenate((c12_c12_1,c12_c12_2),axis=0)
#o16_o16=np.concatenate((o16_o16_1,o16_o16_2),axis=0)
#burn_ne=np.concatenate((burn_ne_1,burn_ne_2),axis=0)
#burn_mg=np.concatenate((burn_mg_1,burn_mg_2),axis=0)

#conv_mx1_top = np.concatenate((conv_mx1_top_1,conv_mx1_top_2),axis=0)
#conv_mx1_bot = np.concatenate((conv_mx1_bot_1,conv_mx1_bot_2),axis=0)
#conv_mx2_top = np.concatenate((conv_mx2_top_1,conv_mx2_top_2),axis=0)
#conv_mx2_bot = np.concatenate((conv_mx2_bot_1,conv_mx2_bot_2),axis=0)
#mx1_top = np.concatenate((mx1_top_1,mx1_top_2),axis=0)
#mx1_bot = np.concatenate((mx1_bot_1,mx1_bot_2),axis=0)
#mx2_top = np.concatenate((mx2_top_1,mx2_top_2),axis=0)
#mx2_bot = np.concatenate((mx2_bot_1,mx2_bot_2),axis=0)

#epsnuc_M_1 = np.concatenate((epsnuc_M_1_1,epsnuc_M_1_2),axis=0)
#epsnuc_M_2 = np.concatenate((epsnuc_M_2_1,epsnuc_M_2_2),axis=0)
#epsnuc_M_3 = np.concatenate((epsnuc_M_3_1,epsnuc_M_3_2),axis=0)
#epsnuc_M_4 = np.concatenate((epsnuc_M_4_1,epsnuc_M_4_2),axis=0)
#epsnuc_M_5 = np.concatenate((epsnuc_M_5_1,epsnuc_M_5_2),axis=0)
#epsnuc_M_6 = np.concatenate((epsnuc_M_6_1,epsnuc_M_6_2),axis=0)
#epsnuc_M_7 = np.concatenate((epsnuc_M_7_1,epsnuc_M_7_2),axis=0)
#epsnuc_M_8 = np.concatenate((epsnuc_M_8_1,epsnuc_M_8_2),axis=0)

#c_core_mass = np.concatenate((c_core_mass_1,c_core_mass_2),axis=0)
#o_core_mass = np.concatenate((o_core_mass_1,o_core_mass_2),axis=0)

tdyn = np.concatenate((tdyn_1,tdyn_2),axis=0)
#tkh = np.concatenate((tkh_1,tkh_2),axis=0)
#tnuc = np.concatenate((tnuc_1,tnuc_2),axis=0)

#log_LC = np.concatenate((log_LC_1,log_LC_2),axis=0)
#log_LZ = np.concatenate((log_LZ_1,log_LZ_2),axis=0)
#log_Lnuc = np.concatenate((log_Lnuc_1,log_Lnuc_2),axis=0)
#log_Lneu = np.concatenate((log_Lneu_1,log_Lneu_2),axis=0)
#log_Lneu_nuc = np.concatenate((log_Lneu_nuc_1,log_Lneu_nuc_2),axis=0)
#log_Lneu_nonnuc = np.concatenate((log_Lneu_nonnuc_1,log_Lneu_nonnuc_2),axis=0)
#mass_loc_of_max_eps_nuc = np.concatenate((mass_loc_of_max_eps_nuc_1,mass_loc_of_max_eps_nuc_2),axis=0)
#log_abs_Lgrav = np.concatenate((log_abs_Lgrav_1,log_abs_Lgrav_2),axis=0)

mass_conv_core = np.concatenate((mass_conv_core_1,mass_conv_core_2),axis=0)
#cz_bot_mass = np.concatenate((cz_bot_mass_1,cz_bot_mass_2),axis=0)
#cz_top_mass = np.concatenate((cz_top_mass_1,cz_top_mass_2),axis=0)
#cz_log_eps_nuc = np.concatenate((cz_log_eps_nuc_1,cz_log_eps_nuc_2),axis=0)
#neutron_rich_core_mass = np.concatenate((neutron_rich_core_mass_1,neutron_rich_core_mass_2),axis=0)
#envelope_mass = np.concatenate((envelope_mass_1,envelope_mass_2),axis=0)

#center_eps_grav = np.concatenate((center_eps_grav_1,center_eps_grav_2),axis=0)
#center_non_nuc_neu = np.concatenate((center_non_nuc_neu_1,center_non_nuc_neu_2),axis=0)
#center_eps_nuc = np.concatenate((center_eps_nuc_1,center_eps_nuc_2),axis=0)

#total_mass_h1 = np.concatenate((total_mass_h1_1,total_mass_h1_2),axis=0)
#total_mass_he4 = np.concatenate((total_mass_he4_1,total_mass_he4_2),axis=0)
#total_mass_c12 = np.concatenate((total_mass_c12_1,total_mass_c12_2),axis=0)
#total_mass_o16 = np.concatenate((total_mass_o16_1,total_mass_o16_2),axis=0)

#total_eps_grav = np.concatenate((total_eps_grav_1,total_eps_grav_2),axis=0)
total_nuclear_heating = np.concatenate((total_nuclear_heating_1,total_nuclear_heating_2),axis=0)
total_non_nuc_neu_cooling = np.concatenate((total_non_nuc_neu_cooling_1,total_non_nuc_neu_cooling_2),axis=0)

#max_eps_nuc = np.concatenate((max_eps_nuc_1,max_eps_nuc_2),axis=0)
#max_eps_nuc_lgT = np.concatenate((max_eps_nuc_lgT_1,max_eps_nuc_lgT_2),axis=0)
#max_eps_nuc_lgRho = np.concatenate((max_eps_nuc_lgRho_1,max_eps_nuc_lgRho_2),axis=0)
#max_eps_nuc_m = np.concatenate((max_eps_nuc_m_1,max_eps_nuc_m_2),axis=0)
#max_eps_nuc_xm = np.concatenate((max_eps_nuc_xm_1,max_eps_nuc_xm_2),axis=0)
#max_eps_nuc_lgP = np.concatenate((max_eps_nuc_lgP_1,max_eps_nuc_lgP_2),axis=0)
#max_eps_nuc_lgR = np.concatenate((max_eps_nuc_lgR_1,max_eps_nuc_lgR_2),axis=0)
#max_eps_nuc_opacity = np.concatenate((max_eps_nuc_opacity_1,max_eps_nuc_opacity_2),axis=0)
#max_eps_nuc_cp = np.concatenate((max_eps_nuc_cp_1,max_eps_nuc_cp_2),axis=0)
max_eps_nuc_t_heat = np.concatenate((max_eps_nuc_t_heat_1,max_eps_nuc_t_heat_2),axis=0)
#max_eps_nuc_csound = np.concatenate((max_eps_nuc_csound_1,max_eps_nuc_csound_2),axis=0)




eta=1.-2.*center_ye
eta_sun=1.4E-3




## Profiles


# Start from the first profile and find the corresponding model_number
prof_index_1=Table.read(dir1+'/profiles.index',format='ascii.fixed_width')
prof_index_2=Table.read(dir2+'/profiles.index',format='ascii.fixed_width')



idx1=[]
idx1a=[]
idx1b=[]
idx1c=[]
index1=[]
for i in range(len(prof_index_1)):
    idx1.append(i)
    idx1a.append(i)
    idx1b.append(i)
    idx1c.append(i)
    index1.append(i)

    idx1[i]=(prof_index_1[i][0][0:])
    idx1a[i],idx1b[i],idx1c[i]=idx1[i].split() # Take the elements from the string
    idx1a[i]=float(idx1a[i]) # Convert them to floating point numbers
    index1[i]=np.where(model_number_1==idx1a[i])[0][0]


idx2=[]
idx2a=[]
idx2b=[]
idx2c=[]
index2=[]
for i in range(len(prof_index_2)):
    idx2.append(i)
    idx2a.append(i)
    idx2b.append(i)
    idx2c.append(i)
    index2.append(i)

    idx2[i]=(prof_index_2[i][0][0:])
    idx2a[i],idx2b[i],idx2c[i]=idx2[i].split() # Take the elements from the string
    idx2a[i]=float(idx2a[i]) # Convert them to floating point numbers
    index2[i]=np.where(model_number_2==idx2a[i])[0][0]

mod_number_1=model_number_1[index1]
st_age_1=star_age_1[index1]

mod_number_2=model_number_2[index2]
st_age_2=star_age_2[index2]

st_age=np.concatenate((st_age_1,st_age_2),axis=0)



prof_1=np.empty((len(mod_number_1),1), dtype=object) #error return without exception set avoided with dtype=object
for i in range(len(mod_number_1)):
    for path, dirs, files in os.walk(dir1): #os.walk avoids [Errno21] Is a directory
        prof_1[i:,]=ms.mesa_profile(dir1,mod_number_1[i])




prof_2=np.empty((len(mod_number_2),1), dtype=object) #error return without exception set avoided with dtype=object
for i in range(len(mod_number_2)):
    for path, dirs, files in os.walk(dir2): #os.walk avoids [Errno21] Is a directory
        prof_2[i:,]=ms.mesa_profile(dir2,mod_number_2[i])




prof=np.concatenate((prof_1,prof_2),axis=0)[:,0] # Remove "dtype=object" from the final profiles to get the information


# Profile variables. Add as many as desired


q=[]
logT=[]
gradr_sub_grada=[]

conv=[]
conv_core=[]

ccsize=[]


for i in range(len(prof)):
    
    q.append(i)
    q[i]=prof[i].get('q')

    logT.append(i)
    logT[i]=prof[i].get('logT')
    
    gradr_sub_grada.append(i)
    gradr_sub_grada[i]=prof[i].get('gradr_sub_grada')

    conv.append(i)
    conv_core.append(i)
    conv[i]=np.where(gradr_sub_grada[i]>0)[0]
    conv_core[i]=np.where((gradr_sub_grada[i]>0)&(q[i]<0.9))[0]

    ccsize.append(i)
    ccsize[i]=conv_core[i].size 
cc=[i for i,x in enumerate(ccsize) if x !=0][0] #Find the first index where convection happens (.size=0 for an empty array)


lenarr=[]
for i in range(len(q)):
    lenarr.append(i)
    lenarr[i]=len(q[i])
maxn=np.maximum.accumulate(lenarr)
maxlength=maxn[-1]


age=[]
for i, item in enumerate(lenarr):
    age.append(i)
    age[i]=[st_age[i]]*lenarr[i] # Replicates every element "i" in st_age "lenarr[i]" times


# Several important definitions

age_ign=age[cc][0] 
c_ign=np.where(mass_conv_core>0)[0][0] # Index where carbon ignition occurs in the history data
#c_ign=np.where(star_age>age_ign)[0][0] # Index where carbon ignition occurs in the history data
eta_ign=eta[c_ign]
mass_ign=star_mass[c_ign]
max_age=np.max(star_age)
#max_age=np.max(st_age)


# Isothermal temperature. It is defined as the one existent in the boundary between the convective
# core and the conductive WD (where the "-1" arises) at the beginning of carbon simmering

logTiso=logT[cc][conv_core[cc][0]-1]




# Final table


# Creates a title which includes some parts of the current work directory. I am removing all the "/" and other characters
title='_'+os.path.dirname(path)[-22:-18]+'_'+os.path.dirname(path)[-12:-8]+'_'+os.path.dirname(path)[-7:] 



Mass = "%.1f"%(star_mass[0])
Rate = 1.E-7#float(os.path.dirname(path)[-22]+'.'+os.path.dirname(path)[-21:-18])
Wind_C_ab = 0.34#float(os.path.dirname(path)[-11]+'.'+os.path.dirname(path)[-10:-8])
Sim_age = "%.0f"%(age_ign)
Sim_logTc = "%.2f"%(log_center_T[c_ign]) 
Sim_logRhoc = "%.2f"%(log_center_Rho[c_ign]) 
Sim_mass = "%.4f"%(mass_ign)
Elapsed_time = "%.0f"%(max_age-age_ign)
Final_age = "%.0f"%(max_age)
Final_mass = "%.4f"%(star_mass[-1])
Final_mass_cc = "%.4f"%(mass_conv_core[-1])
Frac_final_mass_cc = "%.4f"%(mass_conv_core[-1]/star_mass[-1])
Final_eta = "%.5f"%(eta[-1])
Incr_eta_sun = "%.2f"%((eta[-1]-eta_ign)/eta_sun)
th = "%.2f"%(max_eps_nuc_t_heat[-1])
tdyn = "%.2f"%(tdyn[-1])
logTiso = "%.4f"%(logTiso)



val=np.array([Mass, Rate, Wind_C_ab,Sim_age,Sim_logTc,Sim_logRhoc,Sim_mass,Elapsed_time,\
Final_age,Final_mass,Final_mass_cc,Frac_final_mass_cc,Final_eta,Incr_eta_sun,th,tdyn,logTiso])
nam=['Mass','Rate','Wind_C_ab','Sim_age','Sim_logTc','Sim_logRhoc','Sim_mass','Elapsed_time',\
	'Final_age','Final_mass','Final_mass_cc','Frac_final_mass_cc','Final_eta','Incr_eta_sun','th','tdyn','logTiso']
un=['Msun','Msun/yr','','yr','','','Msun','yr','yr','Msun','Msun','','','etasun','s','s','K']


t=Table(val,names=nam)
for i in range(0,len(un)):
	t[nam[i]].unit=un[i]

t.write("thermo_runaway"+title+".txt",format='ascii.fixed_width')  
#t.write(pathw+"thermo_runaway"+title+".txt",format='ascii.fixed_width') 