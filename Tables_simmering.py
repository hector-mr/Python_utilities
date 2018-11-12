
# Information about the thermonuclear runaway in Type Ia supernovae models. Creates a table with all the information and the group of subdirectories


# coding: utf-8

import numpy as np
from astropy.table import Table,vstack,Column
import mesa as ms
import os


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
for i in range(len([x[1] for x in os.walk(pathw)][0])):
	subdirect.append(i)
	subdirect[i] = [x[1] for x in os.walk(pathw)][0][i] 	
	# List of all subdirectories. Change the parameter [0] depending on the needs (i.e., on the depth of the nested grid of subdirectories)
	# It should always be x[1], which corresponds to the directories present in the current work path. x[2] comprises the files present in the current work path
	
subdirect = sorted(subdirect, key=str.lower) 	# Otherwise, capitalized strings go first and the vector becomes a mess


direct = []
for i in range(len(subdirect)):
	direct.append(i)
	direct[i] = pathw + '/' + subdirect[i]



log_center_Rho_1 = []
log_center_T_1 = []
log_center_P_1 = []
star_age_1 = []
star_mass_1 = []
model_number_1 = []
#log_L_1 = []
#log_R_1 = []
#log_Teff_1 = []
#log_g_1 = []
#center_h1_1 = []
#center_he4_1 = []
#center_c12_1 = []
#center_c13_1 = []
#center_n14_1 = []
#center_o16_1 = []
#center_o18_1 = []
#center_ne20_1 = []
#center_ne22_1 = []
#center_ne23_1 = []
#center_na23_1 = []
#center_mg24_1 = []
#center_si28_1 = []
#center_ni56_1 = []
#center_gamma_1 = []
center_ye_1 = []
#surface_h1_1 = []
#surface_he4_1 = []
#surface_c12_1 = []
#surface_o16_1 = []
#log_average_h1_1 = []
#log_average_he4_1 = []
#log_average_c12_1 = []
#log_average_o16_1 = []
#log_average_ne20_1 = []
#burn_c_1 = []
#burn_o_1 = []
#c12_c12_1 = []
#o16_o16_1 = []
#burn_ne_1 = []
#burn_mg_1 = []
#conv_mx1_top_1 = []
#conv_mx1_bot_1 = []
#conv_mx2_top_1 = []
#conv_mx2_bot_1 = []
#mx1_top_1 = []
#mx1_bot_1 = []
#mx2_top_1 = []
#mx2_bot_1 = []
#epsnuc_M_1_1 = []
#epsnuc_M_2 = []
#epsnuc_M_3_1 = []
#epsnuc_M_4_1 = []
#epsnuc_M_5_1 = []
#epsnuc_M_6_1 = []
#epsnuc_M_7_1 = []
#epsnuc_M_8_1 = []
#c_core_mass_1 = []
#o_core_mass_1 = []
tdyn_1 = []
#tkh_1 = []
#tnuc_1 = []
#log_LC_1 = []
#log_LZ_1 = []
#log_Lnuc_1 = []
#log_Lneu_1 = []
#log_Lneu_nuc_1 = []
#log_Lneu_nonnuc_1 = []
#mass_loc_of_max_eps_nuc_1 = []
#log_abs_Lgrav_1 = []
mass_conv_core_1 = []
#cz_bot_mass_1 = []
#cz_top_mass_1 = []
#cz_log_eps_nuc_1 = []
#neutron_rich_core_mass_1 = []
#envelope_mass_1 = []
#center_eps_grav_1 = []
#center_non_nuc_neu_1 = []
#center_eps_nuc_1 = []
#total_mass_h1_1 = []
#total_mass_he4_1 = []
#total_mass_c12_1 = []
#total_mass_o16_1 = []
#total_eps_grav_1 = []
total_nuclear_heating_1 = []
total_non_nuc_neu_cooling_1 = []
#max_eps_nuc_1 = []
#max_eps_nuc_lgT_1 = []
#max_eps_nuc_lgRho_1 = []
#max_eps_nuc_m_1 = []
#max_eps_nuc_xm_1 = []
#max_eps_nuc_lgP_1 = []
#max_eps_nuc_lgR_1 = []
#max_eps_nuc_opacity_1 = []
#max_eps_nuc_cp_1 = []
max_eps_nuc_t_heat_1 = []
#max_eps_nuc_csound_1 = []
log_center_Rho_2 = []
log_center_T_2 = []
log_center_P_2 = []
star_age_2 = []
star_mass_2 = []
model_number_2 = []
#log_L_2 = []
#log_R_2 = []
#log_Teff_2 = []
#log_g_2 = []
#center_h1_2 = []
#center_he4_2 = []
#center_c12_2 = []
#center_c13_2 = []
#center_n14_2 = []
#center_o16_2 = []
#center_o18_2 = []
#center_ne20_2 = []
#center_ne22_2 = []
#center_ne23_2 = []
#center_na23_2 = []
#center_mg24_2 = []
#center_si28_2 = []
#center_ni56_2 = []
#center_gamma_2 = []
center_ye_2 = []
#surface_h1_2 = []
#surface_he4_2 = []
#surface_c12_2 = []
#surface_o16_2 = []
#log_average_h1_2 = []
#log_average_he4_2 = []
#log_average_c12_2 = []
#log_average_o16_2 = []
#log_average_ne20_2 = []
#burn_c_2 = []
#burn_o_2 = []
#c12_c12_2 = []
#o16_o16_2 = []
#burn_ne_2 = []
#burn_mg_2 = []
#conv_mx1_top_2 = []
#conv_mx1_bot_2 = []
#conv_mx2_top_2 = []
#conv_mx2_bot_2 = []
#mx1_top_2 = []
#mx1_bot_2 = []
#mx2_top_2 = []
#mx2_bot_2 = []
#epsnuc_M_1_2 = []
#epsnuc_M_2 = []
#epsnuc_M_3_2 = []
#epsnuc_M_4_2 = []
#epsnuc_M_5_2 = []
#epsnuc_M_6_2 = []
#epsnuc_M_7_2 = []
#epsnuc_M_8_2 = []
#c_core_mass_2 = []
#o_core_mass_2 = []
tdyn_2 = []
#tkh_2 = []
#tnuc_2 = []
#log_LC_2 = []
#log_LZ_2 = []
#log_Lnuc_2 = []
#log_Lneu_2 = []
#log_Lneu_nuc_2 = []
#log_Lneu_nonnuc_2 = []
#mass_loc_of_max_eps_nuc_2 = []
#log_abs_Lgrav_2 = []
mass_conv_core_2 = []
#cz_bot_mass_2 = []
#cz_top_mass_2 = []
#cz_log_eps_nuc_2 = []
#neutron_rich_core_mass_2 = []
#envelope_mass_2 = []
#center_eps_grav_2 = []
#center_non_nuc_neu_2 = []
#center_eps_nuc_2 = []
#total_mass_h1_2 = []
#total_mass_he4_2 = []
#total_mass_c12_2 = []
#total_mass_o16_2 = []
#total_eps_grav_2 = []
total_nuclear_heating_2 = []
total_non_nuc_neu_cooling_2 = []
#max_eps_nuc_2 = []
#max_eps_nuc_lgT_2 = []
#max_eps_nuc_lgRho_2 = []
#max_eps_nuc_m_2 = []
#max_eps_nuc_xm_2 = []
#max_eps_nuc_lgP_2 = []
#max_eps_nuc_lgR_2 = []
#max_eps_nuc_opacity_2 = []
#max_eps_nuc_cp_2 = []
max_eps_nuc_t_heat_2 = []
#max_eps_nuc_csound_2 = []
log_center_Rho = []
log_center_T = []
log_center_P = []
star_age = []
star_mass = []
#model_number = []
#log_L = []
#log_R = []
#log_Teff = []
#log_g = []
#center_h1 = []
#center_he4 = []
#center_c12 = []
#center_c13 = []
#center_n14 = []
#center_o16 = []
#center_o18 = []
#center_ne20 = []
#center_ne22 = []
#center_ne23 = []
#center_na23 = []
#center_mg24 = []
#center_si28 = []
#center_gamma = []
center_ye = []
#surface_h1 = []
#surface_he4 = []
#surface_c12 = []
#surface_o16 = []
#log_average_h1 = []
#log_average_he4 = []
#log_average_c12 = []
#log_average_o16 = []
#log_average_ne20 = []
#burn_c = []
#burn_o = []
#c12_c12 = []
#o16_o16 = []
#burn_ne = []
#burn_mg = []
#conv_mx1_top = []
#conv_mx1_bot = []
#conv_mx2_top = []
#conv_mx2_bot = []
#mx1_top = []
#mx1_bot = []
#mx2_top = []
#mx2_bot = []
#epsnuc_M_1 = []
#epsnuc_M_2 = []
#epsnuc_M_3 = []
#epsnuc_M_4 = []
#epsnuc_M_5 = []
#epsnuc_M_6 = []
#epsnuc_M_7 = []
#epsnuc_M_8 = []
#c_core_mass = []
#o_core_mass = []
tdyn = []
#tkh = []
#tnuc = []
#log_LC = []
#log_LZ = []
#log_Lnuc = []
#log_Lneu = []
#log_Lneu_nuc = []
#log_Lneu_nonnuc = []
#mass_loc_of_max_eps_nuc = []
#log_abs_Lgrav = []
mass_conv_core = []
#cz_bot_mass = []
#cz_top_mass = []
#cz_log_eps_nuc = []
#neutron_rich_core_mass = []
#envelope_mass = []
#center_eps_grav = []
#center_non_nuc_neu = []
#center_eps_nuc = []
#total_mass_h1 = []
#total_mass_he4 = []
#total_mass_c12 = []
#total_mass_o16 = []
#total_eps_grav = []
total_nuclear_heating = []
total_non_nuc_neu_cooling = []
#max_eps_nuc = []
#max_eps_nuc_lgT = []
#max_eps_nuc_lgRho = []
#max_eps_nuc_m = []
#max_eps_nuc_xm = []
#max_eps_nuc_lgP = []
#max_eps_nuc_lgR = []
#max_eps_nuc_opacity = []
#max_eps_nuc_cp = []
max_eps_nuc_t_heat = []
#max_eps_nuc_csound = []
eta = []
prof_index_1 = []
prof_index_2 = []
idx1 = []
idx1a = []
idx1b = []
idx1c = []
index1 = []
idx2 = []
idx2a = []
idx2b = []
idx2c = []
index2 = []
mod_number_1 = []
st_age_1 = []
mod_number_2 = []
st_age_2 = []
st_age = []
prof_1 = []
prof_2 = []
prof = []
#q = []
#logT = []
#gradr_sub_grada = []
#conv = []
#conv_core = []
#ccsize = []
#cc = []
#lenarr = []
#maxn = []
#maxlength = []
age = []
age_ign = []
c_ign = []
eta_ign = []
mass_ign = []
max_age = []
logTiso = []
title = []
Z = []
Mass = []
Rate = []
Sim_age = []
Sim_logTc = []
Sim_logRhoc = []
Sim_mass = []
Elapsed_time = []
Final_age = []
Final_logRhoc = []
Final_mass = []
Final_mass_cc = []
Frac_final_mass_cc = []
Final_eta = []
Diff_eta = []
Incr_eta_sun = []
th = []
tdyn = []
#logTiso = []
val = []
t = []
t_tex = []






dir1 = []
dir2 = []
m1 = []
m2 = []
for j in range(len(direct)):
    dir1.append(j)
    dir2.append(j)
    dir1[j]=direct[j]+'/LOGS1'
    dir2[j]=direct[j]+'/LOGS2'

    m1.append(j)
    m2.append(j)


    for path, dirs, files in os.walk(dir1[j]): #os.walk avoids [Errno21] Is a directory
        m1[j]=ms.history_data(dir1[j],clean_starlog=False)
    for path, dirs, files in os.walk(dir2[j]): 
        m2[j]=ms.history_data(dir2[j],clean_starlog=False)



    log_center_Rho_1.append(j)
    log_center_T_1.append(j)
    log_center_P_1.append(j)
    star_age_1.append(j)
    star_mass_1.append(j)
    model_number_1.append(j)
    #log_L_1.append(j)
    #log_R_1.append(j)
    #log_Teff_1.append(j)
    #log_g_1.append(j)
    #center_h1_1.append(j)
    #center_he4_1.append(j)
    #center_c12_1.append(j)
    #center_c13_1.append(j)
    #center_n14_1.append(j)
    #center_o16_1.append(j)
    #center_o18_1.append(j)
    #center_ne20_1.append(j)
    #center_ne22_1.append(j)
    #center_ne23_1.append(j)
    #center_na23_1.append(j)
    #center_mg24_1.append(j)
    #center_si28_1.append(j)
    #center_ni56_1.append(j)
    #center_gamma_1.append(j)
    center_ye_1.append(j)
    #surface_h1_1.append(j)
    #surface_he4_1.append(j)
    #surface_c12_1.append(j)
    #surface_o16_1.append(j)
    #log_average_h1_1.append(j)
    #log_average_he4_1.append(j)
    #log_average_c12_1.append(j)
    #log_average_o16_1.append(j)
    #log_average_ne20_1.append(j)
    #burn_c_1.append(j)
    #burn_o_1.append(j)
    #c12_c12_1.append(j)
    #o16_o16_1.append(j)
    #burn_ne_1.append(j)
    #burn_mg_1.append(j)
    #conv_mx1_top_1.append(j)
    #conv_mx1_bot_1.append(j)
    #conv_mx2_top_1.append(j)
    #conv_mx2_bot_1.append(j)
    #mx1_top_1.append(j)
    #mx1_bot_1.append(j)
    #mx2_top_1.append(j)
    #mx2_bot_1.append(j)
    #epsnuc_M_1_1.append(j)
    #epsnuc_M_2.append(j)
    #epsnuc_M_3_1.append(j)
    #epsnuc_M_4_1.append(j)
    #epsnuc_M_5_1.append(j)
    #epsnuc_M_6_1.append(j)
    #epsnuc_M_7_1.append(j)
    #epsnuc_M_8_1.append(j)
    #c_core_mass_1.append(j)
    #o_core_mass_1.append(j)
    tdyn_1.append(j)
    #tkh_1.append(j)
    #tnuc_1.append(j)
    #log_LC_1.append(j)
    #log_LZ_1.append(j)
    #log_Lnuc_1.append(j)
    #log_Lneu_1.append(j)
    #log_Lneu_nuc_1.append(j)
    #log_Lneu_nonnuc_1.append(j)
    #mass_loc_of_max_eps_nuc_1.append(j)
    #log_abs_Lgrav_1.append(j)
    mass_conv_core_1.append(j)
    #cz_bot_mass_1.append(j)
    #cz_top_mass_1.append(j)
    #cz_log_eps_nuc_1.append(j)
    #neutron_rich_core_mass_1.append(j)
    #envelope_mass_1.append(j)
    #center_eps_grav_1.append(j)
    #center_non_nuc_neu_1.append(j)
    #center_eps_nuc_1.append(j)
    #total_mass_h1_1.append(j)
    #total_mass_he4_1.append(j)
    #total_mass_c12_1.append(j)
    #total_mass_o16_1.append(j)
    #total_eps_grav_1.append(j)
    total_nuclear_heating_1.append(j)
    total_non_nuc_neu_cooling_1.append(j)
    #max_eps_nuc_1.append(j)
    #max_eps_nuc_lgT_1.append(j)
    #max_eps_nuc_lgRho_1.append(j)
    #max_eps_nuc_m_1.append(j)
    #max_eps_nuc_xm_1.append(j)
    #max_eps_nuc_lgP_1.append(j)
    #max_eps_nuc_lgR_1.append(j)
    #max_eps_nuc_opacity_1.append(j)
    #max_eps_nuc_cp_1.append(j)
    max_eps_nuc_t_heat_1.append(j)
    #max_eps_nuc_csound_1.append(j)
    log_center_Rho_2.append(j)
    log_center_T_2.append(j)
    log_center_P_2.append(j)
    star_age_2.append(j)
    star_mass_2.append(j)
    model_number_2.append(j)
    #log_L_2.append(j)
    #log_R_2.append(j)
    #log_Teff_2.append(j)
    #log_g_2.append(j)
    #center_h1_2.append(j)
    #center_he4_2.append(j)
    #center_c12_2.append(j)
    #center_c13_2.append(j)
    #center_n14_2.append(j)
    #center_o16_2.append(j)
    #center_o18_2.append(j)
    #center_ne20_2.append(j)
    #center_ne22_2.append(j)
    #center_ne23_2.append(j)
    #center_na23_2.append(j)
    #center_mg24_2.append(j)
    #center_si28_2.append(j)
    #center_ni56_2.append(j)
    #center_gamma_2.append(j)
    center_ye_2.append(j)
    #surface_h1_2.append(j)
    #surface_he4_2.append(j)
    #surface_c12_2.append(j)
    #surface_o16_2.append(j)
    #log_average_h1_2.append(j)
    #log_average_he4_2.append(j)
    #log_average_c12_2.append(j)
    #log_average_o16_2.append(j)
    #log_average_ne20_2.append(j)
    #burn_c_2.append(j)
    #burn_o_2.append(j)
    #c12_c12_2.append(j)
    #o16_o16_2.append(j)
    #burn_ne_2.append(j)
    #burn_mg_2.append(j)
    #conv_mx1_top_2.append(j)
    #conv_mx1_bot_2.append(j)
    #conv_mx2_top_2.append(j)
    #conv_mx2_bot_2.append(j)
    #mx1_top_2.append(j)
    #mx1_bot_2.append(j)
    #mx2_top_2.append(j)
    #mx2_bot_2.append(j)
    #epsnuc_M_1_2.append(j)
    #epsnuc_M_2.append(j)
    #epsnuc_M_3_2.append(j)
    #epsnuc_M_4_2.append(j)
    #epsnuc_M_5_2.append(j)
    #epsnuc_M_6_2.append(j)
    #epsnuc_M_7_2.append(j)
    #epsnuc_M_8_2.append(j)
    #c_core_mass_2.append(j)
    #o_core_mass_2.append(j)
    tdyn_2.append(j)
    #tkh_2.append(j)
    #tnuc_2.append(j)
    #log_LC_2.append(j)
    #log_LZ_2.append(j)
    #log_Lnuc_2.append(j)
    #log_Lneu_2.append(j)
    #log_Lneu_nuc_2.append(j)
    #log_Lneu_nonnuc_2.append(j)
    #mass_loc_of_max_eps_nuc_2.append(j)
    #log_abs_Lgrav_2.append(j)
    mass_conv_core_2.append(j)
    #cz_bot_mass_2.append(j)
    #cz_top_mass_2.append(j)
    #cz_log_eps_nuc_2.append(j)
    #neutron_rich_core_mass_2.append(j)
    #envelope_mass_2.append(j)
    #center_eps_grav_2.append(j)
    #center_non_nuc_neu_2.append(j)
    #center_eps_nuc_2.append(j)
    #total_mass_h1_2.append(j)
    #total_mass_he4_2.append(j)
    #total_mass_c12_2.append(j)
    #total_mass_o16_2.append(j)
    #total_eps_grav_2.append(j)
    total_nuclear_heating_2.append(j)
    total_non_nuc_neu_cooling_2.append(j)
    #max_eps_nuc_2.append(j)
    #max_eps_nuc_lgT_2.append(j)
    #max_eps_nuc_lgRho_2.append(j)
    #max_eps_nuc_m_2.append(j)
    #max_eps_nuc_xm_2.append(j)
    #max_eps_nuc_lgP_2.append(j)
    #max_eps_nuc_lgR_2.append(j)
    #max_eps_nuc_opacity_2.append(j)
    #max_eps_nuc_cp_2.append(j)
    max_eps_nuc_t_heat_2.append(j)
    #max_eps_nuc_csound_2.append(j)
    log_center_Rho.append(j)
    log_center_T.append(j)
    log_center_P.append(j)
    star_age.append(j)
    star_mass.append(j)
    #model_number.append(j)
    #log_L.append(j)
    #log_R.append(j)
    #log_Teff.append(j)
    #log_g.append(j)
    #center_h1.append(j)
    #center_he4.append(j)
    #center_c12.append(j)
    #center_c13.append(j)
    #center_n14.append(j)
    #center_o16.append(j)
    #center_o18.append(j)
    #center_ne20.append(j)
    #center_ne22.append(j)
    #center_ne23.append(j)
    #center_na23.append(j)
    #center_mg24.append(j)
    #center_si28.append(j)
    #center_gamma.append(j)
    center_ye.append(j)
    #surface_h1.append(j)
    #surface_he4.append(j)
    #surface_c12.append(j)
    #surface_o16.append(j)
    #log_average_h1.append(j)
    #log_average_he4.append(j)
    #log_average_c12.append(j)
    #log_average_o16.append(j)
    #log_average_ne20.append(j)
    #burn_c.append(j)
    #burn_o.append(j)
    #c12_c12.append(j)
    #o16_o16.append(j)
    #burn_ne.append(j)
    #burn_mg.append(j)
    #conv_mx1_top.append(j)
    #conv_mx1_bot.append(j)
    #conv_mx2_top.append(j)
    #conv_mx2_bot.append(j)
    #mx1_top.append(j)
    #mx1_bot.append(j)
    #mx2_top.append(j)
    #mx2_bot.append(j)
    #epsnuc_M_1.append(j)
    #epsnuc_M_2.append(j)
    #epsnuc_M_3.append(j)
    #epsnuc_M_4.append(j)
    #epsnuc_M_5.append(j)
    #epsnuc_M_6.append(j)
    #epsnuc_M_7.append(j)
    #epsnuc_M_8.append(j)
    #c_core_mass.append(j)
    #o_core_mass.append(j)
    tdyn.append(j)
    #tkh.append(j)
    #tnuc.append(j)
    #log_LC.append(j)
    #log_LZ.append(j)
    #log_Lnuc.append(j)
    #log_Lneu.append(j)
    #log_Lneu_nuc.append(j)
    #log_Lneu_nonnuc.append(j)
    #mass_loc_of_max_eps_nuc.append(j)
    #log_abs_Lgrav.append(j)
    mass_conv_core.append(j)
    #cz_bot_mass.append(j)
    #cz_top_mass.append(j)
    #cz_log_eps_nuc.append(j)
    #neutron_rich_core_mass.append(j)
    #envelope_mass.append(j)
    #center_eps_grav.append(j)
    #center_non_nuc_neu.append(j)
    #center_eps_nuc.append(j)
    #total_mass_h1.append(j)
    #total_mass_he4.append(j)
    #total_mass_c12.append(j)
    #total_mass_o16.append(j)
    #total_eps_grav.append(j)
    total_nuclear_heating.append(j)
    total_non_nuc_neu_cooling.append(j)
    #max_eps_nuc.append(j)
    #max_eps_nuc_lgT.append(j)
    #max_eps_nuc_lgRho.append(j)
    #max_eps_nuc_m.append(j)
    #max_eps_nuc_xm.append(j)
    #max_eps_nuc_lgP.append(j)
    #max_eps_nuc_lgR.append(j)
    #max_eps_nuc_opacity.append(j)
    #max_eps_nuc_cp.append(j)
    max_eps_nuc_t_heat.append(j)
    #max_eps_nuc_csound.append(j)
    eta.append(j)
    prof_index_1.append(j)
    prof_index_2.append(j)
    idx1.append(j)
    idx1a.append(j)
    idx1b.append(j)
    idx1c.append(j)
    index1.append(j)
    idx2.append(j)
    idx2a.append(j)
    idx2b.append(j)
    idx2c.append(j)
    index2.append(j)
    mod_number_1.append(j)
    st_age_1.append(j)
    mod_number_2.append(j)
    st_age_2.append(j)
    st_age.append(j)
    prof_1.append(j)
    prof_2.append(j)
    prof.append(j)
    #q.append(j)
    #logT.append(j)
    #gradr_sub_grada.append(j)
    #conv.append(j)
    #conv_core.append(j)
    #ccsize.append(j)
    #cc.append(j)
    #lenarr.append(j)
    #maxn.append(j)
    #maxlength.append(j)
    #age.append(j)
    age_ign.append(j)
    c_ign.append(j)
    eta_ign.append(j)
    mass_ign.append(j)
    max_age.append(j)
    #logTiso.append(j)
    title.append(j)
    Z.append(j)
    Mass.append(j)
    Rate.append(j)
    Sim_age.append(j)
    Sim_logTc.append(j)
    Sim_logRhoc.append(j)
    Sim_mass.append(j)
    Elapsed_time.append(j)
    Final_age.append(j)
    Final_logRhoc.append(j)
    Final_mass.append(j)
    Final_mass_cc.append(j)
    Frac_final_mass_cc.append(j)
    Final_eta.append(j)
    Diff_eta.append(j)
    Incr_eta_sun.append(j)
    th.append(j)
    tdyn.append(j)
    #logTiso.append(j)
    val.append(j)
    t.append(j)
    t_tex.append(j)



    log_center_Rho_1[j] = m1[j].get('log_center_Rho')
    log_center_T_1[j] = m1[j].get('log_center_T')
    log_center_P_1[j] = m1[j].get('log_center_P')
    star_age_1[j] = m1[j].get('star_age')
    star_mass_1[j] = m1[j].get('star_mass')
    model_number_1[j] = m1[j].get('model_number')
    #log_L_1[j] = m1[j].get('log_L')
    #log_R_1[j] = m1[j].get('log_R')
    #log_Teff_1[j] = m1[j].get('log_Teff')
    #log_g_1[j] = m1[j].get('log_g')

    #center_h1_1[j] = m1[j].get('center_h1')
    #center_he4_1[j] = m1[j].get('center_he4')
    #center_c12_1[j] = m1[j].get('center_c12')
    #center_c13_1[j] = m1[j].get('center_c13')
    #center_n14_1[j] = m1[j].get('center_n14')
    #center_o16_1[j] = m1[j].get('center_o16')
    #center_o18_1[j] = m1[j].get('center_o18')
    #center_ne20_1[j] = m1[j].get('center_ne20')
    #center_ne22_1[j] = m1[j].get('center_ne22')
    #center_ne23_1[j] = m1[j].get('center_ne23')
    #center_na23_1[j] = m1[j].get('center_na23')
    #center_mg24_1[j] = m1[j].get('center_mg24')
    #center_si28_1[j] = m1[j].get('center_si28')
    #center_ni56_1[j] = m1[j].get('center_ni56')
    #center_gamma_1[j] = m1[j].get('center_gamma')
    center_ye_1[j] = m1[j].get('center_ye')

    #surface_h1_1[j] = m1[j].get('surface_h1')
    #surface_he4_1[j] = m1[j].get('surface_he4')
    #surface_c12_1[j] = m1[j].get('surface_c12')
    #surface_o16_1[j] = m1[j].get('surface_o16')
    #log_average_h1_1[j] = m1[j].get('log_average_h1')
    #log_average_he4_1[j] = m1[j].get('log_average_he4')
    #log_average_c12_1[j] = m1[j].get('log_average_c12')
    #log_average_o16_1[j] = m1[j].get('log_average_o16')
    #log_average_ne20_1[j] = m1[j].get('log_average_ne20')

    #burn_c_1[j] = m1[j].get('burn_c')
    #burn_o_1[j] = m1[j].get('burn_o')
    #c12_c12_1[j] = m1[j].get('c12_c12')
    #o16_o16_1[j] = m1[j].get('o16_o16')
    #burn_ne_1[j] = m1[j].get('burn_ne')
    #burn_mg_1[j] = m1[j].get('burn_mg')

    #conv_mx1_top_1[j] = m1[j].get('conv_mx1_top')
    #conv_mx1_bot_1[j] = m1[j].get('conv_mx1_bot')
    #conv_mx2_top_1[j] = m1[j].get('conv_mx2_top')
    #conv_mx2_bot_1[j] = m1[j].get('conv_mx2_bot')
    #mx1_top_1[j] = m1[j].get('mx1_top')
    #mx1_bot_1[j] = m1[j].get('mx1_bot')
    #mx2_top_1[j] = m1[j].get('mx2_top')
    #mx2_bot_1[j] = m1[j].get('mx2_bot')

    #epsnuc_M_1_1[j] = m1[j].get('epsnuc_M_1')
    #epsnuc_M_2[j]_1[j] = m1[j].get('epsnuc_M_2[j]')
    #epsnuc_M_3_1[j] = m1[j].get('epsnuc_M_3')
    #epsnuc_M_4_1[j] = m1[j].get('epsnuc_M_4')
    #epsnuc_M_5_1[j] = m1[j].get('epsnuc_M_5')
    #epsnuc_M_6_1[j] = m1[j].get('epsnuc_M_6')
    #epsnuc_M_7_1[j] = m1[j].get('epsnuc_M_7')
    #epsnuc_M_8_1[j] = m1[j].get('epsnuc_M_8')

    #c_core_mass_1[j] = m1[j].get('c_core_mass')
    #o_core_mass_1[j] = m1[j].get('o_core_mass')

    tdyn_1[j] = m1[j].get('dynamic_timescale')
    #tkh_1[j] = m1[j].get('kh_timescale')
    #tnuc_1[j] = m1[j].get('nuc_timescale')

    #log_LC_1[j] = m1[j].get('log_LC')
    #log_LZ_1[j] = m1[j].get('log_LZ')
    #log_Lnuc_1[j] = m1[j].get('log_Lnuc')
    #log_Lneu_1[j] = m1[j].get('log_Lneu')
    #log_Lneu_nuc_1[j] = m1[j].get('log_Lneu_nuc')
    #log_Lneu_nonnuc_1[j] = m1[j].get('log_Lneu_nonnuc')
    #mass_loc_of_max_eps_nuc_1[j] = m1[j].get('mass_loc_of_max_eps_nuc')
    #log_abs_Lgrav_1[j] = m1[j].get('log_abs_Lgrav')

    mass_conv_core_1[j] = m1[j].get('mass_conv_core')
    #cz_bot_mass_1[j] = m1[j].get('cz_bot_mass')
    #cz_top_mass_1[j] = m1[j].get('cz_top_mass')
    #cz_log_eps_nuc_1[j] = m1[j].get('cz_log_eps_nuc')
    #neutron_rich_core_mass_1[j] = m1[j].get('neutron_rich_core_mass')
    #envelope_mass_1[j] = m1[j].get('envelope_mass')

    #center_eps_grav_1[j] = m1[j].get('center_eps_grav')
    #center_non_nuc_neu_1[j] = m1[j].get('center_non_nuc_neu')
    #center_eps_nuc_1[j] = m1[j].get('center_eps_nuc')

    #total_mass_h1_1[j] = m1[j].get('total_mass_h1')
    #total_mass_he4_1[j] = m1[j].get('total_mass_he4')
    #total_mass_c12_1[j] = m1[j].get('total_mass_c12')
    #total_mass_o16_1[j] = m1[j].get('total_mass_o16')

    #total_eps_grav_1[j] = m1[j].get('total_eps_grav')
    total_nuclear_heating_1[j] = m1[j].get('total_nuclear_heating')
    total_non_nuc_neu_cooling_1[j] = m1[j].get('total_non_nuc_neu_cooling')

    #max_eps_nuc_1[j] = m1[j].get('max_eps_nuc')
    #max_eps_nuc_lgT_1[j] = m1[j].get('max_eps_nuc_lgT')
    #max_eps_nuc_lgRho_1[j] = m1[j].get('max_eps_nuc_lgRho')
    #max_eps_nuc_m_1[j] = m1[j].get('max_eps_nuc_m')
    #max_eps_nuc_xm_1[j] = m1[j].get('max_eps_nuc_xm')
    #max_eps_nuc_lgP_1[j] = m1[j].get('max_eps_nuc_lgP')
    #max_eps_nuc_lgR_1[j] = m1[j].get('max_eps_nuc_lgR')
    #max_eps_nuc_opacity_1[j] = m1[j].get('max_eps_nuc_opacity')
    #max_eps_nuc_cp_1[j] = m1[j].get('max_eps_nuc_cp')
    max_eps_nuc_t_heat_1[j] = m1[j].get('max_eps_nuc_t_heat')
    #max_eps_nuc_csound_1[j] = m1[j].get('max_eps_nuc_csound')



    log_center_Rho_2[j] = m2[j].get('log_center_Rho')
    log_center_T_2[j] = m2[j].get('log_center_T')
    log_center_P_2[j] = m2[j].get('log_center_P')
    star_age_2[j] = m2[j].get('star_age')
    star_mass_2[j] = m2[j].get('star_mass')
    model_number_2[j] = m2[j].get('model_number')
    #log_L_2[j] = m2[j].get('log_L')
    #log_R_2[j] = m2[j].get('log_R')
    #log_Teff_2[j] = m2[j].get('log_Teff')
    #log_g_2[j] = m2[j].get('log_g')

    #center_h1_2[j] = m2[j].get('center_h1')
    #center_he4_2[j] = m2[j].get('center_he4')
    #center_c12_2[j] = m2[j].get('center_c12')
    #center_c13_2[j] = m2[j].get('center_c13')
    #center_n14_2[j] = m2[j].get('center_n14')
    #center_o16_2[j] = m2[j].get('center_o16')
    #center_o18_2[j] = m2[j].get('center_o18')
    #center_ne20_2[j] = m2[j].get('center_ne20')
    #center_ne22_2[j] = m2[j].get('center_ne22')
    #center_ne23_2[j] = m2[j].get('center_ne23')
    #center_na23_2[j] = m2[j].get('center_na23')
    #center_mg24_2[j] = m2[j].get('center_mg24')
    #center_si28_2[j] = m2[j].get('center_si28')
    #center_ni56_2[j] = m2[j].get('center_ni56')
    #center_gamma_2[j] = m2[j].get('center_gamma')
    center_ye_2[j] = m2[j].get('center_ye')

    #surface_h1_2[j] = m2[j].get('surface_h1')
    #surface_he4_2[j] = m2[j].get('surface_he4')
    #surface_c12_2[j] = m2[j].get('surface_c12')
    #surface_o16_2[j] = m2[j].get('surface_o16')
    #log_average_h1_2[j] = m2[j].get('log_average_h1')
    #log_average_he4_2[j] = m2[j].get('log_average_he4')
    #log_average_c12_2[j] = m2[j].get('log_average_c12')
    #log_average_o16_2[j] = m2[j].get('log_average_o16')
    #log_average_ne20_2[j] = m2[j].get('log_average_ne20')

    #burn_c_2[j] = m2[j].get('burn_c')
    #burn_o_2[j] = m2[j].get('burn_o')
    #c12_c12_2[j] = m2[j].get('c12_c12')
    #o16_o16_2[j] = m2[j].get('o16_o16')
    #burn_ne_2[j] = m2[j].get('burn_ne')
    #burn_mg_2[j] = m2[j].get('burn_mg')

    #conv_mx1_top_2[j] = m2[j].get('conv_mx1_top')
    #conv_mx1_bot_2[j] = m2[j].get('conv_mx1_bot')
    #conv_mx2_top_2[j] = m2[j].get('conv_mx2_top')
    #conv_mx2_bot_2[j] = m2[j].get('conv_mx2_bot')
    #mx1_top_2[j] = m2[j].get('mx1_top')
    #mx1_bot_2[j] = m2[j].get('mx1_bot')
    #mx2_top_2[j] = m2[j].get('mx2_top')
    #mx2_bot_2[j] = m2[j].get('mx2_bot')

    #epsnuc_M_1_2[j] = m2[j].get('epsnuc_M_1')
    #epsnuc_M_2[j]_2[j] = m2[j].get('epsnuc_M_2[j]')
    #epsnuc_M_3_2[j] = m2[j].get('epsnuc_M_3')
    #epsnuc_M_4_2[j] = m2[j].get('epsnuc_M_4')
    #epsnuc_M_5_2[j] = m2[j].get('epsnuc_M_5')
    #epsnuc_M_6_2[j] = m2[j].get('epsnuc_M_6')
    #epsnuc_M_7_2[j] = m2[j].get('epsnuc_M_7')
    #epsnuc_M_8_2[j] = m2[j].get('epsnuc_M_8')

    #c_core_mass_2[j] = m2[j].get('c_core_mass')
    #o_core_mass_2[j] = m2[j].get('o_core_mass')

    tdyn_2[j] = m2[j].get('dynamic_timescale')
    #tkh_2[j] = m2[j].get('kh_timescale')
    #tnuc_2[j] = m2[j].get('nuc_timescale')

    #log_LC_2[j] = m2[j].get('log_LC')
    #log_LZ_2[j] = m2[j].get('log_LZ')
    #log_Lnuc_2[j] = m2[j].get('log_Lnuc')
    #log_Lneu_2[j] = m2[j].get('log_Lneu')
    #log_Lneu_nuc_2[j] = m2[j].get('log_Lneu_nuc')
    #log_Lneu_nonnuc_2[j] = m2[j].get('log_Lneu_nonnuc')
    #mass_loc_of_max_eps_nuc_2[j] = m2[j].get('mass_loc_of_max_eps_nuc')
    #log_abs_Lgrav_2[j] = m2[j].get('log_abs_Lgrav')

    mass_conv_core_2[j] = m2[j].get('mass_conv_core')
    #cz_bot_mass_2[j] = m2[j].get('cz_bot_mass')
    #cz_top_mass_2[j] = m2[j].get('cz_top_mass')
    #cz_log_eps_nuc_2[j] = m2[j].get('cz_log_eps_nuc')
    #neutron_rich_core_mass_2[j] = m2[j].get('neutron_rich_core_mass')
    #envelope_mass_2[j] = m2[j].get('envelope_mass')

    #center_eps_grav_2[j] = m2[j].get('center_eps_grav')
    #center_non_nuc_neu_2[j] = m2[j].get('center_non_nuc_neu')
    #center_eps_nuc_2[j] = m2[j].get('center_eps_nuc')

    #total_mass_h1_2[j] = m2[j].get('total_mass_h1')
    #total_mass_he4_2[j] = m2[j].get('total_mass_he4')
    #total_mass_c12_2[j] = m2[j].get('total_mass_c12')
    #total_mass_o16_2[j] = m2[j].get('total_mass_o16')

    #total_eps_grav_2[j] = m2[j].get('total_eps_grav')
    total_nuclear_heating_2[j] = m2[j].get('total_nuclear_heating')
    total_non_nuc_neu_cooling_2[j] = m2[j].get('total_non_nuc_neu_cooling')

    #max_eps_nuc_2[j] = m2[j].get('max_eps_nuc')
    #max_eps_nuc_lgT_2[j] = m2[j].get('max_eps_nuc_lgT')
    #max_eps_nuc_lgRho_2[j] = m2[j].get('max_eps_nuc_lgRho')
    #max_eps_nuc_m_2[j] = m2[j].get('max_eps_nuc_m')
    #max_eps_nuc_xm_2[j] = m2[j].get('max_eps_nuc_xm')
    #max_eps_nuc_lgP_2[j] = m2[j].get('max_eps_nuc_lgP')
    #max_eps_nuc_lgR_2[j] = m2[j].get('max_eps_nuc_lgR')
    #max_eps_nuc_opacity_2[j] = m2[j].get('max_eps_nuc_opacity')
    #max_eps_nuc_cp_2[j] = m2[j].get('max_eps_nuc_cp')
    max_eps_nuc_t_heat_2[j] = m2[j].get('max_eps_nuc_t_heat')
    #max_eps_nuc_csound_2[j] = m2[j].get('max_eps_nuc_csound')



    log_center_Rho[j] = np.concatenate((log_center_Rho_1[j],log_center_Rho_2[j]),axis=0)
    log_center_T[j] = np.concatenate((log_center_T_1[j],log_center_T_2[j]),axis=0)
    log_center_P[j] = np.concatenate((log_center_P_1[j],log_center_P_2[j]),axis=0)
    star_age[j] = np.concatenate((star_age_1[j],star_age_2[j]),axis=0)
    star_mass[j] = np.concatenate((star_mass_1[j],star_mass_2[j]),axis=0)
    #model_number[j] = np.concatenate((model_number_1[j],model_number_2[j]),axis=0)
    #log_L[j] = np.concatenate((log_L_1[j],log_L_2[j]),axis=0)
    #log_R[j] = np.concatenate((log_R_1[j],log_R_2[j]),axis=0)
    #log_Teff[j] = np.concatenate((log_Teff_1[j],log_Teff_2[j]),axis=0)
    #log_g[j] = np.concatenate((log_g_1[j],log_g_2[j]),axis=0)

    #center_h1[j] = np.concatenate((center_h1_1[j],center_h1_2[j]),axis=0)
    #center_he4[j] = np.concatenate((center_he4_1[j],center_he4_2[j]),axis=0)
    #center_c12[j] = np.concatenate((center_c12_1[j],center_c12_2[j]),axis=0)
    #center_c13[j] = np.concatenate((center_c13_1[j],center_c13_2[j]),axis=0)
    #center_n14[j] = np.concatenate((center_n14_1[j],center_n14_2[j]),axis=0)
    #center_o16[j] = np.concatenate((center_o16_1[j],center_o16_2[j]),axis=0)
    #center_o18[j] = np.concatenate((center_o18_1[j],center_o18_2[j]),axis=0)
    #center_ne20[j] = np.concatenate((center_ne20_1[j],center_ne20_2[j]),axis=0)
    #center_ne22[j] = np.concatenate((center_ne22_1[j],center_ne22_2[j]),axis=0)
    #center_ne23[j] = np.concatenate((center_ne23_1[j],center_ne23_2[j]),axis=0)
    #center_na23[j] = np.concatenate((center_na23_1[j],center_na23_2[j]),axis=0)
    #center_mg24[j] = np.concatenate((center_mg24_1[j],center_mg24_2[j]),axis=0)
    #center_si28[j] = np.concatenate((center_si28_1[j],center_si28_2[j]),axis=0)
    #center_gamma[j] = np.concatenate((center_gamma_1[j],center_gamma_2[j]),axis=0)
    center_ye[j] = np.concatenate((center_ye_1[j],center_ye_2[j]),axis=0)

    #surface_h1[j] = np.concatenate((surface_h1_1[j],surface_h1_2[j]),axis=0)
    #surface_he4[j] = np.concatenate((surface_he4_1[j],surface_he4_2[j]),axis=0)
    #surface_c12[j] = np.concatenate((surface_c12_1[j],surface_c12_2[j]),axis=0)
    #surface_o16[j] = np.concatenate((surface_o16_1[j],surface_o16_2[j]),axis=0)
    #log_average_h1[j] = np.concatenate((log_average_h1_1[j],log_average_h1_2[j]),axis=0)
    #log_average_he4[j] = np.concatenate((log_average_he4_1[j],log_average_he4_2[j]),axis=0)
    #log_average_c12[j] = np.concatenate((log_average_c12_1[j],log_average_c12_2[j]),axis=0)
    #log_average_o16[j] = np.concatenate((log_average_o16_1[j],log_average_o16_2[j]),axis=0)
    #log_average_ne20[j] = np.concatenate((log_average_ne20_1[j],log_average_ne20_2[j]),axis=0)

    #burn_c[j] = np.concatenate((burn_c_1[j],burn_c_2[j]),axis=0)
    #burn_o[j] = np.concatenate((burn_o_1[j],burn_o_2[j]),axis=0)
    #c12_c12[j] = np.concatenate((c12_c12_1[j],c12_c12_2[j]),axis=0)
    #o16_o16[j] = np.concatenate((o16_o16_1[j],o16_o16_2[j]),axis=0)
    #burn_ne[j] = np.concatenate((burn_ne_1[j],burn_ne_2[j]),axis=0)
    #burn_mg[j] = np.concatenate((burn_mg_1[j],burn_mg_2[j]),axis=0)

    #conv_mx1_top[j] = np.concatenate((conv_mx1_top_1[j],conv_mx1_top_2[j]),axis=0)
    #conv_mx1_bot[j] = np.concatenate((conv_mx1_bot_1[j],conv_mx1_bot_2[j]),axis=0)
    #conv_mx2_top[j] = np.concatenate((conv_mx2_top_1[j],conv_mx2_top_2[j]),axis=0)
    #conv_mx2_bot[j] = np.concatenate((conv_mx2_bot_1[j],conv_mx2_bot_2[j]),axis=0)
    #mx1_top[j] = np.concatenate((mx1_top_1[j],mx1_top_2[j]),axis=0)
    #mx1_bot[j] = np.concatenate((mx1_bot_1[j],mx1_bot_2[j]),axis=0)
    #mx2_top[j] = np.concatenate((mx2_top_1[j],mx2_top_2[j]),axis=0)
    #mx2_bot[j] = np.concatenate((mx2_bot_1[j],mx2_bot_2[j]),axis=0)

    #epsnuc_M_1[j][j] = np.concatenate((epsnuc_M_1[j]_1[j],epsnuc_M_1[j]_2[j]),axis=0)
    #epsnuc_M_2[j][j] = np.concatenate((epsnuc_M_2[j]_1[j],epsnuc_M_2[j]_2[j]),axis=0)
    #epsnuc_M_3[j] = np.concatenate((epsnuc_M_3_1[j],epsnuc_M_3_2[j]),axis=0)
    #epsnuc_M_4[j] = np.concatenate((epsnuc_M_4_1[j],epsnuc_M_4_2[j]),axis=0)
    #epsnuc_M_5[j] = np.concatenate((epsnuc_M_5_1[j],epsnuc_M_5_2[j]),axis=0)
    #epsnuc_M_6[j] = np.concatenate((epsnuc_M_6_1[j],epsnuc_M_6_2[j]),axis=0)
    #epsnuc_M_7[j] = np.concatenate((epsnuc_M_7_1[j],epsnuc_M_7_2[j]),axis=0)
    #epsnuc_M_8[j] = np.concatenate((epsnuc_M_8_1[j],epsnuc_M_8_2[j]),axis=0)

    #c_core_mass[j] = np.concatenate((c_core_mass_1[j],c_core_mass_2[j]),axis=0)
    #o_core_mass[j] = np.concatenate((o_core_mass_1[j],o_core_mass_2[j]),axis=0)

    tdyn[j] = np.concatenate((tdyn_1[j],tdyn_2[j]),axis=0)
    #tkh[j] = np.concatenate((tkh_1[j],tkh_2[j]),axis=0)
    #tnuc[j] = np.concatenate((tnuc_1[j],tnuc_2[j]),axis=0)

    #log_LC[j] = np.concatenate((log_LC_1[j],log_LC_2[j]),axis=0)
    #log_LZ[j] = np.concatenate((log_LZ_1[j],log_LZ_2[j]),axis=0)
    #log_Lnuc[j] = np.concatenate((log_Lnuc_1[j],log_Lnuc_2[j]),axis=0)
    #log_Lneu[j] = np.concatenate((log_Lneu_1[j],log_Lneu_2[j]),axis=0)
    #log_Lneu_nuc[j] = np.concatenate((log_Lneu_nuc_1[j],log_Lneu_nuc_2[j]),axis=0)
    #log_Lneu_nonnuc[j] = np.concatenate((log_Lneu_nonnuc_1[j],log_Lneu_nonnuc_2[j]),axis=0)
    #mass_loc_of_max_eps_nuc[j] = np.concatenate((mass_loc_of_max_eps_nuc_1[j],mass_loc_of_max_eps_nuc_2[j]),axis=0)
    #log_abs_Lgrav[j] = np.concatenate((log_abs_Lgrav_1[j],log_abs_Lgrav_2[j]),axis=0)

    mass_conv_core[j] = np.concatenate((mass_conv_core_1[j],mass_conv_core_2[j]),axis=0)
    #cz_bot_mass[j] = np.concatenate((cz_bot_mass_1[j],cz_bot_mass_2[j]),axis=0)
    #cz_top_mass[j] = np.concatenate((cz_top_mass_1[j],cz_top_mass_2[j]),axis=0)
    #cz_log_eps_nuc[j] = np.concatenate((cz_log_eps_nuc_1[j],cz_log_eps_nuc_2[j]),axis=0)
    #neutron_rich_core_mass[j] = np.concatenate((neutron_rich_core_mass_1[j],neutron_rich_core_mass_2[j]),axis=0)
    #envelope_mass[j] = np.concatenate((envelope_mass_1[j],envelope_mass_2[j]),axis=0)

    #center_eps_grav[j] = np.concatenate((center_eps_grav_1[j],center_eps_grav_2[j]),axis=0)
    #center_non_nuc_neu[j] = np.concatenate((center_non_nuc_neu_1[j],center_non_nuc_neu_2[j]),axis=0)
    #center_eps_nuc[j] = np.concatenate((center_eps_nuc_1[j],center_eps_nuc_2[j]),axis=0)

    #total_mass_h1[j] = np.concatenate((total_mass_h1_1[j],total_mass_h1_2[j]),axis=0)
    #total_mass_he4[j] = np.concatenate((total_mass_he4_1[j],total_mass_he4_2[j]),axis=0)
    #total_mass_c12[j] = np.concatenate((total_mass_c12_1[j],total_mass_c12_2[j]),axis=0)
    #total_mass_o16[j] = np.concatenate((total_mass_o16_1[j],total_mass_o16_2[j]),axis=0)

    #total_eps_grav[j] = np.concatenate((total_eps_grav_1[j],total_eps_grav_2[j]),axis=0)
    total_nuclear_heating[j] = np.concatenate((total_nuclear_heating_1[j],total_nuclear_heating_2[j]),axis=0)
    total_non_nuc_neu_cooling[j] = np.concatenate((total_non_nuc_neu_cooling_1[j],total_non_nuc_neu_cooling_2[j]),axis=0)

    #max_eps_nuc[j] = np.concatenate((max_eps_nuc_1[j],max_eps_nuc_2[j]),axis=0)
    #max_eps_nuc_lgT[j] = np.concatenate((max_eps_nuc_lgT_1[j],max_eps_nuc_lgT_2[j]),axis=0)
    #max_eps_nuc_lgRho[j] = np.concatenate((max_eps_nuc_lgRho_1[j],max_eps_nuc_lgRho_2[j]),axis=0)
    #max_eps_nuc_m[j] = np.concatenate((max_eps_nuc_m_1[j],max_eps_nuc_m_2[j]),axis=0)
    #max_eps_nuc_xm[j] = np.concatenate((max_eps_nuc_xm_1[j],max_eps_nuc_xm_2[j]),axis=0)
    #max_eps_nuc_lgP[j] = np.concatenate((max_eps_nuc_lgP_1[j],max_eps_nuc_lgP_2[j]),axis=0)
    #max_eps_nuc_lgR[j] = np.concatenate((max_eps_nuc_lgR_1[j],max_eps_nuc_lgR_2[j]),axis=0)
    #max_eps_nuc_opacity[j] = np.concatenate((max_eps_nuc_opacity_1[j],max_eps_nuc_opacity_2[j]),axis=0)
    #max_eps_nuc_cp[j] = np.concatenate((max_eps_nuc_cp_1[j],max_eps_nuc_cp_2[j]),axis=0)
    max_eps_nuc_t_heat[j] = np.concatenate((max_eps_nuc_t_heat_1[j],max_eps_nuc_t_heat_2[j]),axis=0)
    #max_eps_nuc_csound[j] = np.concatenate((max_eps_nuc_csound_1[j],max_eps_nuc_csound_2[j]),axis=0)




    eta[j] = 1.-2.*center_ye[j]
    eta_sun = 1.4E-3




    # Several important definitions

    c_ign[j] = np.where(mass_conv_core[j]>0)[0][0] # Index where carbon ignition occurs in the history data
    age_ign[j] = star_age[j][c_ign[j]]
    #c_ign = np.where(star_age>age_ign)[0][0] # Index where carbon ignition occurs in the history data
    eta_ign[j] = eta[j][c_ign[j]]
    mass_ign[j] = star_mass[j][c_ign[j]]
    max_age[j] = np.max(star_age[j])
    #max_age = np.max(st_age)


    # Final table


    # Creates a title which includes some parts of the current work directory. I am removing all the "/" and other characters
    #title[j] = '_'+os.path.dirname(path)[-22:-18]+'_'+os.path.dirname(path)[-12:-8]+'_'+os.path.dirname(path)[-7:] 
    title[j] = str(j)




    Z[j] = float(os.path.dirname(path)[-22:-18])
    Mass[j] = float(os.path.dirname(path)[-8:-4])
    Rate[j] = float(os.path.dirname(path)[-13]+'.'+os.path.dirname(path)[-12:-9])
    Sim_age[j] = "%.3f"%(1.E-6*age_ign[j])
    Sim_logTc[j] = "%.2f"%(log_center_T[j][c_ign[j]]) 
    Sim_logRhoc[j] = "%.2f"%(log_center_Rho[j][c_ign[j]]) 
    Sim_mass[j] = "%.3f"%(mass_ign[j])
    Elapsed_time[j] = "%.2f"%(1.E-3*(max_age[j]-age_ign[j]))
    Final_age[j] = "%.3f"%(1.E-6*max_age[j])
    Final_logRhoc[j] = "%.2f"%(log_center_Rho[j][-1])
    Final_mass[j] = "%.3f"%(star_mass[j][-1])
    Final_mass_cc[j] = "%.3f"%(mass_conv_core[j][-1])
    Final_eta[j] = "%.2f"%(1.E3*eta[j][-1])
    #Diff_eta[j] = "%.2f"%(1.E3*(eta[j][-1]-0.101*0.014*Z[j]))
    Diff_eta[j] = "%.2f"%(1.E3*(eta[j][-1]-eta[j][0]))
    Incr_eta_sun[j] = "%.2f"%((eta[j][-1]-eta_ign[j])/eta_sun)
    th[j] = "%.2f"%(max_eps_nuc_t_heat[j][-1])
    tdyn[j] = "%.2f"%(tdyn[j][-1])



    val[j]=np.array([Z[j],Mass[j],Rate[j],Sim_age[j],Sim_logTc[j],Sim_logRhoc[j],Sim_mass[j],Elapsed_time[j],Final_age[j],Final_logRhoc[j],Final_mass[j],Final_mass_cc[j],Final_eta[j],Diff_eta[j],Incr_eta_sun[j]])#,th[j],tdyn[j]])

    nam=['Z','Mass','Rate','Sim_age','Sim_logTc','Sim_logRhoc','Sim_mass','Elapsed_time','Final_age','Final_logRhoc','Final_mass','Final_mass_cc','1E3*Final_eta','1E3*Diff_eta_nosim','Incr_eta_sun']#,'th','tdyn']

    nam_tex=[r'$\rm{Z}$',r'$\rm{M^{init}}$',r'$\rm{\dot{M}}$',r'$\rm{t^{sim}}$',r'$\rm{\log T_{c}^{\, sim}}$',r'$\rm{\log \rho_{c}^{\, sim}}$',r'$\rm{M^{sim}}$',r'$\rm{\Delta t^{sim}}$',\
    r'$\rm{t^{end}}$',r'$\rm{\log \rho_{c}^{\, end}}$',r'$\rm{M^{end}}$',r'$\rm{M_{cc}^{\, end}}$',r'$\rm{10^{3} \eta^{end}}$',r'$\rm{10^{3} (\eta-\eta_{0}})$',r'$\rm{\Delta \eta^{sim}}$']#,r'$\rm{t_{h}}$',r'$\rm{t_{dyn}}$']

    un=['Zsun','Msun','Msun/yr','Myr','K','g/cm3','Msun','kyr','Myr','g/cm3','Msun','Msun','','','etasun']#,'s','s']

    un_tex=[r'$[\rm{Z_{\odot}}]$',r'$[\rm{M_{\odot}}]$',r'$[\rm{M_{\odot} \, yr^{-1}}]$',r'$[\rm{Myr}]$',r'$[\rm{K}]$',r'$[\rm{g \, cm^{-3}}]$',r'$[\rm{M_{\odot}}]$',r'[$\rm{kyr}]$',\
    r'$[\rm{Myr}]$',r'$[\rm{g \, cm^{-3}}]$',r'$[\rm{M_{\odot}}]$',r'$[\rm{M_{\odot}}]$','','',r'$[\rm{\eta_{\odot}}]$']#,r'$[\rm{s}]$',r'$[\rm{s}]$']


    t[j]=Table(val[j],names=nam)
    t_tex[j]=Table(val[j],names=nam_tex)
    for i in range(0,len(un)):
        t[j][nam[i]].unit=un[i]
        t_tex[j][nam_tex[i]].unit=un_tex[i]

    #t[j].write("thermo_runaway"+title[j]+".txt",format='ascii.fixed_width')    # Write the individual tables for every subdirectory
    #t[j].write(pathw+"thermo_runaway"+title[j]+".txt",format='ascii.fixed_width') 


Table_end=Table(vstack([x for x in t]))     # Join all the tables in a single one
Table_tex_end=Table(vstack([x for x in t_tex]))
Table_end.write("thermo_runaway_"+os.getcwd()[-4:]+".txt",format='ascii.fixed_width') 
Table_tex_end.write("thermo_runaway_"+os.getcwd()[-4:]+"_tex.txt",format='latex')