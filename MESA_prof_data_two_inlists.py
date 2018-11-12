
# coding: utf-8
# This script can be used to build data structures from MESA profiles. It creates a similar structure for the stellar age in case the user
# wanted to employ it for, e.g., plots. This model consists of two inlists, which profiles are saved in different LOGS directories. The 
# command employed is pickle. Done by Héctor Martínez Rodríguez


import numpy as np
import scipy.stats as stats
import astropy.stats as astats
import numpy.random as random
from astropy.table import Table


## Set the history.data



import mesa as ms

import os

# Start an instance of the history.data


pathwork=os.getcwd()  # Current work directory
dir1=pathwork+'/LOGS1'
dir2=pathwork+'/LOGS2'
#dir1='/../'
#dir2='/../'
#dir1 = raw_input ('Please, insert the name of the directory where the first history.data and profiles are:\n')
#dir2 = raw_input ('Please, insert the name of the directory where the second history.data and profiles are:\n')




for path, dirs, files in os.walk(dir1): #os.walk avoids [Errno21] Is a directory
    m1=ms.history_data(dir1,clean_starlog=False)
for path, dirs, files in os.walk(dir2): 
    m2=ms.history_data(dir2,clean_starlog=False)




star_age_1=m1.get('star_age')
model_number_1=m1.get('model_number')

star_age_2=m2.get('star_age')
model_number_2=m2.get('model_number')

star_age = np.concatenate((star_age_1,star_age_2), axis=0)
model_number = np.concatenate((model_number_1,model_number_2), axis=0)



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


logxq=[]
radius=[]
logT=[]
logRho=[]
logP=[]
lum=[]
q=[]
ye=[]

pp=[]
cno=[]
tri_alpha=[]
c12_c12=[]  
o16_o16=[]
burn_c=[]
burn_ne=[]
burn_mg=[]
gradr_sub_grada=[]
eps_nuc=[]
non_nuc_neu=[]
photo=[]

neut=[]
h1=[]
he3=[]
he4=[]
c12=[]
c13=[]
n13=[]
n14=[]
o16=[]
o18=[]
ne20=[]
ne22=[]
ne23=[]
na23=[]
mg24=[]
si28=[]
s32=[]
ar36=[]
ca40=[]
cr48=[]
ni56=[]
fe56=[]

log_g=[]
log_opacity=[]
gam=[]
net_nuclear_energy=[]  # erg/gm/s from nuclear reactions minus all neutrino losses. The value plotted is 
#net_nuclear_energy = sign(val)*log10(max(1,abs(val))), where val = net nuclear energy minus all neutrino losses.
eps_nuc_neu_total=[] # erg/gm/sec as neutrinos from nuclear reactions
c12_o16=[]
eps_nuc_plus_nuc_neu=[]

mg25=[]
mg26=[]

log_rate_na23_wk_ne23=[]
log_rate_ne23_wk_minus_na23=[]
log_rate_ne23_pn_na23=[]
log_rate_mg23_wk_na23=[]
log_rate_n13_wk_c13=[]
log_rate_h1_h1_wk_h2=[]
log_rate_be7_wk_li7=[]
log_rate_al25_wk_mg25=[]
log_rate_na21_wk_ne21=[]
log_rate_h1_he3_wk_he4=[]
log_rate_na22_wk_ne22=[]
log_rate_f18_wk_minus_ne18=[]
log_rate_b8_wk_he4_he4=[]
log_rate_na23_wk_minus_mg23=[]
log_rate_na24_wk_minus_mg24=[]
log_rate_f18_wk_o18=[]
log_rate_ne18_wk_f18=[]
log_rate_al26_wk_mg26=[]
log_rate_ne19_wk_f19=[]
log_rate_h1_wk_neut=[]
log_rate_ne21_wk_minus_na21=[]
log_rate_ne22_wk_minus_na22=[]
log_rate_mg25_wk_minus_al25=[]
log_rate_mg24_wk_na24=[]
log_rate_neut_wk_minus_h1=[]
log_rate_o14_wk_n14=[]
log_rate_o15_wk_n15=[]
log_rate_f17_wk_o17=[]
log_rate_o17_wk_minus_f17=[]
log_rate_mg26_wk_minus_al26=[]
log_rate_o18_wk_minus_f18=[]
log_rate_f19_wk_minus_ne19=[]

log_rate_n13_np_c13=[]
log_rate_mg23_np_na23=[]
log_rate_na21_np_ne21=[]
log_rate_na22_np_ne22=[]
log_rate_na23_np_ne23=[]
log_rate_mg24_np_na24=[]
log_rate_o16_na_c13=[]
log_rate_mg24_ng_mg25=[]
log_rate_na22_ng_na23=[]
log_rate_na23_ng_na24=[]
log_rate_ne21_ng_ne22=[]
log_rate_ne20_ng_ne21=[]
log_rate_ne19_ng_ne20=[]
log_rate_ne22_ng_ne23=[]
log_rate_c12_ng_c13=[]
log_rate_mg23_ng_mg24=[]
log_rate_n13_ng_n14=[]

mg23=[] 
na21=[]
ne21=[]
f18=[]
ne18=[]
ne19=[]
f19=[]
f17=[]
o17=[]
na24=[]

t_c=[]
pressure_scale_height=[]
conv_vel=[]
t_heat=[]


for i in range(len(prof)):
    
    
    logxq.append(i) #Avoid IndexError: list assignment index out of range with empty lists []
    logxq[i]=prof[i].get('logxq')
    
    radius.append(i)
    radius[i]=prof[i].get('radius')
    
    logT.append(i)
    logT[i]=prof[i].get('logT')
    
    logRho.append(i)
    logRho[i]=prof[i].get('logRho')
    
    logP.append(i)
    logP[i]=prof[i].get('logP')
    
    lum.append(i)
    lum[i]=prof[i].get('luminosity')
    
    q.append(i)
    q[i]=prof[i].get('q')
    
    ye.append(i)
    ye[i]=prof[i].get('ye')

      
        
    pp.append(i)
    pp[i]=prof[i].get('pp')
    
    cno.append(i)
    cno[i]=prof[i].get('cno')
    
    tri_alpha.append(i)
    tri_alpha[i]=prof[i].get('tri_alfa')
    
    c12_c12.append(i)
    c12_c12[i]=prof[i].get('c12_c12')   
    
    o16_o16.append(i)
    o16_o16[i]=prof[i].get('o16_o16')
    
    burn_c.append(i)
    burn_c[i]=prof[i].get('burn_c')
    
    burn_ne.append(i)
    burn_ne[i]=prof[i].get('burn_ne')
    
    burn_mg.append(i)
    burn_mg[i]=prof[i].get('burn_mg')
    
    gradr_sub_grada.append(i)
    gradr_sub_grada[i]=prof[i].get('gradr_sub_grada')
    
    eps_nuc.append(i)
    eps_nuc[i]=prof[i].get('eps_nuc')
    
    non_nuc_neu.append(i)
    non_nuc_neu[i]=prof[i].get('non_nuc_neu')

    photo.append(i)
    photo[i]=prof[i].get('photo')

    
    
    neut.append(i)
    neut[i]=prof[i].get('neut');

    h1.append(i)
    h1[i]=prof[i].get('h1');
    
    he3.append(i)
    he3[i]=prof[i].get('he3');
    
    he4.append(i)
    he4[i]=prof[i].get('he4');
    
    c12.append(i)
    c12[i]=prof[i].get('c12');
    
    c13.append(i)
    c13[i]=prof[i].get('c13');

    n13.append(i)
    n13[i]=prof[i].get('n13');
    
    n14.append(i)
    n14[i]=prof[i].get('n14');
    
    o16.append(i)
    o16[i]=prof[i].get('o16');
    
    o18.append(i)
    o18[i]=prof[i].get('o18');
    
    ne20.append(i)
    ne20[i]=prof[i].get('ne20');
    
    ne22.append(i)
    ne22[i]=prof[i].get('ne22');
    
    ne23.append(i)
    ne23[i]=prof[i].get('ne23');
    
    na23.append(i)
    na23[i]=prof[i].get('na23');
    
    mg24.append(i)
    mg24[i]=prof[i].get('mg24');
    
    si28.append(i)
    si28[i]=prof[i].get('si28');
    
    s32.append(i)
    s32[i]=prof[i].get('s32');
    
    ar36.append(i)
    ar36[i]=prof[i].get('ar36');
    
    ca40.append(i)
    ca40[i]=prof[i].get('ca40');
    
    cr48.append(i)
    cr48[i]=prof[i].get('cr48');
    
    ni56.append(i)
    ni56[i]=prof[i].get('ni56');
    
    fe56.append(i)
    fe56[i]=prof[i].get('fe56');


    
    log_g.append(i)
    log_g[i]=prof[i].get('log_g');

    log_opacity.append(i)
    log_opacity[i]=prof[i].get('log_opacity');

    gam.append(i)
    gam[i]=prof[i].get('gam');

    net_nuclear_energy.append(i)
    net_nuclear_energy[i]=prof[i].get('net_nuclear_energy');

    eps_nuc_neu_total.append(i)
    eps_nuc_neu_total[i]=prof[i].get('eps_nuc_neu_total');

    c12_o16.append(i)
    c12_o16[i]=prof[i].get('c12_o16');

    eps_nuc_plus_nuc_neu.append(i)
    eps_nuc_plus_nuc_neu[i]=prof[i].get('eps_nuc_plus_nuc_neu');



    mg25.append(i)
    mg25[i]=prof[i].get('mg25');

    mg26.append(i)
    mg26[i]=prof[i].get('mg26');


    log_rate_na23_wk_ne23.append(i)
    log_rate_na23_wk_ne23[i]=prof[i].get('log_rate_na23_wk_ne23');
    log_rate_ne23_wk_minus_na23.append(i)
    log_rate_ne23_wk_minus_na23[i]=prof[i].get('log_rate_ne23_wk-minus_na23');
    log_rate_ne23_pn_na23.append(i)
    log_rate_ne23_pn_na23[i]=prof[i].get('log_rate_ne23_pn_na23');
    log_rate_mg23_wk_na23.append(i)
    log_rate_mg23_wk_na23[i]=prof[i].get('log_rate_mg23_wk_na23');
    log_rate_n13_wk_c13.append(i)
    log_rate_n13_wk_c13[i]=prof[i].get('log_rate_n13_wk_c13');
    log_rate_h1_h1_wk_h2.append(i)
    log_rate_h1_h1_wk_h2[i]=prof[i].get('log_rate_h1_h1_wk_h2');
    log_rate_be7_wk_li7.append(i)
    log_rate_be7_wk_li7[i]=prof[i].get('log_rate_be7_wk_li7');
    log_rate_al25_wk_mg25.append(i)
    log_rate_al25_wk_mg25[i]=prof[i].get('log_rate_al25_wk_mg25');
    log_rate_na21_wk_ne21.append(i)
    log_rate_na21_wk_ne21[i]=prof[i].get('log_rate_na21_wk_ne21');
    log_rate_h1_he3_wk_he4.append(i)
    log_rate_h1_he3_wk_he4[i]=prof[i].get('log_rate_h1_he3_wk_he4');
    log_rate_na22_wk_ne22.append(i)
    log_rate_na22_wk_ne22[i]=prof[i].get('log_rate_na22_wk_ne22');
    log_rate_f18_wk_minus_ne18.append(i)
    log_rate_f18_wk_minus_ne18[i]=prof[i].get('log_rate_f18_wk-minus_ne18');
    log_rate_b8_wk_he4_he4.append(i)
    log_rate_b8_wk_he4_he4[i]=prof[i].get('log_rate_b8_wk_he4_he4');
    log_rate_na23_wk_minus_mg23.append(i)
    log_rate_na23_wk_minus_mg23[i]=prof[i].get('log_rate_na23_wk-minus_mg23');
    log_rate_na24_wk_minus_mg24.append(i)
    log_rate_na24_wk_minus_mg24[i]=prof[i].get('log_rate_na24_wk-minus_mg24');
    log_rate_f18_wk_o18.append(i)
    log_rate_f18_wk_o18[i]=prof[i].get('log_rate_f18_wk_o18');
    log_rate_ne18_wk_f18.append(i)
    log_rate_ne18_wk_f18[i]=prof[i].get('log_rate_ne18_wk_f18');
    log_rate_al26_wk_mg26.append(i)
    log_rate_al26_wk_mg26[i]=prof[i].get('log_rate_al26_wk_mg26');
    log_rate_ne19_wk_f19.append(i)
    log_rate_ne19_wk_f19[i]=prof[i].get('log_rate_ne19_wk_f19');
    log_rate_h1_wk_neut.append(i)
    log_rate_h1_wk_neut[i]=prof[i].get('log_rate_h1_wk_neut');
    log_rate_ne21_wk_minus_na21.append(i)
    log_rate_ne21_wk_minus_na21[i]=prof[i].get('log_rate_ne21_wk-minus_na21');
    log_rate_ne22_wk_minus_na22.append(i)
    log_rate_ne22_wk_minus_na22[i]=prof[i].get('log_rate_ne22_wk-minus_na22');
    log_rate_mg25_wk_minus_al25.append(i)
    log_rate_mg25_wk_minus_al25[i]=prof[i].get('log_rate_mg25_wk-minus_al25');
    log_rate_mg24_wk_na24.append(i)
    log_rate_mg24_wk_na24[i]=prof[i].get('log_rate_mg24_wk_na24');
    log_rate_neut_wk_minus_h1.append(i)
    log_rate_neut_wk_minus_h1[i]=prof[i].get('log_rate_neut_wk-minus_h1');
    log_rate_o14_wk_n14.append(i)
    log_rate_o14_wk_n14[i]=prof[i].get('log_rate_o14_wk_n14');
    log_rate_o15_wk_n15.append(i)
    log_rate_o15_wk_n15[i]=prof[i].get('log_rate_o15_wk_n15');
    log_rate_f17_wk_o17.append(i)
    log_rate_f17_wk_o17[i]=prof[i].get('log_rate_f17_wk_o17');
    log_rate_o17_wk_minus_f17.append(i)
    log_rate_o17_wk_minus_f17[i]=prof[i].get('log_rate_o17_wk-minus_f17');
    log_rate_mg26_wk_minus_al26.append(i)
    log_rate_mg26_wk_minus_al26[i]=prof[i].get('log_rate_mg26_wk-minus_al26');
    log_rate_o18_wk_minus_f18.append(i)
    log_rate_o18_wk_minus_f18[i]=prof[i].get('log_rate_o18_wk-minus_f18');
    log_rate_f19_wk_minus_ne19.append(i)
    log_rate_f19_wk_minus_ne19[i]=prof[i].get('log_rate_f19_wk-minus_ne19');

    log_rate_n13_np_c13.append(i)
    log_rate_n13_np_c13[i]=prof[i].get('log_rate_n13_np_c13');
    log_rate_mg23_np_na23.append(i)
    log_rate_mg23_np_na23[i]=prof[i].get('log_rate_mg23_np_na23');
    log_rate_na21_np_ne21.append(i)
    log_rate_na21_np_ne21[i]=prof[i].get('log_rate_na21_np_ne21');
    log_rate_na22_np_ne22.append(i)
    log_rate_na22_np_ne22[i]=prof[i].get('log_rate_na22_np_ne22');
    log_rate_na23_np_ne23.append(i)
    log_rate_na23_np_ne23[i]=prof[i].get('log_rate_na23_np_ne23');
    log_rate_mg24_np_na24.append(i)
    log_rate_mg24_np_na24[i]=prof[i].get('log_rate_mg24_np_na24');
    log_rate_o16_na_c13.append(i)
    log_rate_o16_na_c13[i]=prof[i].get('log_rate_o16_na_c13');
    log_rate_mg24_ng_mg25.append(i)
    log_rate_mg24_ng_mg25[i]=prof[i].get('log_rate_mg24_ng_mg25');
    log_rate_na22_ng_na23.append(i)
    log_rate_na22_ng_na23[i]=prof[i].get('log_rate_na22_ng_na23');
    log_rate_na23_ng_na24.append(i)
    log_rate_na23_ng_na24[i]=prof[i].get('log_rate_na23_ng_na24');
    log_rate_ne21_ng_ne22.append(i)
    log_rate_ne21_ng_ne22[i]=prof[i].get('log_rate_ne21_ng_ne22');
    log_rate_ne20_ng_ne21.append(i)
    log_rate_ne20_ng_ne21[i]=prof[i].get('log_rate_ne20_ng_ne21');
    log_rate_ne19_ng_ne20.append(i)
    log_rate_ne19_ng_ne20[i]=prof[i].get('log_rate_ne19_ng_ne20');
    log_rate_ne22_ng_ne23.append(i)
    log_rate_ne22_ng_ne23[i]=prof[i].get('log_rate_ne22_ng_ne23');
    log_rate_c12_ng_c13.append(i)
    log_rate_c12_ng_c13[i]=prof[i].get('log_rate_c12_ng_c13');
    log_rate_mg23_ng_mg24.append(i)
    log_rate_mg23_ng_mg24[i]=prof[i].get('log_rate_mg23_ng_mg24');
    log_rate_n13_ng_n14.append(i)   
    log_rate_n13_ng_n14[i]=prof[i].get('log_rate_n13_ng_n14');


    mg23.append(i)
    mg23[i]=prof[i].get('mg23');
    na21.append(i)
    na21[i]=prof[i].get('na21');
    ne21.append(i)
    ne21[i]=prof[i].get('ne21');
    f18.append(i)
    f18[i]=prof[i].get('f18');
    ne18.append(i)
    ne18[i]=prof[i].get('ne18');
    ne19.append(i)
    ne19[i]=prof[i].get('ne19');
    f19.append(i)
    f19[i]=prof[i].get('f19');
    f17.append(i)
    f17[i]=prof[i].get('f17');
    o17.append(i)
    o17[i]=prof[i].get('o17');
    na24.append(i)
    na24[i]=prof[i].get('na24');


    t_c.append(i)
    t_c[i]=prof[i].get('t_c');
    pressure_scale_height.append(i)
    pressure_scale_height[i]=prof[i].get('pressure_scale_height')
    conv_vel.append(i)
    conv_vel[i]=prof[i].get('conv_vel')
    t_heat.append(i)
    t_heat[i]=prof[i].get('t_heat');




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




import pickle, pprint # Create the data structures


output_0=open('age.pkl','wb')
pickle.dump(age,output_0)
output_0.close()



output_1=open('logxq.pkl','wb')
pickle.dump(logxq,output_1)
output_1.close()

output_2=open('radius.pkl','wb')
pickle.dump(radius,output_2)
output_2.close()

output_3=open('logT.pkl','wb')
pickle.dump(logT,output_3)
output_3.close()

output_4=open('logRho.pkl','wb')
pickle.dump(logRho,output_4)
output_4.close()

output_5=open('logP.pkl','wb')
pickle.dump(logP,output_5)
output_5.close()

output_6=open('lum.pkl','wb')
pickle.dump(lum,output_6)
output_6.close()

output_7=open('q.pkl','wb')
pickle.dump(q,output_7)
output_7.close()

output_8=open('ye.pkl','wb')
pickle.dump(ye,output_8)
output_8.close()



output_9=open('pp.pkl','wb')
pickle.dump(pp,output_9)
output_9.close()

output_10=open('cno.pkl','wb')
pickle.dump(cno,output_10)
output_10.close()

output_11=open('tri_alpha.pkl','wb')
pickle.dump(tri_alpha,output_11)
output_11.close()

output_12=open('c12_c12.pkl','wb')
pickle.dump(c12_c12,output_12)
output_12.close()

output_13=open('o16_o16.pkl','wb')
pickle.dump(o16_o16,output_13)
output_13.close()

output_14=open('burn_c.pkl','wb')
pickle.dump(burn_c,output_14)
output_14.close()

output_15=open('burn_ne.pkl','wb')
pickle.dump(burn_ne,output_15)
output_15.close()

output_16=open('burn_mg.pkl','wb')
pickle.dump(burn_mg,output_16)
output_16.close()

output_17=open('gradr_sub_grada.pkl','wb')
pickle.dump(gradr_sub_grada,output_17)
output_17.close()

output_18=open('eps_nuc.pkl','wb')
pickle.dump(eps_nuc,output_18)
output_18.close()

output_19=open('non_nuc_neu.pkl','wb')
pickle.dump(non_nuc_neu,output_19)
output_19.close()



output_20=open('h1.pkl','wb')
pickle.dump(h1,output_20)
output_20.close()

output_21=open('he3.pkl','wb')
pickle.dump(he3,output_21)
output_21.close()

output_22=open('he4.pkl','wb')
pickle.dump(he4,output_22)
output_22.close()

output_23=open('c12.pkl','wb')
pickle.dump(c12,output_23)
output_23.close()

output_24=open('c13.pkl','wb')
pickle.dump(c13,output_24)
output_24.close()

output_25=open('n13.pkl','wb')
pickle.dump(n13,output_25)
output_25.close()

output_26=open('n14.pkl','wb')
pickle.dump(n14,output_26)
output_26.close()

output_27=open('o16.pkl','wb')
pickle.dump(o16,output_27)
output_27.close()

output_28=open('o18.pkl','wb')
pickle.dump(o18,output_28)
output_28.close()

output_29=open('ne20.pkl','wb')
pickle.dump(ne20,output_29)
output_29.close()

output_30=open('ne22.pkl','wb')
pickle.dump(ne22,output_30)
output_30.close()

output_31=open('ne23.pkl','wb')
pickle.dump(ne23,output_31)
output_31.close()

output_32=open('na23.pkl','wb')
pickle.dump(na23,output_32)
output_32.close()

output_33=open('mg24.pkl','wb')
pickle.dump(mg24,output_33)
output_33.close()

output_34=open('si28.pkl','wb')
pickle.dump(si28,output_34)
output_34.close()

output_35=open('s32.pkl','wb')
pickle.dump(s32,output_35)
output_35.close()

output_36=open('ar36.pkl','wb')
pickle.dump(ar36,output_36)
output_36.close()

output_37=open('ca40.pkl','wb')
pickle.dump(ca40,output_37)
output_37.close()

output_38=open('cr48.pkl','wb')
pickle.dump(cr48,output_38)
output_38.close()

output_39=open('ni56.pkl','wb')
pickle.dump(ni56,output_39)
output_39.close()

output_40=open('fe56.pkl','wb')
pickle.dump(fe56,output_40)
output_40.close()


output_41=open('log_g.pkl','wb')
pickle.dump(log_g,output_41)
output_41.close()

output_42=open('log_opacity.pkl','wb')
pickle.dump(log_opacity,output_42)
output_42.close()

output_43=open('gam.pkl','wb')
pickle.dump(gam,output_43)
output_43.close()

output_44=open('net_nuclear_energy.pkl','wb')
pickle.dump(net_nuclear_energy,output_44)
output_44.close()

output_45=open('eps_nuc_neu_total.pkl','wb')
pickle.dump(eps_nuc_neu_total,output_45)
output_45.close()

output_46=open('c12_o16.pkl','wb')
pickle.dump(c12_o16,output_46)
output_46.close()

output_47=open('eps_nuc_plus_nuc_neu.pkl','wb')
pickle.dump(eps_nuc_plus_nuc_neu,output_47)
output_47.close()

output_48=open('neut.pkl','wb')
pickle.dump(neut,output_48)
output_48.close()

output_49=open('photo.pkl','wb')
pickle.dump(neut,output_49)
output_49.close()

output_49=open('photo.pkl','wb')
pickle.dump(neut,output_49)
output_49.close()


output_50=open('mg25.pkl','wb')
pickle.dump(mg25,output_50)
output_50.close()


output_51=open('mg26.pkl','wb')
pickle.dump(mg26,output_51)
output_51.close()


output_52=open('log_rate_na23_wk_ne23.pkl','wb')
pickle.dump(log_rate_na23_wk_ne23,output_52)
output_52.close()


output_53=open('log_rate_ne23_wk_minus_na23.pkl','wb')
pickle.dump(log_rate_ne23_wk_minus_na23,output_53)
output_53.close()


output_54=open('log_rate_ne23_pn_na23.pkl','wb')
pickle.dump(log_rate_ne23_pn_na23,output_54)
output_54.close()


output_55=open('log_rate_mg23_wk_na23.pkl','wb')
pickle.dump(log_rate_mg23_wk_na23,output_55)
output_55.close()


output_56=open('log_rate_n13_wk_c13.pkl','wb')
pickle.dump(log_rate_n13_wk_c13,output_56)
output_56.close()


output_57=open('log_rate_h1_h1_wk_h2.pkl','wb')
pickle.dump(log_rate_h1_h1_wk_h2,output_57)
output_57.close()


output_58=open('log_rate_be7_wk_li7.pkl','wb')
pickle.dump(log_rate_be7_wk_li7,output_58)
output_58.close()


output_59=open('log_rate_al25_wk_mg25.pkl','wb')
pickle.dump(log_rate_al25_wk_mg25,output_59)
output_59.close()


output_60=open('log_rate_na21_wk_ne21.pkl','wb')
pickle.dump(log_rate_na21_wk_ne21,output_60)
output_60.close()


output_61=open('log_rate_h1_he3_wk_he4.pkl','wb')
pickle.dump(log_rate_h1_he3_wk_he4,output_61)
output_61.close()


output_62=open('log_rate_na22_wk_ne22.pkl','wb')
pickle.dump(log_rate_na22_wk_ne22,output_62)
output_62.close()


output_63=open('log_rate_f18_wk_minus_ne18.pkl','wb')
pickle.dump(log_rate_f18_wk_minus_ne18,output_63)
output_63.close()


output_64=open('log_rate_b8_wk_he4_he4.pkl','wb')
pickle.dump(log_rate_b8_wk_he4_he4,output_64)
output_64.close()


output_65=open('log_rate_na23_wk_minus_mg23.pkl','wb')
pickle.dump(log_rate_na23_wk_minus_mg23,output_65)
output_65.close()


output_66=open('log_rate_na24_wk_minus_mg24.pkl','wb')
pickle.dump(log_rate_na24_wk_minus_mg24,output_66)
output_66.close()


output_67=open('log_rate_f18_wk_o18.pkl','wb')
pickle.dump(log_rate_f18_wk_o18,output_67)
output_67.close()


output_68=open('log_rate_ne18_wk_f18.pkl','wb')
pickle.dump(log_rate_ne18_wk_f18,output_68)
output_68.close()


output_69=open('log_rate_al26_wk_mg26.pkl','wb')
pickle.dump(log_rate_al26_wk_mg26,output_69)
output_69.close()


output_70=open('log_rate_ne19_wk_f19.pkl','wb')
pickle.dump(log_rate_ne19_wk_f19,output_70)
output_70.close()


output_71=open('log_rate_h1_wk_neut.pkl','wb')
pickle.dump(log_rate_h1_wk_neut,output_71)
output_71.close()


output_72=open('log_rate_ne21_wk_minus_na21.pkl','wb')
pickle.dump(log_rate_ne21_wk_minus_na21,output_72)
output_72.close()


output_73=open('log_rate_ne22_wk_minus_na22.pkl','wb')
pickle.dump(log_rate_ne22_wk_minus_na22,output_73)
output_73.close()


output_74=open('log_rate_mg25_wk_minus_al25.pkl','wb')
pickle.dump(log_rate_mg25_wk_minus_al25,output_74)
output_74.close()


output_75=open('log_rate_mg24_wk_na24.pkl','wb')
pickle.dump(log_rate_mg24_wk_na24,output_75)
output_75.close()


output_76=open('log_rate_neut_wk_minus_h1.pkl','wb')
pickle.dump(log_rate_neut_wk_minus_h1,output_76)
output_76.close()


output_77=open('log_rate_o14_wk_n14.pkl','wb')
pickle.dump(log_rate_o14_wk_n14,output_77)
output_77.close()


output_78=open('log_rate_o15_wk_n15.pkl','wb')
pickle.dump(log_rate_o15_wk_n15,output_78)
output_78.close()


output_79=open('log_rate_f17_wk_o17.pkl','wb')
pickle.dump(log_rate_f17_wk_o17,output_79)
output_79.close()


output_80=open('log_rate_o17_wk_minus_f17.pkl','wb')
pickle.dump(log_rate_o17_wk_minus_f17,output_80)
output_80.close()


output_81=open('log_rate_mg26_wk_minus_al26.pkl','wb')
pickle.dump(log_rate_mg26_wk_minus_al26,output_81)
output_81.close()


output_82=open('log_rate_o18_wk_minus_f18.pkl','wb')
pickle.dump(log_rate_o18_wk_minus_f18,output_82)
output_82.close()


output_83=open('log_rate_f19_wk_minus_ne19.pkl','wb')
pickle.dump(log_rate_f19_wk_minus_ne19,output_83)
output_83.close()


output_84=open('log_rate_n13_np_c13.pkl','wb')
pickle.dump(log_rate_n13_np_c13,output_84)
output_84.close()


output_85=open('log_rate_mg23_np_na23.pkl','wb')
pickle.dump(log_rate_mg23_np_na23,output_85)
output_85.close()


output_86=open('log_rate_na21_np_ne21.pkl','wb')
pickle.dump(log_rate_na21_np_ne21,output_86)
output_86.close()


output_87=open('log_rate_na22_np_ne22.pkl','wb')
pickle.dump(log_rate_na22_np_ne22,output_87)
output_87.close()


output_88=open('log_rate_na23_np_ne23.pkl','wb')
pickle.dump(log_rate_na23_np_ne23,output_88)
output_88.close()


output_89=open('log_rate_mg24_np_na24.pkl','wb')
pickle.dump(log_rate_mg24_np_na24,output_89)
output_89.close()


output_90=open('log_rate_o16_na_c13.pkl','wb')
pickle.dump(log_rate_o16_na_c13,output_90)
output_90.close()


output_91=open('log_rate_mg24_ng_mg25.pkl','wb')
pickle.dump(log_rate_mg24_ng_mg25,output_91)
output_91.close()


output_92=open('log_rate_na22_ng_na23.pkl','wb')
pickle.dump(log_rate_na22_ng_na23,output_92)
output_92.close()


output_93=open('log_rate_na23_ng_na24.pkl','wb')
pickle.dump(log_rate_na23_ng_na24,output_93)
output_93.close()


output_94=open('log_rate_ne21_ng_ne22.pkl','wb')
pickle.dump(log_rate_ne21_ng_ne22,output_94)
output_94.close()


output_95=open('log_rate_ne20_ng_ne21.pkl','wb')
pickle.dump(log_rate_ne20_ng_ne21,output_95)
output_95.close()


output_96=open('log_rate_ne19_ng_ne20.pkl','wb')
pickle.dump(log_rate_ne19_ng_ne20,output_96)
output_96.close()


output_97=open('log_rate_ne22_ng_ne23.pkl','wb')
pickle.dump(log_rate_ne22_ng_ne23,output_97)
output_97.close()


output_98=open('log_rate_c12_ng_c13.pkl','wb')
pickle.dump(log_rate_c12_ng_c13,output_98)
output_98.close()


output_99=open('log_rate_mg23_ng_mg24.pkl','wb')
pickle.dump(log_rate_mg23_ng_mg24,output_99)
output_99.close()


output_100=open('log_rate_n13_ng_n14.pkl','wb')
pickle.dump(log_rate_n13_ng_n14,output_100)
output_100.close()


output_101=open('mg23.pkl','wb')
pickle.dump(mg23,output_101)
output_101.close()


output_102=open('na21.pkl','wb')
pickle.dump(na21,output_102)
output_102.close()


output_103=open('ne21.pkl','wb')
pickle.dump(ne21,output_103)
output_103.close()


output_104=open('f18.pkl','wb')
pickle.dump(f18,output_104)
output_104.close()


output_105=open('ne18.pkl','wb')
pickle.dump(ne18,output_105)
output_105.close()


output_106=open('ne19.pkl','wb')
pickle.dump(ne19,output_106)
output_106.close()


output_107=open('f19.pkl','wb')
pickle.dump(f19,output_107)
output_107.close()


output_108=open('f17.pkl','wb')
pickle.dump(f17,output_108)
output_108.close()


output_109=open('o17.pkl','wb')
pickle.dump(o17,output_109)
output_109.close()


output_110=open('na24.pkl','wb')
pickle.dump(na24,output_110)
output_110.close()


output_111=open('t_c.pkl','wb')
pickle.dump(t_c,output_111)
output_111.close()


output_112=open('pressure_scale_height.pkl','wb')
pickle.dump(pressure_scale_height,output_112)
output_112.close()


output_113=open('conv_vel.pkl','wb')
pickle.dump(conv_vel,output_113)
output_113.close()


output_114=open('t_heat.pkl','wb')
pickle.dump(t_heat,output_114)
output_114.close()



# Import the pickled data in a new script/notebook

# import pickle, pprint

# pkl_file_0=open(path_data + 'age.pkl','rb')
# age=pickle.load(pkl_file_0)

# pkl_file_1=open(path_data + 'logxq.pkl','rb')
# logxq=pickle.load(pkl_file_1)

# pkl_file_2=open(path_data + 'radius.pkl','rb')
# radius=pickle.load(pkl_file_2)

# pkl_file_3=open(path_data + 'logT.pkl','rb')
# logT=pickle.load(pkl_file_3)

# pkl_file_4=open(path_data + 'logRho.pkl','rb')
# logRho=pickle.load(pkl_file_4)

# pkl_file_5=open(path_data + 'logP.pkl','rb')
# logP=pickle.load(pkl_file_5)

# pkl_file_6=open(path_data + 'lum.pkl','rb')
# lum=pickle.load(pkl_file_6)

# pkl_file_7=open(path_data + 'q.pkl','rb')
# q=pickle.load(pkl_file_7)

# pkl_file_8=open(path_data + 'ye.pkl','rb')
# ye=pickle.load(pkl_file_8)

# pkl_file_9=open(path_data + 'pp.pkl','rb')
# pp=pickle.load(pkl_file_9)

# pkl_file_10=open(path_data + 'cno.pkl','rb')
# cno=pickle.load(pkl_file_10)

# pkl_file_11=open(path_data + 'tri_alpha.pkl','rb')
# tri_alpha=pickle.load(pkl_file_11)

# pkl_file_12=open(path_data + 'c12_c12.pkl','rb')
# c12_c12=pickle.load(pkl_file_12)

# pkl_file_13=open(path_data + 'o16_o16.pkl','rb')
# o16_o16=pickle.load(pkl_file_13)

# pkl_file_14=open(path_data + 'burn_c.pkl','rb')
# burn_c=pickle.load(pkl_file_14)

# pkl_file_15=open(path_data + 'burn_ne.pkl','rb')
# burn_ne=pickle.load(pkl_file_15)

# pkl_file_16=open(path_data + 'burn_mg.pkl','rb')
# burn_mg=pickle.load(pkl_file_16)

# pkl_file_17=open(path_data + 'gradr_sub_grada.pkl','rb')
# gradr_sub_grada=pickle.load(pkl_file_17)

# pkl_file_18=open(path_data + 'eps_nuc.pkl','rb')
# eps_nuc=pickle.load(pkl_file_18)

# pkl_file_19=open(path_data + 'non_nuc_neu.pkl','rb')
# non_nuc_neu=pickle.load(pkl_file_19)

# pkl_file_20=open(path_data + 'h1.pkl','rb')
# h1=pickle.load(pkl_file_20)

# pkl_file_21=open(path_data + 'he3.pkl','rb')
# he3=pickle.load(pkl_file_21)

# pkl_file_22=open(path_data + 'he4.pkl','rb')
# he4=pickle.load(pkl_file_22)

# pkl_file_23=open(path_data + 'c12.pkl','rb')
# c12=pickle.load(pkl_file_23)

# pkl_file_24=open(path_data + 'c13.pkl','rb')
# c13=pickle.load(pkl_file_24)

# pkl_file_25=open(path_data + 'n13.pkl','rb')
# n13=pickle.load(pkl_file_25)

# pkl_file_26=open(path_data + 'n14.pkl','rb')
# n14=pickle.load(pkl_file_26)

# pkl_file_27=open(path_data + 'o16.pkl','rb')
# o16=pickle.load(pkl_file_27)

# pkl_file_28=open(path_data + 'o18.pkl','rb')
# o18=pickle.load(pkl_file_28)

# pkl_file_29=open(path_data + 'ne20.pkl','rb')
# ne20=pickle.load(pkl_file_29)

# pkl_file_30=open(path_data + 'ne22.pkl','rb')
# ne22=pickle.load(pkl_file_30)

# pkl_file_31=open(path_data + 'ne23.pkl','rb')
# ne23=pickle.load(pkl_file_31)

# pkl_file_32=open(path_data + 'na23.pkl','rb')
# na23=pickle.load(pkl_file_32)

# pkl_file_33=open(path_data + 'mg24.pkl','rb')
# mg24=pickle.load(pkl_file_33)

# pkl_file_34=open(path_data + 'si28.pkl','rb')
# si28=pickle.load(pkl_file_34)

# pkl_file_35=open(path_data + 's32.pkl','rb')
# s32=pickle.load(pkl_file_35)

# pkl_file_36=open(path_data + 'ar36.pkl','rb')
# ar36=pickle.load(pkl_file_36)

# pkl_file_37=open(path_data + 'ca40.pkl','rb')
# ca40=pickle.load(pkl_file_37)

# pkl_file_38=open(path_data + 'cr48.pkl','rb')
# cr48=pickle.load(pkl_file_38)

# pkl_file_39=open(path_data + 'ni56.pkl','rb')
# ni56=pickle.load(pkl_file_39)

# pkl_file_40=open(path_data + 'fe56.pkl','rb')
# fe56=pickle.load(pkl_file_40)

# pkl_file_41=open(path_data + 'log_g.pkl','rb')
# log_g=pickle.load(pkl_file_41)

# pkl_file_42=open(path_data + 'log_opacity.pkl','rb')
# log_opacity=pickle.load(pkl_file_42)

# pkl_file_43=open(path_data + 'gam.pkl','rb')
# gam=pickle.load(pkl_file_43)

# pkl_file_44=open(path_data + 'net_nuclear_energy.pkl','rb')
# net_nuclear_energy=pickle.load(pkl_file_44)

# pkl_file_45=open(path_data + 'eps_nuc_neu_total.pkl','rb')
# eps_nuc_neu_total=pickle.load(pkl_file_45)

# pkl_file_46=open(path_data + 'cno.pkl','rb')
# cno=pickle.load(pkl_file_46)

# pkl_file_47=open(path_data + 'neut.pkl','rb')
# neut=pickle.load(pkl_file_47)

# pkl_file_48=open(path_data + 'photo.pkl','rb')
# photo=pickle.load(pkl_file_48)

# pkl_file_49=open(path_data + 'c12_o16.pkl','rb')
# c12_o16=pickle.load(pkl_file_49)

# pkl_file_50=open(path_data + 'mg25.pkl','rb')
# mg25=pickle.load(pkl_file_50)

# pkl_file_51=open(path_data + 'mg26.pkl','rb')
# mg26=pickle.load(pkl_file_51)


# pkl_file_52=open(path_data + 'log_rate_na23_wk_ne23.pkl','rb')
# log_rate_na23_wk_ne23=pickle.load(pkl_file_52)

# pkl_file_53=open(path_data + 'log_rate_ne23_wk_minus_na23.pkl','rb')
# log_rate_ne23_wk_minus_na23=pickle.load(pkl_file_53)

# pkl_file_54=open(path_data + 'log_rate_ne23_pn_na23.pkl','rb')
# log_rate_ne23_pn_na23=pickle.load(pkl_file_54)

# pkl_file_55=open(path_data + 'log_rate_mg23_wk_na23.pkl','rb')
# log_rate_mg23_wk_na23=pickle.load(pkl_file_55)

# pkl_file_56=open(path_data + 'log_rate_n13_wk_c13.pkl','rb')
# log_rate_n13_wk_c13=pickle.load(pkl_file_56)

# pkl_file_57=open(path_data + 'log_rate_h1_h1_wk_h2.pkl','rb')
# log_rate_h1_h1_wk_h2=pickle.load(pkl_file_57)

# pkl_file_58=open(path_data + 'log_rate_be7_wk_li7.pkl','rb')
# log_rate_be7_wk_li7=pickle.load(pkl_file_58)

# pkl_file_59=open(path_data + 'log_rate_al25_wk_mg25.pkl','rb')
# log_rate_al25_wk_mg25=pickle.load(pkl_file_59)

# pkl_file_60=open(path_data + 'log_rate_na21_wk_ne21.pkl','rb')
# log_rate_na21_wk_ne21=pickle.load(pkl_file_60)

# pkl_file_61=open(path_data + 'log_rate_h1_he3_wk_he4.pkl','rb')
# log_rate_h1_he3_wk_he4=pickle.load(pkl_file_61)

# pkl_file_62=open(path_data + 'log_rate_na22_wk_ne22.pkl','rb')
# log_rate_na22_wk_ne22=pickle.load(pkl_file_62)

# pkl_file_63=open(path_data + 'log_rate_f18_wk_minus_ne18.pkl','rb')
# log_rate_f18_wk_minus_ne18=pickle.load(pkl_file_63)

# pkl_file_64=open(path_data + 'log_rate_b8_wk_he4_he4.pkl','rb')
# log_rate_b8_wk_he4_he4=pickle.load(pkl_file_64)

# pkl_file_65=open(path_data + 'log_rate_na23_wk_minus_mg23.pkl','rb')
# log_rate_na23_wk_minus_mg23=pickle.load(pkl_file_65)

# pkl_file_66=open(path_data + 'log_rate_na24_wk_minus_mg24.pkl','rb')
# log_rate_na24_wk_minus_mg24=pickle.load(pkl_file_66)

# pkl_file_67=open(path_data + 'log_rate_f18_wk_o18.pkl','rb')
# log_rate_f18_wk_o18=pickle.load(pkl_file_67)

# pkl_file_68=open(path_data + 'log_rate_ne18_wk_f18.pkl','rb')
# log_rate_ne18_wk_f18=pickle.load(pkl_file_68)

# pkl_file_69=open(path_data + 'log_rate_al26_wk_mg26.pkl','rb')
# log_rate_al26_wk_mg26=pickle.load(pkl_file_69)

# pkl_file_70=open(path_data + 'log_rate_ne19_wk_f19.pkl','rb')
# log_rate_ne19_wk_f19=pickle.load(pkl_file_70)

# pkl_file_71=open(path_data + 'log_rate_h1_wk_neut.pkl','rb')
# log_rate_h1_wk_neut=pickle.load(pkl_file_71)

# pkl_file_72=open(path_data + 'log_rate_ne21_wk_minus_na21.pkl','rb')
# log_rate_ne21_wk_minus_na21=pickle.load(pkl_file_72)

# pkl_file_73=open(path_data + 'log_rate_ne22_wk_minus_na22.pkl','rb')
# log_rate_ne22_wk_minus_na22=pickle.load(pkl_file_73)

# pkl_file_74=open(path_data + 'log_rate_mg25_wk_minus_al25.pkl','rb')
# log_rate_mg25_wk_minus_al25=pickle.load(pkl_file_74)

# pkl_file_75=open(path_data + 'log_rate_mg24_wk_na24.pkl','rb')
# log_rate_mg24_wk_na24=pickle.load(pkl_file_75)

# pkl_file_76=open(path_data + 'log_rate_neut_wk_minus_h1.pkl','rb')
# log_rate_neut_wk_minus_h1=pickle.load(pkl_file_76)

# pkl_file_77=open(path_data + 'log_rate_o14_wk_n14.pkl','rb')
# log_rate_o14_wk_n14=pickle.load(pkl_file_77)

# pkl_file_78=open(path_data + 'log_rate_o15_wk_n15.pkl','rb')
# log_rate_o15_wk_n15=pickle.load(pkl_file_78)

# pkl_file_79=open(path_data + 'log_rate_f17_wk_o17.pkl','rb')
# log_rate_f17_wk_o17=pickle.load(pkl_file_79)

# pkl_file_80=open(path_data + 'log_rate_o17_wk_minus_f17.pkl','rb')
# log_rate_o17_wk_minus_f17=pickle.load(pkl_file_80)

# pkl_file_81=open(path_data + 'log_rate_mg26_wk_minus_al26.pkl','rb')
# log_rate_mg26_wk_minus_al26=pickle.load(pkl_file_81)

# pkl_file_82=open(path_data + 'log_rate_o18_wk_minus_f18.pkl','rb')
# log_rate_o18_wk_minus_f18=pickle.load(pkl_file_82)

# pkl_file_83=open(path_data + 'log_rate_f19_wk_minus_ne19.pkl','rb')
# log_rate_f19_wk_minus_ne19=pickle.load(pkl_file_83)

# pkl_file_84=open(path_data + 'log_rate_n13_np_c13.pkl','rb')
# log_rate_n13_np_c13=pickle.load(pkl_file_84)

# pkl_file_85=open(path_data + 'log_rate_mg23_np_na23.pkl','rb')
# log_rate_mg23_np_na23=pickle.load(pkl_file_85)

# pkl_file_86=open(path_data + 'log_rate_na21_np_ne21.pkl','rb')
# log_rate_na21_np_ne21=pickle.load(pkl_file_86)

# pkl_file_87=open(path_data + 'log_rate_na22_np_ne22.pkl','rb')
# log_rate_na22_np_ne22=pickle.load(pkl_file_87)

# pkl_file_88=open(path_data + 'log_rate_na23_np_ne23.pkl','rb')
# log_rate_na23_np_ne23=pickle.load(pkl_file_88)

# pkl_file_89=open(path_data + 'log_rate_mg24_np_na24.pkl','rb')
# log_rate_mg24_np_na24=pickle.load(pkl_file_89)

# pkl_file_90=open(path_data + 'log_rate_o16_na_c13.pkl','rb')
# log_rate_o16_na_c13=pickle.load(pkl_file_90)

# pkl_file_91=open(path_data + 'log_rate_mg24_ng_mg25.pkl','rb')
# log_rate_mg24_ng_mg25=pickle.load(pkl_file_91)

# pkl_file_92=open(path_data + 'log_rate_na22_ng_na23.pkl','rb')
# log_rate_na22_ng_na23=pickle.load(pkl_file_92)

# pkl_file_93=open(path_data + 'log_rate_na23_ng_na24.pkl','rb')
# log_rate_na23_ng_na24=pickle.load(pkl_file_93)

# pkl_file_94=open(path_data + 'log_rate_ne21_ng_ne22.pkl','rb')
# log_rate_ne21_ng_ne22=pickle.load(pkl_file_94)

# pkl_file_95=open(path_data + 'log_rate_ne20_ng_ne21.pkl','rb')
# log_rate_ne20_ng_ne21=pickle.load(pkl_file_95)

# pkl_file_96=open(path_data + 'log_rate_ne19_ng_ne20.pkl','rb')
# log_rate_ne19_ng_ne20=pickle.load(pkl_file_96)

# pkl_file_97=open(path_data + 'log_rate_ne22_ng_ne23.pkl','rb')
# log_rate_ne22_ng_ne23=pickle.load(pkl_file_97)

# pkl_file_98=open(path_data + 'log_rate_c12_ng_c13.pkl','rb')
# log_rate_c12_ng_c13=pickle.load(pkl_file_98)

# pkl_file_99=open(path_data + 'log_rate_mg23_ng_mg24.pkl','rb')
# log_rate_mg23_ng_mg24=pickle.load(pkl_file_99)

# pkl_file_100=open(path_data + 'log_rate_n13_ng_n14.pkl','rb')
# log_rate_n13_ng_n14=pickle.load(pkl_file_100)

# pkl_file_101=open(path_data + 'mg23.pkl','rb')
# mg23=pickle.load(pkl_file_101)

# pkl_file_102=open(path_data + 'na21.pkl','rb')
# na21=pickle.load(pkl_file_102)

# pkl_file_103=open(path_data + 'ne21.pkl','rb')
# ne21=pickle.load(pkl_file_103)

# pkl_file_104=open(path_data + 'f18.pkl','rb')
# f18=pickle.load(pkl_file_104)

# pkl_file_105=open(path_data + 'ne18.pkl','rb')
# ne18=pickle.load(pkl_file_105)

# pkl_file_106=open(path_data + 'ne19.pkl','rb')
# ne19=pickle.load(pkl_file_106)

# pkl_file_107=open(path_data + 'f19.pkl','rb')
# f19=pickle.load(pkl_file_107)

# pkl_file_108=open(path_data + 'f17.pkl','rb')
# f17=pickle.load(pkl_file_108)

# pkl_file_109=open(path_data + 'o17.pkl','rb')
# o17=pickle.load(pkl_file_109)

# pkl_file_110=open(path_data + 'na24.pkl','rb')
# na24=pickle.load(pkl_file_110)

# pkl_file_111=open(path_data + 't_c.pkl','rb')
# t_c=pickle.load(pkl_file_111)

# pkl_file_112=open(path_data + 'conv_vel.pkl','rb')
# conv_vel=pickle.load(pkl_file_112)

# pkl_file_113=open(path_data + 'pressure_scale_height.pkl','rb')
# pressure_scale_height=pickle.load(pkl_file_113)

# pkl_file_114=open(path_data + 't_heat.pkl','rb')
# t_heat=pickle.load(pkl_file_114)