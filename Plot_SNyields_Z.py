import numpy as np                        ##IMPORTS
import os
import glob
import time
from scipy import misc
from scipy import stats
from scipy.interpolate import spline
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from astropy import constants as const
from astropy.io import fits
from astropy.io import ascii
from elements import ELEMENTS

#Read data from files
elemsToRead = ['C','O','Ne','Mg','Si','S','Ar','Ca','Cr','Mn','Fe','Ni']
atNumsToRead = [6,8,10,12,14,16,18,20,24,25,26,28]
massNumsToRead = np.asarray([12.0,18.0,20.18,24.31,28.09,32.07,39.95,40.08,52.0,54.9,55.85,58.69])
nElemsToRead = len(atNumsToRead)

#Model metallicities
zArrBr = np.asarray([0.00025, 0.0025, 0.01, 0.025, 0.075])
nMetallicities = len (zArrBr)
zAsplund = 0.014
zArrSol = zArrBr/zAsplund
zArrBrMod = np.asarray([0.00025, 0.0025, 0.075, 0.025, 0.01]) 
zArrSolMod = zArrBrMod/zAsplund

#####################################
#TYPE IA MODELS
#DDT models
nLayers = 162
dirName = '/Users/badenes/data/explosion_models/'
fileNameList = ['_decay.crank.z0.01.dat','_decay.crank.z0.1.dat','_decay.crank.dat','_decay.crank.z1.dat','_decay.crank.z3.dat']
nModels = len(fileNameList)
expModelList = ['DDTa','DDTc','DDTe','DDTf']
nExpModels = len(expModelList)
yieldArr = np.zeros([nElemsToRead,nExpModels,nModels])
profilesArr = np.zeros([nLayers,nElemsToRead,nExpModels,nModels])
lagMassProfilesArr = np.zeros([nLayers,nExpModels])
massProfilesArr = np.zeros([nLayers,nExpModels])
for expModel in range(nExpModels) :
    for model in range(nModels) :
        print 'Reading '+expModelList[expModel]+fileNameList[model] 
        with open(dirName+expModelList[expModel]+fileNameList[model], 'r') as f:
            for layer in range(nLayers) :
                words = f.readline().split()
                layerMassLag = float(words[1].replace('D','E'))/1.989e33
                if model == 0 : 
                    lagMassProfilesArr[layer,expModel] = layerMassLag
                    
                if layer == 0 : layerMass = layerMassLag
                else : layerMass = layerMassLag-massProfilesArr[layer-1,expModel]

                if model == 0 : 
                    massProfilesArr[layer,expModel] = layerMass
                    
                nNucLayer = int(words[2])
                nNucLinesLayer = int(np.ceil(nNucLayer/4))
                f.readline()
                readNucs = 0
                for line in range(nNucLinesLayer) :
                    words = f.readline().split()
                    nNucToRead = min(nNucLayer-readNucs,4)
                    for nuclide in range(nNucToRead) :
                        atNum = int(words[3*nuclide])
                        if (atNum in atNumsToRead) :
                            index = atNumsToRead.index(atNum)
                            profilesArr[layer,index,expModel,model] = float(words[3*nuclide+2].replace('D','E'))
                            yieldArr[index,expModel,model] += layerMass*float(words[3*nuclide+2].replace('D','E'))

                profilesArr[layer,:,expModel,model] = profilesArr[layer,:,expModel,model]/profilesArr[layer,:,expModel,model].sum()
                f.readline()

#SCH Models
nSChModels = 3
yieldArrSCh = np.zeros([nElemsToRead,nSChModels]) 
fileNameListSCh = ['SCH_decay.crank.dat','SCH3DOP_decay.crank.dat','SCH3DMP_decay.crank.dat']
nLayersArr = [115,130,130]
profilesSCh = np.zeros([130,nElemsToRead,3])
lagMassProfilesSCH =  np.zeros([130,3])
massProfilesSCH =  np.zeros([130,3])
for model in range(nSChModels) :
    print 'Reading '+fileNameListSCh[model] 
    with open(dirName+fileNameListSCh[model], 'r') as f:
        for layer in range(nLayersArr[model]) :
            words = f.readline().split()
            layerMassLag = float(words[1].replace('D','E'))/1.989e33
            lagMassProfilesSCH[layer,model] =  layerMassLag
            if layer == 0 : layerMass = layerMassLag
            else : layerMass = layerMassLag-lagMassProfilesSCH[layer-1,model]
                
            massProfilesSCH[layer,model] =  layerMass
            nNucLayer = int(words[2])
            nNucLinesLayer = int(np.ceil(nNucLayer/4.0))
            f.readline()
            readNucs = 0
            for line in range(nNucLinesLayer) :
                words = f.readline().split()
                nNucToRead = min(nNucLayer-readNucs,4)
                for nuclide in range(nNucToRead) :
                    atNum = int(words[3*nuclide])
                    massNum = int(words[3*nuclide+1])
                    if ((atNum == 26) and (massNum == 55)) : atNum = 25
                    readNucs += 1
                    if (atNum in atNumsToRead) :
                        index = atNumsToRead.index(atNum)
                        profilesSCh[layer,index,model] += float(words[3*nuclide+2].replace('D','E'))

            profilesSCh[layer,:,model] = profilesSCh[layer,:,model]*massNumsToRead
            profilesSCh[layer,:,model] = profilesSCh[layer,:,model]/profilesSCh[layer,:,model].sum()

    for element in range(nElemsToRead) : yieldArrSCh[element,model] = (profilesSCh[:,element,model]*massProfilesSCH[:,model]).sum()
                
#Read new models from Eduardo
yieldArrNew = np.zeros([nElemsToRead,nExpModels,nMetallicities])
yieldArrNewFePeakBR = np.zeros([4,3,nExpModels,nMetallicities])
profilesArrNew = np.zeros([nLayers,nElemsToRead,nExpModels,nMetallicities])
lagMassProfilesArrNew = np.zeros([nLayers,nExpModels])
with open(dirName+'newSNIamodels_mod.dat', 'r') as f:
    for expModel in range(nExpModels) :
        for metallicity in range(nMetallicities) :
            f.readline()
            layerMasses = []
            layerMassesLag = []
            for layer in range(nLayers) :
                words = f.readline().split()
                layerMassLag = float(words[0])
                if layer == 0 :
                    layerMass = layerMassLag
                else :
                    layerMass = layerMassLag-(np.asarray(layerMasses).sum())

                layerMasses.append(layerMass)
                layerMassesLag.append(layerMassLag)
                
                for element in range(nElemsToRead) :
                    yieldArrNew[element,expModel,metallicity] += layerMass*float(words[atNumsToRead[element]+1])
                    profilesArrNew[layer,element,expModel,metallicity] = float(words[atNumsToRead[element]+1])

                #separate yields by burning regime
                ySi = float(words[14+1])
                yCr = float(words[24+1])
                yMn = float(words[25+1])
                yFe = float(words[26+1])
                yNi = float(words[28+1])
                if ((yFe > 0.3) and (yNi > 0.03) and (yMn > 0.003)) :
                    #print expModel, layer, layerMassLag, 'nNSE!'
                    yieldArrNewFePeakBR[0,0,expModel,metallicity] += layerMass*yCr
                    yieldArrNewFePeakBR[1,0,expModel,metallicity] += layerMass*yMn
                    yieldArrNewFePeakBR[2,0,expModel,metallicity] += layerMass*yFe
                    yieldArrNewFePeakBR[3,0,expModel,metallicity] += layerMass*yNi

                elif ((yFe > 0.3) and (yNi > 0.01) and (yMn < 0.003)) :
                    #print expModel, layer, layerMassLag, 'NSE!'
                    yieldArrNewFePeakBR[0,1,expModel,metallicity] += layerMass*yCr
                    yieldArrNewFePeakBR[1,1,expModel,metallicity] += layerMass*yMn
                    yieldArrNewFePeakBR[2,1,expModel,metallicity] += layerMass*yFe
                    yieldArrNewFePeakBR[3,1,expModel,metallicity] += layerMass*yNi

                elif ((yFe > 0.01) and (ySi > 0.00001)) :
                    #print expModel, layer, layerMassLag, 'ESiB!'
                    yieldArrNewFePeakBR[0,2,expModel,metallicity] += layerMass*yCr
                    yieldArrNewFePeakBR[1,2,expModel,metallicity] += layerMass*yMn
                    yieldArrNewFePeakBR[2,2,expModel,metallicity] += layerMass*yFe
                    yieldArrNewFePeakBR[3,2,expModel,metallicity] += layerMass*yNi
                    

            if metallicity==0 : lagMassProfilesArrNew[:,expModel] = np.asarray(layerMassesLag)    

#Read SCH models from Eduardo
nExpModelsSCHSim = 4
yieldArrSCHSim = np.zeros([nElemsToRead,nExpModelsSCHSim,nMetallicities])
profilesArrSCHSim = np.zeros([nLayers,nElemsToRead,nExpModelsSCHSim,nMetallicities])
lagMassProfilesArrSCHSim = np.zeros([nLayers,nExpModelsSCHSim])
fileNameListSCHSim = ['0p88','0p97','1p06','1p15']
metNameListSCHSim = ['Z2p25e-4','Z2p25e-3','Z9e-3','Z2p25e-2','Z6p75e-2']
yieldArrNewFePeakBRSCH = np.zeros([4,3,nExpModelsSCHSim,nMetallicities])
#zArrSCHSim = np.array([0.01,0.1,0.33,1.0,3.0])
zArrSCHSim = np.array([0.00025,0.0025,0.01,0.025,0.075])
for expModel in range(nExpModelsSCHSim) :
    for metallicity in range(nMetallicities) :
        fileNameSCHSim = dirName+fileNameListSCHSim[expModel]+metNameListSCHSim[metallicity]+'.perfiln1.dat'
        print 'Reading '+fileNameSCHSim
        f = open(fileNameSCHSim, 'r')
        layerMasses = []
        layerMassesLag = []
        for layer in range(nLayers) :
            words = f.readline().split()
            layerMassLag = float(words[0])
            if layer == 0 :
                layerMass = layerMassLag
            else :
                layerMass = layerMassLag-(np.asarray(layerMasses).sum())
                
            layerMasses.append(layerMass)
            layerMassesLag.append(layerMassLag)
            
            for element in range(nElemsToRead) :
                yieldArrSCHSim[element,expModel,metallicity] += layerMass*float(words[atNumsToRead[element]+1])
                profilesArrSCHSim[layer,element,expModel,metallicity] = float(words[atNumsToRead[element]+1])

            #separate yields by burning regime
            ySi = float(words[14+1])
            yCr = float(words[24+1])
            yMn = float(words[25+1])
            yFe = float(words[26+1])
            yNi = float(words[28+1])
            if ((yFe > 0.3) and (ySi < 1e-5)) :
                #print expModel, layer, layerMassLag, 'NSE!'
                yieldArrNewFePeakBRSCH[0,1,expModel,metallicity] += layerMass*yCr
                yieldArrNewFePeakBRSCH[1,1,expModel,metallicity] += layerMass*yMn
                yieldArrNewFePeakBRSCH[2,1,expModel,metallicity] += layerMass*yFe
                yieldArrNewFePeakBRSCH[3,1,expModel,metallicity] += layerMass*yNi
                
            elif ((yFe > 0.01) and (ySi > 0.00001)) :
                #print expModel, layer, layerMassLag, 'ESiB!'
                yieldArrNewFePeakBRSCH[0,2,expModel,metallicity] += layerMass*yCr
                yieldArrNewFePeakBRSCH[1,2,expModel,metallicity] += layerMass*yMn
                yieldArrNewFePeakBRSCH[2,2,expModel,metallicity] += layerMass*yFe
                yieldArrNewFePeakBRSCH[3,2,expModel,metallicity] += layerMass*yNi
                
        if metallicity==0 : lagMassProfilesArrSCHSim[:,expModel] = np.asarray(layerMassesLag) 

        f.close()


#Read simmering models from Eduardo
nLayersSimmering = 161
yieldArrSimmering = np.zeros([nElemsToRead,2,nMetallicities])
profilesArrSimmering = np.zeros([nLayersSimmering,nElemsToRead,2,nMetallicities])
lagMassProfilesArrSimmering = np.zeros([nLayersSimmering,2])
with open(dirName+'SNIamodels_simmering_22Ne_mod.dat', 'r') as f:
    for expModel in range(2) :
        for metallicity in range(nMetallicities) :
            words = f.readline().split()
            layerMasses = []
            layerMassesLag = []
            #print words
            for layer in range(nLayersSimmering) :
                words = f.readline().split()
                layerMassLag = float(words[0])
                if layer == 0 :
                    layerMass = layerMassLag
                else :
                    layerMass = layerMassLag-(np.asarray(layerMasses).sum())
                    
                layerMasses.append(layerMass)
                layerMassesLag.append(layerMassLag)
                print layerMassLag
                
                for element in range(nElemsToRead) :
                    yieldArrSimmering[element,expModel,metallicity] += layerMass*float(words[atNumsToRead[element]+1])
                    profilesArrSimmering[layer,element,expModel,metallicity] = float(words[atNumsToRead[element]+1])
                    
        
            lagMassProfilesArrSimmering[:,expModel] = np.asarray(layerMassesLag) 

        
burnRegimeList = ['nNSE','NSE','ESiB']
with open('DDT_BurnRegimes.dat', 'w') as outFile :
    for expModel in range(nExpModels) :
        for metallicity in range(nMetallicities) :
            outFile.write('%4s %5s \n' % (expModelList[expModel], zArrBrMod[metallicity]))
            for burnReg in range(3) :
                outFile.write('%4s   Cr: %5e   Mn: %5e   Fe: %5e   Ni:%5e \n' % 
                              (burnRegimeList[burnReg],
                               yieldArrNewFePeakBR[0,burnReg,expModel,metallicity], 
                               yieldArrNewFePeakBR[1,burnReg,expModel,metallicity],
                               yieldArrNewFePeakBR[2,burnReg,expModel,metallicity],
                               yieldArrNewFePeakBR[3,burnReg,expModel,metallicity]))
            
with open('SCH_BurnRegimes.dat', 'w') as outFile :
    for expModel in range(nExpModelsSCHSim) :
        for metallicity in range(nMetallicities) :
            outFile.write('%4s %5s \n' % (fileNameListSCHSim[expModel], zArrSCHSim[metallicity]))
            for burnReg in range(3) :
                outFile.write('%4s   Cr: %5e   Mn: %5e   Fe: %5e   Ni:%5e \n' % 
                              (burnRegimeList[burnReg],
                               yieldArrNewFePeakBRSCH[0,burnReg,expModel,metallicity], 
                               yieldArrNewFePeakBRSCH[1,burnReg,expModel,metallicity],
                               yieldArrNewFePeakBRSCH[2,burnReg,expModel,metallicity],
                               yieldArrNewFePeakBRSCH[3,burnReg,expModel,metallicity]))


with open('DDT_Yields_Z.dat', 'w') as outFile :
    for expModel in range(nExpModels) :
        for metallicity in range(nMetallicities) :
            outFile.write('%4s %5s \n' % (expModelList[expModel], zArrBrMod[metallicity]))
            for i in range(nElemsToRead) :
                outFile.write('%s: %5e   ' % (elemsToRead[i], yieldArrNew[i,expModel,metallicity]))
            outFile.write('\n')

with open('SCH_Yields_Z.dat', 'w') as outFile :
    for expModel in range(nExpModelsSCHSim) :
        for metallicity in range(nMetallicities) :
            outFile.write('%4s %5s \n' % (fileNameListSCHSim[expModel], zArrSCHSim[metallicity]))
            for i in range(nElemsToRead) :
                outFile.write('%s: %5e   ' % (elemsToRead[i], yieldArrSCHSim[i,expModel,metallicity]))
            outFile.write('\n')


            
                
#Read CC models from Alex Heger
dirCC = '/Users/badenes/data/explosion_models_CC/'
fileNameListCC = ['s15a28c.expl_element','s20a28n.expl_element','s25a28d.expl_element']
nFilesCC = len(fileNameListCC)
yieldArrCC = np.zeros([nElemsToRead,nFilesCC])
for ccfile in range(nFilesCC):
    tab = ascii.read(dirCC+fileNameListCC[ccfile],data_start=1)
    nLines = len(tab)
    for element in range(nElemsToRead) :
        yieldArrCC[element,ccfile] = (tab[nLines-1])[2]*(tab[nLines-1])[atNumsToRead[element]+2]/1.989e33
    
#Read WoWe CC models
modelNamesCCWoWe = ['12A','15A','20A','22A']
nModelsCCWoWe = len(modelNamesCCWoWe)
elemLinesListWoWe = [[10,11],      #C
                     [15,16,17],   #O
                     [19,20,21],   #Ne
                     [24,25,26],   #Mg
                     [29,30,31],   #Si
                     [33,34,35,37],#S 
                     [41,43,44],   #Ar
                     [48,50,51,52],#Ca
                     [74,75,76,78,82],#Cr
                     [81,85,91],   #Mn
                     [84,86,87,98],#Fe
                     [100,102]]    #Ni
nMetallicitiesWoWe = 3
zArrWoWe = 0.9*zAsplund*np.asarray([0.01,0.1,1.0])
yieldArrCCWoWe = np.zeros([nElemsToRead,nModelsCCWoWe,nMetallicitiesWoWe])
tab5a = ascii.read(dirCC+'woosley_table5a.tex',delimiter='&',data_start=0,names=['Nuclide','S11A','S12A','S13A','S15A','S18A','S19A','S20A','S22A','S25A'])
for model in range(nModelsCCWoWe) :
    modelName = 'S'+modelNamesCCWoWe[model]
    for elemIndex in range(len(elemsToRead)) :
        for line in elemLinesListWoWe[elemIndex] :
            yieldArrCCWoWe[elemIndex,model,2] += tab5a[line][modelName] 

tab10a = ascii.read(dirCC+'woosley_table10a.tex',delimiter='&',data_start=0,names=['Nuclide','P12A','P13A','P15A','P18A','P20A','P22A','P25A'])
for model in range(nModelsCCWoWe) :
    modelName = 'P'+modelNamesCCWoWe[model]
    for elemIndex in range(len(elemsToRead)) :
        for line in elemLinesListWoWe[elemIndex] :
            yieldArrCCWoWe[elemIndex,model,1] += tab10a[line][modelName] 

tab12a = ascii.read(dirCC+'woosley_table10a.tex',delimiter='&',data_start=0,names=['Nuclide','T12A','T13A','T15A','T18A','T20A','T22A','T25A'])
for model in range(nModelsCCWoWe) :
    modelName = 'T'+modelNamesCCWoWe[model]
    for elemIndex in range(len(elemsToRead)) :
        for line in elemLinesListWoWe[elemIndex] :
            yieldArrCCWoWe[elemIndex,model,0] += tab12a[line][modelName] 
   
zSun = 0.0143
            
######################
#Plots 
## #plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'weight':'normal'})
plt.rcdefaults()
plt.rcParams.update({'figure.autolayout':'True'})
plt.rcParams.update({'font.size': 8})
plt.rcParams.update({'mathtext.default':'regular'})
plt.rcParams.update({'mathtext.fontset':'stixsans'})
plt.rcParams.update({'axes.linewidth': 1.5})
plt.rcParams.update({'xtick.major.size': 5})
plt.rcParams.update({'xtick.major.width': 1.25})
plt.rcParams.update({'xtick.minor.size': 2.5})
plt.rcParams.update({'xtick.minor.width': 1.25})
plt.rcParams.update({'ytick.major.size': 5})
plt.rcParams.update({'ytick.major.width': 1.25})
plt.rcParams.update({'ytick.minor.size': 2.5})
plt.rcParams.update({'ytick.minor.width': 1.25})

plotFileName = 'yieldsNew.pdf'
plt.figure(1,figsize = [11.0, 8.5])
plt.clf()
f, axArr = plt.subplots(3, 4, sharex = True, sharey = False, squeeze = True)   

for element in range(nElemsToRead) :

    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    #yLabel = symbol+' Yield [M$_{\odot}$]'

    row = int(np.floor(element/4))
    column = element-4*row

    #axArr[row,column].set_xlim([5e-3,5.0])
    axArr[row,column].set_xscale('log')
    #axArr[row,column].set_ylabel(yLabel)
    axArr[row,column].set_yscale('log')
    axArr[row,column].text(0.85, 0.95, symbol, transform=axArr[row,column].transAxes, fontsize=8, verticalalignment='top')
    axArr[row,column].plot(zArrBrMod,yieldArrNew[index,0,:],'bo')
    axArr[row,column].plot(zArrBrMod,yieldArrNew[index,1,:],'bs')
    axArr[row,column].plot(zArrBrMod,yieldArrNew[index,2,:],'b^')
    axArr[row,column].plot(zArrBrMod,yieldArrNew[index,3,:],'b*')

    axArr[row,column].plot(zArrBr,yieldArrSCHSim[index,0,:],'ro')
    axArr[row,column].plot(zArrBr,yieldArrSCHSim[index,1,:],'rs')
    axArr[row,column].plot(zArrBr,yieldArrSCHSim[index,2,:],'r^')
    axArr[row,column].plot(zArrBr,yieldArrSCHSim[index,3,:],'r*')

    axArr[row,column].plot(zArrWoWe,yieldArrCCWoWe[index,0,:],'o',color='grey')
    axArr[row,column].plot(zArrWoWe,yieldArrCCWoWe[index,1,:],'s',color='grey')
    axArr[row,column].plot(zArrWoWe,yieldArrCCWoWe[index,2,:],'^',color='grey')
    axArr[row,column].plot(zArrWoWe,yieldArrCCWoWe[index,3,:],'*',color='grey')

    axArr[row,column].plot(zAsplund*np.ones(nFilesCC),yieldArrCC[index,:],'o',color='k')
    
plt.savefig(plotFileName)

plotFileName = 'yieldsNewRatios.pdf'
ratiosList = ['Mn/Cr','Fe/Cr','Ni/Cr','Ni/Fe','Mn/Fe','Cr/Fe','Cr/Ca','Fe/Ca','Ca/S']
symArr = ['bo','bs','b^','b*']
symArrSCH = ['ro','rs','r^','r*']

plt.figure(1,figsize = [11.0, 8.5])
plt.clf()
f, axArr = plt.subplots(3, 3, sharex = True, squeeze = True)   

for ratio in range(len(ratiosList)) :
    [firstElement, secondElement] = ratiosList[ratio].split('/')
    firstElIndex = elemsToRead.index(firstElement)
    secondElIndex = elemsToRead.index(secondElement)
    row = int(np.floor(ratio/3))
    column = ratio-3*row
    #axArr[row,column].set_xlim([5e-3,5.0])
    axArr[row,column].set_xscale('log')
    axArr[row,column].set_yscale('log')
    axArr[row,column].text(0.75, 0.95, ratiosList[ratio], transform=axArr[row,column].transAxes, fontsize=10, verticalalignment='top')

    for model in range(nExpModels) :
        axArr[row,column].plot(zArrBrMod,(yieldArrNew[firstElIndex,model,:]/yieldArrNew[secondElIndex,model,:]),symArr[model],ms=8,mew=1)
        axArr[row,column].plot(zArrBr,(yieldArrSCHSim[firstElIndex,model,:]/yieldArrSCHSim[secondElIndex,model,:]),symArrSCH[model],ms=8,mew=1)

    axArr[row,column].plot(zAsplund*np.ones(nFilesCC),(yieldArrCC[firstElIndex,:]/yieldArrCC[secondElIndex,:]),'o',color='k',ms=8,mew=1)
    axArr[row,column].plot(zAsplund*np.ones(3),(yieldArrSCh[firstElIndex,:]/yieldArrSCh[secondElIndex,:]),'o',color = 'orange',ms=8,mew=1)

    ## for model in range(nModelsCCWoWe) :
    ##     axArr[row,column].plot(zArrWoWe,(yieldArrCCWoWe[firstElIndex,model,:]/yieldArrCCWoWe[secondElIndex,model,:]),symArr[model],color='grey',ms=8,mew=1)
    

plt.savefig(plotFileName)

plt.rcParams.update({'font.size': 12})
plotFileName = 'profiles_DDT.pdf'

elemsToPlot = ['O','Ne','Mg','Si','S','Ar','Ca','Cr','Mn','Fe','Ni']
colorArr = ['Blue','Green','SteelBlue','Cyan','Orange','Grey','Purple','Olive','OrangeRed','Red','SaddleBrown']

plt.figure(1,figsize = [11.0, 8.5])
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

axArr[0,0].set_xlim([0,1.5])
axArr[0,0].set_xscale('linear')
axArr[0,0].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,0].set_yscale('log')
axArr[0,0].set_ylim([1e-3,2.0])
axArr[0,0].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,0].plot(lagMassProfilesArrNew[:,0],profilesArrNew[:,index,0,2], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,0].plot(lagMassProfilesArrNew[:,0],profilesArrNew[:,index,0,0], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25,linestyle='dotted')
        
bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,0].text(0.05,0.05,'DDTa  ',transform=axArr[0,0].transAxes,bbox=bbox_props)

axArr[0,1].set_xlim([0,1.5])
axArr[0,1].set_xscale('linear')
axArr[0,1].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,1].set_yscale('log')
axArr[0,1].set_ylim([1e-3,2.0])
axArr[0,1].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,1].plot(lagMassProfilesArrNew[:,1],profilesArrNew[:,index,1,2], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,1].plot(lagMassProfilesArrNew[:,1],profilesArrNew[:,index,1,0], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25,linestyle='dotted')

bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,1].text(0.05,0.05,'DDTc  ',transform=axArr[0,1].transAxes,bbox=bbox_props)

axArr[1,0].set_xlim([0,1.5])
axArr[1,0].set_xscale('linear')
axArr[1,0].set_xlabel('Mass [M$_{\odot}$]')
axArr[1,0].set_yscale('log')
axArr[1,0].set_ylim([1e-3,2.0])
axArr[1,0].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[1,0].plot(lagMassProfilesArrNew[:,3],profilesArrNew[:,index,3,2], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[1,0].plot(lagMassProfilesArrNew[:,3],profilesArrNew[:,index,3,0], color = colorArr[elemsToPlot.index(symbol)], linewidth = 1.25,linestyle='dotted')

bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[1,0].text(0.05,0.05,'DDTf  ',transform=axArr[1,0].transAxes,bbox=bbox_props)

axArr[1,0].legend(loc=(1.2,0.1),ncol=2,fontsize=10)


axArr[1,1].axis('off')

plt.savefig(plotFileName)

plotFileName = 'profiles_DDTSimm.pdf'
plt.figure(1,figsize = [11.0, 8.5])
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

axArr[0,0].set_xlim([0,1.5])
axArr[0,0].set_xscale('linear')
axArr[0,0].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,0].set_yscale('log')
axArr[0,0].set_ylim([1e-3,2.0])
axArr[0,0].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,0].plot(lagMassProfilesArrSimmering[:,0],profilesArrSimmering[:,index,0,2], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,0].plot(lagMassProfilesArrSimmering[:,0],profilesArrSimmering[:,index,0,0], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25,linestyle='dotted')
        
bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,0].text(0.05,0.05,'DDTa  ',transform=axArr[0,0].transAxes,bbox=bbox_props)

axArr[0,1].set_xlim([0,1.5])
axArr[0,1].set_xscale('linear')
axArr[0,1].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,1].set_yscale('log')
axArr[0,1].set_ylim([1e-3,2.0])
axArr[0,1].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,1].plot(lagMassProfilesArrSimmering[:,1],profilesArrSimmering[:,index,1,2], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,1].plot(lagMassProfilesArrSimmering[:,1],profilesArrSimmering[:,index,1,0], color = colorArr[elemsToPlot.index(symbol)],linewidth = 1.25,linestyle='dotted')

bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,1].text(0.05,0.05,'DDTf  ',transform=axArr[0,1].transAxes,bbox=bbox_props)

axArr[0,1].legend(loc=(0.2,-0.9),ncol=2,fontsize=10)

axArr[1,0].axis('off')
axArr[1,1].axis('off')

plt.savefig(plotFileName)

plt.rcParams.update({'font.size': 12})
plotFileName = 'profiles_SCH.pdf'

plt.figure(1,figsize = [11.0, 8.5])
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

axArr[0,0].set_xlim([0,1.2])
axArr[0,0].set_xscale('linear')
axArr[0,0].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,0].set_yscale('log')
axArr[0,0].set_ylim([1e-3,2.0])
axArr[0,0].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,0].plot(lagMassProfilesArrSCHSim[:,0],profilesArrSCHSim[:,index,0,4], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,0].plot(lagMassProfilesArrSCHSim[:,0],profilesArrSCHSim[:,index,0,0], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25,linestyle='dotted')
        
bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,0].text(0.05,0.075,'0.88 M$_{\odot}$',transform=axArr[0,0].transAxes,bbox=bbox_props)

axArr[0,1].set_xlim([0,1.2])
axArr[0,1].set_xscale('linear')
axArr[0,1].set_xlabel('Mass [M$_{\odot}$]')
axArr[0,1].set_yscale('log')
axArr[0,1].set_ylim([1e-3,2.0])
axArr[0,1].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[0,1].plot(lagMassProfilesArrSCHSim[:,1],profilesArrSCHSim[:,index,1,4], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[0,1].plot(lagMassProfilesArrSCHSim[:,1],profilesArrSCHSim[:,index,1,0], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25,linestyle='dotted')

bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[0,1].text(0.05,0.075,'0.97  M$_{\odot}$',transform=axArr[0,1].transAxes,bbox=bbox_props)

axArr[1,0].set_xlim([0,1.2])
axArr[1,0].set_xscale('linear')
axArr[1,0].set_xlabel('Mass [M$_{\odot}$]')
axArr[1,0].set_yscale('log')
axArr[1,0].set_ylim([1e-3,2.0])
axArr[1,0].set_ylabel('Mass Fraction')

for element in range(nElemsToRead) :
    index = atNumsToRead.index(atNumsToRead[element])
    ele = ELEMENTS[atNumsToRead[element]]
    symbol = ele.symbol
    if symbol in elemsToPlot :
        axArr[1,0].plot(lagMassProfilesArrSCHSim[:,3],profilesArrSCHSim[:,index,3,4], color = colorArr[elemsToPlot.index(symbol)], label=symbol,linewidth = 1.25)
        axArr[1,0].plot(lagMassProfilesArrSCHSim[:,3],profilesArrSCHSim[:,index,3,0], color = colorArr[elemsToPlot.index(symbol)], linewidth = 1.25,linestyle='dotted')

bbox_props = dict(boxstyle='square,pad=0.4',fc='white')
axArr[1,0].text(0.05,0.075,'1.15 M$_{\odot}$',transform=axArr[1,0].transAxes,bbox=bbox_props)

axArr[1,0].legend(loc=(1.2,0.1),ncol=2,fontsize=10)


axArr[1,1].axis('off')

plt.savefig(plotFileName)

plt.rcParams.update({'font.size': 10})

plotFileName = 'metaratios.pdf'
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

metaRatiosArrOne = ['Ni/Fe','Ni/Fe','Cr/Fe','Fe/Ca'] 
metaRatiosArrTwo = ['Mn/Fe','Mn/Cr','Cr/Ca','Cr/Ca']
nMetaRatios = len(metaRatiosArrOne)

NiToFe3C397 = (58.0/56.0)*0.5*(0.076+0.37)
NiToFe3C397Min = (58.0/56.0)*0.076
NiToFe3C397Max = (58.0/56.0)*0.37
MnToCr3C397 = 1.057*0.5*(0.53+1.9)
MnToCr3C397Min = 1.057*0.53
MnToCr3C397Max = 1.057*1.9
MnToFe3C397 = (55.0/56.0)*0.5*(0.023+0.045)
MnToFe3C397Min = (55.0/56.0)*0.023
MnToFe3C397Max = (55.0/56.0)*0.045

for metaRatio in range(nMetaRatios) :

    [firstElement, secondElement] = metaRatiosArrOne[metaRatio].split('/')
    firstElIndex = elemsToRead.index(firstElement)
    secondElIndex = elemsToRead.index(secondElement)

    [firstElement2, secondElement2] = metaRatiosArrTwo[metaRatio].split('/')
    firstElIndex2 = elemsToRead.index(firstElement2)
    secondElIndex2 = elemsToRead.index(secondElement2)

    row = int(np.floor(metaRatio/2))
    column = metaRatio-2*row
    
    axArr[row,column].set_xscale('log')
    axArr[row,column].set_ylabel(metaRatiosArrTwo[metaRatio])
    axArr[row,column].set_yscale('log')
    axArr[row,column].set_xlabel(metaRatiosArrOne[metaRatio])

    symArr = ['s','o','^','*']

    for model in range(nExpModels) :
        for metallicity in range(nMetallicities) : 
            axArr[row,column].plot((yieldArrNew[firstElIndex,model,metallicity]/yieldArrNew[secondElIndex,model,metallicity])
                                   , (yieldArrNew[firstElIndex2,model,metallicity]/yieldArrNew[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = 'blue', label = 'CH SN Ia Models',ms=8, mew = 1, linestyle = 'None')
            axArr[row,column].plot((yieldArrSCHSim[firstElIndex,model,metallicity]/yieldArrSCHSim[secondElIndex,model,metallicity])
                                   , (yieldArrSCHSim[firstElIndex2,model,metallicity]/yieldArrSCHSim[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = 'red', label = 'SCH SN Ia Models',ms=8, mew = 1, linestyle = 'None')

    axArr[row,column].plot((yieldArrSCh[firstElIndex,:]/yieldArrSCh[secondElIndex,:])
                           ,(yieldArrSCh[firstElIndex2,:]/yieldArrSCh[secondElIndex2,:]),'o', color = 'orange', ms=8, mew = 1, label = 'SCH Models')
    
    for model in range(nModelsCCWoWe) :
        for metallicity in range(nMetallicitiesWoWe) : 
            axArr[row,column].plot((yieldArrCCWoWe[firstElIndex,model,metallicity]/yieldArrCCWoWe[secondElIndex,model,metallicity])
                                   , (yieldArrCCWoWe[firstElIndex2,model,metallicity]/yieldArrCCWoWe[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = 'grey',ms=8, mew = 1, linestyle = 'None')

    axArr[row,column].plot((yieldArrCC[firstElIndex,:]/yieldArrCC[secondElIndex,:])
                           ,(yieldArrCC[firstElIndex2,:]/yieldArrCC[secondElIndex2,:]),'o', color = 'k', ms=8, mew = 1)
            
    if metaRatio == 0 or metaRatio == 1  :
        axArr[row,column].axvline(NiToFe3C397,lw=1.5,color='k')
        axArr[row,column].axvline(NiToFe3C397Min,lw=1.5,ls=':',color='k')
        axArr[row,column].axvline(NiToFe3C397Max,lw=1.5,ls=':',color='k')

    if metaRatio == 0 :
        axArr[row,column].axhline(MnToFe3C397,lw=1.5,color='k')
        axArr[row,column].axhline(MnToFe3C397Min,lw=1.5,ls=':',color='k')
        axArr[row,column].axhline(MnToFe3C397Max,lw=1.5,ls=':',color='k')
    
    if metaRatio == 1 :
        axArr[row,column].axhline(MnToCr3C397,lw=1.5,color='k')
        axArr[row,column].axhline(MnToCr3C397Min,lw=1.5,ls=':',color='k')
        axArr[row,column].axhline(MnToCr3C397Max,lw=1.5,ls=':',color='k')
    
    #axArr[row,column].legend()
    
plt.savefig(plotFileName)


plotFileName = 'metaratios2.pdf'
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

metaRatiosArrOne = ['Ni/Fe','Ni/Fe','Ni/Fe','Ni/Fe'] 
metaRatiosArrTwo = ['Mn/Fe','Mn/Cr','Cr/Ca','Mn/Ca']
nMetaRatios = len(metaRatiosArrOne)

NiToFe3C397 = (58.0/56.0)*0.5*(0.076+0.37)
NiToFe3C397Min = (58.0/56.0)*0.076
NiToFe3C397Max = (58.0/56.0)*0.37
MnToCr3C397 = 1.057*0.5*(0.53+1.9)
MnToCr3C397Min = 1.057*0.53
MnToCr3C397Max = 1.057*1.9
MnToFe3C397 = (55.0/56.0)*0.5*(0.023+0.045)
MnToFe3C397Min = (55.0/56.0)*0.023
MnToFe3C397Max = (55.0/56.0)*0.045

for metaRatio in range(nMetaRatios) :

    [firstElement, secondElement] = metaRatiosArrOne[metaRatio].split('/')
    firstElIndex = elemsToRead.index(firstElement)
    secondElIndex = elemsToRead.index(secondElement)

    [firstElement2, secondElement2] = metaRatiosArrTwo[metaRatio].split('/')
    firstElIndex2 = elemsToRead.index(firstElement2)
    secondElIndex2 = elemsToRead.index(secondElement2)

    row = int(np.floor(metaRatio/2))
    column = metaRatio-2*row
    
    axArr[row,column].set_xscale('log')
    axArr[row,column].set_ylabel(metaRatiosArrTwo[metaRatio])
    axArr[row,column].set_yscale('log')
    axArr[row,column].set_xlabel(metaRatiosArrOne[metaRatio])

    if metaRatio == 0 :
        axArr[row,column].set_ylim(3e-4,1e-1) 
        axArr[row,column].set_xlim(5e-4,1.0)
    
    symArr = ['s','o','^','*']
    
    for metallicity in range(nMetallicities) :
         axArr[row,column].plot((yieldArrSCHSim[firstElIndex,:,metallicity]/yieldArrSCHSim[secondElIndex,:,metallicity])
                                   , (yieldArrSCHSim[firstElIndex2,:,metallicity]/yieldArrSCHSim[secondElIndex2,:,metallicity])
                                   , color = 'red', ls = '--', lw = 1.5)

         lastx = yieldArrSCHSim[firstElIndex,3,metallicity]/yieldArrSCHSim[secondElIndex,3,metallicity]
         lasty = yieldArrSCHSim[firstElIndex2,3,metallicity]/yieldArrSCHSim[secondElIndex2,3,metallicity]
         if metaRatio == 0 :
             metSun = zArrBr[metallicity]/zAsplund
             axArr[row,column].text(lastx/1.5,lasty/2.0,'Z/Z$_{\odot}$ = %2.1f' % metSun, color = 'red',clip_on=True)
         
         ## axArr[row,column].plot([0.5*lastx,lastx]
         ##                        , [lasty,lasty] 
         ##                        , color = 'red', lw = 1.5)

         ## xnew = np.linspace(0.5*yieldArrSCHSim[firstElIndex,0,metallicity]/yieldArrSCHSim[secondElIndex,0,metallicity]
         ##                    , 1.5*yieldArrSCHSim[firstElIndex,3,metallicity]/yieldArrSCHSim[secondElIndex,3,metallicity]
         ##                    , 300)           
         ## ysmooth = spline(yieldArrSCHSim[firstElIndex,:,metallicity]/yieldArrSCHSim[secondElIndex,:,metallicity]
         ##                  , yieldArrSCHSim[firstElIndex2,:,metallicity]/yieldArrSCHSim[secondElIndex2,:,metallicity]
         ##                  , xnew)
         ## axArr[row,column].plot(xnew,ysmooth,color='k')
         
    for model in range(nExpModels) :
        for metallicity in range(nMetallicities) : 
            axArr[row,column].plot((yieldArrNew[firstElIndex,model,metallicity]/yieldArrNew[secondElIndex,model,metallicity])
                                   , (yieldArrNew[firstElIndex2,model,metallicity]/yieldArrNew[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = 'blue', label = 'CH SN Ia Models',ms=8, mew = 1, linestyle = 'None')
            axArr[row,column].plot((yieldArrSCHSim[firstElIndex,model,metallicity]/yieldArrSCHSim[secondElIndex,model,metallicity])
                                   , (yieldArrSCHSim[firstElIndex2,model,metallicity]/yieldArrSCHSim[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = 'red', label = 'SCH SN Ia Models',ms=8, mew = 1, linestyle = 'None')

    for model in range(2) :
        for metallicity in range(nMetallicities) : 
            axArr[row,column].plot((yieldArrSimmering[firstElIndex,model,metallicity]/yieldArrSimmering[secondElIndex,model,metallicity])
                                   , (yieldArrSimmering[firstElIndex2,model,metallicity]/yieldArrSimmering[secondElIndex2,model,metallicity])
                                   , marker = symArr[model*3], color = 'green', label = 'CH SN Ia Models',ms=8, mew = 1, linestyle = 'None')
        

    ## axArr[row,column].plot((yieldArrSCh[firstElIndex,:]/yieldArrSCh[secondElIndex,:])
    ##                        ,(yieldArrSCh[firstElIndex2,:]/yieldArrSCh[secondElIndex2,:]),'o', color = 'orange', ms=8, mew = 1, label = 'SCH Models')

    if metaRatio == 0 :
        #axArr[row,column].errorbar([NiToFe3C397],[MnToFe3C397]
        #                           , xerr = [np.array([NiToFe3C397Min,NiToFe3C397Max])-NiToFe3C397]
        #                           , yerr = [np.array([MnToFe3C397Min,MnToFe3C397Max])-MnToFe3C397]
        #                           , capsize = 0, ecolor = 'k', marker = 'H', mfc='grey', mec='black', ms=10, mew = 1)
        axArr[row,column].axvline(NiToFe3C397Min,color='k',linestyle=':')
        axArr[row,column].axvline(NiToFe3C397Max,color='k',linestyle=':')
        axArr[row,column].axhline(MnToFe3C397Min,color='k',linestyle=':')
        axArr[row,column].axhline(MnToFe3C397Max,color='k',linestyle=':')

        
    if metaRatio == 1 :
        axArr[row,column].errorbar([NiToFe3C397],[MnToCr3C397]
                                   , xerr = [np.array([NiToFe3C397Min,NiToFe3C397Max])-NiToFe3C397]
                                   , yerr = [np.array([MnToCr3C397Min,MnToCr3C397Max])-MnToCr3C397]
                                   , capsize = 0, ecolor = 'k', marker = 'H', mfc='grey', mec='black', ms=10, mew = 1)
    
    ## if metaRatio == 1 :
    ##    axArr[row,column].errorbar()
    
    #axArr[row,column].legend()
    
plt.savefig(plotFileName)

plotFileName = 'ArCa.pdf'
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)   

ratiosList = ['Ar/S','Ca/S']
symArr = ['bo','bs','b^','b*']
symArrSCH = ['ro','rs','r^','r*']


for ratio in range(len(ratiosList)) :
    [firstElement, secondElement] = ratiosList[ratio].split('/')
    firstElIndex = elemsToRead.index(firstElement)
    secondElIndex = elemsToRead.index(secondElement)
    row = int(np.floor(ratio/3))
    column = ratio-3*row
    #axArr[row,column].set_xlim([5e-3,5.0])
    axArr[row,column].set_xscale('log')
    axArr[row,column].set_xlabel('Z')
    axArr[row,column].set_yscale('linear')
    axArr[row,column].text(0.75, 0.95, ratiosList[ratio], transform=axArr[row,column].transAxes, fontsize=10, verticalalignment='top')

    for model in range(nExpModels) :
        axArr[row,column].plot(zArrBrMod,(yieldArrNew[firstElIndex,model,:]/yieldArrNew[secondElIndex,model,:]),symArr[model],ms=8,mew=1)
        axArr[row,column].plot(zArrBr,(yieldArrSCHSim[firstElIndex,model,:]/yieldArrSCHSim[secondElIndex,model,:]),symArrSCH[model],ms=8,mew=1)


metaRatiosArrOne = ['Ca/S','Ar/S'] 
metaRatiosArrTwo = ['Ar/S','Ni/Fe']
nMetaRatios = len(metaRatiosArrOne)

colorArr = ['red','orange','yellow','green','blue']
colorArrMod = ['red','orange','blue','green','yellow']

for metaRatio in range(nMetaRatios) :

    [firstElement, secondElement] = metaRatiosArrOne[metaRatio].split('/')
    firstElIndex = elemsToRead.index(firstElement)
    secondElIndex = elemsToRead.index(secondElement)

    [firstElement2, secondElement2] = metaRatiosArrTwo[metaRatio].split('/')
    firstElIndex2 = elemsToRead.index(firstElement2)
    secondElIndex2 = elemsToRead.index(secondElement2)

    row = int(np.floor(metaRatio/2)) + 1
    column = metaRatio-2*row
    
    axArr[row,column].set_xscale('linear')
    axArr[row,column].set_ylabel(metaRatiosArrTwo[metaRatio])
    axArr[row,column].set_yscale('linear')
    if metaRatio == 1 : axArr[row,column].set_yscale('log')
    axArr[row,column].set_xlabel(metaRatiosArrOne[metaRatio])

    symArr = ['s','o','^','*']
         
    for model in range(nExpModels) :
        for metallicity in range(nMetallicities) : 
            axArr[row,column].plot((yieldArrNew[firstElIndex,model,metallicity]/yieldArrNew[secondElIndex,model,metallicity])
                                   , (yieldArrNew[firstElIndex2,model,metallicity]/yieldArrNew[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = colorArrMod[metallicity], label = 'CH SN Ia Models',ms=8, mew = 1, linestyle = 'None')
            axArr[row,column].plot((yieldArrSCHSim[firstElIndex,model,metallicity]/yieldArrSCHSim[secondElIndex,model,metallicity])
                                   , (yieldArrSCHSim[firstElIndex2,model,metallicity]/yieldArrSCHSim[secondElIndex2,model,metallicity])
                                   , marker = symArr[model], color = colorArr[metallicity], label = 'SCH SN Ia Models',ms=8, mew = 1, linestyle = 'None')
    
plt.savefig(plotFileName)

plotFileName = 'BR.pdf'
plt.clf()
f, axArr = plt.subplots(2, 2, squeeze = True)  


colorList = ['blue','red','orange']
symList = ['o','s']

axArr[0,0].set_xscale('log')
axArr[0,0].set_xlabel('Z/Z$_{\odot}$')
axArr[0,0].set_xlim(1e-2,10.0)
axArr[0,0].set_yscale('log')
axArr[0,0].set_ylabel('Mn/Cr')

for expModel in range(nExpModels) :
    axArr[0,0].plot(zArrSolMod, yieldArrNewFePeakBR[1,0,expModel,:]/yieldArrNewFePeakBR[0,0,expModel,:]
                    , marker = symList[0], color = colorList[0], ms=8, mew = 1, linestyle = 'None')
    axArr[0,0].plot(zArrSolMod, yieldArrNewFePeakBR[1,2,expModel,:]/yieldArrNewFePeakBR[0,2,expModel,:]
                    , marker = symList[0], color = colorList[2], ms=8, mew = 1, linestyle = 'None')
    axArr[0,0].plot(zArrSol, yieldArrNewFePeakBRSCH[1,2,expModel,:]/yieldArrNewFePeakBRSCH[0,2,expModel,:]
                    , marker = symList[1], color = colorList[2], ms=8, mew = 1, linestyle = 'None')

axArr[0,0].plot(zArrSol,
                ((yieldArrNewFePeakBR[1,0,0,2]/yieldArrNewFePeakBR[0,0,0,2])/(5.3*(zArrBr[2]**0.65)))
                *5.3*(zArrBr**0.65),linestyle = ':', color='k',linewidth = 1.5)

axArr[0,1].set_xscale('log')
axArr[0,1].set_xlabel('Z/Z$_{\odot}$')
axArr[0,1].set_xlim(1e-2,10.0)
axArr[0,1].set_yscale('log')
axArr[0,1].set_ylabel('Ni/Fe')

for expModel in range(nExpModels) :
    axArr[0,1].plot(zArrSolMod, yieldArrNewFePeakBR[3,0,expModel,:]/yieldArrNewFePeakBR[2,0,expModel,:]
                    , marker = symList[0], color = colorList[0], ms=8, mew = 1, linestyle = 'None')
    axArr[0,1].plot(zArrSolMod, yieldArrNewFePeakBR[3,1,expModel,:]/yieldArrNewFePeakBR[2,1,expModel,:]
                    , marker = symList[0], color = colorList[1], ms=8, mew = 1, linestyle = 'None')
    axArr[0,1].plot(zArrSolMod, yieldArrNewFePeakBR[3,2,expModel,:]/yieldArrNewFePeakBR[2,2,expModel,:]
                    , marker = symList[0], color = colorList[2], ms=8, mew = 1, linestyle = 'None')
    axArr[0,1].plot(zArrSol, yieldArrNewFePeakBRSCH[3,1,expModel,:]/yieldArrNewFePeakBRSCH[2,1,expModel,:]
                    , marker = symList[1], color = colorList[1], ms=8, mew = 1, linestyle = 'None')
    axArr[0,1].plot(zArrSol, yieldArrNewFePeakBRSCH[3,2,expModel,:]/yieldArrNewFePeakBRSCH[2,2,expModel,:]
                    , marker = symList[1], color = colorList[2], ms=8, mew = 1, linestyle = 'None')

    

#for expModel in range(nExpModelsSCHSim) :
#    for burnReg in range(1,3) :
#        axArr[0,1].plot(zArrSol, yieldArrNewFePeakBRSCH[1,burnReg,expModel,:]/yieldArrNewFePeakBRSCH[0,burnReg,expModel,:]
#                        , marker = symList[burnReg], color = colorList[burnReg], ms=8, mew = 1, linestyle = 'None')
                 
plt.savefig(plotFileName)

