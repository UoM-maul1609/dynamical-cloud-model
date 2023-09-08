#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 21:35:27 2020

@author: mccikpc2

run by importing

import batch_k_distribution
batch_k_distribution.driver()

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from tqdm import tqdm
import time

from multiprocessing import TimeoutError, Pool, Process, cpu_count, Lock, Manager
import subprocess
import itertools
import sys
import os
import scipy.interpolate as scint

import calculation


"""
==============================================================================
for all run variables                                                        =
==============================================================================
"""

h=6.62607015e-34
c=2.99792458e10 # cm/s
k=1.380649e-23
c2=h*c/k
Tref=296.
A=6.02e23
deltanur=0.00005 # cm-1 see Jacobson 2005, page 509
min_resample=5000000
deltanu=0.05
# table 2 of Jacobson, 2005, page 509
weights=[0.227979164257, 0.227979164257, 0.227979164257, 0.227979164257, \
         0.051055694388, 0.021463214949, 0.009022883764, 0.003793114481, \
         0.001594580828, 0.000670343073, 0.000281804364, 0.000118467248, \
         0.000049802241, 0.000020936278, 0.000008801366, 0.000003699991]



# table 3 of Jacobson, 2005, page 511 - surely this is g/cm-2!
mlpath=[0.21274, 0.026143, 2.9051e-5, 2.0783e-5,3.83e-6, 4.1468e-5,10.76,4.4894e-8] 
#,3.8387e-8,6.27e-8,2.1148e-8]


lineByLineFile='/models/mccikpc2/absorption_data/' + \
    '5f13dd32.par.txt'
    
lineByLinePartitionFiles='/models/mccikpc2/absorption_data/' + \
    'partitiondata/'

# lineByLineFile='/Volumes/Macintosh HD/Users/mccikpc2/' + \
#     'Dropbox (The University of Manchester)/absorption_data/' + \
#     '5f13dd32.par.txt'
#     
# lineByLinePartitionFiles='/Volumes/Macintosh HD/Users/mccikpc2/' + \
#     'Dropbox (The University of Manchester)/absorption_data/' + \
#     'partitiondata/'



"""
==============================================================================
define run variables                                                         =
==============================================================================
"""
nc=32
Pref=101325.
"""
page 511, Jacobson, 2005:
    31 pressures between 0.2 and 1050 hPa
    11 temperatures between 150 and 350K
    2 water vapour partial pressures between 0 0.03 atm
    1 partial pressure - for a given total pressure

"""
ntemp =11    # 11 page 511, Jacobson et al 2005
npress=31    # 31 page 511, Jacobson et al 2005
nh2o  =2    # 2 page 511, Jacobson et al 2005
tlow=150.
thigh=350.
plow=0.2e2
phigh=1050.e2
# tlow=270.
# plow=338e2
h2olow=0.0
h2ohigh=0.03

# tlow=270
#plow=1000.e2


Pcalc=np.linspace(plow,phigh,npress)
Tcalc=np.linspace(tlow,thigh,ntemp)
h2ocalc=np.linspace(h2olow,h2ohigh,nh2o)


gen_flag=2

if gen_flag == 1:
	lambdas=np.array([100.e-9,200.e-9,300.e-9,350.e-9,400.e-9,450.e-9,\
			 500.e-9,550.e-9,600.e-9,700.e-9,800.e-9,900.e-9,\
			 1000.e-9,1100.e-9,1200.e-9,1300.e-9,1400.e-9,10500.e-9,\
			16000.e-9,120000.e-9])
	# lambdas=np.logspace(-7,-3,50)
	lambdalow=10.e-9
	lambdahigh=150000.e-9

	lambda_low=np.zeros((len(lambdas),))
	lambda_high=np.zeros((len(lambdas),))
	lambda_low[0]=lambdalow
	lambda_low[1:]=[(lambdas[i+1]+lambdas[i])/2. for i in range(len(lambdas)-1)]
	lambda_high[-1]=lambdahigh
	lambda_high[0:-1]=[(lambdas[i+1]+lambdas[i])/2. for i in range(len(lambdas)-1)]

	lambda_low[10]=705e-9
	lambda_high[9]=705e-9

elif gen_flag == 2:
	lambdas=np.array([0.165, 0.175,0.1825,0.1875,0.1925,0.1975,0.2025,0.2075,0.2125,\
				0.2175,0.2225,0.2275,0.2325,0.2375,0.2425,0.2475,0.2525,0.2575, \
				0.2625,0.2675,0.2725,0.2775,0.2825,0.2875,0.2925,0.2975,0.302,0.305,\
				0.307,0.309,0.311,0.313,0.315,0.3175,0.3225,0.3275,0.3325,0.3375, \
				0.3425,0.3475,0.3525,0.3575,0.3625,0.3675,0.3725,0.3775,0.3825,0.3875,\
				0.3925,0.3975,0.4025,0.4075,0.4125,0.4175,0.4225,0.4275,0.435,0.445, \
				0.455,0.465,0.475,0.485,0.495,0.505,0.515,0.525,0.535,0.545,0.555, \
				0.565,0.575,0.585,0.595,0.605,0.615,0.625,0.635,0.645,0.655,0.665, \
				0.695,0.705,0.730,0.750,0.770,0.790,0.805,0.82,0.85,0.90,0.95,1.00, \
				1.05,1.10,1.15,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0, \
				6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5, \
				14.0,14.5,15.0,15.5,16.0,16.5,17.0,17.5,18.0,18.5,19.0,19.5,20.,22., \
				25.,30.,35.,40.,50.,70.,100.,200.,300,500,700.,1000])*1.e-6
	lambdalow=165.e-9
	lambdahigh=1000000.e-9
	
	lambda_low=lambdas[0:-1]
	lambda_high=lambdas[1:]
	
	# redefine lambdas - which are the mid-points in this version
	lambdas = np.zeros_like(lambda_low)
	lambdas[0] = lambdalow/2.
	for i in range(len(lambdas)):
    	lambdas[i] = 2.*lambda_low[i] - lambdas[i-1]
    
	#lambdas = (lambda_high-lambda_low)/2+lambda_low
	print('Lambda low / high (short / long) ' + str(lambda_low[81]*1e9))



# lambdahigh=1.5e-3


nbands=len(lambda_low)


# to work out main absorber, specify pressure and temp and wv
mtemp=0 # remember n-1 indexing
mpress=0
mh2o=0
# pint=scint.interp1d(Pcalc,np.mgrid[0:npress],kind='nearest')
# mpress=pint(1e5).astype('int') 
# tint=scint.interp1d(Tcalc,np.mgrid[0:ntemp],kind='nearest')
# mtemp=pint(290).astype('int') 

# These are the isotopologue files. 
# see here for details https://hitran.org/docs/iso-meta/
moleculeKey={'H2O':[1,[1,2,3,4,5,6,7],\
                    ['q1.txt','q2.txt','q3.txt','q4.txt',\
                     'q5.txt','q6.txt','q129.txt']],\
              'CO2':[2,[1,2,3,4,5,6,7,8,9,0,11,12],\
                    ['q7.txt','q8.txt','q9.txt','q10.txt',\
                     'q11.txt','q12.txt','q13.txt','q14.txt',\
                         'q121.txt','q15.txt','q120.txt','q122.txt']], \
              'O3':[3,[1,2,3,4,5], \
                    ['q16.txt','q17.txt','q18.txt','q19.txt',\
                     'q20.txt']], \
              'N2O':[4,[1,2,3,4,5],\
                     ['q21.txt','q22.txt','q23.txt','q24.txt',\
                     'q25.txt']], \
              'CO':[5,[1,2,3,4,5,6],\
                    ['q26.txt','q27.txt','q28.txt','q29.txt',\
                     'q30.txt','q31.txt']], \
              'CH4':[6,[1,2,3,4],\
                     ['q32.txt','q33.txt','q34.txt','q35.txt']], \
              'O2':[7,[1,2,3],\
                    ['q36.txt','q37.txt','q38.txt']],\
              'CH3Cl':[24,[1,2],['q73.txt','q74.txt']]}
# moleculeKey={'H2O':[1,[1,2,3,4,5,6,7],\
#                     ['q1.txt','q2.txt','q3.txt','q4.txt',\
#                      'q5.txt','q6.txt','q129.txt']]}              
molecularWeights=[18.01528,44.01,48.,44.013,28.01,16.04,32.,50.49]

moleculeID=[moleculeKey[i][0] for i in moleculeKey]
molecule=[i for i in moleculeKey]


def driver():

    """
    ==========================================================================
    Set up processes                                                         =
    ==========================================================================
    """
    p=[None]*nc
    log=[None]*nc
    submission=[None]*nc
    for i in range(nc):
        p[i] = subprocess.Popen('', shell=True)
    time.sleep(1)
    
    
    """
    ==========================================================================
    Start loop                                                               =
    ==========================================================================
    """
    # up this calls "doRun", creating temporary files for all kabs, and lambdas in bands
    start=time.time()
    
    myelements=[]
    for elements in \
        itertools.product(*[[i for i in range(nh2o)],\
                            [j for j in range(npress)],\
                            [k for k in range(ntemp)]]):
        myelements.append(elements)
        
    i = 0 
    iter1 = 0
    lr = len(myelements)
    while True:
        runsProcessing = False
        for j in range(nc):
            # submit a new run
            if p[j].poll() == 0 and i < lr:
                (ih2o,ipress,itemp)=myelements[i]
                # print("Processing: " + str(myelements[i]))
                runs=sys.executable + ' ' + \
                    os.path.abspath(__file__). \
                        replace('(','\('). \
                            replace(')','\)').replace(' ','\ ') + ' ' + \
                    str(ih2o).zfill(2) + ' ' + \
                    str(ipress).zfill(2) + ' ' + str(itemp).zfill(2) + ' 1'
                # note last 1 is for the flag
                submission[j]=runs
                print(runs)
                p[j]=subprocess.Popen(runs, shell=True)
                i += 1
            
            # resubmit if it didn't work
            if (p[j].poll() != 0) and (p[j].poll() != None):
                print('resubmitting')
                p[j].terminate()
                p[j]=subprocess.Popen(submission[j], shell=True)

            # check if runs are still going
            runsProcessing = runsProcessing or (p[j].poll() != 0)
        
        if not(runsProcessing) and (i>=lr):
            break
        iter1 += 1
        time.sleep(5)
        
    
    
    print("Heavy processing took ", time.time()-start," seconds")
    
    
    
    
    
    """
    ==========================================================================
    Now heavy processing done, use a single processor to 
    find main absorber at specific P and T (and h2o setting?)
    ==========================================================================
    """
    # just reads the temporary files and finds the main absorber for a typical 
    # mlpath 
    import netCDF4 as nc4
    nus=[None]*nbands
    kabs=[[None for x in range(len(moleculeKey))] for y in range(nbands)]
    
    fileName='/tmp/output_' + str(mh2o).zfill(2) + '_' + \
            str(mpress).zfill(2) + '_' + str(mtemp).zfill(2) + '.nc'
    f = nc4.Dataset(fileName,'r')
    absgrp = f.groups['Absorption_data']
    for i in range(nbands):
        nus[i]=absgrp.variables['nus' + str(i).zfill(3)][:]
        for j in range(len(moleculeKey)):
            kabs[i][j]=absgrp.variables['kabs' + str(i).zfill(3)][j,:]
            
    f.close()
    
    sums1=np.zeros((nbands,len(moleculeKey)))
    for i in range(nbands):
        for j in range(len(moleculeKey)):
            sums1[i,j]=sums1[i,j]+np.sum(kabs[i][j]*mlpath[j])
            
    # identify the strongest absorber for each wave length interval
    max1=np.zeros((nbands))
    mainAbsorber=[None]*nbands
    for i in range(nbands):
        ind,=np.where(sums1[i]==np.max(sums1[i],axis=0)) 
        max1[i]=ind[0]
        mainAbsorber[i]=molecule[ind[0]]
    
    print(mainAbsorber)
    
    import csv
    with open('/tmp/mainAbsorber.csv','w') as result_file:
        wr = csv.writer(result_file, dialect='excel')
        wr.writerow(mainAbsorber)
    
    
    
    
    
    """ 
    ==========================================================================
    Now main absorber found, calculate the coefficients on multiple processors
    and output
    ==========================================================================
    """
    # calls "procRun" to calculate absorption coefficient for each band, molecule and weight
    p=[None]*nc
    log=[None]*nc
    submission=[None]*nc
    for i in range(nc):
        p[i] = subprocess.Popen('', shell=True)
    time.sleep(1)
    i=0
    iter1=0
    while True:
        runsProcessing = False
        for j in range(nc):
            # submit a new run
            if p[j].poll() == 0 and i < lr:
                (ih2o,ipress,itemp)=myelements[i]
                runs=sys.executable + ' ' + \
                    os.path.abspath(__file__). \
                        replace('(','\('). \
                            replace(')','\)').replace(' ','\ ') + ' ' + \
                    str(ih2o).zfill(2) + ' ' + \
                    str(ipress).zfill(2) + ' ' + str(itemp).zfill(2) + ' 2'
                # note last 2 is for the flag
                submission[j]=runs
                print(runs)
                p[j]=subprocess.Popen(runs, shell=True)
                i += 1
            
            # resubmit if it didn't work
            if (p[j].poll() != 0) and (p[j].poll() != None):
                print('resubmitting')
                p[j].terminate()
                p[j]=subprocess.Popen(submission[j], shell=True)
            # check if runs are still going
            runsProcessing = runsProcessing or (p[j].poll() != 0)
        
        if not(runsProcessing) and (i>=lr):
            break
        iter1 += 1
        time.sleep(5)
    
    

    """ 
    finally, use a single processor to output to a single file
    """    
    
    f1 = nc4.Dataset('/tmp/bli_data.nc','w',format='NETCDF3_CLASSIC')
    # absgrp = f1.createGroup('Absorption_data')
    f1.createDimension('nh2o', nh2o)
    f1.createDimension('npress', npress)
    f1.createDimension('ntemp', ntemp)
    f1.createDimension('nbands', nbands)
    f1.createDimension('molecule', len(moleculeKey))
    f1.createDimension('weights', len(weights))
    bliout=f1.createVariable('bli','f4', \
            ('nh2o','npress','ntemp','weights','molecule' ,'nbands'))
    cumout=f1.createVariable('cumulat','f4', ('weights'))
    weightout=f1.createVariable('probs','f4', ('weights'))
    lam_lowout=f1.createVariable('lambda_low','f4', ('nbands'))
    lam_highout=f1.createVariable('lambda_high','f4', ('nbands'))
    h2oout=f1.createVariable('h2o','f4', 'nh2o')
    pressout=f1.createVariable('press','f4', 'npress')
    tempout=f1.createVariable('temp','f4', 'ntemp')
    moleculeIDnc=f1.createVariable('moleculeID','i4', 'molecule')
    molecularWeightsnc=f1.createVariable('molecularWeights','f4', 'molecule')
    
    h2oout[:]=h2ocalc[:]
    pressout[:]=Pcalc[:]
    tempout[:]=Tcalc[:]
    moleculeIDnc[:]=moleculeID[:]
    molecularWeightsnc[:]=molecularWeights[:]
    for i in range(nh2o):
        for j in range(npress):
            for k in range(ntemp):
                fileName='/tmp/output_' + str(i).zfill(2) + '_' + \
                    str(j).zfill(2) + '_' + str(k).zfill(2) + '_bli.nc'
    
                f2 = nc4.Dataset(fileName,'r')
                absgrp = f2.groups['Absorption_summary_data']
            
                # create variables
                # ('weights','molecule' ,'nbands')
                bli=absgrp.variables['bli'][:,:,:]
                # weights
                cumulat=absgrp.variables['cumulat'][:]
                
                    
                # put the data in
                lam_lowout[:]=lambda_low
                lam_highout[:]=lambda_high
                cumout[:] = cumulat
                weightout[:] = weights
                bliout[i,j,k,:,:,:]=bli
                
                f2.close()
                
    f1.close()


"""
==============================================================================
calculate average in bin, equation 7, Jacobson 2005                          =
==============================================================================
"""
def calculate_average_in_bin(cumulat,nus,kabs,bli,indStore):
    
    # resample if needed
    if len(kabs)<min_resample:
        resampled=True
        f1=interp1d(nus,kabs)
        nus2=np.linspace(nus[0],nus[-1],min_resample)
        kabs2=f1(nus2)
    else:
        resampled=False
        kabs2=kabs
        nus2=nus
        
    X=np.percentile(kabs2,cumulat[:]*100)
#     print(X)
    for i in range(cumulat.size-1):
        ind,=np.where((kabs2>=X[i]) & (kabs2<X[i+1]))
        indStore[i]=ind
        if len(ind):
            bli[i]=np.mean(kabs2[ind])
#             print(bli[i])
    # sorted1=np.argsort(kabs2) 

    return (bli,indStore,resampled)

"""
==============================================================================
end calculate average in bin, equation 7, Jacobson 2005                      =
==============================================================================
"""





"""
==============================================================================
calculate average in bin, equation 8, Jacobson 2005                          =
==============================================================================
"""
def calculate_2or_average_in_bin(cumulat,nus,kabs,bli,indStore,imol,resampled):
    
    # resample if needed
    if resampled:
        f1=interp1d(nus,kabs)
        nus2=np.linspace(nus[0],nus[-1],min_resample)
        kabs2=f1(nus2)    
    else:
        kabs2=kabs
        nus2=nus
    
    for i in range(cumulat.size-1):
        ind=indStore[i]

        if len(ind):
            bli[i]=-1./mlpath[imol] * \
                np.log(np.mean(np.exp(-kabs2[ind]*mlpath[imol])) )


    return (bli)

"""
==============================================================================
end calculate average in bin, equation 8, Jacobson 2005                      =
==============================================================================
"""



"""
==============================================================================
calculate absorption coefficients from line-by-line data                     =
==============================================================================
"""
def main_loop(MN,MR,ipress,itemp,ih2o,nbands,moleculeID,moleculeKey, \
              IN,TF,LI,EAC,ABW,SBW,LSE,TD,PS,Q1,Q2,MOLW):

    len1=len(MN)
    len2=len(moleculeKey)
    
    nus=[None]*nbands
    kabs=[[None for x in range(len(moleculeKey))] for y in range(nbands)]

    ind,=np.where(MN[:,0]==1) # find the water molecule
    MR[ind,0]=h2ocalc[ih2o] #/Pcalc[ipress]
    Patm1=Pcalc[ipress]/Pref
    for ilam in range(len(lambda_high)):
        print("*********" + str(ilam) + ' ' + str(lambda_low[ilam]) + ' to ' +\
              str(lambda_high[ilam]) + "*********")
        
        for k in range(len(moleculeID)):
            print("MoleculeID is: " + str(moleculeID[k]) + " and temp " \
                + str(Tcalc[itemp])+ " and press " \
                + str(Pcalc[ipress])+ " and h2o " \
                + str(h2ocalc[ih2o]))
                
            (nus[ilam],kabs[ilam][k])= \
                    absorptionCoefficient(len1,len2,moleculeKey, \
                                   MN,IN,TF,LI,EAC,ABW,SBW,LSE,TD,PS,\
                                   Q1,Q2[:,itemp],Tcalc[itemp],Patm1, \
                                   MR,MOLW, \
                             lambda_low[ilam],lambda_high[ilam],moleculeID[k])


    return (nus,kabs,itemp,ipress,ih2o)
"""
==============================================================================
end calculate absorption coefficients from line-by-line data                 =
==============================================================================
"""



"""
==============================================================================
calculate the absorption coefficient at this wavenumber                      =
==============================================================================
"""
def absorptionCoefficient(len1,len2,moleculeKey, \
                           MN,IN,TF,LI,EAC,ABW,SBW,LSE,TD,PS,Q1,Q2,T,Patm, \
                           mr,MOLW,lambda_low,lambda_high,mol_choice):

    # 1. loop through all molecules
    # 2. loop through all isotopologues - find parameters equal to this
    #    including partition functions
    # 3. calculate the sum of all linestrengths at this wavenumber
    gamma_air=ABW
    gamma_self=SBW
    nq=TD
    nuq=TF
    delq=PS
    
    nulow=1./lambda_high/100.
    nuhigh=1./lambda_low/100. # cm-1


    """
    Relevant lines for this band
    """
    ind1,=np.where((LSE[:,0] != 0.) & (TF[:,0]>=nulow*0.9) & \
                   (TF[:,0]<nuhigh*1.1)  &\
                   (MN[:,0]==mol_choice) ) #& (IN[:,0]==1))
    """
    --------------------------------------------------------------------------
    """
    
    
    
    # loop over all molecule  
    """
        loop over the selected molecules / isotopologues++++++++++++++
    """
    
    len3=np.floor((nuhigh-nulow) / deltanu)
    
    # kabs=np.zeros((len3,1),dtype=float)
    print("Calculating absorption coefficient with " + str(len(ind1)) + " lines")
    
    # Calculate gamma_q (equation 9.28, Jacobson, page 296):
#     gamma_q=(Tref/T)**nq[ind1]*(gamma_air[ind1]*(Patm-PP[ind1])+ \
#                                 gamma_self[ind1]*PP[ind1])
    # Line intensity (see https://hitran.org/docs/definitions-and-units/)
#     S=LI[ind1]*Q1[ind1] / \
#         Q2[ind1] \
#         * np.exp(-c2*LSE[ind1]/T)*(1.-np.exp(-c2*nuq[ind1]/T)) / \
#         ( np.exp(-c2*LSE[ind1]/Tref)*(1.-np.exp(-c2*nuq[ind1]/Tref)) )
        
    
    # absorption coeffcient (see 9.27, Jacobson, page 296):
    # kabs=np.sum(A/(np.pi*MOLW[ind1]) * \
    #     (S*gamma_q / \
    #      (gamma_q**2+(nus-(nuq[ind1]-delq[ind1]/Patm))**2)),axis=0)
    if len(ind1)>0:
        nus,kabs=calculation.absorption_calc(len3,A,T,Tref,c2,nq[ind1,0],\
                gamma_air[ind1,0],gamma_self[ind1,0],mr[ind1,0],LI[ind1,0],\
                Q1[ind1,0],Q2[ind1],LSE[ind1,0],MOLW[ind1,0],nuq[ind1,0],\
                delq[ind1,0],Patm,nulow,deltanu,len1=len(ind1))
    else:
        test=np.zeros((1))
        nus=np.mgrid[nulow+deltanu:nuhigh:deltanu]
        kabs=nus.copy()
        kabs[:]=0.
    print("Done")
    """
        --------------------------------------------------------------
    """            
            
            
    return (nus,kabs)
"""
==============================================================================
end calculate the absorption coefficient at this wavenumber                  =
==============================================================================
"""
    



   
"""
==============================================================================
populate partition data                                                      =
==============================================================================
"""
def partitionData(moleculeKey):
    len1=len(moleculeKey)
    for i in moleculeKey:
        moleculeKey[i].append(list())
        moleculeKey[i].append(list())
        len2=len(moleculeKey[i][2])
        for j in range(len2):
            # read the file
            tmp=[]
            with open(lineByLinePartitionFiles + \
                      moleculeKey[i][2][j],'r') as fp:
                for k,line in enumerate(fp):
                    tmp.append([float(item) for item in line.split()])
            moleculeKey[i][3].append(np.array(tmp))
            
        
        for j in range(len2):
            moleculeKey[i][4].append(interp1d(moleculeKey[i][3][j][:,0], \
                                          moleculeKey[i][3][j][:,1]))
            
    return moleculeKey
    
"""
==============================================================================
end populate partition data                                                  =
==============================================================================
"""




"""
==============================================================================
populate partition data at specific temp                                     =
==============================================================================
"""
def partitionDataTemp(moleculeKey,Tcalc,Tref):
    len1=len(moleculeKey)
    for i in moleculeKey:
        moleculeKey[i].append(list()) # at Tref
        moleculeKey[i].append(list()) # at T
        len2=len(moleculeKey[i][2])
        for j in range(len2):
            moleculeKey[i][5].append(moleculeKey[i][4][j](Tref))
            moleculeKey[i][6].append(moleculeKey[i][4][0]([t for t in Tcalc]) )
            
    return moleculeKey
    
"""
==============================================================================
end populate partition data at specific temp                                 =
==============================================================================
"""



def procRun(ih2o,ipress,itemp,fileName):
    import csv
    import netCDF4 as nc4
    from batch_k_distribution import moleculeKey, moleculeID, molecule, \
        molecularWeights, nbands, weights
        
    # read the main absorbers
    with open('/tmp/mainAbsorber.csv', newline='') as f:
        reader = csv.reader(f)
        mainAbsorber = list(reader)[0]

    
    bli=np.zeros((len(weights),len(molecule),nbands))
    sorted1=[[None for x in range(len(molecule))] for y in range(nbands)]    
    cumulat=np.insert(np.cumsum(weights),0,0,axis=0)  
    nus=[None]*nbands
    kabs=[[None for x in range(len(moleculeKey))] for y in range(nbands)]


    
    fileName='/tmp/output_' + str(ih2o).zfill(2) + '_' + \
            str(ipress).zfill(2) + '_' + str(itemp).zfill(2) + '.nc'
    f = nc4.Dataset(fileName,'r')
    absgrp = f.groups['Absorption_data']
    for i in range(nbands):
        nus[i]=absgrp.variables['nus' + str(i).zfill(3)][:]
        for j in range(len(moleculeKey)):
            kabs[i][j]=absgrp.variables['kabs' + str(i).zfill(3)][j,:]
            
    f.close()
    
    

    indStore=[None]*(cumulat.size-1)
    # loop through all bands
    for j in range(len(mainAbsorber)):
        # Loop through all molecules:
        for i in range(len(molecule)):
            # is this molecule a main absorber?
            if molecule[i] == mainAbsorber[j]:
                # for this absorber sort kabs and calculate 
                # the average absorption coefficient
                # in each probability interval using 
                # equation 7, Jacobson, 2005, 510
                # (bli[:,i,j],nus2,kabs2,sorted1[j][i],indStore,resampled) = \
                #     calculate_average_in_bin(\
                #                 cumulat,nus[j],kabs[j][i],bli[:,i,j])
                (bli[:,i,j],indStore,resampled) = \
                    calculate_average_in_bin(\
                                cumulat,nus[j],kabs[j][i],bli[:,i,j],indStore)
                mAbs=i # store the main absorber
        
        # Loop through all molecules again:
        for i in range(len(molecule)):
            # is this molecule a main absorber?
            if molecule[i] != mainAbsorber[j]:
                # for this absorber sort kabs and calculate 
                # the average absorption coefficient
                # in each probability interval using 
                # equation 8, Jacobson, 2005, 510
                (bli[:,i,j]) = calculate_2or_average_in_bin( \
                                    cumulat,nus[j],kabs[j][i],bli[:,i,j], \
                                                indStore,i,resampled)
    

    """
        Now write the processed results to a file
    """
    fileName1=fileName.replace('.nc','_bli.nc')
    print("Creating bli file for: " + fileName1)
    # output to a file here
    from datetime import datetime
    today = datetime.today()
    
    f = nc4.Dataset(fileName1,'w',format='NETCDF4')
    absgrp = f.createGroup('Absorption_summary_data')
    absgrp.createDimension('nbands', nbands)
    absgrp.createDimension('molecule', len(moleculeKey))
    absgrp.createDimension('weights', len(weights))

    # create variables
    bliout=absgrp.createVariable('bli', \
                        'f4', ('weights','molecule' ,'nbands'))
    cumout=absgrp.createVariable('cumulat','f4', ('weights'))
    lam_lowout=absgrp.createVariable('lambda_low','f4', ('nbands'))
    lam_highout=absgrp.createVariable('lambda_high','f4', ('nbands'))
        
    # put the data in
    lam_lowout[:]=lambda_low
    lam_highout[:]=lambda_high
    cumout[:] = cumulat[1:]
    bliout[:,:,:]=bli

            
    f.description = "Temporary absorption dataset"
    f.history = "Created " + today.strftime("%d/%m/%y")
    bliout.units = "cm**2 g**-1"
    cumout.units = ""
    lam_lowout.units = "m"
    lam_highout.units = "m"
    f.close()
    print("Creating bli file done: " + fileName1)    


# create temporary absorption dataset
def doRun(ih2o,ipress,itemp,fileName):
    
    import netCDF4 as nc4
    from batch_k_distribution import moleculeKey, moleculeID, molecule, \
        molecularWeights
    print("Running for: " + fileName)
    print("Reading for: " + fileName)
    str1=list()
    with open(lineByLineFile,'r') as fp:
        for line in fp:
            str1.append(line)
    print("Reading done: " + fileName)



    print("Partition function interpolations for: " + fileName)
    moleculeKey=partitionData(moleculeKey)
    mixing_ratios=[0., 416e-6, 30e-9, 319e-9, 0.2e-6, 1.82e-6, 0.21, 550e-12]
    print("Partition function interpolations done: " + fileName)




    print("Declaring arrays for: " + fileName)
    len1=len(str1)
    MN =np.zeros((len1,1),dtype=int) # molecule number
    IN =np.zeros((len1,1),dtype=int) # isotope number
    TF =np.zeros((len1,1),dtype=float) # transition freq
    LI =np.zeros((len1,1),dtype=float) # line intensity
    EAC=np.zeros((len1,1),dtype=float) # einstein a coefficient
    ABW=np.zeros((len1,1),dtype=float) # air broadened width
    SBW=np.zeros((len1,1),dtype=float) # self "
    LSE=np.zeros((len1,1),dtype=float) # lower state energy
    TD =np.zeros((len1,1),dtype=float) # temperature dependence (of air width)
    PS =np.zeros((len1,1),dtype=float) # pressure shift 
    MOLW =np.zeros((len1,1),dtype=float) # molecular weight
    MR =np.zeros((len1,1),dtype=float) # mixing ratio - to determine main absorber
    # UVQ=np.zeros((len1,1),dtype=float)
    # LVQ=np.zeros((len1,1),dtype=float)
    # ULQ=np.zeros((len1,1),dtype=float)
    # LLQ=np.zeros((len1,1),dtype=float)
    #QS =[None]*len1
    Q1=np.zeros((len1,1),dtype=float)
    Q2=np.zeros((len1,ntemp),dtype=float)
    print("Declaring arrays done: " + fileName)


    print("Extracting information from string for: " + fileName)
    for i in range(len1):
        MN[i]  =int(str1[i][0:2])
        IN[i]  =int(str1[i][2:3])
        TF[i]  =float(str1[i][3:15])
        LI[i]  =float(str1[i][15:25])
        EAC[i] =float(str1[i][25:35])
        ABW[i] =float(str1[i][35:40])
        SBW[i] =float(str1[i][40:45])
        LSE[i] =float(str1[i][45:55])
        TD[i]  =float(str1[i][55:59])
        PS[i]  =float(str1[i][59:67])
        # UVQ[i] =float(str1[i][67:82])
        # LVQ[i] =float(str1[i][82:95])
        # ULQ[i] =float(str1[i][95:110])
        # LLQ[i] =float(str1[i][110:125])
    del str1
    print("Extracting information from string done: " + fileName)





    print("Assigning interpolations of partition functions for: " + fileName)
    len2=len(moleculeKey)
    for i in range(len1):
        # molecule id
        ind1,=np.where(moleculeID==MN[i])
        if len(ind1):
            # isotope id
            ind2,=np.where(moleculeKey[molecule[ind1[0]]][1]==IN[i]) 
            if len(ind2):
                # calculate the value of Q
                #QS[i]=moleculeKey[molecule[ind1[0]]][4][ind2[0]]
                MOLW[i]=molecularWeights[ind1[0]]
                MR[i]=mixing_ratios[ind1[0]]
    # calculate all partition functions you need
    moleculeKey=partitionDataTemp(moleculeKey,Tcalc,Tref)
    print("Assigning interpolations of partition functions done: " + fileName)





    
    print("Assigning partition functions to arrays for: " + fileName)
    for i in moleculeKey: # loop over molecules
        for j in range(len(moleculeKey[i][1])): # loop over isotopologues
            ind,=np.where((MN[:,0]==moleculeKey[i][0]) & \
            (IN[:,0]==moleculeKey[i][1][j]))
            # assign functions to Q1s
            Q1[ind,0]=moleculeKey[i][5][j]
    
            # assign all tcalcs to Q2s
            for k in [itemp]:
                Q2[ind,k]=moleculeKey[i][6][j][k]
    print("Assigning partition functions to arrays done: " + fileName)
    



    print("Calculating absorption lines for: " + fileName)
    (nus,kabs,itemp,ipress,ih2o)=main_loop(MN.copy(),MR.copy(),\
                ipress,itemp,ih2o,nbands,moleculeID,moleculeKey,\
                IN,TF,LI,EAC,ABW,SBW,LSE,TD,PS,Q1,Q2,MOLW)
    print("Calculating absorption lines done: " + fileName)
    
    
    
    
    
    print("Creating temporary file for: " + fileName)
    # output to a file here
    from datetime import datetime
    today = datetime.today()
    
    f = nc4.Dataset(fileName,'w',format='NETCDF4')
    absgrp = f.createGroup('Absorption_data')
    for i in range(nbands):
        absgrp.createDimension('lambda' + str(i).zfill(3), len(kabs[i][0]))
    absgrp.createDimension('molecules', len(moleculeKey))

    # create variables
    kabsout=[]
    nusout=[]
    for i in range(nbands):
        kabsout.append(  \
            absgrp.createVariable('kabs' + str(i).zfill(3), \
                        'f4', ('molecules','lambda' + str(i).zfill(3))))
        nusout.append(  \
            absgrp.createVariable('nus' + str(i).zfill(3), \
                        'f4', 'lambda' + str(i).zfill(3)))
        
    # put the data in
    for i in range(nbands):
        nusout[i][:] = nus[i][:]
        for j in range(len(moleculeKey)):
            kabsout[i][j,:] = kabs[i][j][:]
            
    f.description = "Temporary absorption dataset"
    f.history = "Created " + today.strftime("%d/%m/%y")
    for i in range(nbands):
        kabsout[i].units = "cm**2 g**-1"
        nusout[i].units = "cm**-1"
    f.close()
    print("Creating temporary file done: " + fileName)
    
    
    
    

if __name__ == "__main__":
    """
        Call the mult_job function using commandline arguemnts for testing
        e.g. 
        python3 batch_k_distribution.py 4 10 10 1
    """
    import sys
    ih2o=sys.argv[1]
    ipress=sys.argv[2]
    itemp=sys.argv[3]
    flag=int(sys.argv[4])
    
    fileName='/tmp/output_' + ih2o + '_' + ipress + '_' + itemp + '.nc'

    if flag==1:
        # creates a temporary absorption dataset, with all the kabs and nus for each line
        doRun(int(ih2o),int(ipress),int(itemp),fileName)
    elif flag==2:
        procRun(int(ih2o),int(ipress),int(itemp),fileName)

