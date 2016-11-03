import numpy as np
from const import *
from random import randint

##get the evolution of re
def v_rz(lr_9,lm_sz):
     lr_z = lr_9+0.2*(lm_sz-9)
     return lr_z

## get the initial re for galaxy at logmass=9 
## From van der Wel et al. 2013
def get_lr9(lr9,elr9):
     temp = np.random.normal(lr9,elr9,15000)
     ran = randint(0,14999)
     lr_9 = temp[ran]
     return lr_9

## if rand=True will use the scatter from the paper to generate re randomly
## if rand= False will use the re at the fitting line
def v_lr9(z,rand=False):
     if z>=2.5:
          lr9 = 0.51-0.18*1.699
          elr9 = 0.19
     elif z>=2.0 and z<2.5:
          lr9 = 0.55-0.22*1.699
          elr9 = 0.19
     elif z>=1.5 and z<2.0:
          lr9 = 0.65-0.23*1.699
          elr9 = 0.18
     elif z>=1.0 and z<1.5:
          lr9 = 0.70-0.22*1.699
          elr9 = 0.17
     elif z>=0.5 and z<1.0:
          lr9 = 0.78-0.22*1.699
          elr9 = 0.16
     elif z<0.5:
          lr9 = 0.86-0.25*1.699
          elr9 = 0.16
     if rand==True:
          lr9 = get_lr9(lr9,elr9)
     else:
          lr9 = lr9
     return lr9

##The parameters of sig1, re and mstar
A = 0.25  ##fitting slope of sig1,re and mstar relation
p = 1.0   ##slope of re
q = 1.0   ##slope of mstar

## the log sig1 value from the fitting, unit of sig1 is Msun/kpc^2
def v_sig1(lr_z,lm_sz):
     lsig1 = np.log10(A)+q*lm_sz-p*lr_z
     return lsig1

## log sig1(mass density) vs log Sig1(velocity dispersion) relation fitting parameters
##From Jerome Fang et al. 2013 plot
a = 9.18  ##intercept
b = 1.99  ##slope

## black hole mass vs sig_e (velocity) relation fitting parameters
##From Kormendy and Ho
c = 0.31  ##intercept
d = 4.38  ##slope

##dispersion and Re relation
##From Wache et al. 2012
f = 0.066  ##slope of Re

## the mass of black hole, unit: Msun
def v_mbh(lsig1,lr_z):
     lm_bh = (np.log10(c)+2*d+9-d*np.log10(200.0)-d*a/b)+\
             d*(lsig1/b-f*lr_z)
     return lm_bh

##get the mheat old, unit: Msun/yr
eta = 0.1 ##coefficient before energy input rate Hernques et al 2015 Eq.S25
k_AGN = 6.0e-6 ##Msun/yr  ##coefficient before BH accretion rate Hernques et al. 2015 Eq.S24

##get the black hole accretion rate Msun/yr
def v_mbh_acrate(lm_bh,lVvir):
     mbh_acrate = np.log10(k_AGN)+lm_bh-8+np.log10(fhot)+1+3*(lVvir-np.log10(200))
     return mbh_acrate

def v_mheat_old(lm_bh,lVvir):
     mbh_acrate = v_mbh_acrate(lm_bh,lVvir)
     lm_heat = np.log10(2*eta)+mbh_acrate+2*(np.log10(3.0e5)-lVvir)
     return lm_heat+4

##get the mheat new,from Sandy's new formular, unit: Msun/yr
#def v_mheat_new(lm_sz,lVvir,lre):
#     lm_heat = np.log10(2.36)-14+lm_sz-lre-2*lVvir+2*np.log10(3.0e5)
#     return lm_heat
def v_mheat_new(lm_sz,lVvir,lre):
     lsig1 = v_sig1(lre,lm_sz)
     lm_heat = np.log10(2.36)-14+0.5*lm_sz-lre-2*lVvir+2*np.log10(3.0e5)+(lsig1-4.5)
     return lm_heat

##mheat old from P9 is the same as old
def v_mheat_p9(lm_bh,lm_hz,z):
     lm_heat = np.log10(2.37)-28+np.log10(fhot)+0.33333*lm_hz+lm_bh-np.log10(1+z)+2*np.log10(3e10)
     return lm_heat
 
def v_mheat_p9_2(lm_sz,lr9,z):
     lm_heat = np.log10(1.18)-13+np.log10(fhot)+1.75*lm_sz-2.26*lr9-\
               (11.0*np.log10(1+z)/12.0)
     return lm_heat
