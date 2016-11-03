from scipy.misc import derivative
import numpy as np
from const import *

##calculation the halo and stellar mass
def v_mhalo_p3(lm_sz,z):
     lm_hz = 12.0+(lm_sz+0.42*np.log10(1+z)-10.4)/1.7
     return lm_hz

def v_mstar_p3(lm_hz,z):
     lm_sz = 10.4+1.7*(lm_hz-12)-0.42*np.log10(1+z)
     return lm_sz

def v_mhalo_p17(lm_h0,z):
     lm_hz = -0.35*z+lm_h0-0.11*z*(lm_h0-13)
     return lm_hz

def v_mstar_p18(lm_s0,lm_h0,z):
     lm_sz = lm_s0-0.325*np.log10(1+z)-0.6*z-0.19*z*(lm_h0-13)
     return lm_sz

##Ref. Moster et al. 2013
def v_mstar_moster(lm_hz,redshift):
     temp = redshift/(1.+redshift)
     lgM_1 = 11.590+1.195*temp
     N = 0.0351-0.0247*temp
     beta = 1.376-0.826*temp
     gamma = 0.608+0.329*temp
     temp = lm_hz-lgM_1
     lm_sz = np.log10(2.*N/(np.power(10.0,temp*(-beta))+np.power(10.0,temp*gamma)))+lm_hz
     return lm_sz

     
