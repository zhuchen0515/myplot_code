import numpy as np
from const import *

## formula from Navarro et al. 1997  unit: Kpc
def lomega_z(z):
     lomega = np.log10(omega_l+omega_m*np.power(1+z,3))-3*np.log10(1+z)
     return lomega

def v_Rvir_navarro(lm_hz,z):
     lomega = lomega_z(z)
     lRvir_navarro = 0.33333*(lG-8+lm_hz+np.log10(h0)\
                     -lomega-3.0*np.log10(1+z))-np.log10(h0)
     return lRvir_navarro

def v_Vvir_navarro(lRvir,z):
     lomega = lomega_z(z)
     lVvir_navarro = lRvir+0.5*lomega+1.5*np.log10(1+z)+np.log10(h0)
     return lVvir_navarro

## the evoltion of log H(z) in the unit of cm/g/s
def v_lhz(z):
     lhz = lh0+0.5*np.log10(omega_m*np.power(1+z,3.0)+omega_l)
     return lhz

