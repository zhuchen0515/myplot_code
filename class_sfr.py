from scipy.misc import derivative
import numpy as np
from const import *


class get_Mstar_dot(object):

     def __init__(self):
          self.redshift = 0.0

     def __del__(self):
          pass

     def pass_redshift(self,redshift):
          self.redshift = redshift

     def v_mstar_moster(self,lm_hz):
          temp = self.redshift/(1.+self.redshift)
          lgM_1 = 11.590+1.195*temp
          N = 0.0351-0.0247*temp
          beta = 1.376-0.826*temp
          gamma = 0.608+0.329*temp
          temp = lm_hz-lgM_1
          lm_sz = np.log10(2.*N/(np.power(10.0,temp*(-beta))+np.power(10.0,temp*gamma)))+lm_hz
          return lm_sz

     ## dlogM_star/dlogM_halo
     def v_eta(self,lm_hz):
          eta = derivative(self.v_mstar_moster,lm_hz,dx=1e-4)
          return eta

