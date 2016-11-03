import numpy as np
from astropy.cosmology import FlatLambdaCDM
from class_cosmos import GetCosmos

class SMHM_Sandy(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     def v_mstar_p3(self,lm_hz,z):
          lm_sz = 10.4+1.7*(lm_hz-12)-0.42*np.log10(1+z)
          return lm_sz

     def v_mstar_p18(self,lm_s0,lm_h0,z):
          lm_sz = lm_s0-0.325*np.log10(1+z)-0.6*z-0.19*z*(lm_h0-13)
          return lm_sz

##Ref. Moster et al. 2013
class SMHM_Moster(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     def v_mstar(self,lm_hz,z):          
          temp = z/(1.0+z)
          lgM_1 = 11.590+1.195*temp
          N = 0.0351-0.0247*temp
          beta = 1.376-0.826*temp
          gamma = 0.608+0.329*temp
          temp = lm_hz-lgM_1
          lm_sz = np.log10(2.*N/(np.power(10.0,temp*(-beta))+np.power(10.0,temp*gamma)))+lm_hz
          return lm_sz

class MHGrowth_Sandy(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     def v_mhalo_p3(self,lm_sz,z):
          lm_hz = 12.0+(lm_sz+0.42*np.log10(1+z)-10.4)/1.7
          return lm_hz

     def v_mhalo_p17(self,lm_h0,z):
          lm_hz = -0.35*z+lm_h0-0.11*z*(lm_h0-13)
          return lm_hz

##Ref. Fakhouri et al. 2010
class MHGrowth_Fakhouri(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     def v_mhalo_dot_mean(self,z,lm_hz):
          m_hz = np.power(10.0,lm_hz)     
          unit = self.get_unit(z)
          mhalo_dot_mean = 46.1*np.power((m_hz/10e12),1.1)*(1+1.11*z)*unit
          return mhalo_dot_mean

     def v_mhalo_dot_mid(self,z,lm_hz):
          m_hz = np.power(10.0,lm_hz)
          unit = self.get_unit()
          mhalo_dot_mid = 25.3*np.power((m_hz/10e12),1.1)*(1+1.65*z)*unit
          return mhalo_dot_mid

     def v_mhalo(self,lm_h0,z,form='mean',track='forward'):
          cosmo = FlatLambdaCDM(H0=self.h0*10, Om0=self.omega_m)
          m_halo = np.zeros_like(z)
          m_halo[0] = np.power(10,lm_h0)
          for i in range(np.size(z)-1):
               age_cos_1 = cosmo.age(z[i])
               age_cos_2 = cosmo.age(z[i+1])
               if form=='mean':
                    mhalo_dot = self.v_mhalo_dot_mean(z[i],np.log10(m_halo[i])) 
               elif form=='mid':
                    mhalo_dot = self.v_mhalo_dot_mid(z[i],np.log10(m_halo[i]))
               m_min = mhalo_dot*np.absolute(age_cos_1.value-age_cos_2.value)*10e9
               if track=='forward':
                    m_halo[i+1] = m_halo[i]+m_min
               elif track=='backward':
                    m_halo[i+1] = m_halo[i]-m_min
          return np.log10(m_halo)

##Ref. Aldo et al. 2016
class MHGrowth_Aldo(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     def v_aldo_a0(self,m_h0):
          a0 = 0.592-np.log10(0.113*(np.power(10,15.7)*self.h0/m_h0)+1)
          return a0
   
     def v_aldo_gfunc(self,m_h0,z):
          a0 = self.v_aldo_a0(m_h0)
          a = 1.0/(1.0+z)
          gfunc = 1+np.exp(-3.676*(a-a0))
          return gfunc
    
     def v_aldo_m13func(self,z):
          value = 13.6+np.log10(self.h0)+2.755*np.log10(1+z)-6.351*np.log10(1+0.5*z)+np.log10(np.exp(-0.413*z))
          return np.power(10,value)

     def v_aldo_ffunc(self,m_h0,z):
          m13_0 = self.v_aldo_m13func(0)
          gfunc_0 = self.v_aldo_gfunc(m_h0,0)
          gfunc_z = self.v_aldo_gfunc(m_h0,z)
          ffunc_z = np.log10(m_h0/m13_0)*(gfunc_0/gfunc_z)
          return ffunc_z

     def v_mhalo(self,lm_h0,z):
          m_h0 = np.power(10,lm_h0)
          m13_z = self.v_aldo_m13func(z)
          ffunc_z = self.v_aldo_ffunc(m_h0,z)
          lm_hz = np.log10(m13_z)+ffunc_z
          return lm_hz 
