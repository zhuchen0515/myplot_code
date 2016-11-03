import numpy as np
from const import *
from class_cosmos import GetCosmos

def v_mcool_p7(lm_hz,z):
     lmcool = 10+np.log10(1.9)+1.5*np.log10(fhot)+(5.0*(lm_hz+lsun)/12.0)+\
              (11.0*np.log10(1+z)/8.0)
     lmcool = lmcool-lsun+lyrtosec
     return lmcool

class HaloP(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)

     def __del__(self):
          pass

     ## The virial velocity in the uit of km/s
     def v_lVvir(self,lm_hz,z): 
          Hz = self.v_Hz(z)
          lVvir = 0.33333*(lm_hz+lsun+np.log10(Hz)+lh0_u+lG-7)-5
          return lVvir
     
     ## The virial radius in the unit of kpc
     def v_lRvir(self,lm_hz,z): 
          Hz = self.v_Hz(z)
          lRvir = 0.33333*(lm_hz+lsun-2.0*(np.log10(Hz)+lh0_u)+lG-10)-lkpc
          return lRvir
     

class FuncCool(GetCosmos):
     def __init__(self):
          GetCosmos.__init__(self)
          self.nu = 0.59259
          self.fhot = 0.1
          self.N_eff = 1.0 
          self.lVvir = np.log10(220)
          self.lRvir = np.log10(5)

     def __del__(self):
          pass

     ##mean molecular weight
     ##for Solar surface X=0.70, Y=0.28, Z=0.02
     def pass_nu(self,X=0.75,Y=0.25,Z=0.0):
         nu = 1.0/(2.0*X+0.75*Y+0.5*Z)
         self.nu = nu

     def pass_fhot(self,f_hot):
          self.fhot = f_hot

     ##the index of mass cooling rate
     def pass_Neff(self,N_eff):
          self.N_eff = N_eff  ##1.0 for Herquist15, 0.5 for Croton06
          
     ## unit is km/s
     def pass_lVvir(self,lVvir):
          self.lVvir = lVvir

     ## unit is kpc
     def pass_lRvir(self,lRvir):
          self.lRvir = lRvir

     ## the log virial temperature in the unit of K 
     ## Vvir here in the unit of km/s
     def v_lTemp(self,lm_hz,z):
          return np.log10(35.9)+2.0*(self.lVvir)

     ## the log cooling efficient in the unit of erg/cm^3/s
     def v_lcool(self,lm_hz,z):
          lTemp = self.v_lTemp(lm_hz,z)
          return -18-0.75*lTemp

     ## cooling radius: unit: kpc
     def v_lrcool(self,lm_hz,z,tdyn='Sandy'):
          lTemp = self.v_lTemp(lm_hz,z)
          if tdyn=='Sandy':
                lTdyn = np.log10(self.v_tdyn_sandy(z))-lh0_u
          elif tdyn=='Aldo':
                lTdyn = np.log10(self.v_tdyn(z))-lh0_u
          lcool = self.v_lcool(lm_hz,z)
          lrcool = 0.5*(np.log10(self.fhot)+lm_hz+lcool+lTdyn)-\
                   0.5*np.log10(6.0*PI*k_const*mp*self.nu)-0.5*(lTemp+self.lRvir+lkpc)+lsun*0.5-lkpc
          return lrcool

     ## The mass cooling rate in the unit of Msun/yr
     def v_lmcool(self,lm_hz,z,tdyn='Sandy',break_form=False):
          lmcool = np.zeros_like(lm_hz)
          lrcool = self.v_lrcool(lm_hz,z,tdyn=tdyn)
          if tdyn=='Sandy':
                lTdyn = np.log10(self.v_tdyn_sandy(z))-lh0_u
          elif tdyn=='Aldo':
                lTdyn = np.log10(self.v_tdyn(z))-lh0_u
          Hz = self.v_Hz(z)
          if break_form==False:
               lmcool = np.log10(self.N_eff)+np.log10(self.fhot)+lm_hz+lrcool-self.lRvir-(lTdyn-lyrtosec)
          elif break_form==True:
               for i in range(np.size(lm_hz)):
                    if lrcool[i]<=lRvir[i]:
                         lmcool[i] = np.log10(self.N_eff)+np.log10(self.fhot)+lm_hz[i]+lrcool[i]-self.lRvir[i]-(lTdyn-lyrtosec)
                    else:
                         lmcool[i] = np.log10(self.N_eff)+np.log10(self.fhot)+lm_hz[i]+np.log10(Hz[i])-(lTdyn-lyrtosec)
          return lmcool
     
