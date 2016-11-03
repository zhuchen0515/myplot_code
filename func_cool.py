import numpy as np
from scipy.misc import derivative
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM,z_at_value
from const import *
from func_mass import *

## calculate the delta_vir, the over density
## Bryan et al. 1998 Eq.6
def  v_deltavir(omega_m,omega_l,z):
     omega_z = np.power((1+z),3.0)*omega_m/(np.power((1+z),3.0)*omega_m+omega_l)
     min_x = 1-omega_z
     delta_vir = (18.0*np.pi*np.pi-82.0*min_x-39.0*min_x*min_x)/omega_z
     return delta_vir 

def v_tdyn(omega_m,omega_l,z):
     omega_z = np.power((1+z),3.0)*omega_m/(np.power((1+z),3.0)*omega_m+omega_l)
     delta_vir = v_deltavir(omega_m,omega_l,z)*omega_z
     a_term = np.power(3.0*delta_vir/(8.0*np.pi),-0.5)
     b_term = v_Hz(z,omega_m,omega_l)
     t_dyn = a_term/b_term
     return t_dyn

def v_floss(z,omega_m,omega_l):
     cosmo = FlatLambdaCDM(H0=70, Om0=omega_m)
     age_cos = cosmo.age(z)
     f_loss = 0.05*np.log(1+age_cos.value*10e3/1.4)
     return f_loss

def v_Mvir_dot(z,lm_h0,omega_m,omega_l):
     cosmo = FlatLambdaCDM(H0=70, Om0=omega_m)
     age_cos = cosmo.age(z)
     age_dyn = np.log10(v_tdyn(omega_m,omega_l,z))-lh0_u-lyrtosec-9  ##unit: Gyr
     del_age = age_cos.value-np.power(10.0,age_dyn)
     z2 = np.zeros_like(del_age)
     for i in range(np.size(del_age)):
          z2[i] = z_at_value(cosmo.age,del_age[i]*u.Gyr)
     Mvir_z = v_mhalo_p17(lm_h0,z)
     Mvir_z2 = v_mhalo_p17(lm_h0,z2)
     del_Mvir = np.power(10.0,Mvir_z)-np.power(10.0,Mvir_z2)
     Mvir_dot = np.log10(del_Mvir)-age_dyn-9
     return Mvir_dot

#def v_Mstar_dot(lm_hz):
#     Mstar_dot = derivative(v_mstar_moster,lm_hz,dx=1e-4)
#     return Mstar_dot
     
## the log virial temperature in the unit of K 
## Vvir here in the unit of km/s
def v_lTemp(lVvir):  ##Fun.4
     return np.log10(35.9)+2.0*(lVvir)

## the log cooling efficient in the unit of erg/cm^3/s
def v_lcooling(lTemp):  ##Fun.6
     return -18-0.75*lTemp

## the evolution of H(z) in the unit of km/s/Mpc
def v_Hz(z,omega_m,omega_l):  ##Fun.3
     Hz = H0*np.power(omega_m*np.power(1+z,3.0)+omega_l,0.5)
     return Hz

## The virial velocity in the uit of km/s
def v_Vvir(lm_hz,Hz):  ##Fun.2
     lVvir = 0.33333*(lm_hz+lsun+np.log10(Hz)+lh0_u+lG-7)-5 
     return lVvir

def v_Vvir_2(lm_hz,lRvir):
     lVvir = 0.5*(lm_hz+lsun+lG-8-(lRvir+lkpc))
     return lVvir-5

## The virial radius in the unit of kpc
def v_Rvir(lm_hz,Hz):  ##Fun.1
     lRvir = 0.33333*(lm_hz+lsun-2.0*(np.log10(Hz)+lh0_u)+lG-10)-lkpc
     return lRvir

## The mass cooling rate in the unit of Msun/yr
def v_mcool(lm_hz,lRvir,lrcool,Hz,break_form=False):
     lmcool = np.zeros_like(lm_hz)
     if break_form==False:
          lmcool = np.log10(N_eff)+np.log10(fhot)+lm_hz+lrcool-lRvir+\
                        np.log10(Hz)+np.log10(1.022)-10
     elif break_form==True:
          for i in range(np.size(lm_hz)):
               if lrcool[i]<=lRvir[i]:
                    lmcool[i] = np.log10(N_eff)+np.log10(fhot)+lm_hz[i]+lrcool[i]-lRvir[i]+\
                                np.log10(Hz[i])+np.log10(1.022)-10
               else:
                    lmcool[i] = np.log10(N_eff)+np.log10(fhot)+lm_hz[i]+np.log10(Hz[i])+\
                                np.log10(1.022)-10
     return lmcool

## cooling radius: unit: kpc
def v_rcool(lm_hz,lVvir,lTemp,lcool,nu):
     lrcool = 0.5*(np.log10(fhot)+lm_hz+lcool-lTemp-lVvir)-\
              0.5*np.log10(6.0*PI*k_const*mp*nu)+(lsun-5)*0.5-lkpc
     return lrcool

## The mass cooling rate from Sandy
##p7 is too low
def v_mcool_p7(lm_hz,z):
     lmcool = 10+np.log10(1.9)+1.5*np.log10(fhot)+(5.0*lm_hz/12.0)+\
              (11.0*np.log10(1+z)/8.0)
     return lmcool

def v_mcool_p10(lm_sz,z):
     lmcool = np.log10(63.7)+1.5*np.log10(fhot)+0.245*lm_sz+1.479*np.log10(1+z)
     return lmcool
