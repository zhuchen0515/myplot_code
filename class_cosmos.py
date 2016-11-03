import numpy as np

class GetCosmos(object):
     def __init__(self):
          self.omega_m = 0.3
          self.omega_l = 0.7
          self.omega_b = 0.048
          self.h0 = 0.7

     def __del__(self):
          pass

     def pass_cosmos(self,omega_m,omega_l):
          self.omega_m = omega_m
          self.omega_l = omega_l

     def pass_baryon(self,omega_b):
          self.omega_b = omega_b

     def pass_hubble(self,h):
          self.h0 = h

     def get_unit(self,z):
          a = np.sqrt(self.omega_l+self.omega_m*np.power(1+z,3))
          return a

     ## omega_m at z
     def omega_z(self,z):
          unit = self.get_unit(z)
          omega_z = 1-self.omega_l/(unit*unit)
          return omega_z

     ## the evolution of H(z) in the unit of km/s/Mpc
     def v_Hz(self,z):
          unit = self.get_unit(z)
          Hz = self.h0*100*unit
          return Hz

     def  v_deltavir(self,z):
          omega_z = self.omega_z(z)
          x = omega_z-1
          delta_vir = (18.0*np.pi*np.pi+82.0*x-39.0*x*x)/omega_z
          return delta_vir

     def v_tdyn(self,z):
          omega_z = self.omega_z(z)
          delta_vir = self.v_deltavir(z)*omega_z
          a_term = np.power(3.0*delta_vir/(8.0*np.pi),-0.5)
          b_term = self.v_Hz(z)
          t_dyn = a_term/b_term
          return t_dyn

     def v_tdyn_sandy(self,z):
          Hz = self.v_Hz(z)
          tdyn_sandy = 0.1/Hz
          return tdyn_sandy
