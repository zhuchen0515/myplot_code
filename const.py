##pi value
PI = 3.1415926

#mass omega_m and dark matter omega_b
#omega_m = 0.3  ##matter
#omega_l = 0.7  ##vacuum

##the index of mass cooling rate
N_eff = 1.0 ##1.0 for Herquist15, 0.5 for Croton06

##heating fraction of halo
fhot = 0.1

##light speed km/s
light = 3.0e5

#log grav. constant G (cm/g/s)
#lG = np.log10(6.674)-8 = -7.1756
lG = 0.824393   ## the value of log10(6.6741)

#log hubble constant today (cm/g/s) H0=70
#transfer unit to s^-1
#lh0 = np.log10(2.26)-18
lh0 = -17.6459

##hubble constant 
h0 = 0.7
H0 = 70

#h0 cm/g/s to km/s/Mpc
#lh0_u = np.log10(2.26/70)-18
lh0_u = -19.489351

#mpc to cm
#lmpc = 24+np.log10(3.085678)
lmpc = 24.4894

##kpc to cm
lkpc = 21.4894

#km/s to cm/s
lkmsec = 5

##transfer yr to sec
yrtosec = 3.1536e7
lyrtosec = 7.4988

#solar mass to g
#lsun = np.log10(1.9891)+33
lsun = 33.2986

##mass of proton unit : g
mp = 1.672622e-24
#mp = 1.672622e-24
lmp = 0.2234-24

##the Boltzmann const erg/K (GS system)
k_const = 1.380648e-16
lk_const = 0.140083-16

##mean molecular weight
##for Solar surface X=0.70, Y=0.28, Z=0.02
def get_nu(X=0.75,Y=0.25,Z=0.0):
    nu = 1.0/(2.0*X+0.75*Y+0.5*Z)
    return nu    

