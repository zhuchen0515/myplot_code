from astropy.table import Table
import numpy as np
from func_heat import *
from func_cool import *

def get_data(table,lr9,omega_m=0.3,omega_l=0.7):
     z = table['redshift']
     lm_hz = table['m_halo']
     lm_sz = table['m_star']

     ##cal the cooling rate
     Hz = v_Hz(z,omega_m,omega_l)
     lRvir = v_Rvir(lm_hz,Hz)
     lVvir = v_Vvir(lm_hz,Hz)
     lTemp = v_lTemp(lVvir)
     lcool = v_lcooling(lTemp)
     nu = get_nu(X=0.75,Y=0.25,Z=0.0)
     lrcool = v_rcool(lm_hz,lVvir,lTemp,lcool,nu)
     lm_cool = v_mcool(lm_hz,lRvir,lrcool,Hz)
     lm_cool_p10 = v_mcool_p10(lm_sz,Hz)
     
     lr_z = v_rz(lr9,lm_sz)
     lsig1 = v_sig1(lr_z,lm_sz)
     lm_bh = v_mbh(lsig1,lr_z)
     lm_heat_old = v_mheat_old(lm_bh,lVvir)
     lm_heat_new = v_mheat_new(lm_sz,lVvir,lr_z)

     Rvir = np.power(10,lRvir)
     Vvir= np.power(10,lVvir)
     Rcool = np.power(10,lrcool)
     Hz = np.around(Hz,3)

     data = Table([z,table['m_star'],table['m_halo'],lr_z,
                    lm_cool,lm_heat_old,lm_heat_new,lm_cool_p10,
                    Hz,Rvir,Vvir,lTemp,lcool,Rcool,
                    lsig1,lm_bh],
                   names=('redshift','m_star','m_halo','lre',
                          'lm_cool','lm_heat_old','lm_heat_new','lm_cool_p10',
                          'Hz','R_vir','V_vir','logT','log(cooling)','Rcool',
                          'logsig1','lm_bh'))
     return data

